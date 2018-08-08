SESSIONID = 17; %Session to view examples


%%
sessionList = readtable('../sessionList.csv','FileType','text','Delimiter',',');
areaCols = [      0    0.4470    0.7410; %blue
    0.3010    0.7450    0.9330; %lightblue
    0.4660    0.6740    0.1880; %green
    0.8500    0.3250    0.0980; %orange
    0.4940    0.1840    0.5560; %purple
    0.8 0.8 0.8]; %grey

cDiffCols = [1 0 0 1;
    1 0 0 0.66;
    1 0 0 0.33;
    0 0 0 0.33;
    0 0 1 0.33;
    0 0 1 0.66;
    0 0 1 1];

contrastCols = [1 0 0; 0 0 0; 0 0 1];

responseCols = copper(3);

plotTypes = {'Stim-aligned, split by stimulus CDIFF';
    'Stim-aligned, split by stimulus';
    'Stim-aligned, split by choice';
    'Move-aligned, split by stimulus';
    'Move-aligned, split by choice'};
plotColours = {cDiffCols; contrastCols; responseCols; contrastCols; responseCols};


%% haemodynamic correction & load SVD
num_svd_components = 500;
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    wfFile = [ '../preproc/WF_SVD/' eRef '.mat'];
    
    if ~exist(wfFile,'file')
        corrPath = fullfile(sessionList.widefieldDir{sess}, 'svdTemporalComponents_corr.npy');
        if ~exist(corrPath,'file') %Do haemodynamic correction if it doesnt exist
            quickHemoCorrect(sessionList.widefieldDir{sess},num_svd_components);
        end
        [U,V,wfTime,meanImg]=quickLoadUVt(sessionList.widefieldDir{sess},num_svd_components);
        
        save([ '../preproc/widefield/' eRef '.mat'],'U','V','wfTime','meanImg');
    end
    
end

%% Get activity and behavioural data, align timestamps, compute traces
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    preprocFile = [ '../preproc/WF_ROI/' eRef '.mat'];
    
    if ~exist(preprocFile,'file') & ~isnan(sessionList.row_VISp(sess))
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Load widefield data %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        wf = [ '../preproc/WF_SVD/' eRef '.mat'];
        load(wf,'U','V','wfTime','meanImg');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Compute dF/F %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [U_dff, V_dff] = dffFromSVD(U, V, meanImg);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Compute derivative of dF/F %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        d = designfilt('differentiatorfir','FilterOrder',50, ...
            'PassbandFrequency',8.5,'StopbandFrequency',10, ...
            'SampleRate',1/mean(diff(wfTime)));
        dV = filter(d,V_dff')'; %*Fs;
        delay = mean(grpdelay(d));
        wfTime_dt = wfTime(1:end-delay);
        dV(:,1:delay) = [];
        wfTime_dt(1:delay) = [];
        dV(:,1:delay) = [];
        % fvtool(d,'MagnitudeDisplay','zero-phase','Fs',1/mean(diff(wfTime)))
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Load behavioural data %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [D,meta] = loadData(eRef);
        blockFile = meta.blockFile;
        tlFile = strrep(blockFile,'_Block.mat','_Timeline.mat');
        [mouseName, thisDate, expNum]=dat.parseExpRef(eRef);
        load(blockFile);
        load(tlFile);
        tr = block.trial; tr = tr(1:block.numCompletedTrials);
        trial = struct;
        trial.cond = [tr.condition];
        vcc = [trial.cond.visCueContrast];
        trial.contrastLeft = vcc(1,:);
        trial.contrastRight = vcc(2,:);
        trial.choice = [tr.responseMadeID];
        trial.feedback = [tr.feedbackType];
        trial.repNum = [trial.cond.repeatNum];
        [trial.contrastCondDefs, trial.ia, trial.contrastCondsRaw] = unique(vcc', 'rows');
        timings = struct;
        timings.b_stimOn = [tr.stimulusCueStartedTime];
        timings.b_goCue = [tr.interactiveStartedTime];
        timings.b_responseTime = [tr.responseMadeTime];
        timings.b_trialStarts = [tr.trialStartedTime];
        timings.b_trialEnds = [tr.trialEndedTime];
        timings.b_stimWindowUpdateTimes = block.stimWindowUpdateTimes;
        trial.reactionTime = [tr.responseMadeTime]-timings.b_goCue;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Align block to timeline %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pd = Timeline.rawDAQData(:,2);
        tt = Timeline.rawDAQTimestamps;
        pdFlips = schmittTimes(tt,pd, [3 5]); %Get the times when the photodiode flips
        [blockToTL,blockToTL_weights] = makeCorrection(pdFlips(2:end-1), timings.b_stimWindowUpdateTimes, false);
        timings.t_stimOn = blockToTL(timings.b_stimOn);
        timings.t_goCue = blockToTL(timings.b_goCue);
        timings.t_responseTime = blockToTL(timings.b_responseTime);
        timings.t_trialStarts = blockToTL(timings.b_trialStarts);
        timings.t_trialEnds = blockToTL(timings.b_trialEnds);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Wheel trace movement detection %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tt = Timeline.rawDAQTimestamps;
        wh = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder')); %Wheel encoder trace
        wh = wh(tt<max(timings.t_stimOn)+5); tt = tt(tt<max(timings.t_stimOn)+5); %Trim away wheel trace which extended 5 seconds after the last stimulus onset time
        wh = wheel.correctCounterDiscont(wh); %Correct artifactual values in wheel trace.
        [timings.t_moveOnsets, timings.t_moveOffsets] = wheel.findWheelMoves3(wh, tt, 1000, []);
        resp = trial.choice; hasTurn = resp==1|resp==2; resp = resp(hasTurn);
        mon=nan(size(timings.t_stimOn));
        for tr = 1:length(timings.t_stimOn)
            idx = find(timings.t_stimOn(tr) < timings.t_moveOnsets & timings.t_moveOnsets < timings.t_responseTime(tr),1);
            if trial.choice(tr)<3 && ~isempty(idx)
                mon(tr) = timings.t_moveOnsets( idx );
            end
        end
        rt = [mon-timings.t_stimOn]';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Load ROI definitions %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear px;
        px(1).name = 'VISp'; px(1).xy = [sessionList.row_VISp(sess) sessionList.col_VISp(sess)];
        px(2).name = 'VISlm'; px(2).xy = [sessionList.row_VISlm(sess) sessionList.col_VISlm(sess)];
        px(3).name = 'MOs'; px(3).xy = [sessionList.row_MOs(sess) sessionList.col_MOs(sess)];
        px(4).name = 'MOp'; px(4).xy = [sessionList.row_MOp(sess) sessionList.col_MOp(sess)];
        px(5).name = 'SSp'; px(5).xy = [sessionList.row_S1(sess) sessionList.col_S1(sess)];
        %         px(6).name = 'RSP'; px(6).xy = [sessionList.row_RSP(sess) sessionList.col_RSP(sess)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% INCLUSION CRITERIA FOR TRIALS %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        isDetectionTrial = (trial.contrastRight==0 | trial.contrastLeft==0);
        isGoodPerformance = trial.repNum==1 & trial.feedback==1;
        isGoodWheelMove = trial.choice==3 | isnan(rt) | (rt>0.125 & rt<0.5);
        inclTrials = isDetectionTrial & isGoodWheelMove & isGoodPerformance;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Save %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        save(preprocFile,'eRef','px','wfTime_dt','U_dff','dV','meanImg','mon','trial','timings','inclTrials');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% COMPUTE STIM-ALIGNED TRACES & ONSET LATENCIES %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        aligned_analysis_file = [ '../preproc/WF_aligned/' eRef '.mat'];
        [~, winSamps, periEventdV, label] = ...
            eventLockedAvgSVD(U_dff, dV, wfTime_dt,...
            timings.t_stimOn(inclTrials), ...
            sign(trial.contrastRight(inclTrials)-trial.contrastLeft(inclTrials)), ...
            [-0.3 0.8]);
        numTrials = length(label);
        avgStimAlign = nan(numTrials, length(winSamps), length(px));
        contraLatency = nan(numel(px),1);
        for p =1:numel(px)
            thisU = squeeze(U_dff(px(p).xy(1), px(p).xy(2), :));
            for n = 1:numTrials
                avgStimAlign(n,:,p) = thisU'*squeeze(periEventdV(n,:,:));
            end
            
            %onset latency
            avgcontra = mean( avgStimAlign(label==-1,:,p), 1);
            normcontra = (avgcontra/max(avgcontra)).*(winSamps>0.015);
            contraLatency(p) = winSamps(find(normcontra>0.25,1,'first'));
        end
        avgStimAlign(avgStimAlign<0) = 0; % rectify
        % subtract average pre-stim activity
        avgStimAlign = bsxfun(@minus, avgStimAlign, mean(avgStimAlign(:,winSamps<0,:),2));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Save %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        save(aligned_analysis_file,'winSamps','avgStimAlign','contraLatency','label','px');
        
    end
end

%% 0) Load and prepare data for the example session
eRef = sessionList.expRef{SESSIONID};
fprintf('Session %d %s\n',SESSIONID,eRef);
load([ '../preproc/WF_ROI/' eRef '.mat'],'eRef','px','wfTime_dt','U_dff','dV','meanImg','mon','trial','timings','inclTrials');
load([ '../preproc/WF_mask/' eRef '.mat'],'mask');
load([ '../preproc/bregma/' eRef '_bregma.mat'],'bregma_AP','bregma_ML','lambda_AP','lambda_ML');
load([ '../preproc/WF_aligned/' eRef '.mat'],'winSamps','avgStimAlign','contraLatency','label','px');

disp('Data loaded');


%% 1) plot map aligned to stim onset

figure(101); set(gcf,'color','w');
[avgdV, wS, ~, ~] = ...
    eventLockedAvgSVD(U_dff, dV, wfTime_dt,...
    timings.t_stimOn(inclTrials), ...
    sign(trial.contrastRight(inclTrials)-trial.contrastLeft(inclTrials)), ...
    [0 0.25]);
avgdV = avgdV - avgdV(:,:,1);%remove baseline
sliceIdx = 5:15:length(wS);
ha = tight_subplot(3,length(sliceIdx),0.01,[0.7 0.05],0.01);
for i = 1:length(sliceIdx)
    imL = svdFrameReconstruct(U_dff, permute(avgdV(1,:,sliceIdx(i)),[2 3 1]) ).*mask;
    imagesc( ha(i), imL );
    addAllenCtxOutlines([bregma_AP bregma_ML], [lambda_AP lambda_ML], [1 1 1]*1, ha(i))
    hold(ha(i),'on'); plot(ha(i),bregma_ML,bregma_AP,'k.')
    %overlay mask on the edges
    im2=imagesc(ha(i), zeros(size(imL)));
    im2.AlphaData=~mask;
    
    imR = svdFrameReconstruct(U_dff, permute(avgdV(3,:,sliceIdx(i)),[2 3 1]) ).*mask;
    imagesc( ha(i + length(sliceIdx) ), imR );
    addAllenCtxOutlines([bregma_AP bregma_ML], [lambda_AP lambda_ML], [1 1 1]*1,ha(i + length(sliceIdx) ))
    hold(ha(i + length(sliceIdx) ),'on'); plot(ha(i + length(sliceIdx) ),bregma_ML,bregma_AP,'k.')
    im2=imagesc(ha(i), zeros(size(imL)));
    im2.AlphaData=~mask;
    
    imagesc( ha(i + 2*length(sliceIdx) ), imL-imR );
    addAllenCtxOutlines([bregma_AP bregma_ML], [lambda_AP lambda_ML], [1 1 1]*1,ha(i + 2*length(sliceIdx) ))
    hold(ha(i + 2*length(sliceIdx) ),'on'); plot(ha(i + 2*length(sliceIdx) ),bregma_ML,bregma_AP,'k.')
    im2=imagesc(ha(i), zeros(size(imL)));
    im2.AlphaData=~mask;
    
    title(ha(i), wS(sliceIdx(i)));
end
cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
colormap(flipud(cmap));
set(ha,'clim',[-1 1]*0.008);
set(ha,'xcolor','none','ycolor','none','dataaspectratio',[1 1 1]);

%% 2) Plot traces and ROIs
figure(101);
ha = tight_subplot(5,1,0.01,[0.05 0.35],[0.05 0.66]);
thisTraces_contra = avgStimAlign(label==-1,:,:);
thisTraces_ipsi = avgStimAlign(label==1,:,:);
for p =1:numel(px)
    hold(ha(p),'on');
    plot(ha(p), winSamps, thisTraces_contra(:,:,p), 'k-', 'color', [areaCols(p,:) 0.1] );
    plot(ha(p), winSamps, mean(thisTraces_contra(:,:,p),1), '-', 'linewidth',3, 'color', areaCols(p,:));
    
    %Plot ipsi trace
    %     plot(ha(p), winSamps, -thisTraces_ipsi(:,:,p), 'k-', 'color', [areaCols(p,:) 0.1] );
    plot(ha(p), winSamps, mean(thisTraces_ipsi(:,:,p),1), '-', 'linewidth',1, 'color', areaCols(p,:));
    set(ha(p),'xcolor','none','ycolor','none');
    line(ha(p), [0 0], [-1 1]*10, 'color', [0 0 0 0.1],'LineStyle','-');
    
    %Add text in top left corner
    tx=text(ha(p),-0.3,0.01,px(p).name,'VerticalAlignment','top','color',areaCols(p,:),'FontWeight','bold');
end
line(ha(end), [0 0], [0 1], 'color', [0 0 0 0.1],'LineStyle','-');
line(ha(end), [0 0], [0 1], 'color', [0 0 0 0.1],'LineStyle','-');
linkaxes(ha,'y');
ylim(ha(1),[-0.001 0.01]);
set(ha(end),'xcolor','k','XTickLabelMode','auto','ycolor','k','yTickLabelMode','auto');

linkaxes(ha,'x');
xlim(ha(end),[-0.3 0.8]);
xlabel(ha(end),'Stim onset');
ylabel(ha(end),'d/dt ( dF/F )');

%meanimg and normalised traces
ha = tight_subplot(3,1,0.01,[0.05 0.35],[0.35 0.35]);
imagesc(ha(1),meanImg.*mask); colormap(ha(1),'gray'); hold(ha(1),'on');
addAllenCtxOutlines([bregma_AP bregma_ML], [lambda_AP lambda_ML], [1 1 1]*0.5, ha(1))
plot(ha(1),bregma_ML,bregma_AP,'k.')
for p =1:numel(px)
    h = plot(ha(1),px(p).xy(2), px(p).xy(1), '.', 'Color', areaCols(p,:), 'markersize',25);
end
set(ha(1),'xcolor','none','ycolor','none', 'dataaspectratio', [1 1 1]);

hold(ha(2),'on'); hold(ha(3),'on');
for p =1:numel(px)
    contraNormTrace = mean(thisTraces_contra(:,:,p),1)/max(mean(thisTraces_contra(:,:,p),1));
    hxx=plot(ha(2), winSamps, contraNormTrace, '-','linewidth',3, 'color', areaCols(p,:));
    uistack(hxx,'bottom');
    
    ipsiNormTrace = mean(thisTraces_ipsi(:,:,p),1)/max(mean(thisTraces_contra(:,:,p),1));
    hxx=plot(ha(3), winSamps, ipsiNormTrace, '-','linewidth',1, 'color', areaCols(p,:));
    uistack(hxx,'bottom');
    
end
line(ha(2), [0 0], [0 1], 'color', [0 0 0 0.1],'LineStyle','-');
line(ha(3), [0 0], [0 1], 'color', [0 0 0 0.1],'LineStyle','-');
set(ha(2:3),'ylim',[-0.1 1.1],'xlim',[-0.3 0.8]);

%% 3) Plot latencies over all sessions
mice={'Chain','Radnitz','Cori','Reichstein','Hench'};

figure(101);
ContraOnsetTime = nan(height(sessionList),5);
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    try
    load([ '../preproc/WF_aligned/' eRef '.mat'],'contraLatency');
    ContraOnsetTime(sess,:) = contraLatency';
    catch
    end
end

ha = tight_subplot(1,1,0.01,[0.05 0.35],[0.67 0.01]);
hold(ha(1),'on');
Y_jitter = linspace(-0.2,0.2,length(mice));
allMeans=[];
allMeansAx=gobjects;

for m = 1:length(mice)
    mIdx = strcmp(sessionList.mouseName,mice(m));
    lat = ContraOnsetTime(mIdx,:);
    
    for p =1:5
        m1=mean(lat(:,p));
        stdr1 = std(lat(:,p))/sqrt(length(lat(:,p)));
        
        yline = Y_jitter(m) + [-1 1]*0.05;
        
        line(ha(1),[1 1]*m1, yline,'Color',areaCols(p,:), 'linewidth',3);

        fx=fill(ha(1), [m1-stdr1, m1+stdr1, m1+stdr1, m1-stdr1], [yline(1) yline(1) yline(2) yline(2)], 'k',...
                            'EdgeAlpha',0,'FaceColor',areaCols(p,:),'FaceAlpha',0.5);
        
%         lx=line(ha(1), m1 + [-1 1]*stdr1, [1 1]*Y_jitter(m) + (p-1)*0.001, 'linewidth',3 );
%         lx.Color=areaCols(p,:);
%         
%         allMeansAx(m,p)=line(ha(1),[1 1]*m1, Y_jitter(m)+(p-1)*0.001 + [-1 1]*0.05,'Color',areaCols(p,:), 'linewidth',3);
        %         allMeans(m,p) = m1;
    end
    
end

xlim(ha,[0 0.2]);
xlabel(ha(1),'Stim onset');
set(ha,'ytick',Y_jitter,'YTickLabel',mice,'XTickLabelMode','auto');
ylim(ha,[-0.3 0.3]);

%% MISC: Plot ROI locations
figureROILocations = figure('color','w','name','ROI locations');
haROI = tight_subplot(5,8,0.04,[0.01 0.05],[0.05 0.01]);

for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    preprocFile = [ '../preproc/WF_ROI/' eRef '.mat'];
    
    if exist(preprocFile,'file')
        load(preprocFile,'px');
        load([ '../preproc/bregma/' eRef '_bregma.mat'],'bregma_AP','bregma_ML','lambda_AP','lambda_ML');
        load([ '../preproc/WF_mask/' eRef '.mat'],'mask');
        
        %Plot ROI
        load([ '../preproc/WF_SVD/' eRef '.mat'],'meanImg');
        imagesc(haROI(sess),meanImg); colormap(haROI(sess),'gray');
        addAllenCtxOutlines([bregma_AP bregma_ML], [lambda_AP lambda_ML], [1 1 1]*0.5, haROI(sess))
        
        hold(haROI(sess),'on');
        for p =1:numel(px)
            h = plot(haROI(sess),px(p).xy(2), px(p).xy(1), 'o', 'color', areaCols(p,:));
            set(h, 'MarkerFaceColor', get(h, 'Color'));
        end
        title(haROI(sess),eRef,'interpreter','none');
        set(haROI(sess),'xcolor','none','ycolor','none');
        
        %overlay mask on the edges
        im2=imagesc(haROI(sess), ones(size(meanImg))*max(max(meanImg)));
        im2.AlphaData=~mask;
    end
end

%% MISC: Define imaging mask
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    preprocFile = [ '../preproc/WF_ROI/' eRef '.mat'];
    maskFile = [ '../preproc/WF_mask/' eRef '.mat'];
    if exist(preprocFile,'file')
        
        %Plot ROI
        clear mask meanImg;
        load(preprocFile,'meanImg');
        load(maskFile,'mask');
        
        f=figure;
        imagesc(meanImg);
        
        if ~exist('mask','var')
            mask = roipoly;
            save(preprocFile,'mask','-append')
        end
        imagesc(meanImg.*mask);
        
        save(maskFile,'mask')
    end
end

%% MISC: define bregma and lambda
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    preprocFile = [ '../preproc/WF_ROI/' eRef '.mat'];
    bregmaFile = [ '../preproc/bregma/' eRef '_bregma.mat'];
    
    if exist(preprocFile,'file')
        load(preprocFile,'px');
        
        %Plot ROI
        load([ '../preproc/WF_SVD/' eRef '.mat'],'meanImg');
        
        figure;
        imagesc(meanImg); colormap(gray);
        
        hold on;
        for p =1:numel(px)
            h = plot(px(p).xy(2), px(p).xy(1), 'o', 'color', areaCols(p,:));
            set(h, 'MarkerFaceColor', get(h, 'Color'));
        end
        
       
        [X,Y] = ginput(2);
        bregma_AP = Y(1);
        lambda_AP = Y(2);
        bregma_ML = X(1);
        lambda_ML = X(2);
        addAllenCtxOutlines([bregma_AP bregma_ML], [lambda_AP lambda_ML], [1 1 1]*0.5, gca);
        outlines = findobj(gca,'tag','outline');
        
        go=1;
        while go == 1
            k = waitforbuttonpress;
            % 28 leftarrow
            % 29 rightarrow
            % 30 uparrow
            % 31 downarrow
            value = double(get(gcf,'CurrentCharacter'));
            shiftX = 0; shiftY = 0;
            switch(value)
                case 28
                    shiftX = -1;
                case 29
                    shiftX = +1;
                case 30
                    shiftY = -1;
                case 31
                    shiftY = +1;
                otherwise
                    go = 0; 
                    save(bregmaFile,'bregma_AP','bregma_ML','lambda_AP','lambda_ML');
            end

            for i = 1:length(outlines)
                outlines(i).YData = outlines(i).YData + shiftY;
                outlines(i).XData = outlines(i).XData + shiftX;
            end
            
            bregma_AP = bregma_AP + shiftY;
            lambda_AP = lambda_AP + shiftY;
            bregma_ML = bregma_ML + shiftX;
            lambda_ML = lambda_ML + shiftX;
            
            
        
            
        end
        
    end
end

%% Plot sequence at contralateral contrast
figureContraTrace = figure('color','w','name','Sequence to CONTRALATERAL stimulus');
haContraTrace = tight_subplot(5,8,0.03,[0.01 0.05],[0.05 0.01]);
figureIpsiTrace = figure('color','w','name','Sequence to IPSILATERAL stimulus');
haIpsiTrace = tight_subplot(5,8,0.03,[0.01 0.05],[0.05 0.01]);
ContraOnsetTime = nan(height(sessionList),6);
IpsiOnsetTime = nan(height(sessionList),6);
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    preprocFile = [ '../preproc/WF_ROI/' eRef '.mat'];
    
    if exist(preprocFile,'file')
        %         load(preprocFile,'px','winSamps','thisTraces','avgStimAlign','CR_minus_CL');
        load(preprocFile,'px','winSamps','avgStimAlign','label');
        px(end)=[]; %REMOVE RSP
        
        contraIdx = label{2}==-1;
        ipsiIdx = label{2}==1;
        thisTraces_contra = avgStimAlign{2}(contraIdx,:,:);
        thisTraces_ipsi = avgStimAlign{2}(ipsiIdx,:,:);
        
        %rectify
        thisTraces_contra(thisTraces_contra<0)=0;
        thisTraces_ipsi(thisTraces_ipsi<0)=0;
        
        %subtract baseline
        thisTraces_contra = bsxfun(@minus, thisTraces_contra, mean(thisTraces_contra(:,winSamps{2}<0,:),2) );
        thisTraces_ipsi = bsxfun(@minus, thisTraces_ipsi, mean(thisTraces_ipsi(:,winSamps{2}<0,:),2) );
        
        hold(haContraTrace(sess),'on');
        hold(haIpsiTrace(sess),'on');
        for p = 1:length(px)
            normContraTrace = mean(thisTraces_contra(:,:,p), 1)/max(mean(thisTraces_contra(:,:,p), 1));
            %             normIpsiTrace =  mean(thisTraces_ipsi(:,:,p), 1)/max(mean(thisTraces_contra(:,:,p), 1));
            normIpsiTrace =  mean(thisTraces_ipsi(:,:,p), 1)/max(mean(thisTraces_ipsi(:,:,p), 1));
            
            
            plot(haContraTrace(sess),winSamps{2}, normContraTrace, 'Color',areaCols(p,:),'Tag','trace');
            plot(haIpsiTrace(sess),winSamps{2}, normIpsiTrace, 'Color',areaCols(p,:),'Tag','trace');
            
            
            %Detect onset and plot latency
            normZerod = normContraTrace.*(winSamps{2}>0.015);
            ContraOnsetTime(sess,p) = winSamps{2}(find(normZerod>0.25,1,'first'));
            hx=plot(haContraTrace(sess),ContraOnsetTime(sess,p),1,'.','color',areaCols(p,:),'markersize',10);
            
            %Detect onset and plot latency
            try
                normZerod = normIpsiTrace.*(winSamps{2}>0.015);
                IpsiOnsetTime(sess,p) = winSamps{2}(find(normZerod>0.25,1,'first'));
                hx=plot(haIpsiTrace(sess),IpsiOnsetTime(sess,p),1,'.','color',areaCols(p,:),'markersize',10);
            catch
            end
        end
        title(haContraTrace(sess),eRef,'interpreter','none');
        title(haIpsiTrace(sess),eRef,'interpreter','none');
        plot(haContraTrace(sess),[-0.1 0.4], [0 0], 'k--');
        plot(haContraTrace(sess),[0 0],[-0.1 1],  'k--');
        
        plot(haIpsiTrace(sess),[-0.1 0.4], [0 0], 'k--');
        plot(haIpsiTrace(sess),[0 0],[-0.1 1],  'k--');
    end
end
linkaxes([haContraTrace; haIpsiTrace],'xy');

set(haContraTrace(1),'xlim',[-0.1 0.4],'ylim',[-0.1 1],'xticklabelmode','auto','yticklabelmode','auto');
set(haContraTrace,'xcolor','none','ycolor','none');
legend(flipud(findobj(haContraTrace(1),'Tag','trace')),{px.name},'location','northwest');

figure('color','w','name','latencies');
ha = tight_subplot(2,1,0.1,[0.1 0.05],[0.15 0.05]);
mice={'Chain','Radnitz','Cori','Hench','Reichstein'};
hold(ha(1),'on'); hold(ha(2),'on');
Y_jitter = linspace(-0.2,0.2,length(mice));
allMeans=[];
allMeansAx=gobjects;

for m = 1:length(mice)
    mIdx = strcmp(sessionList.mouseName,mice(m));
    lat = ContraOnsetTime(mIdx,:);
    latIpsi = IpsiOnsetTime(mIdx,:);
    
    for p =1:6
        m1=mean(lat(:,p));
        stdr1 = std(lat(:,p))/sqrt(length(lat(:,p)));
        
        lx=line(ha(1), m1 + [-1 1]*stdr1, [1 1]*Y_jitter(m) + (p-1)*0.001, 'linewidth',3 );
        lx.Color=areaCols(p,:);
        
        allMeansAx(m,p)=line(ha(1),[1 1]*m1, Y_jitter(m)+(p-1)*0.001 + [-1 1]*0.05,'Color',areaCols(p,:), 'linewidth',3);
        allMeans(m,p) = m1;
        
        
        
        m1=mean(latIpsi(:,p));
        stdr1 = std(latIpsi(:,p))/sqrt(length(latIpsi(:,p)));
        
        lx=line(ha(2), m1 + [-1 1]*stdr1, [1 1]*Y_jitter(m) + (p-1)*0.001, 'linewidth',3 );
        lx.Color=areaCols(p,:);
        
    end
    
end

for p =1:6
    %     plot(ha,allMeans(:,p), Y_jitter + (p-1)*0.001, '-','Color',areaCols(p,:));
    plot(ha(1), mean(allMeans(:,p),1), Y_jitter(1)-0.1, 'o', 'color', areaCols(p,:));
end
xlim(ha,[-0.1 0.4]);
xlabel(ha(1),'Stim onset');
set(ha,'ytick',Y_jitter,'YTickLabel',mice,'XTickLabelMode','auto');
ylim(ha,[-0.35 0.25]);
legend(allMeansAx(1,:),{px.name});
