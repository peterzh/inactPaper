SESSIONID = 17; %Session to view examples


%%
sessionList = readtable('./sessionList.csv','FileType','text','Delimiter',',');
areaCols = [      0    0.4470    0.7410; %blue
    0.3010    0.7450    0.9330; %lightblue
    0.4660    0.6740    0.1880; %green
    0.8500    0.3250    0.0980; %orange
    0.4940    0.1840    0.5560; %purple
    0.8 0.8 0.8]; %grey

%% PREPROC: Behavioural data
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    behavFile = [ './preproc/BEHAV/' eRef '.mat'];
    if ~exist(behavFile,'file')
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
        trial.reactionTime_better = rt;
        
        %%%%%%%%%%
        %%% Save %
        %%%%%%%%%%
        save(behavFile,'trial','timings');
    end
end
%% PREPROC: Widefield data
num_svd_components = 500;
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    wfFile = [ './preproc/WF_SVD/' eRef '.mat'];
    if ~exist(wfFile,'file')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Haemodynamic correction & SVD %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        corrPath = fullfile(sessionList.widefieldDir{sess}, 'svdTemporalComponents_corr.npy');
        if ~exist(corrPath,'file') %Do haemodynamic correction if it doesnt exist
            quickHemoCorrect(sessionList.widefieldDir{sess},num_svd_components);
        end
        [U,V,wfTime,meanImg]=quickLoadUVt(sessionList.widefieldDir{sess},num_svd_components);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Compute dF/F %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [U_dff, V_dff] = dffFromSVD(U, V, meanImg);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Align to reference map%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        refImg = load(['./preproc/reference_image/' sessionList.mouseName{sess} '.mat'], 'meanImg');
        [optimizer,metric] = imregconfig('Multimodal');
        tform = imregtform(meanImg, refImg.meanImg, 'similarity' ,optimizer,metric);
        figure('name',eRef);
        ha1=subplot(1,2,1);
        ha2=subplot(1,2,2);
        imshowpair(refImg.meanImg, meanImg,'Parent',ha1);
        meanImg_registered = imwarp(meanImg,tform,'OutputView',imref2d(size(refImg.meanImg)));
        imshowpair(refImg.meanImg, meanImg_registered,'Parent',ha2); drawnow;
        
        %Overwrite the U components with the registered version
        try
            U_dff = imwarp(U_dff,tform,'OutputView',imref2d(size(refImg.meanImg)));
        catch me
            disp(me);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Compute derivative of dF/F %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        %%%%%%%%%%%%%%%%%
        %%% Save %%%%%%%%
        %%%%%%%%%%%%%%%%%
        save(wfFile,'U_dff','dV','wfTime_dt','meanImg_registered');
    end
end
%% PREPROC: Identify bregma and ROIs for each mouse
stack = cell(height(sessionList),1);
meanImgs = cell(height(sessionList),1);
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    wfFile = [ './preproc/WF_SVD/' eRef '.mat'];
    behavFile = [ './preproc/BEHAV/' eRef '.mat'];
    
    wf = load(wfFile);
    b = load(behavFile);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% INCLUSION CRITERIA FOR TRIALS %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    isDetectionTrial = (b.trial.contrastRight==0 | b.trial.contrastLeft==0);
    isGoodPerformance = b.trial.repNum==1 & b.trial.feedback==1;
    isGoodWheelMove = b.trial.choice==3 | isnan(b.trial.reactionTime_better) | (b.trial.reactionTime_better>0.125 & b.trial.reactionTime_better<0.5);
    inclTrials = isDetectionTrial & isGoodWheelMove & isGoodPerformance;
    
    disp(mean(inclTrials))
    %Compute average activity map
    [avgdV, ~, ~, ~] = ...
        eventLockedAvgSVD(wf.U_dff, wf.dV, wf.wfTime_dt,...
        b.timings.t_stimOn(inclTrials), ...
        sign(b.trial.contrastRight(inclTrials)-b.trial.contrastLeft(inclTrials)), ...
        [0 0.25]);
    avgdV = avgdV - avgdV(:,:,1);%remove baseline
    
    stack{sess} = svdFrameReconstruct(wf.U_dff, permute(avgdV(1,:,:),[2 3 1]) );
    meanImgs{sess} = wf.meanImg_registered;
end

%Average sessions for each mouse
names = unique(sessionList.mouseName);
sliceIdx = 5:15:126;
for mouse = 1:length(names)
    sessIDs = find(contains(sessionList.mouseName, names{mouse}));
    
    imL = mean( cat(4,stack{sessIDs}), 4);
    
    figure;
    for i = 1:9
        subplot(3,3,i);
        imagesc(imL(:,:,sliceIdx(i)));
    end
    set(get(gcf,'children'), 'clim', [-1 1]*quantile(abs(imL(:)),0.99))
        cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
    colormap(flipud(cmap));
    
    px = struct('name',[],'row',[],'col',[]);
    ROIs = {'Bregma','VISp','VISlm','MOs','MOp','SSp'};
    for r = 1:length(ROIs)
        px(r).name = ROIs{r};
        disp(px(r).name);
        [col,row] = ginput(1); px(r).col=round(col); px(r).row=round(row);
        hold on; plot(col,row,'k+');
    end
    
    %Now overlay points on the meanImg as well as the allen atlas to check
    %the bregma position made sense
    figure;
    imagesc(meanImgs{sessIDs(1)}); hold on;
    for r = 1:length(ROIs)
        plot(px(r).col,px(r).row,'k+');
    end

    pxFile = [ './preproc/WF_ROI/' names{mouse} '.mat'];
    
    addAllenCtxOutlines([px(1).row px(1).col], [px(1).row-1 px(1).col], [1 1 1], gca);
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
                save(pxFile,'px');
        end
        
        for i = 1:length(outlines)
            outlines(i).YData = outlines(i).YData + shiftY;
            outlines(i).XData = outlines(i).XData + shiftX;
        end
        
        px(1).row = px(1).row + shiftY;
        px(1).col = px(1).col + shiftX;
    end
end
%% PREPROC: Extract epoch-aligned activity maps
%Get stimulus-aligned and movement-aligned activity
%Averages over repeats of that condition
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    wfFile = [ './preproc/WF_SVD/' eRef '.mat'];
    behavFile = [ './preproc/BEHAV/' eRef '.mat'];
    pxFile = [ './preproc/WF_ROI/' sessionList.mouseName{sess} '.mat'];
 
    wfAlignedFile = [ './preproc/WF_aligned/' eRef '.mat'];
    if ~exist(wfAlignedFile,'file')
        wf = load(wfFile);
        b = load(behavFile);
        roi = load(pxFile);
        
        %Define mask around the outermost point of the allen atlas
        bregma=roi.px(1);
        f=figure;
        ha = axes;
        addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col],[1 1 1]*0,ha);
        contours = get(ha,'children');
        mask = zeros(size(wf.meanImg_registered));
        for q = 1:length(contours)
            mask = mask | poly2mask(contours(q).XData, contours(q).YData, size(wf.meanImg_registered,1), size(wf.meanImg_registered,2));
        end
        mask = imgaussfilt(double(mask),3);
        f.delete;
        
        %Identify good trials to include in averaging
        isGoodPerformance = b.trial.repNum==1;
        isGoodWheelMove = b.trial.choice==3 | isnan(b.trial.reactionTime_better) | (b.trial.reactionTime_better>0.125 & b.trial.reactionTime_better<0.5);
        inclTrials = isGoodWheelMove & isGoodPerformance;

        %Mark stim/choice combination for every trial
        stim_resp_set = [b.trial.contrastLeft(inclTrials)'>0 b.trial.contrastRight(inclTrials)'>0 b.trial.choice(inclTrials)'];
        [stim_resp,~,stim_resp_id] = unique(stim_resp_set,'rows');
        tab = tabulate(stim_resp_id); 
        stim_resp_counts = tab(:,2);
        
        %Get Stimulus-aligned activity for every trial
        stimulus_times = b.timings.t_stimOn(inclTrials);
        assert(all(diff(stimulus_times)>0),'Timestamps not monotonically increasing');
        [avgperiStimV, periStimT, ~, ~] = ...
            eventLockedAvgSVD(wf.U_dff, wf.dV, wf.wfTime_dt,...
            stimulus_times, stim_resp_id, [0 0.25]);
        
        %Get movement-aligned activity for every trial
        movement_times = stimulus_times + b.trial.reactionTime_better(inclTrials)';
        movement_times(isnan(movement_times)) = stimulus_times(isnan(movement_times)) + nanmedian(b.trial.reactionTime_better(inclTrials)); %emulate nogo "movement" time
        [avgperiMoveV, periMoveT, ~, ~] = ...
            eventLockedAvgSVD(wf.U_dff, wf.dV, wf.wfTime_dt,...
            movement_times, stim_resp_id, [-0.2 0.2]);
        
        %Reconstruct the full map at each timepoint for each alignment
        MAPstim = nan(size(wf.U_dff,1), size(wf.U_dff,2), length(periStimT), 2, 2, 3);
        MAPmove = nan(size(wf.U_dff,1), size(wf.U_dff,2), length(periMoveT), 2, 2, 3);
        stim_resp_counts = nan(2,2,3);
        f = waitbar(0,'Please wait...');
        for cond = 1:size(stim_resp,1)
            waitbar(cond/size(stim_resp,1),f,'Loading your data');
            
            CL = stim_resp(cond,1)+1;
            CR = stim_resp(cond,2)+1;
            R = stim_resp(cond,3);
            
            MAPstim(:,:,:,CL,CR,R) = svdFrameReconstruct(wf.U_dff, permute(avgperiStimV(cond,:,:),[2 3 1]) ).*mask;
            MAPmove(:,:,:,CL,CR,R) = svdFrameReconstruct(wf.U_dff, permute(avgperiMoveV(cond,:,:),[2 3 1]) ).*mask;
            stim_resp_counts(CL,CR,R) = sum(stim_resp_id==cond);
        end
        close(f);
        
        %Save maps to HDFS
        mapFile = strrep(wfAlignedFile,'.mat','.h5');        
        h5create(mapFile,'/MAPstim',size(MAPstim));
        h5create(mapFile,'/MAPmove',size(MAPmove));
        h5write(mapFile,'/MAPstim',MAPstim);
        h5write(mapFile,'/MAPmove',MAPmove);
        
                %Save smaller files to .mat
        save(wfAlignedFile,'stim_resp_counts','periStimT','periMoveT','mask');
        
    end
end


%% OLD: Get activity and behavioural data, align timestamps, compute traces
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    preprocFile = [ './preproc/WF_ROI/' eRef '.mat'];
    
    if ~exist(preprocFile,'file') & ~isnan(sessionList.row_VISp(sess))
       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Load ROI definitions %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%% INCLUSION CRITERIA FOR TRIALS %%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         isDetectionTrial = (trial.contrastRight==0 | trial.contrastLeft==0);
%         isGoodPerformance = trial.repNum==1 & trial.feedback==1;
%         isGoodWheelMove = trial.choice==3 | isnan(rt) | (rt>0.125 & rt<0.5);
%         inclTrials = isDetectionTrial & isGoodWheelMove & isGoodPerformance;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Save %%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         save(preprocFile,'eRef','px','wfTime_dt','U_dff','dV','meanImg','mon','trial','timings','inclTrials');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% COMPUTE STIM-ALIGNED TRACES & ONSET LATENCIES %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
    end
end
%% 0) Load and prepare data for the example session
eRef = sessionList.expRef{SESSIONID};
fprintf('Session %d %s\n',SESSIONID,eRef);
wfFile = [ './preproc/WF_SVD/' eRef '.mat'];
behavFile = [ './preproc/BEHAV/' eRef '.mat'];
roiFile = [ './preproc/WF_ROI/' sessionList.mouseName{SESSIONID} '.mat'];
wfAlignedFile = [ './preproc/WF_aligned/' eRef '.mat'];
wfAlignedMapFile = [ './preproc/WF_aligned/' eRef '.h5'];

wf = load(wfFile);
b = load(behavFile);
roi = load(roiFile);
align = load(wfAlignedFile);
MAPstim = h5read(wfAlignedMapFile,'/MAPstim');
MAPmove = h5read(wfAlignedMapFile,'/MAPmove');

bregma=roi.px(1);
roi.px(1)=[];

disp('Data loaded');
%% 1) plot map aligned to stim onset
figure(101); set(gcf,'color','w');

%Rectify and baseline stim-aligned maps
MAPstim(MAPstim<0)=0;
MAPstim = MAPstim - MAPstim(:,:,1,:);


sliceTimes = [0 0.07 0.12 0.2 0.25];
sliceTimes = [0.008 0.038 0.068 0.098 0.128 0.158 0.188 0.218 0.248];
ha = tight_subplot(3,length(sliceTimes),0,[0.7 0.02],0);
for i = 1:length(sliceTimes)
    %get nearest time index for this sliceTime
    idx = find( align.periStimT >= sliceTimes(i),1,'first');
    
    imL = MAPstim(:,:,idx,2,1,1);
    imagesc( ha(i), imL );
    addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*1, ha(i));
    hold(ha(i),'on'); plot(ha(i),bregma.col,bregma.row,'k.')
    
    imR = MAPstim(:,:,idx,1,2,2);
    imagesc( ha(i + length(sliceTimes) ), imR );
    addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*1,ha(i + length(sliceTimes) ));
    hold(ha(i + length(sliceTimes) ),'on'); plot(ha(i + length(sliceTimes) ),bregma.col,bregma.row,'k.')
    
    imagesc( ha(i + 2*length(sliceTimes) ), imL-imR );
    addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*1,ha(i + 2*length(sliceTimes) ));
    hold(ha(i + 2*length(sliceTimes) ),'on'); plot(ha(i + 2*length(sliceTimes) ),bregma.col,bregma.row,'k.')
    
    title(ha(i), sliceTimes(i) );
end
cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
colormap(flipud(cmap));
set(ha,'clim',[-1 1]*0.008);
set(ha,'xcolor','none','ycolor','none','dataaspectratio',[1 1 1]);
% set(ha,'xlim',XLIM,'ylim',YLIM);
%% 2) Plot traces and ROIs

[~, winSamps, periEventdV, label] = ...
    eventLockedAvgSVD(wf.U_dff, wf.dV, wf.wfTime_dt,...
    b.timings.t_stimOn(inclTrials), ...
    sign(b.trial.contrastRight(inclTrials)-b.trial.contrastLeft(inclTrials)), ...
    [-0.3 0.8]);
numTrials = length(label);
avgStimAlign = nan(numTrials, length(winSamps), length(roi.px));
contraLatency = nan(numel(roi.px),1);
for p =1:numel(roi.px)
    thisU = squeeze(U_dff(roi.px(p).row, roi.px(p).col, :));
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

figure(101);
ha = tight_subplot(5,1,0.01,[0.05 0.35],[0.05 0.66]);
thisTraces_contra = avgStimAlign(label==-1,:,:);
thisTraces_ipsi = avgStimAlign(label==1,:,:);
for p =1:numel(roi.px)
    hold(ha(p),'on');
    plot(ha(p), winSamps, thisTraces_contra(:,:,p), 'k-', 'color', [areaCols(p,:) 0.1] );
    plot(ha(p), winSamps, mean(thisTraces_contra(:,:,p),1), '-', 'linewidth',3, 'color', areaCols(p,:));
    
    %Plot ipsi trace
    %     plot(ha(p), winSamps, -thisTraces_ipsi(:,:,p), 'k-', 'color', [areaCols(p,:) 0.1] );
    plot(ha(p), winSamps, mean(thisTraces_ipsi(:,:,p),1), '-', 'linewidth',1, 'color', areaCols(p,:));
    set(ha(p),'xcolor','none','ycolor','none');
    line(ha(p), [0 0], [-1 1]*10, 'color', [0 0 0 0.1],'LineStyle','-');
    
    %Add text in top left corner
    tx=text(ha(p),-0.3,0.01,roi.px(p).name,'VerticalAlignment','top','color',areaCols(p,:),'FontWeight','bold');
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
imagesc(ha(1),wf.meanImg_registered); colormap(ha(1),'gray'); hold(ha(1),'on');
addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*0.5, ha(1))
plot(ha(1),bregma.col,bregma.row,'k.')
for p =1:numel(roi.px)
    h = plot(ha(1),roi.px(p).col, roi.px(p).row, '.', 'Color', areaCols(p,:), 'markersize',25);
end
set(ha(1),'xcolor','none','ycolor','none', 'dataaspectratio', [1 1 1]);

hold(ha(2),'on'); hold(ha(3),'on');
for p =1:numel(roi.px)
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
%% 3) Plot latencies over all sessions TODO
mice={'Chain','Radnitz','Cori','Reichstein','Hench'};

figure(101);
ContraOnsetTime = nan(height(sessionList),5);
rt_vals = cell(height(sessionList),1);
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    try
        load([ './preproc/WF_aligned/' eRef '.mat'],'contraLatency');
        load([ './preproc/WF_ROI/' eRef '.mat'],'mon','trial','timings','inclTrials');
        ContraOnsetTime(sess,:) = contraLatency';
        
        rt = [mon-timings.t_stimOn];
        rt_vals{sess} = rt(inclTrials);
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
    rts_subj = cat(1,rt_vals{mIdx});
    
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
    
    %now add reaction time distribution
    rt_stats = quantile(rts_subj, [0.25 0.5 0.75]);
    
    line(ha(1), rt_stats([1,3]), [1 1]*Y_jitter(m), 'Color', 'k');
    plot(ha(1), rt_stats(2), Y_jitter(m), 'k.', 'markersize',15);
    
end

xlim(ha,[0 0.3]);
xlabel(ha(1),'Stim onset');
set(ha,'ytick',Y_jitter,'YTickLabel',mice,'XTickLabelMode','auto');
ylim(ha,[-0.3 0.3]);
%% MISC: Plot ROI locations
figureROILocations = figure('color','w','name','ROI locations');
haROI = tight_subplot(5,8,0.04,[0.01 0.05],[0.05 0.01]);

for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    preprocFile = [ './preproc/WF_ROI/' eRef '.mat'];
    
    if exist(preprocFile,'file')
        load(preprocFile,'px');
        load([ './preproc/bregma/' eRef '_bregma.mat'],'bregma_AP','bregma_ML','lambda_AP','lambda_ML');
        
        %Plot ROI
        load([ './preproc/WF_SVD/' eRef '.mat'],'meanImg');
        imagesc(haROI(sess),meanImg); colormap(haROI(sess),'gray');
        addAllenCtxOutlines([bregma_AP bregma_ML], [lambda_AP lambda_ML], [1 1 1]*0.5, haROI(sess))
        
        hold(haROI(sess),'on');
        for p =1:numel(px)
            h = plot(haROI(sess),px(p).xy(2), px(p).xy(1), 'o', 'color', areaCols(p,:));
            set(h, 'MarkerFaceColor', get(h, 'Color'));
        end
        title(haROI(sess),eRef,'interpreter','none');
        set(haROI(sess),'xcolor','none','ycolor','none');
        
    end
end
%% OLD: Plot sequence at contralateral contrast
figureContraTrace = figure('color','w','name','Sequence to CONTRALATERAL stimulus');
haContraTrace = tight_subplot(5,8,0.03,[0.01 0.05],[0.05 0.01]);
figureIpsiTrace = figure('color','w','name','Sequence to IPSILATERAL stimulus');
haIpsiTrace = tight_subplot(5,8,0.03,[0.01 0.05],[0.05 0.01]);
ContraOnsetTime = nan(height(sessionList),6);
IpsiOnsetTime = nan(height(sessionList),6);
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    preprocFile = [ './preproc/WF_ROI/' eRef '.mat'];
    
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
%% MISC: Plot RT histograms

figureROILocations = figure('color','w','name','RT');
ha = tight_subplot(5,8,0.04,[0.01 0.05],[0.05 0.01]);

for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    preprocFile = [ './preproc/WF_ROI/' eRef '.mat'];
    if exist(preprocFile,'file')
        load([ './preproc/WF_ROI/' eRef '.mat'],'mon','trial','timings','inclTrials');
        rt = [mon-timings.t_stimOn]';
        
        histogram(ha(sess),rt(inclTrials), 100)
    end
end

%% GROUP LEVEL: average maps
mouse = 'Hench';
refImg = load(['./preproc/reference_image/' mouse '.mat'], 'meanImg');
sessIDs = find(contains(sessionList.mouseName,mouse));

stimSliceTimes = [0 0.07 0.12 0.2 0.25];
stimSliceTimes = linspace(0,0.25,10);
[~,stimSliceIdx]=find(abs(align.periStimT - stimSliceTimes')==min(abs(align.periStimT - stimSliceTimes'),[],2));

moveSliceTimes = [-0.2 -0.1 0 0.1 0.2];
moveSliceTimes = linspace(-0.1,0.2,10);
[~,moveSliceIdx]=find(abs(align.periMoveT - moveSliceTimes')==min(abs(align.periMoveT - moveSliceTimes'),[],2));

MAPstimStack = nan(size(refImg.meanImg,1), size(refImg.meanImg,2), length(stimSliceIdx),2,2,3,length(sessIDs));
MAPmoveStack = nan(size(refImg.meanImg,1), size(refImg.meanImg,2), length(moveSliceIdx),2,2,3,length(sessIDs));
for sess = 1:length(sessIDs)
    eRef = sessionList.expRef{sessIDs(sess)};
    fprintf('Session %d %s\n',sessIDs(sess),eRef);
    wfFile = [ './preproc/WF_SVD/' eRef '.mat'];
    behavFile = [ './preproc/BEHAV/' eRef '.mat'];
    roiFile = [ './preproc/WF_ROI/' sessionList.mouseName{SESSIONID} '.mat'];
    wfAlignedFile = [ './preproc/WF_aligned/' eRef '.mat'];
    wfAlignedMapFile = [ './preproc/WF_aligned/' eRef '.h5'];
    
    wf = load(wfFile);
    b = load(behavFile);
    roi = load(roiFile);
    align = load(wfAlignedFile);
    MAPstim = h5read(wfAlignedMapFile,'/MAPstim');
    MAPmove = h5read(wfAlignedMapFile,'/MAPmove');
    
    bregma=roi.px(1);
    roi.px(1)=[];
    
    %Rectify and baseline stim-aligned maps
    MAPstim(MAPstim<0)=0;
    MAPstim = MAPstim - MAPstim(:,:,1,:,:,:);
    
%     %baseline mov-aligned maps
%     MAPstim = MAPmove - MAPmove(:,:,1,:,:,:);
    
    MAPstimStack(:,:,:,:,:,:,sess) = MAPstim(:,:,stimSliceIdx,:,:,:);
    MAPmoveStack(:,:,:,:,:,:,sess) = MAPmove(:,:,moveSliceIdx,:,:,:);
end
avgMAPstimStack = nanmean(MAPstimStack,7);
avgMAPmoveStack = nanmean(MAPmoveStack,7);

%Every session: left choice on CL
f=figure('name','Every session: Left choice on CL');
haL = tight_subplot( length(sessIDs)+1, length(stimSliceTimes), 0, [0 0.05], 0);
for sess = 1:length(sessIDs)
    for i = 1:length(stimSliceTimes)
        imagesc(haL( length(stimSliceTimes)*(sess-1) + i), MAPstimStack(:,:,i,2,1,1,sess) );        
    end
end
for i = 1:length(stimSliceTimes)
    imagesc(haL( length(stimSliceTimes)*(length(sessIDs)) + i), avgMAPstimStack(:,:,i,2,1,1) );
%     imagesc(haR( length(sliceIdx)*(length(sessIDs)) + i), avgMap(:,:,2, i) );
end
% set(haL,'dataaspectratio',[ 1 1 1],'clim',[-1 1]*0.008,'xcolor','none','ycolor','none');
set(haL,'dataaspectratio',[ 1 1 1],'clim',[-1 1]*quantile(abs(avgMAPstimStack(:)),0.99),'xcolor','none','ycolor','none');
cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
colormap(flipud(cmap));

%Average only: all set
f=figure('name','Average across sessions for this mouse STIM ALIGNED');
haTop = tight_subplot( 2, length(stimSliceTimes)/2, 0.01, [0 0.05], 0);
for i = 1:length(stimSliceTimes)
    bounds = get(haTop(i),'position');
    set(haTop(i),'xcolor','none','ycolor','none');
    title(haTop(i),stimSliceTimes(i));
        
     %replace with many subplots within the bounds of the original
     ha = tight_subplot(3,3,0, [bounds(2) 1-bounds(2)-bounds(4)], [bounds(1) 1-bounds(1)-bounds(3)]);
     
     %Plot a 3x3 grid of respxstim condition maps
     
     imagesc(ha(1), avgMAPstimStack(:,:,i,2,1,1)); %Left choice & Left stim
     imagesc(ha(2), avgMAPstimStack(:,:,i,1,1,1)); %Left choice & Zero stim
     imagesc(ha(3), avgMAPstimStack(:,:,i,1,2,1)); %Left choice & Right stim
     
     imagesc(ha(4), avgMAPstimStack(:,:,i,2,1,3)); %NoGo choice & Left stim
     imagesc(ha(5), avgMAPstimStack(:,:,i,1,1,3)); %NoGo choice & Zero stim
     imagesc(ha(6), avgMAPstimStack(:,:,i,1,2,3)); %NoGo choice & Right stim
     
     imagesc(ha(7), avgMAPstimStack(:,:,i,2,1,2)); %Right choice & Left stim
     imagesc(ha(8), avgMAPstimStack(:,:,i,1,1,2)); %Right choice & Zero stim
     imagesc(ha(9), avgMAPstimStack(:,:,i,1,2,2)); %Right choice & Right stim
     set(ha,'dataaspectratio',[ 1 1 1],'clim',[-1 1]*quantile(abs(avgMAPstimStack(:)),0.999),'xcolor','none','ycolor','none');
end
colormap(flipud(cmap));

%Average only: all set
f=figure('name','Average across sessions for this mouse MOVE ALIGNED');
haTop = tight_subplot( 2, length(moveSliceTimes)/2, 0.01, [0 0.05], 0);
for i = 1:length(moveSliceTimes)
    bounds = get(haTop(i),'position');
    set(haTop(i),'xcolor','none','ycolor','none');
    title(haTop(i),moveSliceTimes(i));
        
     %replace with many subplots within the bounds of the original
     ha = tight_subplot(3,3,0, [bounds(2) 1-bounds(2)-bounds(4)], [bounds(1) 1-bounds(1)-bounds(3)]);
     
     %Plot a 3x3 grid of respxstim condition maps
     
     imagesc(ha(1), avgMAPmoveStack(:,:,i,2,1,1)); %Left choice & Left stim
     imagesc(ha(2), avgMAPmoveStack(:,:,i,1,1,1)); %Left choice & Zero stim
     imagesc(ha(3), avgMAPmoveStack(:,:,i,1,2,1)); %Left choice & Right stim
     
     imagesc(ha(4), avgMAPmoveStack(:,:,i,2,1,3)); %NoGo choice & Left stim
     imagesc(ha(5), avgMAPmoveStack(:,:,i,1,1,3)); %NoGo choice & Zero stim
     imagesc(ha(6), avgMAPmoveStack(:,:,i,1,2,3)); %NoGo choice & Right stim
     
     imagesc(ha(7), avgMAPmoveStack(:,:,i,2,1,2)); %Right choice & Left stim
     imagesc(ha(8), avgMAPmoveStack(:,:,i,1,1,2)); %Right choice & Zero stim
     imagesc(ha(9), avgMAPmoveStack(:,:,i,1,2,2)); %Right choice & Right stim
     set(ha,'dataaspectratio',[ 1 1 1],'clim',[-1 1]*quantile(abs(avgMAPmoveStack(:)),0.999),'xcolor','none','ycolor','none');
end
colormap(flipud(cmap));





