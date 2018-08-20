SESSIONID = 17; %Session to view examples


%%
sessionList = readtable('./sessionList.csv','FileType','text','Delimiter',',');
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
        
        %Mark stim/choice combination for every trial
        stim_resp_set = [b.trial.contrastLeft'>0 b.trial.contrastRight'>0 b.trial.choice'];
        [stim_resp,~,stim_resp_id] = unique(stim_resp_set,'rows');
        
        %Get Stimulus-aligned activity for every trial
        stimulus_times = b.timings.t_stimOn;
        assert(all(diff(stimulus_times)>0),'Timestamps not monotonically increasing');
        [avgperiStimV, periStimT, ~, ~] = ...
            eventLockedAvgSVD(wf.U_dff, wf.dV, wf.wfTime_dt,...
            stimulus_times, stim_resp_id, [0 0.25]);
        
        %Get movement-aligned activity for every trial
        movement_times = stimulus_times + b.trial.reactionTime_better';
        movement_times(isnan(movement_times)) = stimulus_times(isnan(movement_times)) + nanmedian(b.trial.reactionTime_better); %emulate nogo "movement" time
        [avgperiMoveV, periMoveT, ~, ~] = ...
            eventLockedAvgSVD(wf.U_dff, wf.dV, wf.wfTime_dt,...
            movement_times, stim_resp_id, [-0.2 0.2]);
        
        %Reconstruct the full map at each timepoint for each alignment
        MAPstim = nan(size(wf.U_dff,1), size(wf.U_dff,2), length(periStimT),size(stim_resp,1));
        MAPmove = nan(size(wf.U_dff,1), size(wf.U_dff,2), length(periMoveT),size(stim_resp,1));
        f = waitbar(0,'Please wait...');
        for cond = 1:size(stim_resp,1)
            waitbar(cond/size(stim_resp,1),f,'Loading your data');
            MAPstim(:,:,:,cond) = svdFrameReconstruct(wf.U_dff, permute(avgperiStimV(cond,:,:),[2 3 1]) ).*mask;
            MAPmove(:,:,:,cond) = svdFrameReconstruct(wf.U_dff, permute(avgperiMoveV(cond,:,:),[2 3 1]) ).*mask;
        end
        close(f);
        
        %     %Extract activity at each ROI
        %     ROIstim = nan(numTrials,length(periStimT),numel(roi.px));
        %     ROImove = nan(numTrials,length(periMoveT),numel(roi.px));
        %     for p =1:numel(roi.px)
        %         thisU = squeeze(wf.U_dff(roi.px(p).row, roi.px(p).col, :));
        %         for n = 1:numTrials
        %             ROIstim(n,:,p) = thisU'*squeeze(periStimV(n,:,:));
        %             ROImove(n,:,p) = thisU'*squeeze(periMoveV(n,:,:));
        %         end
        %     end
        %
        
        %Save smaller files to .mat
        save(wfAlignedFile,'stim_resp','periStimT','periMoveT','mask');
        
        %Save maps to HDFS
        mapFile = strrep(wfAlignedFile,'.mat','.h5');        
        h5create(mapFile,'/MAPstim',size(MAPstim));
        h5create(mapFile,'/MAPmove',size(MAPmove));
        h5write(mapFile,'/MAPstim',MAPstim);
        h5write(mapFile,'/MAPmove',MAPmove);
    end
end

%% PLOT 1: SESSION-AVERAGED MAP FOR ONE MOUSE
mouse = 'Hench';




%% Get activity and behavioural data, align timestamps, compute traces
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
wfFile = [ './preproc/WF_SVD/' eRef '.mat'];
behavFile = [ './preproc/BEHAV/' eRef '.mat'];
roiFile = [ './preproc/WF_ROI/' sessionList.mouseName{SESSIONID} '.mat'];
    
    wf = load(wfFile);
    b = load(behavFile);
    roi = load(roiFile);
    
    bregma=roi.px(1);
    roi.px(1)=[];
   

%Define mask around the outermost point of the allen atlas
load('D:\ctxOutlines.mat');
f=figure;
ha = axes;
addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col],[1 1 1]*0,ha);
contours = get(ha,'children');
mask = zeros(size(meanImg_registered));
for q = 1:length(contours)
    mask = mask | poly2mask(contours(q).XData, contours(q).YData, size(meanImg_registered,1), size(meanImg_registered,2));
end

mask = imgaussfilt(double(mask),3);
XLIM = [find(sum(mask,1)>0,1,'first') find(sum(mask,1)>0,1,'last')];
YLIM = [find(sum(mask,2)>0,1,'first') find(sum(mask,2)>0,1,'last')];
close(f);
disp('Data loaded');
%% 1) plot map aligned to stim onset
isDetectionTrial = (b.trial.contrastRight==0 | b.trial.contrastLeft==0);
isGoodPerformance = b.trial.repNum==1 & b.trial.feedback==1;
isGoodWheelMove = b.trial.choice==3 | isnan(b.trial.reactionTime_better) | (b.trial.reactionTime_better>0.125 & b.trial.reactionTime_better<0.5);
inclTrials = isDetectionTrial & isGoodWheelMove & isGoodPerformance;

figure(101); set(gcf,'color','w');
[avgdV, wS, ~, ~] = ...
    eventLockedAvgSVD(wf.U_dff, wf.dV, wf.wfTime_dt,...
    b.timings.t_stimOn(inclTrials), ...
    sign(b.trial.contrastRight(inclTrials)-b.trial.contrastLeft(inclTrials)), ...
    [0 0.25]);
avgdV = avgdV - avgdV(:,:,1);%remove baseline
sliceIdx = 5:15:length(wS);
ha = tight_subplot(3,length(sliceIdx),0.01,[0.7 0.05],0.01);

for i = 1:length(sliceIdx)
    imL = svdFrameReconstruct(U_dff, permute(avgdV(1,:,sliceIdx(i)),[2 3 1]) ).*mask;
    imagesc( ha(i), imL );
    addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*1, ha(i));
    hold(ha(i),'on'); plot(ha(i),bregma.col,bregma.row,'k.')
    
    imR = svdFrameReconstruct(U_dff, permute(avgdV(3,:,sliceIdx(i)),[2 3 1]) ).*mask;
    imagesc( ha(i + length(sliceIdx) ), imR );
    addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*1,ha(i + length(sliceIdx) ));
    hold(ha(i + length(sliceIdx) ),'on'); plot(ha(i + length(sliceIdx) ),bregma.col,bregma.row,'k.')
    
    imagesc( ha(i + 2*length(sliceIdx) ), imL-imR );
    addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*1,ha(i + 2*length(sliceIdx) ));
    hold(ha(i + 2*length(sliceIdx) ),'on'); plot(ha(i + 2*length(sliceIdx) ),bregma.col,bregma.row,'k.')
    
    title(ha(i), wS(sliceIdx(i)));
%     
%     imL = svdFrameReconstruct(U_dff, permute(avgdV(1,:,sliceIdx(i)),[2 3 1]) ).*mask;
%     imR = svdFrameReconstruct(U_dff, permute(avgdV(3,:,sliceIdx(i)),[2 3 1]) ).*mask;
%    
%     imagesc( ha,xcoords+ size(meanImg,2)*i,ycoords,imL );
%     addAllenCtxOutlines([bregma_AP bregma_ML+size(meanImg,2)*i], [lambda_AP lambda_ML+ size(meanImg,2)*i], [1 1 1]*1, ha);
%     
%     imagesc( ha,xcoords+ size(meanImg,2)*i,ycoords + size(meanImg,1),imR );
%     addAllenCtxOutlines([bregma_AP+size(meanImg,1) bregma_ML+size(meanImg,2)*i], [lambda_AP+size(meanImg,1) lambda_ML+ size(meanImg,2)*i], [1 1 1]*1, ha);
% 
%     imagesc( ha,xcoords+ size(meanImg,2)*i,ycoords + 2*size(meanImg,1),imL-imR );
%     addAllenCtxOutlines([bregma_AP+2*size(meanImg,1) bregma_ML+size(meanImg,2)*i], [lambda_AP+2*size(meanImg,1) lambda_ML+ size(meanImg,2)*i], [1 1 1]*1, ha);
end
cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
colormap(flipud(cmap));
set(ha,'clim',[-1 1]*0.008);
set(ha,'xcolor','none','ycolor','none','dataaspectratio',[1 1 1]);
set(ha,'xlim',XLIM,'ylim',YLIM);
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
sessIDs = find(contains(sessionList.mouseName,'Hench'));
eRefs = sessionList.expRef(sessIDs);
%Load first session meanImg 

meanImg_1=load([ './preproc/WF_ROI/' eRefs{1} '.mat'],'meanImg'); meanImg_1 = meanImg_1.meanImg;
[optimizer,metric] = imregconfig('Multimodal');

mapStack = nan(size(meanImg_1,1), size(meanImg_1,2), 2,  9,length(sessIDs));
for sess = 1:length(sessIDs)
    thisSID = sessIDs(sess);
    load([ './preproc/WF_ROI/' eRefs{sess} '.mat'],'eRef','px','wfTime_dt','U_dff','dV','meanImg','mon','trial','timings','inclTrials');
    
    %register this session's meanImage to the standard 
    transform = imregtform(meanImg, meanImg_1, 'similarity' ,optimizer,metric);
    
%     %Test that the registration worked
    figure;
    ha1=subplot(1,2,1);
    ha2=subplot(1,2,2);
    imshowpair(meanImg_1, meanImg,'Parent',ha1);
    meanImg_registered = imwarp(meanImg,transform,'OutputView',imref2d(size(meanImg_1)));
    imshowpair(meanImg_1, meanImg_registered,'Parent',ha2); drawnow;
    
    %Register the U components
    U_dff_registered = nan(size(meanImg_1,1), size(meanImg_1,2), size(U_dff,3));
    for c = 1:size(U_dff,3)
        U_dff_registered(:,:,c) = imwarp(U_dff(:,:,c),transform,'OutputView',imref2d(size(meanImg_1)));
    end
    
    [avgdV, wS, ~, ~] = ...
        eventLockedAvgSVD([], dV, wfTime_dt,...
        timings.t_stimOn(inclTrials), ...
        sign(trial.contrastRight(inclTrials)-trial.contrastLeft(inclTrials)), ...
        [0 0.25]);
    avgdV = avgdV - avgdV(:,:,1);%remove baseline
    
    sliceIdx = 5:15:length(wS);
    for i = 1:length(sliceIdx)
        mapStack(:,:,1,i,sess) = svdFrameReconstruct(U_dff_registered, permute(avgdV(1,:,sliceIdx(i)),[2 3 1]) );
        mapStack(:,:,2,i,sess) = svdFrameReconstruct(U_dff_registered, permute(avgdV(3,:,sliceIdx(i)),[2 3 1]) );
    end
end



fL=figure('name','CL');
haL = tight_subplot( length(sessIDs)+1, length(sliceIdx), 0, [0 0.05], 0);
fR=figure('name','CR');
haR = tight_subplot( length(sessIDs)+1, length(sliceIdx), 0, [0 0.05], 0);
for sess = 1:length(sessIDs)

    sliceIdx = 5:15:length(wS);
    for i = 1:length(sliceIdx)
        imagesc(haL( length(sliceIdx)*(sess-1) + i), mapStack(:,:,1, i, sess) );
        imagesc(haR( length(sliceIdx)*(sess-1) + i), mapStack(:,:,2, i, sess) );
        
        if sess == 1
            title( haL( length(sliceIdx)*(sess-1) + i), 1000*wS(sliceIdx(i)));
            title( haR( length(sliceIdx)*(sess-1) + i), 1000*wS(sliceIdx(i)));
        end
    end
end

avgMap = mean(mapStack,5);
for i = 1:length(sliceIdx)
    imagesc(haL( length(sliceIdx)*(length(sessIDs)) + i), avgMap(:,:,1, i) );
    imagesc(haR( length(sliceIdx)*(length(sessIDs)) + i), avgMap(:,:,2, i) );
end


set([haL; haR],'dataaspectratio',[ 1 1 1],'clim',[-1 1]*0.008,'xcolor','none','ycolor','none');

cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];

figure(fL);
colormap(flipud(cmap));
figure(fR);
colormap(flipud(cmap));

