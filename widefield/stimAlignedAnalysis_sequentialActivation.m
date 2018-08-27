%%
sessionList = readtable('./sessionList.csv','FileType','text','Delimiter',',');
areaCols = [      0    0.4470    0.7410; %blue
    0.3010    0.7450    0.9330; %lightblue
    0.4660    0.6740    0.1880; %green
    0.8500    0.3250    0.0980; %orange
    0.4940    0.1840    0.5560; %purple
    0.8 0.8 0.8]; %grey

%% PREPROC: Behavioural data
clearvars -except sessionList
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
%% PREPROC: Exclude behavioural trials which have mismatched RTs between left and right
warning('not deterministic! uses random sampling, therefore re-running will exclude different trials!');
allBehav = struct;
[names,~,subjID]=unique(sessionList.mouseName);
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    %Load behavioural data
    behavFile = [ './preproc/BEHAV/' eRef '.mat'];
    b = load(behavFile);
    
    %Get good trials based on good performance/wheel movement
    isGoodPerformance = b.trial.repNum==1;
    isGoodWheelMove = b.trial.choice==3 | isnan(b.trial.reactionTime_better) | (b.trial.reactionTime_better>0.125 & b.trial.reactionTime_better<0.5);
    inclTrials = isGoodWheelMove & isGoodPerformance;
    
    row = struct;
    row.contrastLeft = b.trial.contrastLeft(inclTrials)';
    row.contrastRight = b.trial.contrastRight(inclTrials)';
    row.choice = b.trial.choice(inclTrials)';
    row.RT = b.trial.reactionTime_better(inclTrials)';
    row.t_stimOn = b.timings.t_stimOn(inclTrials);
    row.sessionID = ones(size(row.choice))*sess;
    row.subjectID = ones(size(row.choice))*subjID(sess);
    allBehav = addstruct(allBehav, row);
end
allBehav.keep = ones(size(allBehav.choice));

%Go through each mouse and exclude trials in order to match the RT between
%them
allBehavMatchedRT = struct;
figure;
for m = 1:length(names)
    subjData = getrow(allBehav, allBehav.subjectID==m);
    [N,EDGES,BIN]=histcounts(subjData.RT);
    
    subplot(length(names),2,2*(m-1) + 1);
    histogram(subjData.RT(subjData.choice==1),EDGES); hold on;
    histogram(subjData.RT(subjData.choice==2),EDGES); 
    ylabel(names{m}); title('uncorrected');
    
    Lchoice = histcounts(subjData.RT(subjData.choice==1), EDGES);
    Rchoice = histcounts(subjData.RT(subjData.choice==2), EDGES);
    for bin = 1:length(Lchoice)
        %For each bin, identify the trials to remove
        numberToRemove = abs(Lchoice(bin) - Rchoice(bin));
        removeIdx=[];
        if Lchoice(bin) > Rchoice(bin) %trim left choices
            removeIdx = randsample( find(BIN==bin & subjData.choice==1), numberToRemove);
        elseif Lchoice(bin) < Rchoice(bin) %trim right choices
            removeIdx = randsample( find(BIN==bin & subjData.choice==2), numberToRemove);
        end
        subjData.keep(removeIdx) = 0;
    end
    subjData = getrow(subjData, subjData.keep==1);
    
    subplot(length(names),2,2*(m-1) + 2);
    histogram(subjData.RT(subjData.choice==1),EDGES); hold on;
    histogram(subjData.RT(subjData.choice==2),EDGES);  title('corrected');
    
    allBehavMatchedRT = addstruct(allBehavMatchedRT, subjData);
end

%Save this new behavioural data
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    b = getrow(allBehavMatchedRT, allBehavMatchedRT.sessionID==sess);
    behavFile = [ './preproc/BEHAV_RT_Corrected/' eRef '.mat'];
    save(behavFile,'-struct','b');
end
%% PREPROC: Widefield SVD and filtering
clearvars -except sessionList
num_svd_components = 500;
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    wfSVDfile = [ './preproc/WF_SVD/' eRef '.mat'];
    if ~exist(wfSVDfile,'file')
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
        falign=figure('name',eRef);
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
        
        save(wfSVDfile, 'U_dff', 'dV', 'wfTime_dt','meanImg_registered');
    end
end
%% PREPROC: Widefield extract aligned activity
clearvars -except sessionList
window_stim_aligned = [-0.2 0.8];
window_move_aligned = [-0.3 0.3];

MAPstimT = linspace(0.01,0.25,9); %intended times
MAPmoveT = linspace(-0.2,0.2,9); %intended times
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    wfAlignedFile = [ './preproc/WF_aligned/' eRef '.h5'];
    if ~exist(wfAlignedFile,'file')
        wfSVDfile = [ './preproc/WF_SVD/' eRef '.mat'];
        svd = load(wfSVDfile);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Define stim+response combination groups %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        behavFile = [ './preproc/BEHAV_RT_Corrected/' eRef '.mat'];
        b = load(behavFile);
        stim_resp_set = [b.contrastLeft>0 b.contrastRight>0 b.choice];
        [stim_resp,~,stim_resp_id] = unique(stim_resp_set,'rows');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Epoch-aligned fluorescence %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        stimulus_times = b.t_stimOn;
        assert(all(diff(stimulus_times)>0),'Timestamps not monotonically increasing');
        [~, stimT, periStimV, ~] = ...
            eventLockedAvgSVD(svd.U_dff, svd.dV, svd.wfTime_dt,...
            stimulus_times, stim_resp_id, window_stim_aligned);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Get average Movement-aligned activity  %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        movement_times = stimulus_times + b.RT;
        movement_times(isnan(movement_times)) = stimulus_times(isnan(movement_times)) + nanmedian(b.RT); %emulate nogo "movement" time
        [~, moveT, periMoveV, ~] = ...
            eventLockedAvgSVD(svd.U_dff, svd.dV, svd.wfTime_dt,...
            movement_times, stim_resp_id, window_move_aligned);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Reconstruct the trial averaged map at each sliceTime for each alignment
        %%%% Only use "good" trials for this %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         isGoodPerformance = b.trial.repNum==1;
%         isGoodWheelMove = b.trial.choice==3 | isnan(b.trial.reactionTime_better) | (b.trial.reactionTime_better>0.125 & b.trial.reactionTime_better<0.5);
%         inclTrials = isGoodWheelMove & isGoodPerformance;
%         fprintf('Keeping %d%% of trials in averaging\n', round(100*mean(inclTrials)));
%         
        [~,stimsliceIdx]=max(stimT' >= MAPstimT,[],1);
        [~,movesliceIdx]=max(moveT' >= MAPmoveT,[],1);
        MAPstimT = stimT(stimsliceIdx); %actual slice times
        MAPmoveT = moveT(movesliceIdx); %actual slice times
        
        MAPstim = nan(size(svd.U_dff,1), size(svd.U_dff,2), length(MAPstimT), 2, 2, 3);
        MAPmove = nan(size(svd.U_dff,1), size(svd.U_dff,2), length(MAPmoveT), 2, 2, 3);
        for cond = 1:size(stim_resp,1)           
            CL = stim_resp(cond,1)+1;
            CR = stim_resp(cond,2)+1;
            R = stim_resp(cond,3);
            
            %compute average periStimV and periMoveV
            condIdx = stim_resp_id==cond;
            
            %Get average over trials
            avgperiStimV = mean(periStimV(condIdx,:,:),1);
            avgperiMoveV = mean(periMoveV(condIdx,:,:),1);
            
            %Baseline
            avgperiStimV = avgperiStimV - mean( avgperiStimV(:,:,stimT<0),3); 
            
            %Take time slices
            avgperiStimV = avgperiStimV(1,:,stimsliceIdx);
            avgperiMoveV = avgperiMoveV(1,:,movesliceIdx);
            
            
            MAPstim(:,:,:,CL,CR,R) = svdFrameReconstruct(svd.U_dff, permute(avgperiStimV,[2 3 1]) );
            MAPmove(:,:,:,CL,CR,R) = svdFrameReconstruct(svd.U_dff, permute(avgperiMoveV,[2 3 1]) );
        end
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Write aligned data to HDFS %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h5create(wfAlignedFile,'/MAPstim',size(MAPstim));
        h5write(wfAlignedFile,'/MAPstim',MAPstim);
        h5create(wfAlignedFile,'/MAPstimT',size(MAPstimT));
        h5write(wfAlignedFile,'/MAPstimT',MAPstimT);
        
        h5create(wfAlignedFile,'/MAPmove',size(MAPmove));
        h5write(wfAlignedFile,'/MAPmove',MAPmove);
        h5create(wfAlignedFile,'/MAPmoveT',size(MAPmoveT));
        h5write(wfAlignedFile,'/MAPmoveT',MAPmoveT);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% If ROI locations available, get epoch-aligned values for EVERY TRIAL %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %This section is badly designed. As it requires the ROI locations
        %which were only defined AFTER running this code previously. So if
        %you don't have the ROIs defined yet, you'll need to run this
        %widefield preproc, then run the ROI section, then rerun this
        %section
        roiFile = [ './preproc/WF_ROI/' sessionList.mouseName{sess} '.mat'];
        if exist(roiFile,'file')
            roi = load(roiFile); roi.px(1)=[];
            
            numTrials = length(stim_resp_id);
            ROIstim = nan(numTrials,length(roi.px),length(stimT)); 
            ROImove = nan(numTrials,length(roi.px),length(moveT)); 
            for r = 1:length(roi.px)
                thisU = squeeze(svd.U_dff(roi.px(r).row, roi.px(r).col, :));
                for n = 1:numTrials
                    ROIstim(n,r,:) = thisU'*squeeze(periStimV(n,:,:));
                    ROImove(n,r,:) = thisU'*squeeze(periMoveV(n,:,:));
                end
            end
            
            %Write
            h5create(wfAlignedFile,'/ROIstim',size(ROIstim));
            h5write(wfAlignedFile,'/ROIstim',ROIstim);
            h5create(wfAlignedFile,'/ROIstimT',size(stimT));
            h5write(wfAlignedFile,'/ROIstimT',stimT);
            
            
            h5create(wfAlignedFile,'/ROImove',size(ROImove));
            h5write(wfAlignedFile,'/ROImove',ROImove);
            h5create(wfAlignedFile,'/ROImoveT',size(moveT));
            h5write(wfAlignedFile,'/ROImoveT',moveT);
        end
    end
end

%% PREPROC: Identify bregma and ROIs for each mouse
clearvars -except sessionList
stack = cell(height(sessionList),1);
meanImgs = cell(height(sessionList),1);
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    wfSVDfile = [ './preproc/WF_SVD/' eRef '.mat'];
    svd = load(wfSVDfile,'meanImg_registered');

    wfAlignedFile = [ './preproc/WF_aligned/' eRef '.h5'];
    MAPstim = h5read(wfAlignedFile,'/MAPstim');
    MAPstimT = h5read(wfAlignedFile,'/MAPstimT');

    %rectify and baseline
    MAPstim(MAPstim<0)=0;
    MAPstim = MAPstim - MAPstim(:,:,MAPstimT==0,:,:,:);
    
    stack{sess} = MAPstim(:,:, :, 2,1,1); %Left choice on CL
    meanImgs{sess} = svd.meanImg_registered;
end

%Average sessions for each mouse
names = unique(sessionList.mouseName);
for mouse = 1:length(names)
    sessIDs = find(contains(sessionList.mouseName, names{mouse}));
    
    imL = mean( cat(4,stack{sessIDs}), 4);
    
    figure;
    for i = 1:9
        subplot(3,3,i);
        imagesc(imL(:,:,i));
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
    
    f1=figure;
    for i = 1:9
        subplot(3,3,i); hold on;
        imagesc(imL(:,:,i));
        addAllenCtxOutlines([px(1).row px(1).col], [px(1).row-1 px(1).col], [1 1 1], gca);
        for r = 1:length(ROIs)
            plot(px(r).col,px(r).row,'k+');
        end
    end
    set(get(gcf,'children'), 'clim', [-1 1]*quantile(abs(imL(:)),0.99), 'ydir','reverse')
        cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
    colormap(flipud(cmap));
    
%     addAllenCtxOutlines([px(1).row px(1).col], [px(1).row-1 px(1).col], [1 1 1], gca);
    outlines = findobj(f1,'tag','outline');
    
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
%% 0a) Load and prepare data for ONE session
clearvars -except sessionList
SESSIONID = 17;
eRef = sessionList.expRef{SESSIONID};
fprintf('Session %d %s\n',SESSIONID,eRef);

%Load meanImage
wfSVDfile = [ './preproc/WF_SVD/' eRef '.mat'];
svd = load(wfSVDfile,'meanImg_registered');

%Load behavioural data
behavFile = [ './preproc/BEHAV_RT_Corrected/' eRef '.mat'];
b = load(behavFile);

%Load ROI locations
roiFile = [ './preproc/WF_ROI/' sessionList.mouseName{SESSIONID} '.mat'];
roi = load(roiFile);

%Load epoch-aligned data
wfAlignedFile = [ './preproc/WF_aligned/' eRef '.h5'];
MAPstim = h5read(wfAlignedFile,'/MAPstim');
ROIstim = h5read(wfAlignedFile,'/ROIstim');
MAPstimT = h5read(wfAlignedFile,'/MAPstimT');
ROIstimT = h5read(wfAlignedFile,'/ROIstimT');

%Define mask around the outermost point of the allen atlas
bregma=roi.px(1);roi.px(1)=[];

f=figure;
ha = axes;
addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col],[1 1 1]*0,ha);
contours = get(ha,'children');
mask = zeros(size(svd.meanImg_registered));
for q = 1:length(contours)
    mask = mask | poly2mask(contours(q).XData, contours(q).YData, size(svd.meanImg_registered,1), size(svd.meanImg_registered,2));
end
mask = imgaussfilt(double(mask),3);
f.delete;

%Get inclTrials for ROI traces
contraIdx = b.contrastLeft>0 & b.contrastRight==0 & b.choice==1;
ipsiIdx = b.contrastLeft==0 & b.contrastRight>0 & b.choice==2;

%Normalise by peak contra stim firing
ROIstim_norm = ROIstim ./ max( mean(ROIstim(contraIdx,:,:),1) ,[], 3 );

disp('Data loaded');


%% 0b) Load and prepare data for all sessions of one mouse (session average)
clearvars -except sessionList
mouse = 'Hench';
sessIDs = find(contains(sessionList.mouseName,mouse));

%todo: concatenate across sessions and match RTs (exclude from map average
%sessions with significantly different left and right RTs).
MAPstimStack = cell(length(sessIDs),1); 
ROIstimStack = cell(length(sessIDs),1); 
ROIstimNormStack = cell(length(sessIDs),1); 
contraIpsiIdxStack = cell(length(sessIDs),1); 
for sess = 1:length(sessIDs)
    eRef = sessionList.expRef{sessIDs(sess)};
    fprintf('Session %d %s\n',sessIDs(sess),eRef);
    
    %Load meanImage
    wfSVDfile = [ './preproc/WF_SVD/' eRef '.mat'];
    svd = load(wfSVDfile,'meanImg_registered');
    
    %Load ROI locations
    roiFile = [ './preproc/WF_ROI/' mouse '.mat'];
    roi = load(roiFile);
    
    %Load epoch-aligned data
    wfAlignedFile = [ './preproc/WF_aligned/' eRef '.h5'];
    MAPstim = h5read(wfAlignedFile,'/MAPstim');
    ROIstim = h5read(wfAlignedFile,'/ROIstim');
    MAPstimT = h5read(wfAlignedFile,'/MAPstimT');
    ROIstimT = h5read(wfAlignedFile,'/ROIstimT');
    

    %Get inclTrials for ROI traces
    %Load behavioural data
    behavFile = [ './preproc/BEHAV_RT_Corrected/' eRef '.mat'];
    b = load(behavFile);
    contraIdx = b.contrastLeft>0 & b.contrastRight==0 & b.choice==1;
    ipsiIdx = b.contrastLeft==0 & b.contrastRight>0 & b.choice==2;
    ROIstim_norm = ROIstim ./ max( mean(ROIstim(contraIdx,:,:),1) ,[], 3 );

    bregma=roi.px(1);roi.px(1)=[];
   
    ROIstimStack{sess} = ROIstim;
    MAPstimStack{sess} = MAPstim;
     contraIpsiIdxStack{sess} = (contraIdx + 2*ipsiIdx);
end

%Stack ROIs
ROIstimStack = cat(1,ROIstimStack{:});
contraIpsiIdxStack = cat(1,contraIpsiIdxStack{:});

%Average map
MAPstim = nanmean(cat(7,MAPstimStack{:}),7);
clear MAPstimStack;

%Define mask
f=figure;
ha = axes;
addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col],[1 1 1]*0,ha);
contours = get(ha,'children');
mask = zeros(size(svd.meanImg_registered));
for q = 1:length(contours)
    mask = mask | poly2mask(contours(q).XData, contours(q).YData, size(svd.meanImg_registered,1), size(svd.meanImg_registered,2));
end
mask = imgaussfilt(double(mask),3);
f.delete;

disp('Data loaded');

%% 1) plot map aligned to stim onset
figure(101); 
set(gcf,'color','w');

ha = tight_subplot(3,length(MAPstimT),0,[0.7 0.02],0);
for i = 1:length(MAPstimT)
    imL = MAPstim(:,:,i,2,1,1).*mask;
    imagesc( ha(i), imL );
    addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*1, ha(i));
    hold(ha(i),'on'); plot(ha(i),bregma.col,bregma.row,'k.')
    
    imR = MAPstim(:,:,i,1,2,2).*mask;
    imagesc( ha(i + length(MAPstimT) ), imR );
    addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*1,ha(i + length(MAPstimT) ));
    hold(ha(i + length(MAPstimT) ),'on'); plot(ha(i + length(MAPstimT) ),bregma.col,bregma.row,'k.')
    
    imagesc( ha(i + 2*length(MAPstimT) ), imL-imR );
    addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*1,ha(i + 2*length(MAPstimT) ));
    hold(ha(i + 2*length(MAPstimT) ),'on'); plot(ha(i + 2*length(MAPstimT) ),bregma.col,bregma.row,'k.')
    
    title(ha(i), MAPstimT(i) );
end
cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
colormap(flipud(cmap));
set(ha,'clim',[-1 1]*quantile(abs(MAPstim(:)),0.97));
set(ha,'xcolor','none','ycolor','none','dataaspectratio',[1 1 1]);
%% 2) Plot traces and ROIs
areaCols = [      0    0.4470    0.7410; %blue
    0.3010    0.7450    0.9330; %lightblue
    0.4660    0.6740    0.1880; %green
    0.8500    0.3250    0.0980; %orange
    0.4940    0.1840    0.5560; %purple
    0.8 0.8 0.8]; %grey

figure(101);
ha = tight_subplot(5,1,0.01,[0.05 0.35],[0.05 0.66]);
for p =1:numel(roi.px)
    hold(ha(p),'on');
    
    if exist('ROIstimStack','var')
        avg = squeeze( mean( ROIstimStack(contraIpsiIdxStack==1,p,:) ,1) )';
%         err = mad( ROIstimStack(p,:,2,1,1,:),1,6);
        err = squeeze( std( ROIstimStack(contraIpsiIdxStack==1,p,:) ,[],1) )';
        
        fill(ha(p),[ROIstimT fliplr(ROIstimT)],[avg-err fliplr(avg+err)],'k','FaceColor',areaCols(p,:),'FaceAlpha',0.3,'EdgeAlpha',0 )
        plot(ha(p), ROIstimT, avg, 'k-', 'color', areaCols(p,:),'linewidth',3 );
        
        avg = squeeze( mean( ROIstimStack(contraIpsiIdxStack==2,p,:) ,1) )';
%         err = mad( ROIstimStack(p,:,2,1,1,:),1,6);
        err = squeeze( std( ROIstimStack(contraIpsiIdxStack==2,p,:) ,[],1) )';
        fill(ha(p),[ROIstimT fliplr(ROIstimT)],[avg-err fliplr(avg+err)],'k','FaceColor',areaCols(p,:),'FaceAlpha',0.3,'EdgeAlpha',0 )
        plot(ha(p), ROIstimT, avg, 'k-', 'color', areaCols(p,:),'linewidth',1 );
    else
        plot(ha(p), ROIstimT, squeeze(mean(ROIstim(contraIdx,p,:),1)), 'k-', 'color', areaCols(p,:),'linewidth',3 );
        plot(ha(p), ROIstimT, squeeze(mean(ROIstim(ipsiIdx,p,:),1)), 'k-', 'color', areaCols(p,:),'linewidth',1 );
    end
    set(ha(p),'xcolor','none','ycolor','none');
    line(ha(p), [0 0], [-1 1]*10, 'color', [0 0 0 0.1],'LineStyle','-');
 
    %Add text in top left corner
    tx=text(ha(p),-0.2,0.01,roi.px(p).name,'VerticalAlignment','top','color',areaCols(p,:),'FontWeight','bold');
end
line(ha(end), [0 0], [0 1], 'color', [0 0 0 0.1],'LineStyle','-');
line(ha(end), [0 0], [0 1], 'color', [0 0 0 0.1],'LineStyle','-');
linkaxes(ha,'y');
ylim(ha(1),[-0.001 0.011]);
set(ha(end),'xcolor','k','XTickLabelMode','auto','ycolor','k','yTickLabelMode','auto');

linkaxes(ha,'x');
xlim(ha(end),[-0.2 0.8]);
xlabel(ha(end),'Stim onset');
ylabel(ha(end),'d/dt ( dF/F )');

%meanimg and normalised traces
ha = tight_subplot(3,1,0.01,[0.05 0.35],[0.35 0.35]);
imagesc(ha(1),svd.meanImg_registered); colormap(ha(1),'gray'); hold(ha(1),'on');
addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*0.5, ha(1));
plot(ha(1),bregma.col,bregma.row,'k.')
for p =1:numel(roi.px)
    h = plot(ha(1),roi.px(p).col, roi.px(p).row, '.', 'Color', areaCols(p,:), 'markersize',25);
end
set(ha(1),'xcolor','none','ycolor','none', 'dataaspectratio', [1 1 1]);

hold(ha(2),'on'); hold(ha(3),'on');
for p =1:numel(roi.px)
    if exist('ROIstimStack_norm','var')
        avg = mean(ROIstimStack_norm(p,:,2,1,1,:),6);
         err = std( ROIstimStack_norm(p,:,2,1,1,:),[],6);
         
        fill(ha(2),[wf.stimT fliplr(wf.stimT)],[avg-err fliplr(avg+err)],'k','FaceColor',areaCols(p,:),'FaceAlpha',0.3,'EdgeAlpha',0 )
        hxx=plot(ha(2), wf.stimT, avg, '-','linewidth',3, 'color', areaCols(p,:));
        uistack(hxx,'bottom');
        
        avg = mean(ROIstimStack_norm(p,:,1,2,2,:),6);
        err = std( ROIstimStack_norm(p,:,1,2,2,:),[],6);
         
        fill(ha(3),[wf.stimT fliplr(wf.stimT)],[avg-err fliplr(avg+err)],'k','FaceColor',areaCols(p,:),'FaceAlpha',0.3,'EdgeAlpha',0 )
        hxx=plot(ha(3), wf.stimT, avg, '-','linewidth',1, 'color', areaCols(p,:));
        uistack(hxx,'bottom');
    else
        hxx=plot(ha(2), ROIstimT, squeeze(mean(ROIstim_norm(contraIdx,p,:),1)), '-','linewidth',3, 'color', areaCols(p,:));
        uistack(hxx,'bottom');
        hxx=plot(ha(3), ROIstimT, squeeze(mean(ROIstim_norm(ipsiIdx,p,:),1)), '-','linewidth',1, 'color', areaCols(p,:));
        uistack(hxx,'bottom');
    end
end
line(ha(2), [0 0], [0 1], 'color', [0 0 0 0.1],'LineStyle','-');
line(ha(3), [0 0], [0 1], 'color', [0 0 0 0.1],'LineStyle','-');
set(ha(2:3),'ylim',[-0.1 1.1],'xlim',[-0.2 0.8]);

%% 3) Plot latencies over all sessions TODO
mice={'Chain','Radnitz','Cori','Reichstein','Hench'};

figure(101);
onsetTime = nan(height(sessionList),5,2);
rt_vals = cell(height(sessionList),1);

%Get all onset times
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    %Load meanImage
    wfSVDfile = [ './preproc/WF_SVD/' eRef '.mat'];
    svd = load(wfSVDfile,'meanImg_registered');
    
    %Load behavioural data
    behavFile = [ './preproc/BEHAV_RT_Corrected/' eRef '.mat'];
    b = load(behavFile);
    
    %Load epoch-aligned data
    wfAlignedFile = [ './preproc/WF_aligned/' eRef '.h5'];
    ROIstim = h5read(wfAlignedFile,'/ROIstim');
    ROIstimT = h5read(wfAlignedFile,'/ROIstimT');

    %Get inclTrials for ROI traces
    contraIdx = b.contrastLeft>0 & b.contrastRight==0 & b.choice==1;
    ipsiIdx = b.contrastLeft==0 & b.contrastRight>0 & b.choice==2;
    ROIstim_norm = ROIstim ./ max( mean(ROIstim(contraIdx,:,:),1) ,[], 3 );
    
    rt_vals{sess} = b.RT;
 
    [~,onsetIdx]=max( mean(ROIstim_norm(contraIdx,:,:),1) >=0.25,[],3);
    onsetTime(sess,:,1) = ROIstimT(onsetIdx);
    
    [~,onsetIdx]=max( mean(ROIstim_norm(ipsiIdx,:,:),1) >=0.25,[],3);
    onsetTime(sess,:,2) = ROIstimT(onsetIdx);
end

ha = tight_subplot(1,1,0.01,[0.05 0.35],[0.67 0.01]);
hold(ha(1),'on');
Y_jitter = linspace(-0.2,0.2,length(mice));
allAve = nan(length(mice),5,2);
for m = 1:length(mice)
    mIdx = strcmp(sessionList.mouseName,mice(m));
%     lat = IpsiOnsetTime(mIdx,:);
    lat = onsetTime(mIdx,:,:);
    
    avg = median(lat,1);
    err = mad(lat,1);
    
    allAve(m,:,:)=avg;
    
    rts_subj = cat(1,rt_vals{mIdx});
    
    for p =1:5
        yline = Y_jitter(m) + [-1 1]*0.05;
        
        line(ha(1),[1 1]*avg(1,p,1), yline,'Color',areaCols(p,:), 'linewidth',3);

        fx=fill(ha(1), [avg(1,p,1)-err(1,p,1), avg(1,p,1)+err(1,p,1), avg(1,p,1)+err(1,p,1), avg(1,p,1)-err(1,p,1)], [yline(1) yline(1) yline(2) yline(2)], 'k',...
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

%add median and MAD across mice
avg=median(allAve,1);
err=mad(allAve,1,1);
for p = 1:5
    plot(ha(1), avg(1,p,1), -0.3, '.', 'markersize', 20, 'Color', areaCols(p,:));
    line(ha(1), [avg(1,p,1)-err(1,p,1) avg(1,p,1)+err(1,p,1)], [1 1]*-0.3, 'Color', areaCols(p,:), 'linewidth',2);

    plot(ha(1), avg(1,p,2), -0.32, 'o', 'markersize', 7, 'Color', areaCols(p,:));
    line(ha(1), [avg(1,p,2)-err(1,p,2) avg(1,p,2)+err(1,p,2)], [1 1]*-0.32, 'Color', areaCols(p,:), 'linewidth',2);
end
        
xlim(ha,[0 0.32]);
xlabel(ha(1),'Stim onset');
set(ha,'ytick',Y_jitter,'YTickLabel',mice,'XTickLabelMode','auto');
ylim(ha,[-0.35 0.3]);

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
clearvars -except sessionList

mouse = 'Hench';
refImg = load(['./preproc/reference_image/' mouse '.mat'], 'meanImg');
sessIDs = find(contains(sessionList.mouseName,mouse));

stimSliceTimes = [0 0.07 0.12 0.2 0.25];
stimSliceTimes = linspace(0,0.25,6);

moveSliceTimes = [-0.2 -0.1 0 0.1 0.2];
moveSliceTimes = linspace(-0.1,0.2,6);

MAPstimStack = nan(size(refImg.meanImg,1), size(refImg.meanImg,2), length(stimSliceTimes),2,2,3,length(sessIDs));
MAPmoveStack = nan(size(refImg.meanImg,1), size(refImg.meanImg,2), length(moveSliceTimes),2,2,3,length(sessIDs));
for sess = 1:length(sessIDs)
    eRef = sessionList.expRef{sessIDs(sess)};
    fprintf('Session %d %s\n',sessIDs(sess),eRef);
    wfFile = [ './preproc/WF_SVD/' eRef '.mat'];
    behavFile = [ './preproc/BEHAV/' eRef '.mat'];
    roiFile = [ './preproc/WF_ROI/' mouse '.mat'];
    wfAlignedFile = [ './preproc/WF_aligned/' eRef '.mat'];
    wfAlignedFile = [ './preproc/WF_aligned/' eRef '.h5'];
    
    wf = load(wfFile);
    b = load(behavFile);
    roi = load(roiFile);
    align = load(wfAlignedFile);
    MAPstim = h5read(wfAlignedFile,'/MAPstim');
    MAPmove = h5read(wfAlignedFile,'/MAPmove');
    
    bregma=roi.px(1);
    roi.px(1)=[];
    
    %Rectify and baseline stim-aligned maps
    MAPstim(MAPstim<0)=0;
    MAPstim = MAPstim - MAPstim(:,:,1,:,:,:);
    
%     %baseline mov-aligned maps
%     MAPstim = MAPmove - MAPmove(:,:,1,:,:,:);
    [~,stimSliceIdx]=find(abs(align.stimT - stimSliceTimes')==min(abs(align.stimT - stimSliceTimes'),[],2));
    [~,moveSliceIdx]=find(abs(align.moveT - moveSliceTimes')==min(abs(align.moveT - moveSliceTimes'),[],2));

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





