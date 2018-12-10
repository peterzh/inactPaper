%%
sessionList = readtable('./sessionList.csv','FileType','text','Delimiter',',');
areaCols = [      0    0.4470    0.7410; %blue
    0.3010    0.7450    0.9330; %lightblue
    0.4660    0.6740    0.1880; %green
    0.8500    0.3250    0.0980; %orange
    0.4940    0.1840    0.5560; %purple
    0.8 0.8 0.8]; %grey
% 
% for sess = 1:height(sessionList)
%     eRef = sessionList.expRef{sess};
%     [D,meta] = loadData(eRef);
% %     p = load(strrep(meta.blockFile, 'Block', 'parameters'));
% %     fprintf('Inter-trial delay: %0.1f\n', p.parameters.interTrialDelay);
% %     fprintf('Pre-stim quiesence: %0.1f - %0.1f \n', p.parameters.preStimQuiescentPeriod(1),  p.parameters.preStimQuiescentPeriod(2));
% %     fprintf('Post-stim delay to go cue: %0.1f - %0.1f \n', p.parameters.cueInteractiveDelay(1),  p.parameters.cueInteractiveDelay(2));
% %     fprintf('Negative feedback sound duration: %0.1f  \n', p.parameters.negFeedbackSoundDuration);
% %     fprintf('\n');
% end

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
        
        inter_tsteps = linspace(0, 3, 1000);
        for tr = 1:length(timings.t_stimOn)
            idx = find(timings.t_stimOn(tr) < timings.t_moveOnsets & timings.t_moveOnsets < timings.t_responseTime(tr),1);
            if trial.choice(tr)<3 && ~isempty(idx)
                mon(tr) = timings.t_moveOnsets( idx );
            end
            
            %If nogo trial, see if there was a flinch or incomplete
            %movement
            if trial.choice(tr)==3
                idx = timings.t_stimOn(tr) < tt & tt < timings.t_responseTime(tr);
                tt_temp = tt(idx); tt_temp = tt_temp-tt_temp(1);
                wh_temp = wh(idx); wh_temp = wh_temp-wh_temp(1);
                plot( tt_temp, abs(wh_temp))

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

MAPstimT = [0 30 60 90 120 150 180 210]/1000; %intended times
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
%         stim_resp_set = [b.contrastLeft>0 b.contrastRight>0 b.choice];
        stim_resp_set = [b.contrastLeft b.contrastRight b.choice];
        CLS = unique(b.contrastLeft);
        CRS = unique(b.contrastRight);
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
        [~,stimsliceIdx]=max(stimT' >= MAPstimT,[],1);
        [~,movesliceIdx]=max(moveT' >= MAPmoveT,[],1);
        MAPstimT = stimT(stimsliceIdx); %actual slice times
        MAPmoveT = moveT(movesliceIdx); %actual slice times
        
%         MAPstim = nan(size(svd.U_dff,1), size(svd.U_dff,2), length(MAPstimT), 2, 2, 3);
%         MAPmove = nan(size(svd.U_dff,1), size(svd.U_dff,2), length(MAPmoveT), 2, 2, 3);
        MAPstim = nan(size(svd.U_dff,1), size(svd.U_dff,2), length(MAPstimT), length(CLS), length(CRS), 3);
        MAPmove = nan(size(svd.U_dff,1), size(svd.U_dff,2), length(MAPmoveT), length(CLS), length(CRS), 3);
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
            
            
%             MAPstim(:,:,:,CL,CR,R) = svdFrameReconstruct(svd.U_dff, permute(avgperiStimV,[2 3 1]) );
%             MAPmove(:,:,:,CL,CR,R) = svdFrameReconstruct(svd.U_dff, permute(avgperiMoveV,[2 3 1]) );
            MAPstim(:,:,:,stim_resp(cond,1)==CLS,stim_resp(cond,2)==CRS,stim_resp(cond,3)) = svdFrameReconstruct(svd.U_dff, permute(avgperiStimV,[2 3 1]) );
            MAPmove(:,:,:,stim_resp(cond,1)==CLS,stim_resp(cond,2)==CRS,stim_resp(cond,3)) = svdFrameReconstruct(svd.U_dff, permute(avgperiMoveV,[2 3 1]) );
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

%     stack{sess} = MAPstim(:,:, :, 2,1,1); %Left choice on CL
    stack{sess} = cat(3,MAPstim(:,:, :, 2,1,1),MAPstim(:,:, :, 1,2,2)); 
    meanImgs{sess} = svd.meanImg_registered;
end

%Average sessions for each mouse
names = unique(sessionList.mouseName);
for mouse = 5:length(names)
    sessIDs = find(contains(sessionList.mouseName, names{mouse}));
    
    imL = nanmean( cat(4,stack{sessIDs}), 4);
    
    figure;
    for i = 1:size(imL,3)
        subplot(4,4,i);
        imagesc(imL(:,:,i));
    end
    set(get(gcf,'children'), 'clim', [-1 1]*quantile(abs(imL(:)),0.99))
    colormap(BlueWhiteRed(100,1));
    
    px = struct('name',[],'row',[],'col',[]);
    ROIs = {'Bregma','RVISp','RVISal','RMOs','RMOp','RSSp','LVISp','LVISal','LMOs','LMOp','LSSp'};
    for r = 1:length(ROIs)
        px(r).name = ROIs{r};
        disp(px(r).name);
        set(gcf,'name',ROIs{r});
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
    colormap(BlueWhiteRed(100,1));
    
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
MAPmoveStack = cell(length(sessIDs),1); 
ROIstimStack = cell(length(sessIDs),1); 
ROIstimNormStack = cell(length(sessIDs),1); 
contraIpsiIdxStack = cell(length(sessIDs),1); 
allBehav = struct; D = struct;
for sess = 1:length(sessIDs)
    eRef = sessionList.expRef{sessIDs(sess)};
    fprintf('Session %d %s\n',sessIDs(sess),eRef);
    
    %Load meanImage
    wfSVDfile = [ './preproc/WF_SVD/' eRef '.mat'];
    svd = load(wfSVDfile,'meanImg_registered');
    
    %Load ROI locations
    roiFile = [ './preproc/WF_ROI/' mouse '.mat'];
    roi = load(roiFile);
    bregma=roi.px(1);roi.px(1)=[];
    
    %Load epoch-aligned data
    wfAlignedFile = [ './preproc/WF_aligned/' eRef '.h5'];
    MAPstim = h5read(wfAlignedFile,'/MAPstim');
    MAPmove = h5read(wfAlignedFile,'/MAPmove');
    ROIstim = h5read(wfAlignedFile,'/ROIstim');
    MAPstimT = h5read(wfAlignedFile,'/MAPstimT');
    MAPmoveT = h5read(wfAlignedFile,'/MAPmoveT');
    ROIstimT = h5read(wfAlignedFile,'/ROIstimT');
    
    %Get inclTrials for ROI traces
    %Load behavioural data
    behavFile = [ './preproc/BEHAV_RT_Corrected/' eRef '.mat'];
    b = load(behavFile);
    allBehav = addstruct(allBehav, b);
    
    contraIdx = b.contrastLeft>0 & b.contrastRight==0 & b.choice==1;
    ipsiIdx = b.contrastLeft==0 & b.contrastRight>0 & b.choice==2;

   
    ROIstimStack{sess} = ROIstim;
    MAPstimStack{sess} = MAPstim;
    MAPmoveStack{sess} = MAPmove;   
    contraIpsiIdxStack{sess} = (contraIdx + 2*ipsiIdx);
     
    D = addstruct(D, loadData(eRef));
end

%Stack ROIs
ROIstimStack = cat(1,ROIstimStack{:});
contraIpsiIdxStack = cat(1,contraIpsiIdxStack{:});

%Average map
MAPstim = nanmean(cat(7,MAPstimStack{:}),7);
MAPmove = nanmean(cat(7,MAPmoveStack{:}),7);

clear MAPstimStack MAPmoveStack;

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
pixSize = 0.0217; mmBar = 1/pixSize;
line(400 + [0 mmBar], [50 50]);
set(gca,'dataaspectratio',[1 1 1]);
f.delete;

disp('Data loaded');

%% 0) Plot misc maps
figure('color','w');
% ha = tight_subplot(4,length(MAPstimT),0,0.05,0);
ha = tight_subplot(length(MAPstimT),4,0,0.05,0.05);
numConds = 4;
for i = 1:length(MAPstimT)
    
    %C0 and NOGO
    idx = numConds*(i-1) + 1;
    im = MAPstim(:,:,i,1,1,3).*mask; 
    imx=imagesc( ha(idx), im ); alpha(imx,mask);
    addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*0.5, ha(idx));
    hold(ha(idx),'on'); plot(ha(idx),bregma.col,bregma.row,'k.');
    title(ha(idx), MAPstimT(i) );
        
    %CL and NOGO
    idx = numConds*(i-1) + 2;
    im = nanmean(MAPstim(:,:,i,2:end,1,3),4).*mask;
    imx=imagesc( ha(idx), im ); alpha(imx,mask);
    addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*0.5, ha(idx));
    hold(ha(idx),'on'); plot(ha(idx),bregma.col,bregma.row,'k.');

    %CL and LEFT
    idx = numConds*(i-1) + 3;
    im = nanmean(MAPstim(:,:,i,2:end,1,1),4).*mask;
    imx=imagesc( ha(idx), im ); alpha(imx,mask);
    addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*0.5, ha(idx));
    hold(ha(idx),'on'); plot(ha(idx),bregma.col,bregma.row,'k.');
    
    
    %CL and Left - CR and Right
    idx = numConds*(i-1) + 4;
    im1 = nanmean(MAPstim(:,:,i,2:end,1,1),4);
    im2 = nanmean(MAPstim(:,:,i,1,2:end,2),5);
    im = (im1 - im2).*mask;
    imx=imagesc( ha(idx), im ); alpha(imx,mask);
    addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*0.5, ha(idx));
    hold(ha(idx),'on'); plot(ha(idx),bregma.col,bregma.row,'k.');
    
    
%     %CL and NOGO - %C0 and NOGO
%     idx = i + 3*length(MAPstimT);
%     im1 = nanmean(MAPstim(:,:,i,2:end,1,3),4); 
%     im2 = MAPstim(:,:,i,1,1,3);
%     im = (im1 - im2).*mask;
%     imx=imagesc( ha(idx), im ); alpha(imx,mask);
%     addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*0.5, ha(idx));
%     hold(ha(idx),'on'); plot(ha(idx),bregma.col,bregma.row,'k.');
%     
%     
%     %CL and LEFT - %CL and NOGO
%     idx = i + 4*length(MAPstimT);
%     im1 = nanmean(MAPstim(:,:,i,2:end,1,1),4);
%     im2 = nanmean(MAPstim(:,:,i,2:end,1,3),4);
%     im = (im1 - im2).*mask;
%     imx=imagesc( ha(idx), im ); alpha(imx,mask);
%     addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*0.5, ha(idx));
%     hold(ha(idx),'on'); plot(ha(idx),bregma.col,bregma.row,'k.');
    
%     %CR and NOGO
%     idx = i + 3*length(MAPstimT);
%     im = nanmean(MAPstim(:,:,i,1,2:end,3),5).*mask;
%     imx=imagesc( ha(idx), im ); alpha(imx,mask);
%     addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*0.5, ha(idx));
%     hold(ha(idx),'on'); plot(ha(idx),bregma.col,bregma.row,'k.');
% 
%     %CR and RIGHT
%     idx = i + 4*length(MAPstimT);
%     im = nanmean(MAPstim(:,:,i,1,2:end,2),5).*mask;
%     imx=imagesc( ha(idx), im ); alpha(imx,mask);
%     addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*0.5, ha(idx));
%     hold(ha(idx),'on'); plot(ha(idx),bregma.col,bregma.row,'k.');
%     
    

%     %CL and LEFT - CL and NOGO
%     im1 = nanmean(MAPstim(:,:,i,2:end,1,1),4);
%     im2 = nanmean(MAPstim(:,:,i,2:end,1,3),4);
%     im3L = (im1 - im2).*mask;
%     idx = i + 3*length(MAPstimT);
%     imx= imagesc( ha(idx), im3L );  alpha(imx,mask);
%     addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*0.5,ha(idx));
%     hold(ha(idx),'on'); plot(ha(idx),bregma.col,bregma.row,'k.')
%     
%     %CR and RIGHT - CR and NOGO
%     im1 = nanmean(MAPstim(:,:,i,1,2:end,2),5);
%     im2 = nanmean(MAPstim(:,:,i,1,2:end,3),5);
%     im3R = (im1 - im2).*mask;
%     idx = i + 4*length(MAPstimT);
%     imx= imagesc( ha(idx), im3R );  alpha(imx,mask);
%     addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*0.5,ha(idx));
%     hold(ha(idx),'on'); plot(ha(idx),bregma.col,bregma.row,'k.')
% 
%     im = im3R - im3L;
%     idx = i + 5*length(MAPstimT);
%     imx= imagesc( ha(idx), im );  alpha(imx,mask);
%     addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*0.5,ha(idx));
%     hold(ha(idx),'on'); plot(ha(idx),bregma.col,bregma.row,'k.')


end
colormap(BlueWhiteRed(100,1));
set(ha,'clim',[-1 1]*quantile(abs(MAPstim(:)),0.98));
set(ha,'xcolor','none','ycolor','none','dataaspectratio',[1 1 1]);
set(ha,'color',[1 1 1]*0.8);





figure('color','w');
% ha = tight_subplot(4,length(MAPmoveT),0,0.05,0);
ha = tight_subplot(length(MAPmoveT),4,0,0.05,0.05);
numConds = 4;
for i = 1:length(MAPmoveT)
    
    %C0 and NOGO
    idx = numConds*(i-1) + 1;
    im = MAPmove(:,:,i,1,1,3).*mask; 
    imx=imagesc( ha(idx), im ); alpha(imx,mask);
    addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*0.5, ha(idx));
    hold(ha(idx),'on'); plot(ha(idx),bregma.col,bregma.row,'k.');
    title(ha(idx), MAPmoveT(i) );

        
    %CL and NOGO
    idx = numConds*(i-1) + 2;
    im = nanmean(MAPmove(:,:,i,2:end,1,3),4).*mask;
    imx=imagesc( ha(idx), im ); alpha(imx,mask);
    addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*0.5, ha(idx));
    hold(ha(idx),'on'); plot(ha(idx),bregma.col,bregma.row,'k.');

    %CL and LEFT
    idx = numConds*(i-1) + 3;
    im = nanmean(MAPmove(:,:,i,2:end,1,1),4).*mask;
    imx=imagesc( ha(idx), im ); alpha(imx,mask);
    addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*0.5, ha(idx));
    hold(ha(idx),'on'); plot(ha(idx),bregma.col,bregma.row,'k.');
    
    
    %CL and Left - CR and Right
    idx = numConds*(i-1) + 4;
    im1 = nanmean(MAPmove(:,:,i,2:end,1,1),4);
    im2 = nanmean(MAPmove(:,:,i,1,2:end,2),5);
    im = (im1 - im2).*mask;
    imx=imagesc( ha(idx), im ); alpha(imx,mask);
    addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*0.5, ha(idx));
    hold(ha(idx),'on'); plot(ha(idx),bregma.col,bregma.row,'k.');

end
colormap(BlueWhiteRed(100,1));
set(ha,'clim',[-1 1]*quantile(abs(MAPmove(:)),0.98));
set(ha,'xcolor','none','ycolor','none','dataaspectratio',[1 1 1]);
set(ha,'color',[1 1 1]*0.8);



%% 1) plot map aligned to stim onset

figure(101); 
set(gcf,'color','w');
ha = tight_subplot(2,length(MAPstimT),0,[0.7 0.02],0);
for i = 1:length(MAPstimT)
    im = nanmean(MAPstim(:,:,i,2:end,1,3),4).*mask;
    imx=imagesc( ha(i), im ); alpha(imx,mask);
    addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*0.5, ha(i));
    hold(ha(i),'on'); plot(ha(i),bregma.col,bregma.row,'k.');
    
    im = nanmean(MAPstim(:,:,i,2:end,1,1),4).*mask;
    imx= imagesc( ha(i + length(MAPstimT) ), im );  alpha(imx,mask);
    addAllenCtxOutlines([bregma.row bregma.col], [bregma.row-1 bregma.col], [1 1 1]*0.5,ha(i + length(MAPstimT) ));
    hold(ha(i + length(MAPstimT) ),'on'); plot(ha(i + length(MAPstimT) ),bregma.col,bregma.row,'k.')

    title(ha(i), MAPstimT(i) );
end
colormap(BlueWhiteRed(100,1));
set(ha,'clim',[-1 1]*quantile(abs(MAPstim(:)),0.98));
set(ha,'xcolor','none','ycolor','none','dataaspectratio',[1 1 1]);
set(ha,'color',[1 1 1]*0.8)
%% 2) Plot traces and ROIs
areaCols = [      0    0.4470    0.7410; %blue
    0.3010    0.7450    0.9330; %lightblue
    0.4660    0.6740    0.1880; %green
    0.8500    0.3250    0.0980; %orange
    0.4940    0.1840    0.5560; %purple
    0.8 0.8 0.8]; %grey

figure(101);
ha = tight_subplot(5,1,0.01,[0.05 0.35],[0.05 0.66]);
for p =1:5
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
for p =1:5
    h = plot(ha(1),roi.px(p).col, roi.px(p).row, '.', 'Color', areaCols(p,:), 'markersize',25);
end
set(ha(1),'xcolor','none','ycolor','none', 'dataaspectratio', [1 1 1]);

hold(ha(2),'on'); hold(ha(3),'on');
for p =1:5
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
areaCols = [      0    0.4470    0.7410; %blue
    0.3010    0.7450    0.9330; %lightblue
    0.4660    0.6740    0.1880; %green
    0.8500    0.3250    0.0980; %orange
    0.4940    0.1840    0.5560; %purple
    0.8 0.8 0.8]; %grey
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
    ROIstim = ROIstim(:,1:5,:); %get only RIGHT 5 ROIs
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

%save onset times
save('onsetTime.mat','onsetTime');

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
    
%     line(ha(1), rt_stats([1,3]), [1 1]*Y_jitter(m), 'Color', 'k');
%     plot(ha(1), rt_stats(2), Y_jitter(m), 'k.', 'markersize',15);
    
    plot(ha(1), nanmedian(rt_stats) ,Y_jitter(m), 'k.', 'markersize',15);
     line(ha(1),nanmedian(rt_stats) + [-1 1]*mad(rt_stats,1), [1 1]*Y_jitter(m), 'Color', 'k');
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

%plot ipsi vs contra statistical test results
figure;
axes; hold on;
for region = 1:5
    plot( region + [-0.1 0.1], squeeze(avg(1,region,:)), '.-', 'markersize', 20, 'Color', areaCols(region,:));
    lx=line( [1 1]*(region-0.1), avg(1,region,1)+[-1 1]*err(1,region,1), 'color', areaCols(region,:), 'linewidth',2);
    lx=line( [1 1]*(region+0.1), avg(1,region,2)+[-1 1]*err(1,region,2), 'color', areaCols(region,:), 'linewidth',2);
    
    disp(['region', num2str(region)]);
    [H,P,CI,STATS]=ttest( allAve(:,region,1), allAve(:,region,2) );
    text(region, 0.1, num2str(P));
end
set(gca,'xtick',sort([(1:5)-0.1 (1:5)+0.1]));
ylabel('Latency');

figure;
axes; hold on;
%another kind of latency plot
avg=median(onsetTime(:,:,1),1);
err=mad(onsetTime(:,:,1),1,1);
for region = 1:5
    
    fx=fill([avg(region)-err(region) avg(region)+err(region) avg(region)+err(region) avg(region)-err(region) ],[0 0 length(mice)+1 length(mice)+1],'k');
    fx.FaceAlpha=0.3;
    fx.EdgeAlpha=0;
    fx.FaceColor = areaCols(region,:);

    for m = 1:length(mice)
        mIdx = strcmp(sessionList.mouseName,mice(m));
        lat = onsetTime(mIdx,region,1);
        avg_mouse = median(lat,1);
        err_mouse = mad(lat,1,1);
        
        plot( avg_mouse, m, '.', 'markersize',20,'color',areaCols(region,:));
        line( avg_mouse + [-1 1]*err_mouse, [1 1]*m, 'Color',areaCols(region,:) )
    end
end


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
colormap(BlueWhiteRed(100,1));

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
colormap(BlueWhiteRed(100,1));

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
colormap(BlueWhiteRed(100,1));


%% Fit cross-predict model
 
%1) convert widefield to ephys firing using translation

%load dual dataset
fr = load("C:\Users\Peter\Documents\MATLAB\inactPaper\widefield\mechanistic_fit\visMOs_allConds.mat");

VIS_window = [0.075 0.125];
MOS_window = [0.125 0.175];
wf_window_delay = 0.03;

%get ephys data at critical timewindow
vis = mean( mean( fr.vis(:,:,:,VIS_window(1) < fr.binsS & fr.binsS < VIS_window(2)), 4), 3);
m2 = mean( mean( fr.mos(:,:,:,MOS_window(1) < fr.binsS & fr.binsS < MOS_window(2)), 4), 3);

figure('color','w');
ha = tight_subplot(1,2,0.1,0.1,0.1);
imagesc(ha(1), vis); 
imagesc(ha(2), m2);title(ha(2),'LM2 ephys');
set(ha,'dataaspectratio',[1 1 1]);
xlabel(ha(1),'CL'); ylabel(ha(1),'CR'); title(ha(1),'LVIS ephys');
set(ha,'xtick', 1:4, 'ytick', 1:4, 'xticklabel', fr.ucl, 'yticklabel', fr.ucl, 'ydir','normal');

%get widefield data at critical timewindow
dat = struct;
[names,~,mouseID]=unique(sessionList.mouseName);
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    %Load behavioural data
    behavFile = [ './preproc/BEHAV_RT_Corrected/' eRef '.mat'];
    b = load(behavFile);
    
    %Load meanImage
    wfSVDfile = [ './preproc/WF_SVD/' eRef '.mat'];
    svd = load(wfSVDfile,'meanImg_registered');
    
    %Load ROI locations
    roiFile = [ './preproc/WF_ROI/' sessionList.mouseName{sess} '.mat'];
    roi = load(roiFile); roi.px(1)=[];
    
    areaIdx = [find(contains({roi.px.name},'LVISp')),...
            find(contains({roi.px.name},'RVISp')),...
            find(contains({roi.px.name},'LMOs')),...
            find(contains({roi.px.name},'RMOs'))];   
    
    %Load epoch-aligned data
    wfAlignedFile = [ './preproc/WF_aligned/' eRef '.h5'];
    ROIstim = h5read(wfAlignedFile,'/ROIstim');
    ROIstimT = h5read(wfAlignedFile,'/ROIstimT');
    
    ROIstim(ROIstim<0) = 0; %Rectify
    
    %Get ROI only at L and R VIS and MOs at the times relevant
    avgF = [];
    avgF(:,1:2) = mean( ROIstim(:, areaIdx(1:2), (VIS_window(1)+wf_window_delay) <= ROIstimT & ROIstimT <= (VIS_window(2)+wf_window_delay)), 3);
    avgF(:,3:4) = mean( ROIstim(:, areaIdx(3:4), (MOS_window(1)+wf_window_delay) <= ROIstimT & ROIstimT <= (MOS_window(2)+wf_window_delay)), 3);

%     avgF(avgF<0)=0; % rectify
    
    row=struct;
    row.firing_rate = avgF;
    row.choice = b.choice;
    row.subjectID = ones(size(row.choice))*mouseID(sess);
    row.sessionID = ones(size(row.choice))*sess;
    row.contrastLeft = b.contrastLeft;
    row.contrastRight = b.contrastRight;
    row.RT = b.RT;
    dat = addstruct(dat, row);
end

%create average widefield firing for the 16 contrast conditions
vis_wf = nan(size(vis));
m2_wf = nan(size(m2));
for cl = 1:length(fr.ucl)
    for cr = 1:length(fr.ucr)
        vis_wf(cr,cl) = mean( dat.firing_rate( dat.contrastLeft==fr.ucl(cl) & dat.contrastRight==fr.ucr(cr) ,1) ,1); %Left VIS
        m2_wf(cr,cl) = mean( dat.firing_rate( dat.contrastLeft==fr.ucl(cl) & dat.contrastRight==fr.ucr(cr) ,3) ,1);%Left MOS
    end
end
figure('color','w');
subplot(1,2,1); plot(vis_wf(:),vis(:),'ko'); xlabel('VIS widefield (dF/F)'); ylabel('VIS ephys (sp/sec)'); lsline; axis square;
subplot(1,2,2); plot(m2_wf(:),m2(:),'ko'); xlabel('M2 widefield (dF/F)'); ylabel('M2 ephys (sp/sec)');lsline; axis square;
set(get(gcf,'children'),'box','off');

%Fit translation
vis_fit = polyfit( vis_wf(:), vis(:), 1);
m2_fit = polyfit(m2_wf(:), m2(:), 1 );

%convert widefield firing to ephys firing using fits
dat.firing_rate(:,1:2) = vis_fit(2) + vis_fit(1)*dat.firing_rate(:,1:2);
dat.firing_rate(:,3:4) = m2_fit(2) + m2_fit(1)*dat.firing_rate(:,3:4);

% %rectify
dat.firing_rate(dat.firing_rate<0) = 0;

%2) Fit mechanistic model
fit_mech = bfit.fitModel('Two-Level-Mechanistic',dat);

% Fit original phenomenological model:
fit_phen = bfit.fitModel('Two-Level',dat);
%need to fix model code...

%% Plot model log2 likelihood for each contrast condition

likelihood = [];
likelihood(1,:) = exp(mean(fit_mech.posterior.log_lik,1));
likelihood(2,:) = exp(mean(fit_phen.posterior.log_lik,1));


CL = unique(fit_mech.data.contrastLeft);
CR = unique(fit_mech.data.contrastRight);

CL = [0 0.25 0.5 1];
CR = [0 0.25 0.5 1];

figure;
lik = nan(length(CL),length(CR),2);

for i = 1:2
    for cl = 1:length(CL)
        for cr = 1:length(CR)
            idx = fit_mech.data.contrastLeft == CL(cl) & fit_mech.data.contrastRight == CR(cr);
            lik(cl,cr,i) = mean( log2(likelihood(i,idx)) );
        end
    end
    
    subplot(1,2,i);
    imagesc(lik(:,:,i)); 
    set(gca,'ydir','normal','xtick',1:4,'xticklabel',CL,'ytick',1:4,'yticklabel',CR);
end

figure;
imagesc( lik(:,:,1) - lik(:,:,2) );
set(gca,'ydir','normal','xtick',1:4,'xticklabel',CL,'ytick',1:4,'yticklabel',CR);
caxis([-1 1]*0.3)

%predicted probability for every trial based on the grand average
ZL = mean(fit_phen.posterior.logOdds(:,:,1),1);
ZR = mean(fit_phen.posterior.logOdds(:,:,2),1);
pL = exp(ZL)./(1+exp(ZL) + exp(ZR));
pR = exp(ZR)./(1+exp(ZL) + exp(ZR));
pNG = 1./(1+exp(ZL) + exp(ZR));
p = [pL; pR; pNG];
[~,resp]=max(p,[],1);

pCorr = nan(length(CL),length(CR));
for cl = 1:length(CL)
    for cr = 1:length(CR)
        idx = fit_mech.data.contrastLeft == CL(cl) & fit_mech.data.contrastRight == CR(cr);
        pCorr(cl,cr) = mean( fit_mech.data.choice(idx) == resp(idx)' );
    end
end

set(gca,'ydir','normal','xtick',1:4,'xticklabel',CL,'ytick',1:4,'yticklabel',CR);
caxis([-1 1]*0.3)

%% Plot cross-prediction model on widefield data
% Model with old timings
fit_mech = load("C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\tp03f24ce4_5470_4181_b95f_cfaf8abe5b33.mat");

%Model with earlier VIS time window
% fit_mech = load("C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\tp711b3be9_2238_439a_8bf5_1e3a1ce710fa.mat");

fit_mech.data.pedestal = min(fit_mech.data.contrastLeft, fit_mech.data.contrastRight);
fit_mech.data.cDiff = fit_mech.data.contrastRight - fit_mech.data.contrastLeft;
%%%%%%%%%%%%%%%%%%%%%%%% PLOT CURVES
[counts,~,~,labels] = crosstab(fit_mech.data.pedestal,...
    fit_mech.data.cDiff,...
    fit_mech.data.choice,...
    fit_mech.data.sessionID);
prob = counts./sum(counts,3);%Convert to probability over choices
prob = nanmean(prob,4);
cdU = unique(fit_mech.data.cDiff);

colours = [ 0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.4940    0.1840    0.5560];


%Plot simulated model fit
CL = [linspace(1,0,200), zeros(1,200)];
CR = [zeros(1,200), linspace(0,1,200)];

%Interpolate the ephys data to get the firing rate for each of the contrast
%values
F=[];
F(1,:) = interp2(fr.ucl, fr.ucr, vis, CL , CR, 'linear' );
F(2,:) = interp2(fr.ucl, fr.ucr, vis', CL , CR, 'linear' );
F(3,:) = interp2(fr.ucl, fr.ucr, m2, CL , CR, 'linear' );
F(4,:) = interp2(fr.ucl, fr.ucr, m2', CL , CR, 'linear' );

NL_ph=bplot.MECH_SYMETRICAL(fit_mech.posterior.w, F);
mean_ph = mean(NL_ph,1);
q_ph = quantile(NL_ph,[0.025 0.975],1);

figure('color','w');
ha = tight_subplot(1,3,0.01,0.01,0.1);
for r = 1:3
    hold( ha( r ), 'on');
    plot(ha( r ), CR-CL, mean_ph(1,:,r) , 'k-');
    fx=fill( ha( r ), [CR-CL fliplr(CR-CL)],[q_ph(1,:,r) fliplr(q_ph(2,:,r))],'k'); fx.FaceAlpha=0.2; fx.EdgeAlpha=0;
    
    plot( ha( r ), cdU, prob(1,:,r), 'k.','markersize',10);
end
set(ha,'ylim',[0 1],'dataaspectratio',[1 1 1]);
set(ha(1),'XTickLabelMode','auto','YTickLabelMode','auto')
xlabel(ha(1),'Contrast'); ylabel(ha(1),'Probability of choice');

%%%%%%%%%%%%%%%%%%%%%%%% PLOT WEIGHTS
areaCols = [ 109 183 229;
                77 131 158;
                160 206 87;
             119 147 62]/255;
         
%Plot posterior of Weights, copying over symmetrical values
WL = fit_mech.posterior.w(:,2+[1:4]);
WR = fit_mech.posterior.w(:,2+[2 1 4 3]);
figure('color','w');
axes; hold on;
for region = 1:4
    mu = mean([WL(:,region), WR(:,region)]);
    Sigma = cov([WL(:,region), WR(:,region)]);
    x1 = (mu(1) - 100*max(diag(Sigma))):0.001:(mu(1) + 100*max(diag(Sigma)));
    x2 = (mu(2) - 100*max(diag(Sigma))):0.001:(mu(2) + 100*max(diag(Sigma)));
    [X1,X2] = meshgrid(x1,x2);
    F = mvnpdf([X1(:) X2(:)],mu,Sigma);
    F = reshape(F,length(x2),length(x1));
    [c1,c2]=contour(x1,x2,F,8);
    c2.LineColor = areaCols(region,:);
    c2.LineWidth=1;
end
set(gca,'dataaspectratio',[1 1 1]);
line([0 0],[-1 1]*0.5); line([-1 1]*0.5,[0 0]);
hh=ezplot('y=x'); hh.LineStyle='--';
set(gca,'xlim',[-1 1]*0.25,'ylim',[-1 1]*0.25);
xlabel('W_L'); ylabel('W_R'); title('Posterior dist of weights');

%% Use mechanistic model to predict trial-by-trial choices
%Look within a single contrast condition. Does spontaneous variation in ZL
%and ZR correlate with choice?
fit_mech.data.pedestal = min(fit_mech.data.contrastLeft, fit_mech.data.contrastRight);
fit_mech.data.cDiff = fit_mech.data.contrastRight - fit_mech.data.contrastLeft;

stim_cond = [0.5 0;
    0 0.5];



visBins = 20;
% 
% %For each possible contrast conditions
% CLs = unique(fit_mech.data.contrastLeft);
% CRs = unique(fit_mech.data.contrastRight);
% pCorr = nan(length(CLs),length(CRs));
% for cl = 1:length(CLs)
%     for cr = 1:length(CRs)
%         idx = fit_mech.data.contrastLeft==CLs(cl) & fit_mech.data.contrastRight==CRs(cr);
%         if sum(idx)>0
%         choices = fit_mech.data.choice(idx);
%         firing = fit_mech.data.firing_rate(idx,:);
%         ph=bplot.MECH_SYMETRICAL(fit_mech.posterior.w, firing');
%         ZL = log(ph(:,:,1)./ph(:,:,3));
%         mean_ZL = mean(ZL,1);
%         ZR = log(ph(:,:,2)./ph(:,:,3));
%         mean_ZR = mean(ZR,1);
%         B = mnrfit([mean_ZL' mean_ZR'], choices);
%         phat = mnrval(B, [mean_ZL' mean_ZR']);
%         [~,classified_choice]=max(phat,[],2);
%         pCorr(cl,cr) = mean(choices == classified_choice);
%         end
%     end
% end
figure('color','w');
subplot(1,2,1); hold on;

%Medium contrast on left
idx = fit_mech.data.contrastLeft==0.5 & fit_mech.data.contrastRight==0 & fit_mech.data.choice~=2;
choices = fit_mech.data.choice(idx);
firing = fit_mech.data.firing_rate(idx,:);
ph=bplot.MECH_SYMETRICAL(fit_mech.posterior.w, firing');
ZL = log(ph(:,:,1)./ph(:,:,3));
mean_ZL = mean(ZL,1);
ZR = log(ph(:,:,2)./ph(:,:,3));
mean_ZR = mean(ZR,1);
hx=histogram(mean_ZL(choices==1),visBins); hx.EdgeAlpha=0; 
plot( mean(mean_ZL(choices==1)), 0, 'k.', 'markersize',20);
hx=histogram(mean_ZL(choices==3),visBins);hx.EdgeAlpha=0; 
plot( mean(mean_ZL(choices==3)), 0, 'k.', 'markersize',20);
mdl = fitglm(mean_ZL,choices==1,'Distribution','binomial','Link','logit');
scores = mdl.Fitted.Probability;
[X,Y,T,AUC] = perfcurve(choices==1,scores,1);
title(AUC);

% B = mnrfit([mean_ZL' mean_ZR'], choices);
% phat = mnrval(B, [mean_ZL' mean_ZR']);
% [~,classified_choice]=max(phat,[],2);

%Medium contrast on right
subplot(1,2,2); hold on;

idx = fit_mech.data.contrastLeft==0 & fit_mech.data.contrastRight==0.5 & fit_mech.data.choice~=1;
choices = fit_mech.data.choice(idx);
firing = fit_mech.data.firing_rate(idx,:);
ph=bplot.MECH_SYMETRICAL(fit_mech.posterior.w, firing');
ZR = log(ph(:,:,2)./ph(:,:,3));
mean_ZR = mean(ZR,1);
hx=histogram(mean_ZR(choices==2),visBins); hx.EdgeAlpha=0; 
plot( mean(mean_ZR(choices==2)), 0, 'k.', 'markersize',20);
hx=histogram(mean_ZR(choices==3),visBins);hx.EdgeAlpha=0; 
plot( mean(mean_ZR(choices==3)), 0, 'k.', 'markersize',20);
mdl = fitglm(mean_ZR,choices==2,'Distribution','binomial','Link','logit');
scores = mdl.Fitted.Probability;
[X,Y,T,AUC] = perfcurve(choices==2,scores,1);
title(AUC);


% all conditions CP
contrastLeft = fit_mech.data.contrastLeft;
contrastRight = fit_mech.data.contrastRight;
choice = fit_mech.data.choice;
firing = fit_mech.data.firing_rate;
ph=bplot.MECH_SYMETRICAL(fit_mech.posterior.w, firing');
ZR = log(ph(:,:,2)./ph(:,:,3));
mean_ZR = mean(ZR,1)';
ZL = log(ph(:,:,1)./ph(:,:,3));
mean_ZL = mean(ZL,1)';

% group low med high
contrastLeft(contrastLeft==0.2)=0.25; contrastLeft(contrastLeft==0.4)=0.5; contrastLeft(contrastLeft==0.8)=1;  
contrastRight(contrastRight==0.2)=0.25; contrastRight(contrastRight==0.4)=0.5; contrastRight(contrastRight==0.8)=1;  

uCL = unique(contrastLeft);
uCR = unique(contrastRight);

figure;
ha = tight_subplot(4,4,0.08,0.1,0.1); for i=1:16; hold(ha(i),'on'); end;
for L = 1:4
    for R = 1:4
        idx = contrastLeft==uCL(L) & contrastRight==uCR(R) & fit_mech.data.choice<3;
        ph=bplot.MECH_SYMETRICAL(fit_mech.posterior.w, firing(idx,:)');
        ZR = log(ph(:,:,2)./ph(:,:,3));
        mean_ZR = mean(ZR,1);
        ZL = log(ph(:,:,1)./ph(:,:,3));
        mean_ZL = mean(ZL,1);
        dZ = mean_ZL-mean_ZR;
        
        hx = ha(4*(L-1) + R);
        histogram(hx, dZ(ch(idx)==1),visBins,'EdgeAlpha',0); 
        plot(hx, mean(dZ(ch(idx)==1)), 0, 'b.', 'markersize',20);
        histogram(hx, dZ(ch(idx)==2),visBins,'EdgeAlpha',0);
        plot(hx, mean(dZ(ch(idx)==2)), 0, 'r.', 'markersize',20);
        mdl = fitglm(dZ,ch(idx)==2,'Distribution','binomial','Link','logit');
        scores = mdl.Fitted.Probability;
        [X,Y,T,AUC] = perfcurve(ch(idx)==2,scores,1);
        pval=ranksum(dZ(ch(idx)==1), dZ(ch(idx)==2));
        title(hx,{100*AUC,pval});
        xlabel(hx,sprintf('CL=%0.2f CR=%0.2f', uCL(L), uCR(R)))
    end
end

figure('color','w');
haX = tight_subplot(size(stim_cond,1),2,0.1,0.1,0.1);
for c = 1:size(stim_cond,1)
    idx = fit_mech.data.contrastLeft==stim_cond(c,1) & fit_mech.data.contrastRight==stim_cond(c,2);
    choices = fit_mech.data.choice(idx);
    firing = fit_mech.data.firing_rate(idx,:);
    ph=bplot.MECH_SYMETRICAL(fit_mech.posterior.w, firing');
    ZL = log(ph(:,:,1)./ph(:,:,3));
    mean_ZL = mean(ZL,1);
    ZR = log(ph(:,:,2)./ph(:,:,3));
    mean_ZR = mean(ZR,1);
    
    ha = [haX( 2*(c-1) + 1 ) haX(2*(c-1) + 2)];
    hold(ha(1),'on'); hold(ha(2),'on');
    for r=1:3
        bar(ha(1),r,mean(mean_ZL(choices==r)));
        errorbar(ha(1),r, mean(mean_ZL(choices==r)), std(mean_ZL(choices==r))/sqrt(sum(choices==r)) ,'color', 'k');
        
        bar(ha(2),r,mean(mean_ZR(choices==r)));
        errorbar(ha(2),r, mean(mean_ZR(choices==r)), std(mean_ZR(choices==r))/sqrt(sum(choices==r)) ,'color', 'k');
    end
    linkaxes(ha,'xy');
    ylabel(ha(1),'Mechanistic model ZL');
    ylabel(ha(2),'Mechanistic model ZR');
    set(ha,'xtick',1:3,'xticklabel',{'Left','Right','NoGo'});
    title(ha(1),sprintf('CL=%0.1f, CR=%0.1f',stim_cond(c,1),stim_cond(c,2)));
end
set(haX,'yticklabelmode','auto');

%% Cross-validate auROC

%Fit mechanistic model to half the data
c = cvpartition(length(dat.choice), 'kfold',2);
dat1 = getrow(dat, c.training(1));
dat2 = getrow(dat, c.training(2));
fit_mech = bfit.fitModel('Two-Level-Mechanistic-SYMMETRICAL',dat1);

%Use weights fit from 1st half to predict trial by trial choice variability
%on 2nd half
fit_mech.data.pedestal = min(fit_mech.data.contrastLeft, fit_mech.data.contrastRight);
fit_mech.data.cDiff = fit_mech.data.contrastRight - fit_mech.data.contrastLeft;


figure('color','w');
subplot(1,2,1); hold on;
visBins = 30;

%Medium contrast on left
idx = dat2.contrastLeft==0.5 & dat2.contrastRight==0 & dat2.choice~=2;
choices = dat2.choice(idx);
firing = dat2.firing_rate(idx,:);
ph=bplot.MECH_SYMETRICAL(fit_mech.posterior.w, firing');
ZL = log(ph(:,:,1)./ph(:,:,3));
mean_ZL = mean(ZL,1);
histogram(mean_ZL(choices==1),visBins);
histogram(mean_ZL(choices==3),visBins);
mdl = fitglm(mean_ZL,choices==1,'Distribution','binomial','Link','logit');
scores = mdl.Fitted.Probability;
[X,Y,T,AUC] = perfcurve(choices==1,scores,1);
title({AUC, mean( scores>0.5 == (choices==1) )});

%Medium contrast on right
subplot(1,2,2); hold on;

idx = dat2.contrastLeft==0 & dat2.contrastRight==0.5 & dat2.choice~=1;
choices = dat2.choice(idx);
firing = dat2.firing_rate(idx,:);
ph=bplot.MECH_SYMETRICAL(fit_mech.posterior.w, firing');
ZR = log(ph(:,:,2)./ph(:,:,3));
mean_ZR = mean(ZR,1);
histogram(mean_ZR(choices==2),visBins);
histogram(mean_ZR(choices==3),visBins);
mdl = fitglm(mean_ZR,choices==2,'Distribution','binomial','Link','logit');
scores = mdl.Fitted.Probability;
[X,Y,T,AUC] = perfcurve(choices==2,scores,1);
% title(AUC);
title({AUC, mean( scores>0.5 == (choices==2) )});


%% OLD CODE
% %% Fit widefield data to mechanistic model (4 ROIs at once)
% 
% %1) Go through all sessions and compile data 
% clearvars -except sessionList
% 
% tsteps = linspace(0, 0.25, 20);
% avg_window_width = 0.05;
% [names,~,mouseID]=unique(sessionList.mouseName);
% 
% mega_w = nan(2000,8,20);
% for t=1:length(tsteps)
%     time_window = tsteps(t) + [-0.5 0.5]*avg_window_width;
%     dat = struct;
%     
%     for sess = 1:height(sessionList)
%         eRef = sessionList.expRef{sess};
%         fprintf('Session %d %s\n',sess,eRef);
%         
%         %Load meanImage
%         wfSVDfile = [ './preproc/WF_SVD/' eRef '.mat'];
%         svd = load(wfSVDfile,'meanImg_registered');
%         
%         %Load ROI locations
%         roiFile = [ './preproc/WF_ROI/' sessionList.mouseName{sess} '.mat'];
%         roi = load(roiFile); roi.px(1)=[];
%         
%         
%         areaIdx = [find(contains({roi.px.name},'LVISp')),...
%             find(contains({roi.px.name},'RVISp')),...
%             find(contains({roi.px.name},'LMOs')),...
%             find(contains({roi.px.name},'RMOs'))];   
%         
%         %Load epoch-aligned data
%         wfAlignedFile = [ './preproc/WF_aligned/' eRef '.h5'];
%         ROIstim = h5read(wfAlignedFile,'/ROIstim');
%         ROIstimT = h5read(wfAlignedFile,'/ROIstimT');
%         
%         %Get ROI only at L and R VIS and MOs at the times relevant
%         time_idx = find(time_window(1) <= ROIstimT & ROIstimT <= time_window(2));
%         avgF = mean( ROIstim(:, areaIdx, time_idx), 3);
%         
%         %Load behavioural data
%         behavFile = [ './preproc/BEHAV_RT_Corrected/' eRef '.mat'];
%         b = load(behavFile);
%         
%         row=struct;
%         row.firing_rate = avgF;
%         row.choice = b.choice;
%         row.subjectID = ones(size(row.choice))*mouseID(sess);
%         row.sessionID = ones(size(row.choice))*sess;
%         row.contrastLeft = b.contrastLeft;
%         row.contrastRight = b.contrastRight;
%         dat = addstruct(dat, row);
%     end
%     dat.firing_rate = dat.firing_rate*1000;
% 
%     if any(isnan(dat.firing_rate(:)))
%         warning('Some firing rates were NaNs, replacing with zeros');
%         disp(sum(isnan(dat.firing_rate(:))));
%         dat.firing_rate( isnan(dat.firing_rate) ) = 0;
%     end
% 
%       fit_mech = bfit.fitModel('Two-Level-Mechanistic',dat);
%       mega_w(:,:,t) = fit_mech.posterior.w(:,3:end);
% end
% save('./mechanistic_fit/mega_w_4Together.mat','mega_w');
% 
% %% Plot widefield mechanistic model
% mega_w = load('./mechanistic_fit/mega_w_4Together.mat');
% mega_w = load('./mechanistic_fit/mega_w_4Separate.mat');
% mega_w = mega_w.mega_w;
% 
% areaCols = [ 109 183 229;
%                 77 131 158;
%                 160 206 87;
%              119 147 62]/255;
% figure;
% ha = tight_subplot(2,size(mega_w,3)/2,0.01,0.01,0.01);
% for t = 1:length(tsteps)
%     hold(ha(t),'on');
%     for region = 1:4
%         WL = mega_w(:,2*(region-1)+1,t);
%         WR = mega_w(:,2*(region-1)+2,t);
%         mu = mean([WL, WR],1);
%         Sigma = cov([WL, WR]);
%         x1 = (mu(1) - 100*max(diag(Sigma))):0.01:(mu(1) + 100*max(diag(Sigma)));
%         x2 = (mu(2) - 100*max(diag(Sigma))):0.01:(mu(2) + 100*max(diag(Sigma)));
%         [X1,X2] = meshgrid(x1,x2);
%         F = mvnpdf([X1(:) X2(:)],mu,Sigma);
%         F = reshape(F,length(x2),length(x1));
%         [c1,c2]=contour(ha(t),x1,x2,F,8);
%         c2.LineColor = areaCols(region,:);
%         c2.LineWidth=1;
%     end
%     set(ha(t),'dataaspectratio',[1 1 1]);
%     line(ha(t),[0 0],ylim(ha(t))); line(ha(t),xlim(ha(t)),[0 0]);
%     xlabel(ha(t),'w_L'); ylabel(ha(t),'w_R');
%     set(ha(t),'XTickLabelMode','auto','YTickLabelMode','auto');
%     title(ha(t), tsteps(t));
%     
% end
% linkaxes(ha,'xy');
% 
% 
% 
% 
% 
% areaCols = [ 109 183 229;
%                 77 131 158;
%                 160 206 87;
%              119 147 62]/255;
%          
% figure;
% ha = tight_subplot(2,length(tsteps)/2,0.01,0.05,0.05);
% for t = 1:length(tsteps)
%     
%     hold(ha(t),'on');
%     axes(ha(t));
%     for region = 1:4
%         WL = mega_w(:,2*(region-1)+1,t);
%         WR = mega_w(:,2*(region-1)+2,t);
%         
%         plot(WL,WR,'k+');
%         mu = mean([WL, WR],1);
%         Sigma = cov([WL, WR]);
%         
% %         cx=error_ellipse(Sigma,mu,'conf',0.95);
% %         %cx.Parent=ha(t);
% %         x1 = (mu(1) - 100*max(diag(Sigma))):0.1:(mu(1) + 100*max(diag(Sigma)));
% %         x2 = (mu(2) - 100*max(diag(Sigma))):0.1:(mu(2) + 100*max(diag(Sigma)));
% %         [X1,X2] = meshgrid(x1,x2);
% %         F = mvnpdf([X1(:) X2(:)],mu,Sigma);
% %         F = reshape(F,length(x2),length(x1));
% %         [c1,c2]=contour(ha(t),x1,x2,F,8);
% %         cx.Color = areaCols(region,:);
% %         cx.LineWidth=1;
%     end
%     set(ha(t),'dataaspectratio',[1 1 1]);
%     line(ha(t),[0 0],ylim(ha(t))); line(ha(t),xlim(ha(t)),[0 0]);
%     xlabel(ha(t),'w_L'); ylabel(ha(t),'w_R');
%     set(ha(t),'XTickLabelMode','auto','YTickLabelMode','auto');
%     title(ha(t), tsteps(t));
% end
% linkaxes(ha,'xy');
% 
% areaCols = [ 109 183 229;
%                 77 131 158;
%                 160 206 87;
%              119 147 62]/255;
%          
% figure;
% ha = tight_subplot(4,1,0.01,0.05,0.05);
% m = mean(mega_w,1);
% q = quantile(mega_w,[0.025 0.975],1);
% for region = 1:4
%     hx=plot(ha(region), tsteps, squeeze(m(1,region,:)), '-',...
%                     tsteps, squeeze(m(1,4+region,:)), '--');
%     set(hx,'Color',areaCols(region,:));
%     
%     hold(ha(region),'on');
%     
%     fx=fill(ha(region), [tsteps fliplr(tsteps)], [squeeze(q(1,region,:))' fliplr(squeeze(q(2,region,:))') ] ,'k');
%     fx.EdgeAlpha=0; fx.FaceAlpha=0.3; 
%     fx.FaceColor = areaCols(region,:);
%     
%         
%     fx=fill(ha(region), [tsteps fliplr(tsteps)], [squeeze(q(1,4+region,:))' fliplr(squeeze(q(2,4+region,:))') ] ,'k');
%     fx.EdgeAlpha=0; fx.FaceAlpha=0.3; 
%     fx.FaceColor = areaCols(region,:);
% end

%% Concatenate ROI fluoresence for all sessions

ROIstimStack = cell(height(sessionList),1); 
allBehav = struct; D = struct;
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    %Load meanImage
    wfSVDfile = [ './preproc/WF_SVD/' eRef '.mat'];
    svd = load(wfSVDfile,'meanImg_registered');
    
    %Load ROI locations
    roiFile = [ './preproc/WF_ROI/' sessionList.mouseName{sess} '.mat'];
    roi = load(roiFile);
    bregma=roi.px(1);roi.px(1)=[];
    
    %Load epoch-aligned data
    wfAlignedFile = [ './preproc/WF_aligned/' eRef '.h5'];
    ROIstim = h5read(wfAlignedFile,'/ROIstim');
    ROIstimT = h5read(wfAlignedFile,'/ROIstimT');
    ROIstim(ROIstim<0) = 0;
    
    %Get inclTrials for ROI traces
    %Load behavioural data
    behavFile = [ './preproc/BEHAV_RT_Corrected/' eRef '.mat'];
    b = load(behavFile);
    allBehav = addstruct(allBehav, b);
   
    ROIstimStack{sess} = ROIstim;
end

%Stack ROIs
ROIstimStack = cat(1,ROIstimStack{:});

%% Plot ROI activity across all trials pooling subjects and sessions

%Left VIS and MOs
CR = [0 0.25 0.5 1];

figure;
ha = tight_subplot(1,2,0.1,0.1,0.1); for i=1:2; hold(ha(i),'on'); end;
for cr = 1:4
    idx = allBehav.contrastLeft==0 & allBehav.contrastRight==CR(cr);
    
    v = ROIstimStack(idx, 6,:);
    mn = squeeze(nanmean(v,1));
    err = squeeze(nanstd(v,1)/sqrt(size(v,1)));
    plot(ha(1), ROIstimT, mn );
    fill(ha(1), [ROIstimT fliplr(ROIstimT)], [mn'-err' fliplr(mn'+err')], 'k', 'facealpha',0.3);
    
    
    v = ROIstimStack(idx, 8,:);
    mn = squeeze(nanmean(v,1));
    err = squeeze(nanstd(v,1)/sqrt(size(v,1)));
    plot(ha(2), ROIstimT, mn );
    fill(ha(2), [ROIstimT fliplr(ROIstimT)], [mn'-err' fliplr(mn'+err')], 'k', 'facealpha',0.3)
   
end
set(ha,'xlim',[-0.05 0.2],'xticklabelmode','auto','yticklabelmode','auto');

% ROIidx = 6; %Left VIS
% %Zero contrast, different actions
% 
% figure; axes; hold on;
% for r = 1:3
%     idx = allBehav.contrastLeft>0 & allBehav.contrastRight==0 & allBehav.choice==r;
% 
%     dat = permute( ROIstimStack(idx,ROIidx,:), [1 3 2] );
%     
%     m = nanmean(dat,1);
%     se = nanstd(dat,[],1)/sqrt(sum(idx));
%     
%     plot( ROIstimT, m );
%     fill([ROIstimT fliplr(ROIstimT)],[m-se fliplr(m+se)],'k','Facealpha',0.3);
%     
% end
