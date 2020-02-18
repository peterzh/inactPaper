%Dependencies: 
% https://github.com/cortex-lab/widefield
% https://github.com/cortex-lab/wheelAnalysis

%% Load behavioural data and match reaction time between Left and Right choices
clear all; close all;
sessionList = readtable('../data/widefield/sessionList.csv','FileType','text','Delimiter',',');
[mice,~,miceID] = unique(sessionList.mouseName);

E = table;
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Load behavioural data %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,meta] = loadBehaviouralData(eRef);
    blockFile = meta.blockFile;
    tlFile = strrep(blockFile,'_Block.mat','_Timeline.mat');
    load(blockFile);
    load(tlFile);
    tr = block.trial; tr = tr(1:block.numCompletedTrials);
    cond = [tr.condition];
    vcc = [cond.visCueContrast];
    
    dd = table;
    dd.contrastLeft = vcc(1,:)';
    dd.contrastRight = vcc(2,:)';
    dd.choice = [tr.responseMadeID]';
    dd.feedbackType = [tr.feedbackType]';
    dd.repeatNum = [cond.repeatNum]';
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Align behavioural timings to timeline %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    photodiode = Timeline.rawDAQData(:,2);
    tt = Timeline.rawDAQTimestamps;
    pdFlips = schmittTimes(tt,photodiode, [3 5]); %Get the times when the photodiode flips
    [blockToTL,blockToTL_weights] = makeCorrection(pdFlips(2:end-1), block.stimWindowUpdateTimes, false);
    dd.time_stimOn = blockToTL([tr.stimulusCueStartedTime]');
    dd.time_goCue = blockToTL([tr.interactiveStartedTime]');
    dd.time_responseMade = blockToTL([tr.responseMadeTime]');
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Use wheel trace to detect movement pre-goCue %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tt = Timeline.rawDAQTimestamps;
    wheelPos = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder')); %Wheel encoder trace
    wheelPos = wheelPos(tt<max(dd.time_stimOn)+5); tt = tt(tt<max(dd.time_stimOn)+5); %Trim away wheel trace which extended 5 seconds after the last stimulus onset time
    wheelPos = wheel.correctCounterDiscont(wheelPos); %Correct artifactual values in wheel trace.
    [moveOnsets, moveOffsets] = wheel.findWheelMoves3(wheelPos, tt, 1000, []);
    moveType = wheel.classifyWheelMoves(tt, wheelPos, moveOnsets, moveOffsets, dd.time_stimOn, dd.time_responseMade, dd.choice);
    
    dd.firstMovement = nan(size(dd.choice)); %Choice L/R during in the open-loop period
    dd.firstMovementTime = nan(size(dd.choice));
    for tr = 1:length(dd.time_stimOn)
        idx = find(dd.time_stimOn(tr) < moveOnsets & moveOnsets < dd.time_responseMade(tr),1);
        if ~isempty(idx)
            assert(length(idx)==1,'Too many movements detected');
            
            dd.firstMovement(tr) = moveType(idx);
            dd.firstMovementTime(tr) = moveOnsets(idx);
        end
    end
    
    %Add wheel trace
    [ee,~] = loadBehaviouralData(eRef);
    dd.wheel_stimOn = ee.wheel_stimulusOn;
    
    %add session and subject information
    dd.sessionID = ones(size(dd.choice,1),1)*sess;
    dd.subjectID = ones(size(dd.choice,1),1)*miceID(sess);
    
    wheel_stimulusOn_timesteps = ee.wheel_stimulusOn_timesteps(1,:);
    E = vertcat(E, dd);
end

%Get good trials based on good recent performance
E = E(E.repeatNum==1,:);
E.repeatNum = [];

%Get good trials based on good wheel movements
warning('Experimental exclusion based on wheel movement. See if it helps');
RT = E.firstMovementTime - E.time_stimOn;
% isGoodNoGo = E.choice==3 & isnan(E.firstMovement); %No detected movement for whole trial
% isGoodLeft = E.firstMovement==1 & (RT>0.125 & RT<0.5);
% isGoodRight = E.firstMovement==2 & (RT>0.125 & RT<0.5);
% E = E(isGoodNoGo | isGoodLeft | isGoodRight,:);
isGoodWheelMove = E.choice==3 | isnan(RT) | (RT>0.125 & RT<0.5);
E = E(isGoodWheelMove,:);

%Go through data from each mouse and exclude trials until Left and Right
%choices have similar reaction time histograms
warning('not deterministic! uses random sampling, therefore re-running will exclude different trials!');
D = table;
figure;
for m = 1:length(mice)
    subjData = E(E.subjectID==m,:);
    keepIdx = ones(size(subjData.choice));
    
    RT = subjData.firstMovementTime - subjData.time_stimOn;
    [N,EDGES,BIN]=histcounts(RT);
    
    subplot(length(mice),2,2*(m-1) + 1);
    histogram(RT(subjData.choice==1),EDGES); hold on;
    histogram(RT(subjData.choice==2),EDGES); 
    ylabel(mice{m}); title('uncorrected');
    
    Lchoice = histcounts(RT(subjData.choice==1), EDGES);
    Rchoice = histcounts(RT(subjData.choice==2), EDGES);
    for bin = 1:length(Lchoice)
        %For each bin, identify the trials to remove
        numberToRemove = abs(Lchoice(bin) - Rchoice(bin));
        removeIdx=[];
        if Lchoice(bin) > Rchoice(bin) %trim left choices
            removeIdx = randsample( find(BIN==bin & subjData.choice==1), numberToRemove);
        elseif Lchoice(bin) < Rchoice(bin) %trim right choices
            removeIdx = randsample( find(BIN==bin & subjData.choice==2), numberToRemove);
        end
        keepIdx(removeIdx) = 0;
    end
    subjData = subjData(logical(keepIdx),:);
    
    subplot(length(mice),2,2*(m-1) + 2);
    histogram(RT(subjData.choice==1),EDGES); hold on;
    histogram(RT(subjData.choice==2),EDGES);  title('corrected');
    
    D = vertcat(D, subjData);
end

save('../data/widefield/behaviouralData.mat','D','wheel_stimulusOn_timesteps');

%% Load widefield imaging data: image alignment, downsampling
clear all; close all;
sessionList = readtable('../data/widefield/sessionList.csv','FileType','text','Delimiter',',');
load('../data/widefield/behaviouralData.mat','D');

%Parameters 
window_stim_aligned = [-0.2 0.8];
window_move_aligned = [-0.3 0.3];

%Grid for downsampling
pixSize = 0.0217; % mm/pix. Thris is for PCO edge 5.5 with 0.6x mag (as kilotrode)
% [x1,y1]=meshgrid(-4:0.2:4, -5:0.2:3.5); x1=x1(:); y1=y1(:);
% coordSet = [x1, y1];
% coordSet_pix = round(coordSet/pixSize);

%Extract activity for each session
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    activityFile = [ '../data/widefield/eventAlignedActivity/' eRef '.mat'];
    if ~exist(activityFile,'file')

        %Load SVD of raw data, computed earlier
        load([ '../data/widefield/preproc_SVD_dFF/' eRef '.mat'],'meanImg','U','V','wfTime');
        
        %%%%%%%
        %%% If first session for this mouse, then save meanImg as the
        %%% reference image, and select bregma to define origin of
        %%% coordinates
        %%%%%%
        thisMouse = sessionList.mouseName{sess};
        thisMouseSessIDs = sessionList.sessionID(contains(sessionList.mouseName,thisMouse));
        meanImg = meanImg/quantile(meanImg(:),0.99);%Normalise meanImage 
        if sessionList.sessionID(sess) == min(thisMouseSessIDs) %If first session of subject
            f=figure; imshow(meanImg); title('Select bregma');
            bregma = round(ginput(1)); 
            
            %Create mask 
            title('Draw boundary');
            mask = drawpolygon;
            % Set activity outside mask to NaN
            maskImg = nan(size(U,1),size(U,2));
            for row = 1:size(U,1)
                for col = 1:size(U,2)
                    if inpolygon(col,row,mask.Position(:,1),mask.Position(:,2))
                        maskImg(row,col) = 1;
                    end
                end
            end
            U = U.*maskImg;
            
            RA = imref2d(size(meanImg),pixSize,pixSize); %Coordinates of image in mm
            RA.XWorldLimits = RA.XWorldLimits - bregma(1)*pixSize;
            RA.YWorldLimits = RA.YWorldLimits - bregma(2)*pixSize;

            save(['../data/widefield/reference_image/' thisMouse],'meanImg','RA','maskImg');
            meanImg_registered = meanImg;
            RA_registered = RA;
            U_registered = U;
            
        else %Otherwise, align to the first session
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%% Register image to first session of the same subject %%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            refImg = load(['../data/widefield/reference_image/' thisMouse], 'meanImg','RA','maskImg');
            [optimizer,metric] = imregconfig('Multimodal');
            tform = imregtform(meanImg, refImg.meanImg, 'similarity' ,optimizer,metric);
            falign=figure('name',eRef);
            ha1=subplot(1,2,1);
            ha2=subplot(1,2,2);
            imshowpair(refImg.meanImg, meanImg,'Parent',ha1);
            meanImg_registered = imwarp(meanImg,tform,'OutputView',imref2d(size(refImg.meanImg)));
            imshowpair(refImg.meanImg, meanImg_registered,'Parent',ha2); drawnow;
            RA_registered = refImg.RA;
            
            %Overwrite the U components with the registered version
            U_registered = imwarp(U,tform,'OutputView',imref2d(size(refImg.meanImg)));
            
            %Mask
            U_registered = U_registered.*refImg.maskImg;
        end        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Align & crop image to the first mouse %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~contains(thisMouse,'Chain')
            %Align the meanImg for each mouse to the meanImg of the
            %first mouse (Hench)
            chain = load(['../data/widefield/reference_image/chain.mat']);
            
            %Show alongside
            figure;imshowpair(chain.meanImg,chain.RA,meanImg_registered,RA_registered);

            %Get all world coordinates of Chain's image
            
            meanImg_cropped = nan(size(chain.meanImg,1),size(chain.meanImg,2));
            for row = 1:size(chain.meanImg,1)
                for col = 1:size(chain.meanImg,2)
                    
                    [chain_world_x,chain_world_y]=chain.RA.intrinsicToWorld(row,col);
                    [this_row,this_col]=RA_registered.worldToIntrinsic(chain_world_x,chain_world_y);
                    this_row = round(this_row); this_col=round(this_col);
                    
                    try
                        meanImg_cropped(row,col) = meanImg_registered(this_row,this_col);
                    catch
                    end
                end
            end
            %Get parts of this session's image which correspond to Chain's
            %image world coordinates
            
            error('cropping does not work properly. fix');
            meanImg_cropped=imcrop(RA_registered.XWorldLimits,RA_registered.YWorldLimits,meanImg_registered,...
                [chain.RA.XWorldLimits(1) chain.RA.YWorldLimits(1) chain.RA.ImageExtentInWorldX chain.RA.ImageExtentInWorldY]);
            
            %Crop U component
            U_cropped = nan(size(meanImg_cropped,1),size(meanImg_cropped,2), size(U,3));
            for i = 1:size(U_cropped,3)
                U_cropped(:,:,i) = ...
                    imcrop(RA.XWorldLimits,RA.YWorldLimits,U_registered(:,:,i),...
                    [chain.RA.XWorldLimits(1) chain.RA.YWorldLimits(1) chain.RA.ImageExtentInWorldX chain.RA.ImageExtentInWorldY]);
            end

        else
            meanImg_cropped = meanImg_registered;
            U_cropped = U_registered;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Downsample in space %%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        target_size = [50, 50]; %number of [row,col] in downsampled image
        U_downsampled = nan(target_size(1),target_size(2),size(U,3));
        for i = 1:size(U,3)
            U_downsampled(:,:,i) = imresize(U_cropped(:,:,i), target_size,'Antialiasing',false);
        end
        meanImg_downsampled = imresize(meanImg_cropped, target_size,'Antialiasing',false);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Filter with derivative to approximate spikes %%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        d = designfilt('differentiatorfir','FilterOrder',50, ...
            'PassbandFrequency',8.5,'StopbandFrequency',10, ...
            'SampleRate',1/mean(diff(wfTime)));
        dV = filter(d,V')'; %*Fs;
        delay = mean(grpdelay(d));
        wfTime_dt = wfTime(1:end-delay);
        dV(:,1:delay) = [];
        wfTime_dt(1:delay) = [];
        dV(:,1:delay) = [];
        V = dV;
        wfTime = wfTime_dt;
%         % fvtool(d,'MagnitudeDisplay','zero-phase','Fs',1/mean(diff(wfTime)))

%         %Filter with Andy's deconvolution
%         V=AP_deconv_wf(V);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Get stimulus-aligned activity%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        b = D(D.sessionID==sess,:);
        stim_resp_set = [b.contrastLeft b.contrastRight b.choice];
        [stim_resp,~,stim_resp_id] = unique(stim_resp_set,'rows');
        numTrials = length(stim_resp_id);
        
        assert(all(diff(b.time_stimOn)>0),'Timestamps not monotonically increasing');
        [~, stimT, periStimV, ~] = ...
            eventLockedAvgSVD(U_downsampled, V, wfTime,...
            b.time_stimOn, stim_resp_id, window_stim_aligned);
        
        stimAlignedActivity = nan(size(meanImg_downsampled,1), size(meanImg_downsampled,2), length(stimT), numTrials);
        psV = permute(periStimV, [2 3 1]);
        for n = 1:numTrials
            stimAlignedActivity(:,:,:,n) = svdFrameReconstruct(U_downsampled, psV(:,:,n));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Decode task variables from activity in each pixel %%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        labels = {'Go/NoGo decoding','Left Stim decoding','Choice decoding'};
        stimAlignedDecoding = nan(size(meanImg_downsampled,1), size(meanImg_downsampled,2), length(stimT), length(labels));
        decoded_variable = {b.choice<3,
                            b.contrastLeft>0,
                            b.choice(b.choice<3)==1};
        splitting_condition = {[b.contrastLeft b.contrastRight],...
                               [b.choice b.contrastRight],...
                               [b.contrastLeft(b.choice<3,:) b.contrastRight(b.choice<3,:)]};
        for p = 1:3
            fprintf('\tComputing %s\n',labels{p});
            [conds,~,trial_condition]=unique(splitting_condition{p},'rows');
            dec = decoded_variable{p};
            
            %ensure a large enough number of trials to compute the statistic
            q = crosstab(trial_condition, dec);
            nComp = sum(q(:,1).*q(:,2));
            if nComp>10
                %create shuffle labels
                shufLabels = cell(1,max(trial_condition));
                for c = 1:max(trial_condition)
                    idx = trial_condition==c;
                    chA = dec & idx;
                    nA = sum(chA);
                    shufLabels{c} = (1:nA)';
                end
                
                for row = 1:size(stimAlignedActivity,1)
                    for col = 1:size(stimAlignedActivity,2)
                        for t = 1:length(stimT)
                            switch(labels{p})
                                case 'Choice decoding'
                                    activity = permute(stimAlignedActivity(row,col,t,b.choice<3),[4 1 2 3]);
                                otherwise
                                    activity = permute(stimAlignedActivity(row,col,t,:),[4 1 2 3]);
                            end
                            cpSummary = choiceProbShuf(activity, dec, trial_condition, shufLabels);
                            stimAlignedDecoding(row,col,t,p) =  cpSummary(1);
                        end
                    end
                end
            else
                warning('%s has too few trials to compute %s',eRef,labels{p});
            end
        end
        
        save(activityFile,'stimT','stimAlignedActivity','stimAlignedDecoding');
    end
end

%% Define ROIs for each subject

%% Fit neurometric model to VIS and MOs ROIs

