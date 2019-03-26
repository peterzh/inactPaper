%% Get sparse_unilateral data
clear all;
expRefs = readtable('./sessionList_sparse_unilateral.csv','FileType','text','Delimiter',',');
expRefs = expRefs.expRef;
mouseName = cellfun(@(e)e{3},cellfun(@(s) strsplit(s,'_'), expRefs,'uni',0),'uni',0);
[mice,~,subjID] = unique(mouseName);

%Load standard coordset
load('26CoordSet.mat','coordSet');
coordSet = [coordSet; -coordSet(:,1), coordSet(:,2)];

D = struct;
ii=1;
for session = 1:length(expRefs)
    [dd,meta] = loadData(expRefs{session});
%     
%     p = load(strrep(meta.blockFile, 'Block', 'parameters'));
%     fprintf('Inter-trial delay: %0.1f\n', p.parameters.interTrialDelay);
%     fprintf('Pre-stim quiesence: %0.1f - %0.1f \n', p.parameters.preStimQuiescentPeriod(1),  p.parameters.preStimQuiescentPeriod(2));
%     fprintf('Post-stim delay to go cue: %0.1f \n', p.parameters.cueInteractiveDelay);
%         fprintf('Negative feedback sound duration: %0.1f  \n', p.parameters.negFeedbackSoundDuration);
% 
%     fprintf('\n');

    dd = structfun(@(x)(x(6:(end-14),:)),dd,'uni',0); %trim first 5 trials and last 15
    dd.sessionID = ones(length(dd.response),1)*ii;
    dd.subjectID = ones(length(dd.response),1)*subjID(session);
    
    %Add extra inactivation labels
    dd.laserRegion(dd.laserCoord(:,1)==0.5 & dd.laserCoord(:,2)==-2) = 'RightRSP';
    dd.laserRegion(dd.laserCoord(:,1)==-0.5 & dd.laserCoord(:,2)==-2) = 'LeftRSP';
    dd.laserRegion(dd.laserCoord(:,1)==-3.5 & dd.laserCoord(:,2)==-0.5) = 'LeftS1';
    dd.laserRegion(dd.laserCoord(:,1)==3.5 & dd.laserCoord(:,2)==-0.5) = 'RightS1';
    dd.laserRegion(dd.laserCoord(:,1)==1 & dd.laserCoord(:,2)==-0.5) = 'RightM1';
    dd.laserRegion(dd.laserCoord(:,1)==-1 & dd.laserCoord(:,2)==-0.5) = 'LeftM1';
    
    %Remove trials with inactivation on the non-standard
    %coordinate set
    keep = ones(size(dd.response));
    for tr = 1:length(dd.response)
        if dd.laserType(tr)>0
            if any(sum(dd.laserCoord(tr,:) == coordSet,2)==2)
                keep(tr)=1;
            else
                keep(tr)=0;
            end
        end
    end
    dd = getrow(dd, find(keep));
    
    %If there aren't any laser trials
    if ~any(dd.laserCoord(~isnan(dd.laserCoord(:,1)),1) ~= 0)
        %                         keyboard;
    end
    
    if any(dd.laserType>0)
        D = addstruct(D,dd);
        ii=ii+1;
    end
end
D = getrow(D,D.repeatNum==1);

D.laserIdx = zeros(size(D.response));
%Define laser Idx
for i = 1:size(coordSet,1)
    id = D.laserCoord(:,1) == coordSet(i,1) & D.laserCoord(:,2) == coordSet(i,2);
    D.laserIdx(id) = i;
end

%Save the data somewhere
save('.\data\sparse_unilateral.mat','D','mice');

% fit = load("C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\sparse_unilateral.mat");

%% Plot map on CL=CR condition
D = getrow(D, D.stimulus(:,1) == D.stimulus(:,2) & D.stimulus(:,1)>0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Compute empirical change in choices across coordinates  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the session-averaged change in choices from inactiavtion
[counts,~,~,labels] = crosstab(D.response,...
    D.laserIdx,...
    D.sessionID);
prob = counts./sum(counts,1);%Convert to probability over choices
%Compute delta from the non-laser condition
deltaProb = prob(:,2:end,:) - prob(:,1,:);
deltaProb = nanmean(deltaProb,3); %Average over sessions
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Plot map of choice effects %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('ctxOutlines.mat','coords');
figure('color','w');
for r = 1:3
    subplot(1,3,r); hold on;
    for q = 1:numel(coords) % coords is from ctxOutlines.mat
        cx = coords(q).x/100 - 5.7;
        cy = coords(q).y/100 - 5.4;
        plot(cx,-cy, 'LineWidth', 0.5, 'Color', [1 1 1]*0.5, 'Tag', 'outline');
    end
    
    %plot NS points as small dots
    h=scatter(coordSet(:,1),coordSet(:,2),100,'k','o','filled'); axis equal;  drawnow;
    h.MarkerEdgeColor=[1 1 1]*0.75;
    h.CData = deltaProb(r,:);
    
    caxis([-1 1]*0.8);
    set(gca,'xtick','','ytick','','xcolor','w','ycolor','w','xlim',[-1 1]*5.2,'ylim',[-5 4]);
end
cmap = [ linspace(0,1,100)' linspace(0,1,100)' ones(100,1);
    ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
colormap(cmap);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Plot map of contraversive and ipsiversive effects %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Contraversive: L hemi, change in R choices. R hemi, change in L choices
hemi = sign(coordSet(:,1)); %-1 Left, +1 Right
dContra = mean([deltaProb(2, hemi==-1); deltaProb(1, hemi==+1)],1);
dIpsi = mean([deltaProb(1, hemi==-1); deltaProb(2, hemi==+1)],1);

plotVals = [dContra dContra; dIpsi dIpsi];
figure('color','w');
for r = 1:2
    subplot(1,2,r); hold on;
    for q = 1:numel(coords) % coords is from ctxOutlines.mat
        cx = coords(q).x/100 - 5.7;
        cy = coords(q).y/100 - 5.4;
        plot(cx,-cy, 'LineWidth', 0.5, 'Color', [1 1 1]*0.5, 'Tag', 'outline');
    end
    
    %plot NS points as small dots
    h=scatter(coordSet(:,1),coordSet(:,2),100,'k','o','filled'); axis equal;  drawnow;
    h.MarkerEdgeColor=[1 1 1]*0.75;
    h.CData = plotVals(r,:)';
    
    caxis([-1 1]*0.8);
    set(gca,'xtick','','ytick','','xcolor','w','ycolor','w','xlim',[-1 1]*5.2,'ylim',[-5 4]);
end
cmap = [ linspace(0,1,100)' linspace(0,1,100)' ones(100,1);
    ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
colormap(cmap);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Compute sampling distribution of null: no change in choices %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numIter = 10000;
shufDeltaProb = nan(3,52,numIter);
f = waitbar(0,'Please wait...');
for i = 1:numIter
    waitbar(i/numIter,f);
    
    idx = 1:length(D.laserIdx);
    for sess = 1:max(D.sessionID)
        sessIdx = idx(D.sessionID==sess);
        sessIdx = sessIdx(randperm(length(sessIdx)));
        idx(D.sessionID==sess) = sessIdx;
    end
    counts = crosstab(D.response,...
        D.laserIdx( idx ),...
        D.sessionID);
    prob = counts./sum(counts,1);
    dP = prob(:,2:end,:) - prob(:,1,:);
    shufDeltaProb(:,:,i) = nanmean(dP,3);
end
f.delete;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Plot map of choice effects, showing signifiance from shuffle test %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('ctxOutlines.mat','coords');
q1 = quantile(shufDeltaProb,[0.025 0.975],3);
q2 = quantile(shufDeltaProb,[0.005 0.995],3);
q3 = quantile(shufDeltaProb,[0.0005 0.9995],3);
figure('color','w');
for r = 1:3
    subplot(1,3,r); hold on;
    for q = 1:numel(coords) % coords is from ctxOutlines.mat
        cx = coords(q).x/100 - 5.7;
        cy = coords(q).y/100 - 5.4;
        plot(cx,-cy, 'LineWidth', 0.5, 'Color', [1 1 1]*0.5, 'Tag', 'outline');
    end
    
    %plot NS points as small dots
    idx = q1(r,:,1) < deltaProb(r,:) & deltaProb(r,:) < q1(r,:,2);
    h=scatter(coordSet(idx,1),coordSet(idx,2),20,'k','o','filled'); axis equal;  drawnow;
    h.MarkerEdgeColor=[1 1 1]*0.75;
    h.CData = deltaProb(r,idx);
    
    %Plot *
    idx = deltaProb(r,:) < q1(r,:,1) |  q1(r,:,2) < deltaProb(r,:);
    h=scatter(coordSet(idx,1),coordSet(idx,2),80,'k','o','filled'); axis equal;  drawnow;
    h.MarkerEdgeColor=[1 1 1]*0.75;
    h.CData = deltaProb(r,idx);
    
    %Plot **
    idx = deltaProb(r,:) < q2(r,:,1) |  q2(r,:,2) < deltaProb(r,:);
    h=scatter(coordSet(idx,1),coordSet(idx,2),150,'k','o','filled'); axis equal;  drawnow;
    h.MarkerEdgeColor=[1 1 1]*0.75;
    h.CData = deltaProb(r,idx);
    
    %Plot ***
    idx = deltaProb(r,:) < q3(r,:,1) |  q3(r,:,2) < deltaProb(r,:);
    h=scatter(coordSet(idx,1),coordSet(idx,2),300,'k','o','filled'); axis equal;  drawnow;
    h.MarkerEdgeColor=[1 1 1]*0.75;
    h.CData = deltaProb(r,idx);
    
    caxis([-1 1]*0.8);
    set(gca,'xtick','','ytick','','xcolor','w','ycolor','w','xlim',[-1 1]*5.2,'ylim',[-5 4]);
end
cmap = [ linspace(0,1,100)' linspace(0,1,100)' ones(100,1);
    ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
colormap(cmap);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Group coordinates into regions, and perform shuffle test again %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
perturbationRegion = {'LeftVIS','RightVIS','LeftM2','RightM2','LeftS1','RightS1'};
D1 = getrow(D, any(D.laserRegion == perturbationRegion,2) | D.laserType==0);
D1.perturbation = zeros(size(D1.response));
for p = 1:length(perturbationRegion)
    D1.perturbation(D1.laserRegion == perturbationRegion{p}) = p;
end

shufdContra = nan(numIter,3);
shufdIpsi = nan(numIter,3);
shufdContraRT = nan(numIter,3);
shufdIpsiRT = nan(numIter,3);
f = waitbar(0,'Please wait...');
dRT = nan(2,length(perturbationRegion),max(D.sessionID));
for i = 1:numIter
    waitbar(i/numIter,f);
    idx = 1:length(D1.laserIdx);
    for sess = 1:max(D1.sessionID)
        sessIdx = idx(D1.sessionID==sess);
        sessIdx = sessIdx(randperm(length(sessIdx)));
        idx(D1.sessionID==sess) = sessIdx;
    end
    
    %Compute dProb
    counts = crosstab(D1.response,...
        D1.perturbation( idx ),...
        D1.sessionID);
    prob = counts./sum(counts,1);
    dProb = prob(:,2:end,:) - prob(:,1,:);
    dProb = nanmean(dProb,3);
    
    shufdContra(i,:) = mean([dProb(2,1:2:end); dProb(1,2:2:end)],1);
    shufdIpsi(i,:) = mean([dProb(1,1:2:end); dProb(2,2:2:end)],1);
    
    %compute dRT
    for r = 1:2
        for laser = 1:length(perturbationRegion)
            for sess = 1:max(D1.sessionID)
                dRT(r,laser,sess) = median( D1.RT(D1.response==r & D1.perturbation(idx)==laser & D1.sessionID==sess) ) - ...
                    median( D1.RT(D1.response==r & D1.perturbation(idx)==0 & D1.sessionID==sess) );
            end
        end
    end
    dRT = nanmean(dRT,3);
    
    shufdContraRT(i,:) = mean([dRT(2,1:2:end); dRT(1,2:2:end)],1);
    shufdIpsiRT(i,:) = mean([dRT(1,1:2:end); dRT(2,2:2:end)],1);
end
f.delete;

counts = crosstab(D1.response,...
                D1.perturbation,...
                D1.sessionID);
prob = counts./sum(counts,1);
dProb = prob(:,2:end,:) - prob(:,1,:);
dProb = nanmean(dProb,3);
dContra = mean([dProb(2,1:2:end); dProb(1,2:2:end)],1); %average across hemispheres
dIpsi = mean([dProb(1,1:2:end); dProb(2,2:2:end)],1); %average across hemispheres

dRT = nan(2,length(perturbationRegion),max(D.sessionID));
for r = 1:2
    for laser = 1:length(perturbationRegion)
        for sess = 1:max(D1.sessionID)
            dRT(r,laser,sess) = median( D1.RT(D1.response==r & D1.perturbation==laser & D1.sessionID==sess) ) - ...
                median( D1.RT(D1.response==r & D1.perturbation==0 & D1.sessionID==sess) );
        end
    end
end
dRT = nanmean(dRT,3);
dContraRT = mean([dRT(2,1:2:end); dRT(1,2:2:end)],1);
dIpsiRT = mean([dRT(1,1:2:end); dRT(2,2:2:end)],1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Print results of grouped-by-region shuffle test %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
regions = {'VIS','MOs','SSp'};
for region = 1:length(regions)
    
    %CONTRA choice
    null = shufdContra(:,region);
    P = min( mean(dContra(region) >= null), mean(dContra(region) <= null) );
    fprintf('Inactivate %s affects contra choices: effect=%0.5g , p=%0.5g\n',regions{region},dContra(region),P);
    
    %IPSI choice
    null = shufdIpsi(:,region);
    P = min( mean(dIpsi(region) >= null), mean(dIpsi(region) <= null) );
    fprintf('Inactivate %s affects ipsi choices: effect=%0.5g , p=%0.5g\n',regions{region},dIpsi(region),P);
    
    %CONTRA RT
    null = shufdContraRT(:,region);
    P = min( mean(dContraRT(region) >= null), mean(dContraRT(region) <= null) );
    fprintf('Inactivate %s affects contra rt: effect=%0.5g , p=%0.5g\n',regions{region},dContraRT(region),P);
    
    %IPSI RT
    null = shufdIpsiRT(:,region);
    P = min( mean(dIpsiRT(region) >= null), mean(dIpsiRT(region) <= null) );
    fprintf('Inactivate %s affects ipsi rt: effect=%0.5g , p=%0.5g\n',regions{region},dIpsiRT(region),P);
    
end

%% Plot wheel 

%smooth wheel position
dt = D.wheel_stimulusOn_timesteps(1,2)-D.wheel_stimulusOn_timesteps(1,1);
smooth_t = 30/1000; %10msec
D.wheel_stimulusOn_smoothed = smoothdata(D.wheel_stimulusOn,2,'movmean', 10);

%Compute wheel angular velocity
D.wheel_vel = diff(D.wheel_stimulusOn_smoothed')'./dt;
D.wheel_vel_t = D.wheel_stimulusOn_timesteps(:,1:end-1);

figure;
ha = tight_subplot(1, max(D.subjectID), 0.01,0.01,0.01);
for subj = 1:max(D.subjectID)
    %non laser case,
    idx = D.laserType==0 & D.subjectID==subj;
    
    hold(ha(subj),'on');
    plot(ha(subj),D.wheel_vel_t(idx & D.response==1,:)',D.wheel_vel(idx & D.response==1,:)','Color',[0 0 1 0.1])
    plot(ha(subj),D.wheel_vel_t(idx & D.response==2,:)',D.wheel_vel(idx & D.response==2,:)','Color',[1 0 0 0.1])
    plot(ha(subj),D.wheel_vel_t(idx & D.response==1,:)',mean( D.wheel_vel(idx & D.response==1,:),1)','Color',[0 0 1 1],'linewidth',4)
    plot(ha(subj),D.wheel_vel_t(idx & D.response==2,:)',mean( D.wheel_vel(idx & D.response==2,:),1)','Color',[1 0 0 1],'linewidth',4)
end

% Overlay subj average and compute grand avg
figure; subplot(1,2,1); hold on;
subplot(1,2,2); hold on;
subjAvg = nan(max(D.subjectID), size(D.wheel_vel_t,2), 2);
for subj = 1:max(D.subjectID)
    idx = D.laserType==0 & D.subjectID==subj;

    subjAvg(subj,:,1) = mean( D.wheel_vel(idx & D.response==1,:),1);
    subjAvg(subj,:,2) = mean( D.wheel_vel(idx & D.response==2,:),1);
    
    subplot(1,2,1);
    plot(D.wheel_vel_t(1,:),subjAvg(subj,:,1),'Color',[0 0 1],'linewidth',4)
    
    subplot(1,2,2); 
    plot(D.wheel_vel_t(1,:),-subjAvg(subj,:,2),'Color',[1 0 0],'linewidth',4)
end
%add pooled avg
subplot(1,2,1); 
plot(D.wheel_vel_t(1,:),mean(D.wheel_vel(D.laserType==0 & D.response==1,:),1),'Color',[0 0 0],'linewidth',4);
subplot(1,2,2); 
plot(D.wheel_vel_t(1,:),-mean(D.wheel_vel(D.laserType==0 & D.response==2,:),1),'Color',[0 0 0],'linewidth',4);

% Now look at the case with inactivation


figure('color','w');
% regions = {'LeftVIS','RightVIS','LeftM2','RightM2','LeftS1','RightS1'};
regions = {'LeftVIS','RightVIS'};
areaCols = [ 109 183 229;
    77 131 158;
    160 206 87;
    119 147 62;
    50 50 50;
    100 100 100]/255;

for r = 1:2
    hold on;
    plot(D.wheel_vel_t(1,:),mean(D.wheel_vel(D.laserType==0 & D.response==r,:),1),'Color',[0 0 0],'linewidth',4);

    for loc = 1:length(regions)
        plot(D.wheel_vel_t(1,:),mean(D.wheel_vel(D.laserRegion==regions{loc} & D.response==r,:),1),'Color',areaCols(loc,:),'linewidth',4);
    end
end
legend(['OFF' regions])



figure('color','w');
ha = tight_subplot( length(regions), 2, 0.01,0.05,0.05);
timestamps = D.wheel_vel_t(1,:);
for r = 1:2
    
    vel_off = D.wheel_vel(D.laserType==0 & D.response==r,:);
    avg_off = mean(vel_off,1);
    err_off = std(vel_off,[],1);

    for loc = 1:length(regions)
        hold(ha(2*(loc-1) +r),'on');
        
        %Plot off
        
        plot( ha(2*(loc-1) +r), D.wheel_vel_t(1,:), avg_off, 'k-','linewidth',4);
        fill( ha(2*(loc-1) +r), [timestamps fliplr(timestamps)], [avg_off-err_off fliplr(avg_off+err_off)],'k','FaceAlpha',0.2)

        %Plot on
        vel_on = D.wheel_vel(D.laserRegion==regions{loc} & D.response==r,:);

        avg_on = mean(vel_on,1);
        err_on = std(vel_on,[],1);
        plot( ha(2*(loc-1) +r), D.wheel_vel_t(1,:), avg_on, 'r-','linewidth',4);
        fill( ha(2*(loc-1) +r), [timestamps fliplr(timestamps)], [avg_on-err_on fliplr(avg_on+err_on)],'r','FaceAlpha',0.2);
        
        if r==1
            ylabel(ha(2*(loc-1) +r), regions{loc});
        end
    end
    
end
title(ha(1),'Left choice');
title(ha(2),'Right choice');


%% Get pulse experiment data
clear all;

expRefs = readtable('./sessionList_pulse.csv','FileType','text','Delimiter',',');
expRefs = expRefs.expRef;
mouseName = cellfun(@(e)e{3},cellfun(@(s) strsplit(s,'_'), expRefs,'uni',0),'uni',0);
[mice,~,subjID] = unique(mouseName);

D = struct;
ii=1;
for session = 1:length(expRefs)
    [dd,meta] = loadData(expRefs{session});
    
        p = load(strrep(meta.blockFile, 'Block', 'parameters'));
%     fprintf('Pre-stim quiesence: %0.1f \n', p.parameters.preStimQuiescentPeriod);
    fprintf('Post-stim delay to go cue: %0.1f \n', p.parameters.interactiveDelay);
    fprintf('\n');
    
    dd = structfun(@(x)(x(6:(end-14),:)),dd,'uni',0); %trim first 5 trials and last 15
    dd.sessionID = ones(length(dd.response),1)*ii;
    dd.subjectID = ones(length(dd.response),1)*subjID(session);
    
    %Add extra inactivation labels
    dd.laserRegion(dd.laserCoord(:,1)==0.5 & dd.laserCoord(:,2)==-2) = 'RightRSP';
    dd.laserRegion(dd.laserCoord(:,1)==-0.5 & dd.laserCoord(:,2)==-2) = 'LeftRSP';
    dd.laserRegion(dd.laserCoord(:,1)==-3.5 & dd.laserCoord(:,2)==-0.5) = 'LeftS1';
    dd.laserRegion(dd.laserCoord(:,1)==3.5 & dd.laserCoord(:,2)==-0.5) = 'RightS1';
    dd.laserRegion(dd.laserCoord(:,1)==1 & dd.laserCoord(:,2)==-0.5) = 'RightM1';
    dd.laserRegion(dd.laserCoord(:,1)==-1 & dd.laserCoord(:,2)==-0.5) = 'LeftM1';
    
    %add marker for whether the choice is contra or ipsi to the inactivated side is contra or ipsi to stimulus side
    
    dd.ipsiversive = ( contains(string(dd.laserRegion),'Left') & (dd.response==1) ) | ...
               ( contains(string(dd.laserRegion),'Right') & (dd.response==2) );
    
    
    if any(dd.laserType>0)
        D = addstruct(D,dd);
        ii=ii+1;
    end
end
D = getrow(D,D.repeatNum==1);

%define stimulus conditions where one side is higher than the other
D.CL_gt_CR = D.stimulus(:,1) > D.stimulus(:,2);
D.CR_gt_CL = D.stimulus(:,2) > D.stimulus(:,1);

%prepare laser onset times for plotting later
D.laserOnset = D.laserOnset + randn(length(D.laserOnset),1)/100000;
[~,sortIdx]=sort(D.laserOnset);
D = getrow(D, sortIdx);

%Remove trials with pre-stimulus wheel movement
D.keepIdx = zeros(size(D.response));
for sess = 1:max(D.sessionID)
    times = D.wheel_stimulusOn_timesteps(D.sessionID==sess,:);
    times = times(1,:);
    wheel = D.wheel_stimulusOn(D.sessionID==sess,:);
    choice = D.response(D.sessionID==sess);
    
    %Whiten
    wheel_whiten =  wheel/std(wheel(:,end));
    
    %Threshold
    wheel_whiten_window = mean( abs(wheel_whiten(:,-0.150 < times & times < 0.05 & times ~= 0)) , 2);
    D.keepIdx(D.sessionID==sess) = wheel_whiten_window < 0.05;
end
D = getrow(D, D.keepIdx==1);

% %Exclude trials where they move before the laser
% D = getrow(D, D.laserOnset < D.RT | isnan(D.RT));

%Save the data somewhere
% save('.\data\galvo_pulse.mat','D','mice');

%get non-laser RT and performance
NL =  getrow(D, (D.CR_gt_CL | D.CL_gt_CR) & D.laserType==0 );
VIS = getrow(D, (D.laserRegion == 'LeftVIS' & D.CR_gt_CL) | (D.laserRegion == 'RightVIS' & D.CL_gt_CR) );
M2 = getrow(D, (D.laserRegion == 'LeftM2' & D.CR_gt_CL) | (D.laserRegion == 'RightM2' & D.CL_gt_CR) );
M1 = getrow(D, (D.laserRegion == 'LeftM1' & D.CR_gt_CL) | (D.laserRegion == 'RightM1' & D.CL_gt_CR) );

%% Sliding window performance and RT

window_sec = 0.1;
alpha = 0.0001;

% M2.RT(M2.laserOnset>0.1)=0.4;

% figure(301);
figure;
subplot(3,2,2); hold on; line([-1 1]*0.3, [1 1]*mean(NL.feedbackType),'Color',[0 0 0])
subplot(3,2,4); hold on; line([-1 1]*0.3, [1 1]*mean(NL.response==3),'Color',[0 0 0])
subplot(3,2,6); hold on; line([-1 1]*0.3, [1 1]*median(NL.RT(NL.feedbackType==1)),'Color',[0 0 0])

tsteps = linspace(-0.3,0.3,1000);
perf = nan(2,1000);
perf_ci = nan(2,1000,2);
nogo = nan(2,1000);
nogo_ci = nan(2,1000,2);
rt = nan(2,1000);
rt_ci = zeros(2,1000,2);
for t = 1:length(tsteps)
    %VIS
    idx = WithinRanges(VIS.laserOnset, tsteps(t) + [-0.5 0.5]*window_sec);
    [counts,~,pVal] = crosstab( [NL.feedbackType; VIS.feedbackType(idx==1)], [NL.laserType; VIS.laserType(idx==1)]);
    [perf(1,t),perf_ci(1,t,:)] = binofit(counts(2,2),sum(counts(:,2)));
    if pVal<alpha
        subplot(3,2,2);
        plot(tsteps(t), 0.8,'bo');
    end
    
    [counts,~,pVal] = crosstab( [NL.response==3; VIS.response(idx==1)==3], [NL.laserType; VIS.laserType(idx==1)]);
    [nogo(1,t),nogo_ci(1,t,:)] = binofit(counts(2,2),sum(counts(:,2)));
    if pVal<alpha
        subplot(3,2,4);
        plot(tsteps(t), 0,'bo');
    end
    
    pVal = ranksum( NL.RT(NL.feedbackType==1), VIS.RT(idx & VIS.feedbackType==1), 'tail', 'left');
    rt_vals = VIS.RT(idx & VIS.feedbackType==1);
    rt(1,t) = median( rt_vals );
    rt_ci(1,t,:) = median( rt_vals ) + tinv([0.025  0.975],length(rt_vals)-1)*std(rt_vals)/sqrt(length(rt_vals));
    if pVal<alpha
        subplot(3,2,6);
        plot(tsteps(t), 0.31,'bo');
    end
    
    %M2
    idx = WithinRanges(M2.laserOnset, tsteps(t) + [-0.5 0.5]*window_sec);
    [counts,~,pVal] = crosstab( [NL.feedbackType; M2.feedbackType(idx==1)], [NL.laserType; M2.laserType(idx==1)]);
    [perf(2,t),perf_ci(2,t,:)] = binofit(counts(2,2),sum(counts(:,2)));
    if pVal<alpha
        subplot(3,2,2);
        plot(tsteps(t), 0.7,'ro');
    end
    
    [counts,~,pVal] = crosstab( [NL.response==3; M2.response(idx==1)==3], [NL.laserType; M2.laserType(idx==1)]);
    [nogo(2,t),nogo_ci(2,t,:)] = binofit(counts(2,2),sum(counts(:,2)));
    if pVal<alpha
        subplot(3,2,4);
        plot(tsteps(t), 0.05,'ro');
    end
    
    pVal = ranksum( NL.RT(NL.feedbackType==1), M2.RT(idx & M2.feedbackType==1), 'tail', 'left' );
    rt_vals = M2.RT(idx & M2.feedbackType==1);
    rt(2,t) = median( rt_vals );
    rt_ci(2,t,:) = median( rt_vals ) + tinv([0.025  0.975],length(rt_vals)-1)*std(rt_vals)/sqrt(length(rt_vals));
    if pVal<alpha
        subplot(3,2,6);
        plot(tsteps(t), 0.3,'ro');
    end
end

subplot(3,2,2);
plot(tsteps, perf(1,:), 'b'); 
fill([tsteps fliplr(tsteps)],[perf_ci(1,:,1) fliplr(perf_ci(1,:,2))],'k','FaceColor','b', 'EdgeAlpha', 0, 'FaceAlpha', 0.3);
plot(tsteps, perf(2,:), 'r'); 
fill([tsteps fliplr(tsteps)],[perf_ci(2,:,1) fliplr(perf_ci(2,:,2))],'k','FaceColor','r', 'EdgeAlpha', 0, 'FaceAlpha', 0.3);
xlim([-1 1]*0.3); xlabel('Stim onset'); ylabel('Proportion correct');


subplot(3,2,4);
plot(tsteps, nogo(1,:), 'b'); 
fill([tsteps fliplr(tsteps)],[nogo_ci(1,:,1) fliplr(nogo_ci(1,:,2))],'k','FaceColor','b', 'EdgeAlpha', 0, 'FaceAlpha', 0.3);
plot(tsteps, nogo(2,:), 'r'); 
fill([tsteps fliplr(tsteps)],[nogo_ci(2,:,1) fliplr(nogo_ci(2,:,2))],'k','FaceColor','r', 'EdgeAlpha', 0, 'FaceAlpha', 0.3);
xlim([-1 1]*0.3); xlabel('Stim onset'); ylabel('Proportion NoGo');

subplot(3,2,6);
plot(tsteps,rt(1,:),'b');
fill([tsteps fliplr(tsteps)],[rt_ci(1,:,1) fliplr(rt_ci(1,:,2))],'k','FaceColor','b', 'EdgeAlpha', 0, 'FaceAlpha', 0.3);
plot(tsteps,rt(2,:),'r');
fill([tsteps fliplr(tsteps)],[rt_ci(2,:,1) fliplr(rt_ci(2,:,2))],'k','FaceColor','r', 'EdgeAlpha', 0, 'FaceAlpha', 0.3);
xlim([-1 1]*0.3); xlabel('Stim onset'); ylabel('Median RT');
 
%% Plot pulse for other regions
A = {VIS, M2, M1};

figure;
% ha = tight_subplot(3,2,0.1,0.1,0.1);
% for i = 1:length(ha); hold(ha(i),'on'); end;
% 


subplot(3,2,2); hold on; line([-1 1]*0.3, [1 1]*mean(NL.feedbackType),'Color',[0 0 0])
subplot(3,2,4); hold on; line([-1 1]*0.3, [1 1]*median(NL.RT(NL.feedbackType==1)),'Color',[0 0 0])

tsteps = linspace(-0.1,0.3,1000);
perf = nan(1,1000);
perf_ci = nan(1,1000,2);
nogo = nan(1,1000);
nogo_ci = nan(1,1000,2);
rt = nan(1,1000);
rt_ci = zeros(1,1000,2);
for t = 1:length(tsteps)
    %M1
    idx = WithinRanges(M1.laserOnset, tsteps(t) + [-0.5 0.5]*window_sec);
    [counts,~,pVal] = crosstab( [NL.feedbackType; M1.feedbackType(idx==1)], [NL.laserType; M1.laserType(idx==1)]);
    [perf(1,t),perf_ci(1,t,:)] = binofit(counts(2,2),sum(counts(:,2)));
    if pVal<alpha
        subplot(3,2,2);
        plot(tsteps(t), 0.8,'bo');
    end
    
    pVal = ranksum( NL.RT(NL.feedbackType==1), M1.RT(idx & M1.feedbackType==1), 'tail', 'left');
    rt_vals = M1.RT(idx & M1.feedbackType==1);
    rt(1,t) = median( rt_vals );
    rt_ci(1,t,:) = median( rt_vals ) + tinv([0.025  0.975],length(rt_vals)-1)*std(rt_vals)/sqrt(length(rt_vals));
    if pVal<alpha
        subplot(3,2,4);
        plot(tsteps(t), 0.31,'bo');
    end
  
end

subplot(3,2,2);
plot(tsteps, perf(1,:), 'k'); 
fill([tsteps fliplr(tsteps)],[perf_ci(1,:,1) fliplr(perf_ci(1,:,2))],'k','FaceColor','k', 'EdgeAlpha', 0, 'FaceAlpha', 0.3);
xlim([min(tsteps) max(tsteps)]); xlabel('Stim onset'); ylabel('Proportion correct');

subplot(3,2,4);
plot(tsteps,rt(1,:),'k');
fill([tsteps fliplr(tsteps)],[rt_ci(1,:,1) fliplr(rt_ci(1,:,2))],'k','FaceColor','k', 'EdgeAlpha', 0, 'FaceAlpha', 0.3);
xlim([min(tsteps) max(tsteps)]); xlabel('Stim onset'); ylabel('Median RT');
 
%% Fit model to pulse data at bins
perturbationRegion = {'LeftVIS','RightVIS','LeftM2','RightM2'};
D1 = getrow(D, any(D.laserRegion == perturbationRegion,2) | D.laserType==0);
D1.perturbation = zeros(size(D1.response));
for p = 1:length(perturbationRegion)
    D1.perturbation(D1.laserRegion == perturbationRegion{p}) = p;
end

window_sec = 0.2;
tsteps = linspace(-0.150,0.250,10);

deltas = nan(2000,4,length(perturbationRegion),length(tsteps));
for t = 1:length(tsteps)
	idx = WithinRanges(D1.laserOnset, tsteps(t) + [-0.5 0.5]*window_sec);
    D2 = getrow(D1, idx==1);
    
    %Reset session and subject IDs
    [~,~,D2.sessionID]=unique(D2.sessionID);
    [~,~,D2.subjectID]=unique(D2.subjectID);
    
    dat = struct('contrastLeft', D2.stimulus(:,1),...
              'contrastRight', D2.stimulus(:,2),...
              'choice', D2.response,...
              'sessionID', D2.sessionID,...
              'subjectID', D2.subjectID,...
              'perturbation', D2.perturbation);
    fit = bfit.fitModel('Basic_CN-Perturbation',dat);
    deltas(:,:,:,t) = fit.posterior.delta;
end

%Plot deltas over time
pLabel = {'BL','BR','SL','SR'};
figure;
ha = tight_subplot(4,4, 0.05,0.1,0.1);
for p = 1:4
    for roi = 1:4
        dM = squeeze( mean( deltas(:,p,roi,:), 1) );
        dQ = squeeze( quantile(deltas(:,p,roi,:),[0.025 0.975],1) );
        
        hold(ha(4*(p-1) + roi),'on');
        plot(ha(4*(p-1) + roi),  tsteps, dM);
        fill(ha(4*(p-1) + roi),[tsteps fliplr(tsteps)],[dQ(1,:) fliplr(dQ(2,:))],'k','FaceColor','b', 'EdgeAlpha', 0, 'FaceAlpha', 0.3);
        xlim(ha(4*(p-1) + roi),[min(tsteps) max(tsteps)]);
        line(ha(4*(p-1) + roi),get(ha(4*(p-1) + roi),'xlim'),[0 0], 'color','k');
        
        if p == 1
            title(ha(4*(p-1) + roi),perturbationRegion{roi});
        end
        
        if roi == 1
            ylabel(ha(4*(p-1) + roi),pLabel{p});
        end
    end
end
linkaxes(ha(1:8),'y');
linkaxes(ha(9:end),'y');
set(ha(13),'xticklabelmode','auto');



 %% Get multi-power inactivation data
clear all;
expRefs = readtable('./sessionList_unilateral_multiple_powers.csv','FileType','text','Delimiter',',');
expRefs = expRefs.expRef;
mouseName = cellfun(@(e)e{3},cellfun(@(s) strsplit(s,'_'), expRefs,'uni',0),'uni',0);
[mice,~,subjID] = unique(mouseName);

D = struct;
ii=1;
for session = 1:length(expRefs)
    [dd,meta] = loadData(expRefs{session});
           p = load(strrep(meta.blockFile, 'Block', 'parameters'));
%     fprintf('Inter-trial delay: %0.1f\n', p.parameters.interTrialDelay);
%     fprintf('Pre-stim quiesence: %0.1f \n', p.parameters.preStimQuiescentPeriod);
    fprintf('Post-stim delay to go cue: %0.1f \n', p.parameters.interactiveDelay);
    fprintf('\n');
    dd = structfun(@(x)(x(6:(end-14),:)),dd,'uni',0); %trim first 5 trials and last 15
    dd.sessionID = ones(length(dd.response),1)*ii;
    dd.subjectID = ones(length(dd.response),1)*subjID(session);
    
    %Add extra inactivation labels
    dd.laserRegion(dd.laserCoord(:,1)==0.5 & dd.laserCoord(:,2)==-2) = 'RightRSP';
    dd.laserRegion(dd.laserCoord(:,1)==-0.5 & dd.laserCoord(:,2)==-2) = 'LeftRSP';
    dd.laserRegion(dd.laserCoord(:,1)==-3.5 & dd.laserCoord(:,2)==-0.5) = 'LeftS1';
    dd.laserRegion(dd.laserCoord(:,1)==3.5 & dd.laserCoord(:,2)==-0.5) = 'RightS1';
    dd.laserRegion(dd.laserCoord(:,1)==1 & dd.laserCoord(:,2)==-0.5) = 'RightM1';
    dd.laserRegion(dd.laserCoord(:,1)==-1 & dd.laserCoord(:,2)==-0.5) = 'LeftM1';
    
    if any(dd.laserType>0)
        D = addstruct(D,dd);
        ii=ii+1;
    end
end
D = getrow(D,D.repeatNum==1);

coordSet = unique(D.laserCoord(D.laserType==1,:),'rows');
D.laserIdx = zeros(size(D.response));
%Define laser Idx
for i = 1:size(coordSet,1)
    id = D.laserCoord(:,1) == coordSet(i,1) & D.laserCoord(:,2) == coordSet(i,2);
    D.laserIdx(id) = i;
end

%Remove trials with pre-stimulus wheel movement
D.keepIdx = zeros(size(D.response));
for sess = 1:max(D.sessionID)
    times = D.wheel_stimulusOn_timesteps(D.sessionID==sess,:);
    times = times(1,:);
    wheel = D.wheel_stimulusOn(D.sessionID==sess,:);
    choice = D.response(D.sessionID==sess);
    
    %Whiten
    wheel_whiten =  wheel/std(wheel(:,end));
    
    %Threshold
    wheel_whiten_window = mean( abs(wheel_whiten(:,-0.150 < times & times < 0.05 & times ~= 0)) , 2);
    D.keepIdx(D.sessionID==sess) = wheel_whiten_window < 0.05;
end
D = getrow(D, D.keepIdx==1);

save('.\data\galvo_multipower.mat','D','mice');
% %% Multi-power: model-free statistical testing of laser effects
% 
% %%%%%%%%%%% GET DATA %%%%%%%%%%%%%%
% perturbationRegion = {'LeftVIS','RightVIS','LeftM2','RightM2','LeftM1','RightM1'};
% regionLabels = {'VIS','MOs','MOp'};
% areaCols = [ 73 148 208;
%              119 172 66;
%              100 100 100]/255;
% D2 = getrow(D, any(D.laserRegion == perturbationRegion,2) | D.laserType==0);
% D2.perturbation = zeros(size(D2.response));
% for p = 1:length(perturbationRegion)
%     D2.perturbation(D2.laserRegion == perturbationRegion{p}) = p;
% end

% D2 = getrow(D2,D2.stimulus(:,1)>=0.1 | D2.stimulus(:,2)>=0.1);
% 
% %%%%%%%%%%% EFFECT ON FRACTION OF CHOICES %%%%%%%%%%%%%%
% 
% %%%%%%%%% Compute empirical change in probability of moving contraversive
% %%%%%%%%% or ipsiversive for contralateral stimuli only 
[counts,~,~,labels] = crosstab(D2.stimulus(:,1),...
    D2.stimulus(:,2),...
    D2.response,...
    D2.perturbation,...
    D2.laserPower,...
    D2.sessionID);
% prob = counts./sum(counts,1);
% dP = prob(:,2:end,:,:,:,2:end) - prob(:,1,:,:,:,1); %difference from no-laser condition
% dP = nanmean(dP,3); %avg over sessions
% dContraversive = nanmean([dP(1,2:2:end,:,2,1,:);... %change in Left choices on CL, when inactivating RIGHT hemisphere
%     dP(2,1:2:end,:,1,2,:)],1); %change in Right choices on CR, when inactivating LEFT hemisphere
% dIpsiversive = nanmean([dP(2,2:2:end,:,2,1,:);... %change in Right choices on CL, when inactivating RIGHT hemisphere
%     dP(1,1:2:end,:,1,2,:)],1); %change in Left choices on CR, when inactivating LEFT hemisphere
% dContraversive = permute(dContraversive,[2 6 1 3 4 5]); %rows regions, cols powers
% dIpsiversive = permute(dIpsiversive,[2 6 1 3 4 5]);
% 
% %%%%%%%Compute null sampling distribution for no difference
% numIter = 10000;
% dContraversive_null = nan(3,3,numIter);
% dIpsiversive_null = nan(3,3,numIter);
% f = waitbar(0,'Please wait...');
% for i = 1:numIter
%     waitbar(i/numIter,f);
%     idx = 1:length(D2.laserIdx);
%     for sess = 1:max(D2.sessionID)
%         sessIdx = idx(D2.sessionID==sess);
%         sessIdx = sessIdx(randperm(length(sessIdx)));
%         idx(D2.sessionID==sess) = sessIdx;
%     end
%     
%     counts = crosstab(D2.response,...
%         D2.perturbation( idx ),...
%         D2.sessionID,...
%         D2.stimulus(:,1)>0,...
%         D2.stimulus(:,2)>0,...
%         D2.laserPower( idx )); %%%%<<<<< DO I ALSO WANT TO BE SHUFFLING LASER POWER?
%     prob = counts./sum(counts,1);
%     dP = prob(:,2:end,:,:,:,2:end) - prob(:,1,:,:,:,1); %difference from no-laser condition
%     dP = nanmean(dP,3); %avg over sessions
%     dC = nanmean([dP(1,2:2:end,:,2,1,:);... %change in Left choices on CL, when inactivating RIGHT hemisphere
%         dP(2,1:2:end,:,1,2,:)],1); %change in Right choices on CR, when inactivating LEFT hemisphere
%     dI = nanmean([dP(2,2:2:end,:,2,1,:);... %change in Right choices on CL, when inactivating RIGHT hemisphere
%         dP(1,1:2:end,:,1,2,:)],1); %change in Left choices on CR, when inactivating LEFT hemisphere
% %     dContraversive_null(:,1,i) = permute(dC,[2 6 1 3 4 5]); %rows regions, cols powers
% %     dIpsiversive_null(:,1,i) = permute(dI,[2 6 1 3 4 5]);
%     
%     dContraversive_null(:,:,i) = permute(dC,[2 6 1 3 4 5]); %rows regions, cols powers
%     dIpsiversive_null(:,:,i) = permute(dI,[2 6 1 3 4 5]);
% end
% f.delete;
% % 
% % %%%%%%Plot change in probability, overlaying pvalues
% % figure('color','w');
% % for region = 1:3
% %     subplot(1,2,1); hold on;
% %     plot(region + [-1 0 1]*0.2, dContraversive(region,:),'.-','markersize',20);
% %     null = squeeze(dContraversive_null(region,:,:))';
% %     pVal = min(mean(dContraversive(region,:) <= null,1), mean(dContraversive(region,:) >= null,1));
% %     text(region + [-1 0 1]*0.2, dContraversive(region,:),arrayfun(@(n) num2str(n),pVal,'uni',0) );
% %     
% %     subplot(1,2,2); hold on;
% %     plot(region + [-1 0 1]*0.2, dIpsiversive(region,:),'.-','markersize',20);
% %     null = squeeze(dIpsiversive_null(region,:,:))';
% %     pVal = min(mean(dIpsiversive(region,:) <= null,1), mean(dIpsiversive(region,:) >= null,1));
% %     text(region + [-1 0 1]*0.2, dIpsiversive(region,:),arrayfun(@(n) num2str(n),pVal,'uni',0) );
% % end
% % set(get(gcf,'children'),'xtick',1:3,'xticklabel',regionLabels);
% % 
% % 
% 
% %%%%%%Plot absolute probabilities, overlaying pvalues
% counts = crosstab(D2.response,...
%     D2.perturbation,...
%     D2.sessionID,...
%     D2.stimulus(:,1)>0,...
%     D2.stimulus(:,2)>0,...
%     D2.laserPower);
% prob = counts./sum(counts,1);
% prob_ave = nanmean(prob,3);
% pContraNL = nanmean([prob_ave(1,1,1,2,1,1);... % Left choices on CL
%                      prob_ave(2,1,1,1,2,1)],1); %Right choices on CR
% pIpsiNL = nanmean([prob_ave(2,1,1,2,1,1);... % Right choices on CL
%                      prob_ave(1,1,1,1,2,1)],1); %Left choices on CR
% 
% prob_ave = prob_ave(:,2:end,:,:,:,2:end);
% pContraversive = nanmean([prob_ave(1,2:2:end,:,2,1,:);... % Left choices on CL, when inactivating RIGHT hemisphere
%     prob_ave(2,1:2:end,:,1,2,:)],1); %Right choices on CR, when inactivating LEFT hemisphere
% pIpsiversive = nanmean([prob_ave(2,2:2:end,:,2,1,:);... %Right choices on CL, when inactivating RIGHT hemisphere
%     prob_ave(1,1:2:end,:,1,2,:)],1); %Left choices on CR, when inactivating LEFT hemisphere
% pContraversive = permute(pContraversive,[2 6 1 3 4 5]); %rows regions, cols powers
% pIpsiversive = permute(pIpsiversive,[2 6 1 3 4 5]);
% 
% figure('color','w');
% subplot(1,2,1); hold on; line([1 3],[1 1]*pContraNL);
% subplot(1,2,2); hold on; line([1 3],[1 1]*pIpsiNL);
% for region = 1:3
% 	subplot(1,2,1); hold on;
%     h=plot( region + [-1 0 1]*0.2, pContraversive(region,:),'.-','markersize',20)
%     null = squeeze(dContraversive_null(region,:,:))';
%     pVal = min(mean(dContraversive(region,:) <= null,1), mean(dContraversive(region,:) >= null,1));
%     text(region + [-1 0 1]*0.2, pContraversive(region,:),arrayfun(@(n) num2str(n),pVal,'uni',0) );
%     
%     subplot(1,2,2); hold on;
%     plot( region + [-1 0 1]*0.2, pIpsiversive(region,:),'.-','markersize',20)
%     null = squeeze(dIpsiversive_null(region,:,:))';
%     pVal = min(mean(dIpsiversive(region,:) <= null,1), mean(dIpsiversive(region,:) >= null,1));
%     text(region + [-1 0 1]*0.2, pIpsiversive(region,:),arrayfun(@(n) num2str(n),pVal,'uni',0) );
% end
% set(get(gcf,'children'),'ylim',[0 1]);
% set(get(gcf,'children'),'xtick',1:3,'xticklabel',regionLabels);
% 
% %%%%%%%%%%% PAIRED TEST BETWEEN VIS AND M2 on IPSIVERSIVE CHOICES %%%%%%%%%%%%%%
% counts = crosstab(D2.response,...
%     D2.perturbation,...
%     D2.sessionID,...
%     D2.stimulus(:,1)>0,...
%     D2.stimulus(:,2)>0);
% prob = counts./sum(counts,1);
% prob = nanmean(prob,3); %avg over sessions
% prob(:,1,:,:,:) = []; %Delete nolaser 
% pIpsiversive = nanmean([prob(2,2:2:end,:,2,1);... %Right choices on CL, when inactivating RIGHT hemisphere
%     prob(1,1:2:end,:,1,2)],1); %Left choices on CR, when inactivating LEFT hemisphere
% d = pIpsiversive(1) - pIpsiversive(2); %Difference in the ipsiversive rate between VIS and M2 inactivation
% 
% numIter = 10000;
% d_null = nan(1,numIter);
% f = waitbar(0,'Please wait...');
% for i = 1:numIter
%     waitbar(i/numIter,f);
%     idx = 1:length(D2.laserIdx);
%     for sess = 1:max(D2.sessionID)
%         sessIdx = idx(D2.sessionID==sess);
%         sessIdx = sessIdx(randperm(length(sessIdx)));
%         idx(D2.sessionID==sess) = sessIdx;
%     end
%     
%     counts = crosstab(D2.response,...
%         D2.perturbation( idx ),...
%         D2.sessionID,...
%         D2.stimulus(:,1)>0,...
%         D2.stimulus(:,2)>0);
%     prob = counts./sum(counts,1);
%     prob = nanmean(prob,3); %avg over sessions
%     prob(:,1,:,:,:) = []; %Delete nolaser
%     pIpsiversive = nanmean([prob(2,2:2:end,:,2,1);... %Right choices on CL, when inactivating RIGHT hemisphere
%         prob(1,1:2:end,:,1,2)],1); %Left choices on CR, when inactivating LEFT hemisphere
%     d_null(i) = pIpsiversive(1) - pIpsiversive(2); %Difference in the ipsiversive rate between VIS and M2 inactivation
% end
% pVal = min( mean(d >= d_null), mean(d <= d_null) );
% disp('Test between VIS and M2 pooled over powers:');
% disp(pVal);
% f.delete;
% 
% 
% 
% 
% % %%%%%%%%%%% EFFECT ON REACTION TIME %%%%%%%%%%%%%%
% % [G,resp,pert,sess,cl,cr,pow] = findgroups(D2.response,...
% %                                           D2.perturbation,...
% %                                           D2.sessionID,...
% %                                           D2.stimulus(:,1)>0,...
% %                                           D2.stimulus(:,2)>0,...
% %                                           D2.laserPower);
% % rt = splitapply(@nanmedian,D2.RT,G); %Calculate median RT for each combination
% % %reshape into useful array
% % vars = {resp,pert,sess,cl,cr,pow};
% % [~,~,vid] = cellfun(@(v) unique(v), vars, 'uni', 0);
% % dim = cellfun(@(v) length(unique(v)), vars);
% % RT = nan(dim); RT( sub2ind(dim,vid{1},vid{2},vid{3},vid{4},vid{5},vid{6}) ) = rt;
% % RT = nanmean(RT, 3);%average over sessions
% % RT = RT(:,2:end,:,:,:,2:end) - RT(:,1,:,:,:,1);% convert to difference from nolaser condition
% % RTContra = nanmean([RT(1,2:2:end,:,2,1,:);... %change in Left choice RT on CL, when inactivating RIGHT hemisphere
% %                     RT(2,1:2:end,:,1,2,:)],1); %change in Right choice RT on CR, when inactivating LEFT hemisphere
% % RTIpsi = nanmean([RT(2,2:2:end,:,2,1,:);... %change in Right choice RT on CL, when inactivating RIGHT hemisphere
% %                     RT(1,1:2:end,:,1,2,:)],1); %change in Left choice RT on CR, when inactivating LEFT hemisphere
% % RTContra = permute(RTContra,[2 6 1 3 4 5]); %rows regions, cols powers
% % RTIpsi = permute(RTIpsi,[2 6 1 3 4 5]);        
% % 
% % 
% % 
% % %Overlay P values by shuffling
% % numIter = 5000;
% % RTContra_null = nan(3,3,numIter);
% % RTIpsi_null = nan(3,3,numIter);
% % f = waitbar(0,'Please wait...');
% % for i = 1:numIter
% %     waitbar(i/numIter,f);
% %     idx = 1:length(D2.laserIdx);
% %     for sess = 1:max(D2.sessionID)
% %         sessIdx = idx(D2.sessionID==sess);
% %         sessIdx = sessIdx(randperm(length(sessIdx)));
% %         idx(D2.sessionID==sess) = sessIdx;
% %     end
% %     
% %     %%%%%%%%%% Compute empirical reaction times
% %     [G,resp,pert,sess,cl,cr,pow] = findgroups(D2.response,...
% %         D2.perturbation( idx ),...
% %         D2.sessionID,...
% %         D2.stimulus(:,1)>0,...
% %         D2.stimulus(:,2)>0,...
% %         D2.laserPower( idx ));
% %     rt = splitapply(@nanmedian,D2.RT,G); %Calculate median RT for each combination
% %     %reshape into useful array
% %     vars = {resp,pert,sess,cl,cr,pow};
% %     [~,~,vid] = cellfun(@(v) unique(v), vars, 'uni', 0);
% %     dim = cellfun(@(v) length(unique(v)), vars);
% %     RT = nan(dim); RT( sub2ind(dim,vid{1},vid{2},vid{3},vid{4},vid{5},vid{6}) ) = rt;
% %     RT = nanmean(RT, 3);%average over sessions
% %     RT = RT(:,2:end,:,:,:,2:end) - RT(:,1,:,:,:,1);% convert to difference from nolaser condition
% %     RTC = nanmean([RT(1,2:2:end,:,2,1,:);... %change in Left choice RT on CL, when inactivating RIGHT hemisphere
% %         RT(2,1:2:end,:,1,2,:)],1); %change in Right choice RT on CR, when inactivating LEFT hemisphere
% %     RTI = nanmean([RT(2,2:2:end,:,2,1,:);... %change in Right choice RT on CL, when inactivating RIGHT hemisphere
% %         RT(1,1:2:end,:,1,2,:)],1); %change in Left choice RT on CR, when inactivating LEFT hemisphere
% %     RTContra_null(:,:,i) = permute(RTC,[2 6 1 3 4 5]); %rows regions, cols powers
% %     RTIpsi_null(:,:,i) = permute(RTI,[2 6 1 3 4 5]);
% % end
% % f.delete;
% % 
% % figure('color','w');
% % for region = 1:3
% %     subplot(1,2,1); hold on;
% %     plot(region + [-1 0 1]*0.2, RTContra(region,:),'.-','markersize',20);
% %     null = squeeze(RTContra_null(region,:,:))';
% %     pVal = min(mean(RTContra(region,:) <= null,1), mean(RTContra(region,:) >= null,1));
% %     text(region + [-1 0 1]*0.2, RTContra(region,:),arrayfun(@(n) num2str(n),pVal,'uni',0) );
% %     
% %     subplot(1,2,2); hold on;
% %     plot(region + [-1 0 1]*0.2, RTIpsi(region,:),'.-','markersize',20);
% %     null = squeeze(RTIpsi_null(region,:,:))';
% %     pVal = min(mean(RTIpsi(region,:) <= null,1), mean(RTIpsi(region,:) >= null,1));
% %     text(region + [-1 0 1]*0.2, RTIpsi(region,:),arrayfun(@(n) num2str(n),pVal,'uni',0) );
% % end
% % set(get(gcf,'children'),'xtick',1:3,'xticklabel',regionLabels);
% % 
% %% Multi-power: model-free statistical testing of laser effects COMBINING POWERS
% 
% %%%%%%%%%%% GET DATA %%%%%%%%%%%%%%
% perturbationRegion = {'LeftVIS','RightVIS','LeftM2','RightM2','LeftM1','RightM1'};
% regionLabels = {'VIS','MOs','MOp'};
% areaCols = [ 73 148 208;
%              119 172 66;
%              100 100 100]/255;
% D2 = getrow(D, any(D.laserRegion == perturbationRegion,2) | D.laserType==0);
% D2.perturbation = zeros(size(D2.response));
% for p = 1:length(perturbationRegion)
%     D2.perturbation(D2.laserRegion == perturbationRegion{p}) = p;
% end
% 
% D2 = getrow(D2,D2.stimulus(:,1)>=0.1 | D2.stimulus(:,2)>=0.1);
% 
% %%%%%%%%%%% EFFECT ON FRACTION OF CHOICES %%%%%%%%%%%%%%
% counts = crosstab(D2.response,...
%     D2.perturbation,...
%     D2.sessionID,...
%     D2.stimulus(:,1)>0,...
%     D2.stimulus(:,2)>0);
% prob = counts./sum(counts,1);
% pContraNL = nanmean([prob(1,1,:,2,1);... % Left choices on CL
%                      prob(2,1,:,1,2)],1); %Right choices on CR
% pIpsiNL = nanmean([prob(2,1,:,2,1);... % Right choices on CL
%                      prob(1,1,:,1,2)],1); %Left choices on CR
% pContraNL = squeeze(pContraNL);
% pIpsiNL = squeeze(pIpsiNL);
% 
% prob = prob(:,2:end,:,:,:);
% pContraversive_sess = nanmean([prob(1,2:2:end,:,2,1);... % Left choices on CL, when inactivating RIGHT hemisphere
%     prob(2,1:2:end,:,1,2)],1); %Right choices on CR, when inactivating LEFT hemisphere
% pIpsiversive_sess = nanmean([prob(2,2:2:end,:,2,1);... %Right choices on CL, when inactivating RIGHT hemisphere
%     prob(1,1:2:end,:,1,2)],1); %Left choices on CR, when inactivating LEFT hemisphere
% pContraversive_sess = squeeze(pContraversive_sess)';
% pIpsiversive_sess = squeeze(pIpsiversive_sess)';
% 
% figure('color','w');
% subplot(1,2,1); hold on;
% bar(nanmean(pContraversive_sess,1),'BaseValue',nanmean(pContraNL),'EdgeAlpha',0); 
% se = nanstd(pContraversive_sess,[],1)./sqrt(sum(~isnan(pContraversive_sess),1));
% h=errorbar(1:3, nanmean(pContraversive_sess,1), se, 'LineStyle' ,'none');
% for i = 1:3
%     [H,Pval] = ttest(pContraNL,pContraversive_sess(:,i));
%     text(i, nanmean(pContraversive_sess(:,i)), num2str(Pval) );
% end
% 
% subplot(1,2,2); hold on;
% bar(nanmean(pIpsiversive_sess,1),'BaseValue',nanmean(pIpsiNL),'EdgeAlpha',0); 
% se = nanstd(pIpsiversive_sess,[],1)./sqrt(sum(~isnan(pIpsiversive_sess),1));
% h=errorbar(1:3, nanmean(pIpsiversive_sess,1), se, 'LineStyle' ,'none');
% for i = 1:3
%     [H,Pval] = ttest(pIpsiNL,pIpsiversive_sess(:,i));
%     text(i, nanmean(pIpsiversive_sess(:,i)), num2str(Pval) );
% end
% set(get(gcf,'children'),'ylim',[0 1]);
% set(get(gcf,'children'),'xtick',1:3,'xticklabel',regionLabels);
% 
% %pairwise test between VIS and MOs ipsiversive
% [~,Pval]=ttest( pIpsiversive_sess(:,1), pIpsiversive_sess(:,2));
%     text(1.5, 0.5, num2str(Pval) );
%% Multi-power: fisher exact test: per session + pooling across power
 
%Test that inactivation significantly decreases contraversive choices
region = {'VIS','M2','M1'};
for reg = 1:3
    pvals = nan(max(D.sessionID),1);
    pChoice = nan(max(D.sessionID),2);
    for sess = 1:max(D.sessionID)
        d3 = getrow(D,D.sessionID==sess);
        X = nan(2,2);
        
        X(1,1) =  sum( (d3.laserType==0 & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==2)) | ...
            (d3.laserType==0 & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==1)) );
        
        X(2,1) =  sum( (d3.laserType==0 & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=2)) | ...
            (d3.laserType==0 & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=1)) );
        
        X(1,2) =  sum( (d3.laserRegion==['Left' region{reg}] & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==2)) | ...
            (d3.laserRegion==['Right' region{reg}] & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==1)) );
        
        X(2,2) =  sum( (d3.laserRegion==['Left' region{reg}] & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=2)) | ...
            (d3.laserRegion==['Right' region{reg}] & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=1)) );
        
        [~,pvals(sess)] = fishertest(X); %Testing DECREASE
        pChoice(sess,:) = X(1,:)./sum(X,1);
    end
    avg_pchoice = nanmean(pChoice,1);
    chi_vals = -2.*log(pvals);
    group_pval = 1 - chi2cdf(sum(chi_vals),2*length(pvals));
    fprintf('Contraversive %s: %0.2f vs %0.2f  %0.10f\n',region{reg},avg_pchoice(1), avg_pchoice(2), group_pval);
end



%Test that inactivation significantly changes ipsiversive choices
region = {'VIS','M2','M1'};
for reg = 1:3
    pvals = nan(max(D.sessionID),1);
    pChoice = nan(max(D.sessionID),2);
    for sess = 1:max(D.sessionID)
        d3 = getrow(D,D.sessionID==sess);
        X = nan(2,2);
        
        X(1,1) =  sum( (d3.laserType==0 & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==1)) | ...
            (d3.laserType==0 & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==2)) );
        
        X(2,1) =  sum( (d3.laserType==0 & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=1)) | ...
            (d3.laserType==0 & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=2)) );
        
        X(1,2) =  sum( (d3.laserRegion==['Left' region{reg}] & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==1)) | ...
            (d3.laserRegion==['Right' region{reg}] & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==2)) );
        
        X(2,2) =  sum( (d3.laserRegion==['Left' region{reg}] & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=1)) | ...
            (d3.laserRegion==['Right' region{reg}] & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=2)) );
        
        [~,pvals(sess)] = fishertest(X, 'tail', 'left'); %Testing INCREASE
        pChoice(sess,:) = X(1,:)./sum(X,1);
    end
    avg_pchoice = nanmean(pChoice,1);
    chi_vals = -2.*log(pvals);
    group_pval = 1 - chi2cdf(sum(chi_vals),2*length(pvals));
    fprintf('Ipsiversive %s: %0.2f vs %0.2f  %0.10f\n',region{reg},avg_pchoice(1), avg_pchoice(2), group_pval);
end




%Test that the rate of ipsiversive choices for VIS vs M2 inactivation
pvals = nan(max(D.sessionID),1);
pChoice = nan(max(D.sessionID),2);
for sess = 1:max(D.sessionID)
    d3 = getrow(D,D.sessionID==sess);
    X = nan(2,2);
      
    X(1,1) =  sum( (d3.laserRegion=='LeftVIS' & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==1)) | ...
                   (d3.laserRegion=='RightVIS' & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==2)) );
    
    X(2,1) =  sum( (d3.laserRegion=='LeftVIS' & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=1)) | ...
                   (d3.laserRegion=='RightVIS' & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=2)) );
   
    X(1,2) =  sum( (d3.laserRegion=='LeftM2' & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==1)) | ...
                   (d3.laserRegion=='RightM2' & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==2)) );
    
    X(2,2) =  sum( (d3.laserRegion=='LeftM2' & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=1)) | ...
                   (d3.laserRegion=='RightM2' & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=2)) );
       
    [~,pvals(sess)] = fishertest(X,'tail','right');
    pChoice(sess,:) = X(1,:)./sum(X,1);
end
avg_pchoice = nanmean(pChoice,1);
chi_vals = -2.*log(pvals);
group_pval = 1 - chi2cdf(sum(chi_vals),2*length(pvals));
fprintf('Ipsiversive VIS vs MOs: %0.2f vs %0.2f %0.10f\n',avg_pchoice(1), avg_pchoice(2), group_pval);





%Plot
region = {'VIS','M2','M1'};
pChoice = nan(max(D.sessionID),4,2);
for sess = 1:max(D.sessionID)
    d3 = getrow(D,D.sessionID==sess);
    
    a1 =  sum( (d3.laserType==0 & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==2)) | ...
        (d3.laserType==0 & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==1)) );
    
    a2 =  sum( (d3.laserType==0 & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=2)) | ...
        (d3.laserType==0 & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=1)) );
    
    pChoice(sess,1,1) = a1./(a1+a2);
    
    
    b1 =  sum( (d3.laserType==0 & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==1)) | ...
            (d3.laserType==0 & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==2)) );
        
    b2 =  sum( (d3.laserType==0 & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=1)) | ...
            (d3.laserType==0 & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=2)) );
        
    pChoice(sess,1,2) = b1./(b1+b2);


    
    for reg = 1:3
        a1 =  sum( (d3.laserRegion==['Left' region{reg}] & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==2)) | ...
            (d3.laserRegion==['Right' region{reg}] & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==1)) );
        
        a2 =  sum( (d3.laserRegion==['Left' region{reg}] & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=2)) | ...
            (d3.laserRegion==['Right' region{reg}] & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=1)) );
        pChoice(sess,1+reg,1) = a1./(a1+a2);
        
        
        b1 =  sum( (d3.laserRegion==['Left' region{reg}] & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==1)) | ...
            (d3.laserRegion==['Right' region{reg}] & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==2)) );
        
        b2 =  sum( (d3.laserRegion==['Left' region{reg}] & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=1)) | ...
            (d3.laserRegion==['Right' region{reg}] & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=2)) );
         pChoice(sess,1+reg,2) = b1./(b1+b2);
        
    end
end

pChoice(:,:,3) = 1 - sum(pChoice,3);
avg = nanmean(pChoice,1);
se = nanstd(pChoice,[],1)./sqrt(sum(~isnan(pChoice),1));
lab = {'pContraversive','pIpsiversive','pNG'};
figure;
x=repmat(1:4,max(D.sessionID),1);
x=x+ randn(size(x))/30;
for i = 1:3
    subplot(1,3,i); hold on;

    plot(x', pChoice(:,:,i)', '.-', 'color', [ 1 1 1]*0.75, 'markersize',10);
    plot(1:4, avg(:,:,i), 'k.','markersize',25);
    for s = 1:4
        line([1 1]*s, avg(:,s,i) + [-1 1]*se(:,s,i), 'Color','k')
    end
    xlim([0 5]);
    set(gca, 'xticklabel', ['off' region], 'xtick', 1:4);
    ylabel(lab{i});
end
% 
% figure;
% for i = 1:2
%     subplot(1,2,i); hold on;
%     
%     boxplot( pChoice(:,:,i), 'notch', 'on');% , 'plotstyle', 'compact');
%     
%     xlim([0 5]);
%     set(gca, 'xticklabel', ['off' region], 'xtick', 1:4)
%     ylabel(lab{i});
% end



%% Multi-power: fisher exact test: per subject + pooling across power
 
%Test that inactivation significantly decreases contraversive choices
region = {'VIS','M2','M1'};
for reg = 1:3
    pvals = nan(max(D.subjectID),1);
    pChoice = nan(max(D.subjectID),2);
    for subj = 1:max(D.subjectID)
        d3 = getrow(D,D.subjectID==subj);
        X = nan(2,2);
        
        X(1,1) =  sum( (d3.laserType==0 & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==2)) | ...
            (d3.laserType==0 & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==1)) );
        
        X(2,1) =  sum( (d3.laserType==0 & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=2)) | ...
            (d3.laserType==0 & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=1)) );
        
        X(1,2) =  sum( (d3.laserRegion==['Left' region{reg}] & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==2)) | ...
            (d3.laserRegion==['Right' region{reg}] & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==1)) );
        
        X(2,2) =  sum( (d3.laserRegion==['Left' region{reg}] & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=2)) | ...
            (d3.laserRegion==['Right' region{reg}] & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=1)) );
        
        [~,pvals(subj)] = fishertest(X);
        pChoice(subj,:) = X(1,:)./sum(X,1);
    end
    avg_pchoice = nanmean(pChoice,1);
    chi_vals = -2.*log(pvals);
    group_pval = 1 - chi2cdf(sum(chi_vals),2*length(pvals));
    fprintf('Contraversive %s: %0.2f vs %0.2f  %0.10f\n',region{reg},avg_pchoice(1), avg_pchoice(2), group_pval);
end



%Test that inactivation significantly changes ipsiversive choices
region = {'VIS','M2','M1'};
for reg = 1:3
    pvals = nan(max(D.subjectID),1);
    pChoice = nan(max(D.subjectID),2);
    for subj = 1:max(D.subjectID)
        d3 = getrow(D,D.subjectID==subj);
        X = nan(2,2);
        
        X(1,1) =  sum( (d3.laserType==0 & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==1)) | ...
            (d3.laserType==0 & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==2)) );
        
        X(2,1) =  sum( (d3.laserType==0 & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=1)) | ...
            (d3.laserType==0 & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=2)) );
        
        X(1,2) =  sum( (d3.laserRegion==['Left' region{reg}] & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==1)) | ...
            (d3.laserRegion==['Right' region{reg}] & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==2)) );
        
        X(2,2) =  sum( (d3.laserRegion==['Left' region{reg}] & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=1)) | ...
            (d3.laserRegion==['Right' region{reg}] & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=2)) );
        
        [~,pvals(subj)] = fishertest(X);
        pChoice(subj,:) = X(1,:)./sum(X,1);
    end
    avg_pchoice = nanmean(pChoice,1);
    chi_vals = -2.*log(pvals);
    group_pval = 1 - chi2cdf(sum(chi_vals),2*length(pvals));
    fprintf('Ipsiversive %s: %0.2f vs %0.2f  %0.10f\n',region{reg},avg_pchoice(1), avg_pchoice(2), group_pval);
end






%Test that the rate of ipsiversive choices for VIS vs M2 inactivation
pvals = nan(max(D.subjectID),1);
pChoice = nan(max(D.subjectID),2);
for subj = 1:max(D.subjectID)
    d3 = getrow(D,D.subjectID==subj);
    X = nan(2,2);
      
    X(1,1) =  sum( (d3.laserRegion=='LeftVIS' & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==1)) | ...
                   (d3.laserRegion=='RightVIS' & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==2)) );
    
    X(2,1) =  sum( (d3.laserRegion=='LeftVIS' & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=1)) | ...
                   (d3.laserRegion=='RightVIS' & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=2)) );
   
    X(1,2) =  sum( (d3.laserRegion=='LeftM2' & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==1)) | ...
                   (d3.laserRegion=='RightM2' & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==2)) );
    
    X(2,2) =  sum( (d3.laserRegion=='LeftM2' & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=1)) | ...
                   (d3.laserRegion=='RightM2' & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=2)) );
       
    [~,pvals(subj)] = fishertest(X);
    pChoice(subj,:) = X(1,:)./sum(X,1);
end
avg_pchoice = nanmean(pChoice,1);
chi_vals = -2.*log(pvals);
group_pval = 1 - chi2cdf(sum(chi_vals),2*length(pvals));
fprintf('Ipsiversive VIS vs MOs: %0.2f vs %0.2f %0.10f\n',avg_pchoice(1), avg_pchoice(2), group_pval);
%% Multi-power: fisher exact test: per subject + per power level
 
%Test that inactivation significantly decreases contraversive choices
region = {'VIS','M2','M1'};
powers = [1.5 2.9 4.25];
for reg = 1:3
    for pow = 1:3
        pvals = nan(max(D.subjectID),1);
        pChoice = nan(max(D.subjectID),2);
        for subj = 1:max(D.subjectID)
            d3 = getrow(D,D.subjectID==subj & (D.laserPower==powers(pow) | D.laserType==0));
            X = nan(2,2);
            
            X(1,1) =  sum( (d3.laserType==0 & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==2)) | ...
                (d3.laserType==0 & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==1)) );
            
            X(2,1) =  sum( (d3.laserType==0 & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=2)) | ...
                (d3.laserType==0 & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=1)) );
            
            X(1,2) =  sum( (d3.laserRegion==['Left' region{reg}] & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==2)) | ...
                (d3.laserRegion==['Right' region{reg}] & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==1)) );
            
            X(2,2) =  sum( (d3.laserRegion==['Left' region{reg}] & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=2)) | ...
                (d3.laserRegion==['Right' region{reg}] & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=1)) );
            
            [~,pvals(subj)] = fishertest(X);
            pChoice(subj,:) = X(1,:)./sum(X,1);
        end
        avg_pchoice = nanmean(pChoice,1);
        chi_vals = -2.*log(pvals);
        group_pval = 1 - chi2cdf(sum(chi_vals),2*length(pvals));
        fprintf('Contraversive %s (%0.2f): %0.2f vs %0.2f  %0.10f\n',region{reg},powers(pow),avg_pchoice(1), avg_pchoice(2), group_pval);
    end
end


%Test that inactivation significantly changes ipsiversive choices
region = {'VIS','M2','M1'};
powers = [1.5 2.9 4.25];
for reg = 1:3
    for pow = 1:3
        pvals = nan(max(D.subjectID),1);
        pChoice = nan(max(D.subjectID),2);
        for subj = 1:max(D.subjectID)
            d3 = getrow(D,D.subjectID==subj & (D.laserPower==powers(pow) | D.laserType==0));
            X = nan(2,2);
            
            X(1,1) =  sum( (d3.laserType==0 & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==1)) | ...
                (d3.laserType==0 & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==2)) );
            
            X(2,1) =  sum( (d3.laserType==0 & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=1)) | ...
                (d3.laserType==0 & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=2)) );
            
            X(1,2) =  sum( (d3.laserRegion==['Left' region{reg}] & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==1)) | ...
                (d3.laserRegion==['Right' region{reg}] & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==2)) );
            
            X(2,2) =  sum( (d3.laserRegion==['Left' region{reg}] & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=1)) | ...
                (d3.laserRegion==['Right' region{reg}] & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=2)) );
            
            [~,pvals(subj)] = fishertest(X);
            pChoice(subj,:) = X(1,:)./sum(X,1);
        end
        avg_pchoice = nanmean(pChoice,1);
        chi_vals = -2.*log(pvals);
        group_pval = 1 - chi2cdf(sum(chi_vals),2*length(pvals));
        fprintf('Ipsiversive %s (%0.2f): %0.2f vs %0.2f  %0.10f\n',region{reg},powers(pow),avg_pchoice(1), avg_pchoice(2), group_pval);
    end
end



%Test that the rate of ipsiversive choices for VIS vs M2 inactivation
powers = [1.5 2.9 4.25];
for pow = 1:3
    pvals = nan(max(D.subjectID),1);
    pChoice = nan(max(D.subjectID),2);
    for subj = 1:max(D.subjectID)
        d3 = getrow(D,D.subjectID==subj & (D.laserPower==powers(pow) | D.laserType==0));

        X = nan(2,2);
        
        X(1,1) =  sum( (d3.laserRegion=='LeftVIS' & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==1)) | ...
            (d3.laserRegion=='RightVIS' & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==2)) );
        
        X(2,1) =  sum( (d3.laserRegion=='LeftVIS' & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=1)) | ...
            (d3.laserRegion=='RightVIS' & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=2)) );
        
        X(1,2) =  sum( (d3.laserRegion=='LeftM2' & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response==1)) | ...
            (d3.laserRegion=='RightM2' & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response==2)) );
        
        X(2,2) =  sum( (d3.laserRegion=='LeftM2' & (d3.stimulus(:,1)==0 & d3.stimulus(:,2)>0 & d3.response~=1)) | ...
            (d3.laserRegion=='RightM2' & (d3.stimulus(:,1)>0 & d3.stimulus(:,2)==0 & d3.response~=2)) );
        
        [~,pvals(subj)] = fishertest(X);
        pChoice(subj,:) = X(1,:)./sum(X,1);
    end
    avg_pchoice = nanmean(pChoice,1);
    chi_vals = -2.*log(pvals);
    group_pval = 1 - chi2cdf(sum(chi_vals),2*length(pvals));
    fprintf('Ipsiversive VIS vs MOs (%0.2f): %0.2f vs %0.2f %0.10f\n',powers(pow),avg_pchoice(1), avg_pchoice(2), group_pval);
end
%% Multi-power: fit mechanistic model 

%1) Load fit from widefield to get weight posterior, to use as prior for next
%fit
fit_mech = load("C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\tp03f24ce4_5470_4181_b95f_cfaf8abe5b33.mat");

% % % Wrong one?=
% % % fit_mech = load("C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\tp2ade4d28_eb4c_4b67_9974_d8ab0391ee9c.mat");


w_mu = mean( fit_mech.posterior.w, 1);
w_S = cov( fit_mech.posterior.w);

%2) Compute firing_rate values on no-laser trials
perturbationRegion = {'LeftVIS','RightVIS','LeftM2','RightM2'};
D1 = getrow(D, any(D.laserRegion == perturbationRegion,2) | D.laserType==0);
D1.perturbation = zeros(size(D1.response));
for p = 1:length(perturbationRegion)
    D1.perturbation(D1.laserRegion == perturbationRegion{p}) = p;
end
NL_idx = D1.laserType==0;

fr = load("C:\Users\Peter\Documents\MATLAB\inactPaper\widefield\mechanistic_fit\visMOs_allConds.mat");

VIS_window = [0.075 0.125];
MOS_window = [0.125 0.175];

%get ephys data at critical timewindow
vis = mean( mean( fr.vis(:,:,:,VIS_window(1) < fr.binsS & fr.binsS < VIS_window(2)), 4), 3);
m2 = mean( mean( fr.mos(:,:,:,MOS_window(1) < fr.binsS & fr.binsS < MOS_window(2)), 4), 3);

F=[];
F(:,1) = interp2(fr.ucl, fr.ucr, vis, D1.stimulus(:,1) , D1.stimulus(:,2), 'linear' );
F(:,2) = interp2(fr.ucl, fr.ucr, vis', D1.stimulus(:,1) , D1.stimulus(:,2), 'linear' );
F(:,3) = interp2(fr.ucl, fr.ucr, m2, D1.stimulus(:,1) , D1.stimulus(:,2), 'linear' );
F(:,4) = interp2(fr.ucl, fr.ucr, m2', D1.stimulus(:,1) , D1.stimulus(:,2), 'linear' );

dat=struct;
dat.choice = D1.response(NL_idx);
dat.sessionID = D1.sessionID(NL_idx);
dat.subjectID = D1.subjectID(NL_idx);
dat.firing_rate = F(NL_idx,:);
dat.PRIOR_MU = w_mu;
dat.PRIOR_S = w_S;
dat.PRIOR_S(1,1) = 1; %BROADEN PRIOR ON ALPHAS
dat.PRIOR_S(2,2) = 1; %BROADEN PRIOR ON ALPHAS
dat.PRIOR_S = dat.PRIOR_S + diag(diag(w_S)); %BROADEN PRIOR ON ALPHAS
fit_inact = bfit.fitModel('Two-Level-Mechanistic-SYMMETRICAL-PRESETW',dat);
%% Multi-power: plot mechanistic model (cross-prediction from widefield)
fit_inact = load('C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\tpad8496b2_6089_416e_941c_8dee1fb09ce0.mat');

%Others?
% fit_inact_A = load('C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\tp930c9721_46be_4191_8dbe_1e9d8f544eb7.mat');
fit_inact = load('C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\tp40284ff8_d2c3_4249_974e_7ace3b387de9.mat');

figure('color','w');

%Plot posterior of Weights, copying over symmetrical values
areaCols = [ 109 183 229;
    77 131 158;
    160 206 87;
    119 147 62]/255;
subplot(1,2,1); hold on;
WL = fit_mech.posterior.w(:,2+[1:4]);
WR = fit_mech.posterior.w(:,2+[2 1 4 3]);
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
xlabel('W_L'); ylabel('W_R'); title('Mech fit from widefield');


%Plot posterior of Weights, copying over symmetrical values
subplot(1,2,2); hold on;
WL = fit_inact.posterior.w(:,2+[1:4]);
WR = fit_inact.posterior.w(:,2+[2 1 4 3]);
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
xlabel('W_L'); ylabel('W_R'); title('Refit on inactivation sessions');


%Plot model on behavioural data on detection contrast stimuli
D1.pedestal = min(D1.stimulus(:,1), D1.stimulus(:,2));
D1.cDiff = D1.stimulus(:,2) - D1.stimulus(:,1);
counts = crosstab(D1.pedestal,...
                  D1.cDiff,...
                  D1.response,...
                  D1.sessionID,...
                  D1.perturbation);
prob = counts./sum(counts,3); 
prob = nanmean(prob, 4); %avg over sessions
prob = prob(1,:,:,:,:); %only care about zero pedestal
prob = permute(prob,[2 3 5 1 4]); % cD x Resp x Perturb
cdU = unique(D1.cDiff);


CL = [linspace(1,0,500), zeros(1,500)];
CR = [zeros(1,500), linspace(0,1,500)];
F=[];
F(1,:) = interp2(fr.ucl, fr.ucr, vis, CL , CR, 'linear' );
F(2,:) = interp2(fr.ucl, fr.ucr, vis', CL , CR, 'linear' );
F(3,:) = interp2(fr.ucl, fr.ucr, m2, CL , CR, 'linear' );
F(4,:) = interp2(fr.ucl, fr.ucr, m2', CL , CR, 'linear' );

NL_ph=bplot.MECH_SYMETRICAL(fit_inact.posterior.w, F);
mean_ph = mean(NL_ph,1);
q_ph = quantile(NL_ph,[0.025 0.975],1);

figure('color','w');
ha = tight_subplot(4,3,0.05,0.1,0.1);
for f = 1:4
    
    F_inact = F;
    F_inact(f,:) = 0;
    Lph=bplot.MECH_SYMETRICAL(fit_inact.posterior.w, F_inact);
    mean_Lph = mean(Lph,1);
    q_Lph = quantile(Lph,[0.025 0.975],1);
    
    for r = 1:3
        hold( ha( 3*(f-1) + r ), 'on');
        plot(ha( 3*(f-1) + r ), CR-CL, mean_ph(1,:,r) , 'k-');
        fx=fill( ha( 3*(f-1) + r ), [CR-CL fliplr(CR-CL)],[q_ph(1,:,r) fliplr(q_ph(2,:,r))],'k'); fx.FaceAlpha=0.2; fx.EdgeAlpha=0;

        plot( ha( 3*(f-1) + r ), cdU, prob(:,r,1), 'k.', 'markersize',15);
        
        plot( ha( 3*(f-1) + r ), CR-CL, mean_Lph(1,:,r), '-', 'color', areaCols(f,:));
        fx=fill( ha( 3*(f-1) + r ), [CR-CL fliplr(CR-CL)],[q_Lph(1,:,r) fliplr(q_Lph(2,:,r))],'r'); fx.FaceAlpha=0.2; fx.EdgeAlpha=0;
        fx.FaceColor = areaCols(f,:);
        plot( ha( 3*(f-1) + r ), cdU, prob(:,r,1+f), '.','color',areaCols(f,:), 'markersize',15);

    end

end
set(ha,'ylim',[0 1],'dataaspectratio',[1 1 1]);

% Plot change to contraversive and ipsiversive on trials with stimuli on
% the contralateral side only
% Since model is symmetrical i can just look at one side: high CL only,
% inactivating RIGHT hemisphere and looking at changes in Left and Right
% choics
CL = unique(D1.stimulus(:,1)); CL(1) = [];
CR = zeros(size(CL));
proportionLeftContrast = tabulate(D1.stimulus(D1.stimulus(:,2)==0,1));
proportionLeftContrast = proportionLeftContrast(:,2)'; proportionLeftContrast(1)=[];
proportionLeftContrast = proportionLeftContrast/sum(proportionLeftContrast);
proportionLeftContrast = [1 1 1]/3;

F=[];
F(1,:) = interp2(fr.ucl, fr.ucr, vis, CL , CR, 'linear' );
F(2,:) = interp2(fr.ucl, fr.ucr, vis', CL , CR, 'linear' );
F(3,:) = interp2(fr.ucl, fr.ucr, m2, CL , CR, 'linear' );
F(4,:) = interp2(fr.ucl, fr.ucr, m2', CL , CR, 'linear' );

NL_ph=bplot.MECH_SYMETRICAL(fit_inact.posterior.w, F);
mean_ph = mean(NL_ph,1);
pContraNL = sum( mean_ph(:,:,1).*proportionLeftContrast );
pIpsiNL = sum( mean_ph(:,:,2).*proportionLeftContrast );

%Now look at proportion of contra and ipsi when simulating inactivation of
%right hemisphere

%R VIS:
F_inact=F;
F_inact(2,:)=0;
ph=bplot.MECH_SYMETRICAL(fit_inact.posterior.w, F_inact);
mean_ph = mean(ph,1);
pContraVIS = sum( mean_ph(:,:,1).*proportionLeftContrast );
pIpsiVIS = sum( mean_ph(:,:,2).*proportionLeftContrast );

%R MOs:
F_inact=F;
F_inact(4,:)=0;
ph=bplot.MECH_SYMETRICAL(fit_inact.posterior.w, F_inact);
mean_ph = mean(ph,1);
pContraM2 = sum( mean_ph(:,:,1).*proportionLeftContrast );
pIpsiM2 = sum( mean_ph(:,:,2).*proportionLeftContrast );

figure;
subplot(1,2,1);
hx=bar([pContraVIS pContraM2],'BaseValue',pContraNL); 
ylim([0 1]);
subplot(1,2,2);
hx=bar([pIpsiVIS pIpsiM2],'BaseValue',pIpsiNL); 
ylim([0 0.5]);


% 
% %%%TODO: UPDATE THIS CODE
% %Create a contraversive and ipsiversive plot for stimuli present only on
% %the contralateral side
% %Since model is symmetrical i can just look at one side: high CL only
% % CL = unique(fit_mech.data.contrastLeft); CL(1) = [];
% CL = 1;
% CR = zeros(size(CL));
% 
% F=[];
% F(1,:) = interp2(fr.ucl, fr.ucr, vis, CL , CR, 'linear' );
% F(2,:) = interp2(fr.ucl, fr.ucr, vis', CL , CR, 'linear' );
% F(3,:) = interp2(fr.ucl, fr.ucr, m2, CL , CR, 'linear' );
% F(4,:) = interp2(fr.ucl, fr.ucr, m2', CL , CR, 'linear' );
% NL_ph=bplot.MECH_SYMETRICAL(fit_inact.posterior.w, F);
% mean_ph = mean(NL_ph,1);
% pContraNL = mean_ph(:,:,1);
% pIpsiNL = mean_ph(:,:,2);
% 
% %Inactivate right hemisphere VIS
% F_inact = F;
% F_inact(2,:) = 0;
% ph=bplot.MECH_SYMETRICAL(fit_inact.posterior.w, F_inact);
% mean_ph_v = mean(ph,1);
% q_ph_v = quantile(ph,[0.025 0.975],1);
% 
% %inactivate right hemisphere MOs
% F_inact = F;
% F_inact(4,:) = 0;
% ph=bplot.MECH_SYMETRICAL(fit_inact.posterior.w, F_inact);
% mean_ph_m = mean(ph,1);
% q_ph_m = quantile(ph,[0.025 0.975],1);
% 
% %Plot:
% figure('color','w');
% subplot(1,2,1); hold on;
% bar([mean_ph_v(:,:,1) mean_ph_m(:,:,1)],'BaseValue',pContraNL,'EdgeAlpha',0); 
% line([1 1],q_ph_v(:,:,1),'color',[1 1 1]*0.5,'linewidth',8);
% line([2 2],q_ph_m(:,:,1),'color',[1 1 1]*0.5,'linewidth',8);
% ylabel('pContraversive'); ylim([0 1]); xlim([0.5 2.5]);
% 
% subplot(1,2,2); hold on;
% bar([mean_ph_v(:,:,2) mean_ph_m(:,:,2)],'BaseValue',pIpsiNL,'EdgeAlpha',0); 
% line([1 1],q_ph_v(:,:,2),'color',[1 1 1]*0.5,'linewidth',8);
% line([2 2],q_ph_m(:,:,2),'color',[1 1 1]*0.5,'linewidth',8);
% ylabel('pIpsiversive');  ylim([0 0.5]); xlim([0.5 2.5]);
% 
% regionLabels = {'VIS','MOs'};
% set(get(gcf,'children'),'xtick',1:2,'xticklabel',regionLabels);

%% Search for the right model file for the contraversive/ipsiversive plot
files = dir('C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\*.mat');

files = {'tp1c5203a6_f177_4efe_974c_3459eaa9c6fc.mat',
'tp40284ff8_d2c3_4249_974e_7ace3b387de9.mat',
'tp5241d093_6606_4ebe_8672_8c4a7d200a6f.mat',
'tp576ac406_847a_40e9_8d40_7a3a01b9bf07.mat',
'tp5e943931_553a_42e9_a3c4_4b36c299ee2f.mat',
'tp930c9721_46be_4191_8dbe_1e9d8f544eb7.mat',
'tpaa5d2776_4c38_4620_87e6_dcae019a97e6.mat',
'tpad8496b2_6089_416e_941c_8dee1fb09ce0.mat',
'tpb3a49ae1_6438_410b_a11a_4d644ca18dc9.mat',
'tpb5e0c635_7825_4749_b17c_c81db9b7ee03.mat',
'tpc5dd3f7a_2bcc_4d7d_bc98_1f6cb7093f42.mat',
'tpe1ae1bc7_676e_442a_96c8_c3a6164d2ee7.mat',
'tpeebd643f_21bc_476c_836c_a1ff8b275594.mat'
};

CL = unique(D1.stimulus(:,1)); CL(1) = [];
CR = zeros(size(CL));
proportionLeftContrast = tabulate(D1.stimulus(D1.stimulus(:,2)==0,1));
proportionLeftContrast = proportionLeftContrast(:,2)'; proportionLeftContrast(1)=[];
proportionLeftContrast = proportionLeftContrast/sum(proportionLeftContrast);
% proportionLeftContrast = [1 1 1]/3;

F=[];
F(1,:) = interp2(fr.ucl, fr.ucr, vis, CL , CR, 'linear' );
F(2,:) = interp2(fr.ucl, fr.ucr, vis', CL , CR, 'linear' );
F(3,:) = interp2(fr.ucl, fr.ucr, m2, CL , CR, 'linear' );
F(4,:) = interp2(fr.ucl, fr.ucr, m2', CL , CR, 'linear' );


for f = 1:length(files)
    fit = load(fullfile('C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\', files{f}));
    if contains(fit.modelName, 'Two-Level-Mechanistic-SYMMETRICAL-PRESETW')
        disp(files{f});
        
        
        NL_ph=bplot.MECH_SYMETRICAL(fit.posterior.w, F);
        mean_ph = mean(NL_ph,1);
        pContraNL = sum( mean_ph(:,:,1).*proportionLeftContrast );
        pIpsiNL = sum( mean_ph(:,:,2).*proportionLeftContrast );

        %Now look at proportion of contra and ipsi when simulating inactivation of
        %right hemisphere

        %R VIS:
        F_inact=F;
        F_inact(2,:)=0;
        ph=bplot.MECH_SYMETRICAL(fit.posterior.w, F_inact);
        mean_ph = mean(ph,1);
        pContraVIS = sum( mean_ph(:,:,1).*proportionLeftContrast );
        pIpsiVIS = sum( mean_ph(:,:,2).*proportionLeftContrast );

        %R MOs:
        F_inact=F;
        F_inact(4,:)=0;
        ph=bplot.MECH_SYMETRICAL(fit.posterior.w, F_inact);
        mean_ph = mean(ph,1);
        pContraM2 = sum( mean_ph(:,:,1).*proportionLeftContrast );
        pIpsiM2 = sum( mean_ph(:,:,2).*proportionLeftContrast );

        figure;
        subplot(1,2,1);
        hx=bar([pContraVIS pContraM2],'BaseValue',pContraNL);
        ylim([0 1]);
        subplot(1,2,2);
        hx=bar([pIpsiVIS pIpsiM2],'BaseValue',pIpsiNL);
        ylim([0 0.5]);

        set(gcf,'name', files{f});
        drawnow;
    end
end

%% Multi-power: Fit phenomenological model
perturbationRegion = {'LeftVIS','RightVIS','LeftM2','RightM2'};
D1 = getrow(D, any(D.laserRegion == perturbationRegion,2) | D.laserType==0);

%Hard-coded firing rates
fr_vis_contra = @(c_contra) 1 + 7.83*c_contra.^0.74;
fr_mos_contra = @(c_contra) 1 + 2.66*c_contra.^0.41;
fr_vis_ipsi   = @(c_ipsi) 1 + 0*c_ipsi.^0;
fr_mos_ipsi = @(c_ipsi) 1 + 0.33*c_ipsi.^0.56;

% %Interpolated firing rates
% fr = load("C:\Users\Peter\Documents\MATLAB\inactPaper\widefield\mechanistic_fit\visMOs_allConds.mat");
% vis = mean( mean( fr.vis(:,:,:,0.05 < fr.binsS & fr.binsS < 0.1), 4), 3);
% m2 = mean( mean( fr.mos(:,:,:,0.075 < fr.binsS & fr.binsS < 0.125), 4), 3);
% fr_vis_contra = @(c_contra) interp1(fr.ucr, vis(:,1), c_contra);
% fr_mos_contra = @(c_contra) interp1(fr.ucr, m2(:,1), c_contra);

%1) only contra-lateral stimulus has firing
%Generate firing rate
D1.firing_rate = nan(length(D1.response), 4);
%Left VIS
D1.firing_rate(:,1) = fr_vis_contra(D1.stimulus(:,2)) .* ~(D1.laserRegion=='LeftVIS');
%Right VIS
D1.firing_rate(:,2) = fr_vis_contra(D1.stimulus(:,1)) .* ~(D1.laserRegion=='RightVIS');
%Left M2
D1.firing_rate(:,3) = fr_mos_contra(D1.stimulus(:,2)) .* ~(D1.laserRegion=='LeftM2');
%Right M2
D1.firing_rate(:,4) = fr_mos_contra(D1.stimulus(:,1)) .* ~(D1.laserRegion=='RightM2');

dat = struct('contrastLeft', D1.stimulus(:,1),...
              'contrastRight', D1.stimulus(:,2),...
              'choice', D1.response,...
              'sessionID', D1.sessionID,...
              'subjectID', D1.subjectID,...
              'firing_rate',D1.firing_rate);
dat.perturbation = zeros(size(dat.choice));
for p = 1:length(perturbationRegion)
    dat.perturbation(D1.laserRegion == perturbationRegion{p}) = p;
end

fit = bfit.fitModel('Two-Level-Perturbation',dat);
fit=load("C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\Two-level_removed_pre_stim_movement_trials.mat");
%% Multi-power: Plot phenomenological model
% areaCols = [ 0 0 1;
%              0 0 0.5;
%              0 1 0;
%              0 0.5 0;
%              1 0 0;
%              0.5 0 0;
%              1 1 0;
%              0.5 0.5 0];

areaCols = [ 73 148 208;
                141 179 206
                119 172 66;
             147 170 119]/255;


figure(401);

ha_B=subplot(2,4,7); %Bias perturbation 
ha_S=subplot(2,4,8); %Bias perturbation 
hold(ha_B,'on'); hold(ha_S,'on');

for p = 1:4
%     plot(fit.posterior.delta(:,1,p), fit.posterior.delta(:,2,p), 'k.');
    
    mu = mean([fit.posterior.delta(:,1,p), fit.posterior.delta(:,2,p)],1);
    Sigma = cov([fit.posterior.delta(:,1,p), fit.posterior.delta(:,2,p)]);
    x1 = (mu(1) - 10*max(diag(Sigma))):0.01:(mu(1) + 10*max(diag(Sigma)));
    x2 = (mu(2) - 10*max(diag(Sigma))):0.01:(mu(2) + 10*max(diag(Sigma)));
    [X1,X2] = meshgrid(x1,x2);
    F = mvnpdf([X1(:) X2(:)],mu,Sigma);
    F = reshape(F,length(x2),length(x1));
    [c1,c2]=contour(ha_B,x1,x2,F,8);
    c2.LineColor = areaCols(p,:);
    c2.LineWidth=1;
end
set(ha_B,'dataaspectratio',[1 1 1]);
line(ha_B,[0 0],ylim(ha_B)); line(ha_B,xlim(ha_B),[0 0]);
xlabel(ha_B,'\Delta b_L'); ylabel(ha_B,'\Delta b_R');


for p = 1:4
%     plot(fit.posterior.delta(:,3,p), fit.posterior.delta(:,4,p), 'k.');
    
    mu = mean([fit.posterior.delta(:,3,p), fit.posterior.delta(:,4,p)],1);
    Sigma = cov([fit.posterior.delta(:,3,p), fit.posterior.delta(:,4,p)]);
    x1 = (mu(1) - 10*max(diag(Sigma))):0.01:(mu(1) + 10*max(diag(Sigma)));
    x2 = (mu(2) - 10*max(diag(Sigma))):0.01:(mu(2) + 10*max(diag(Sigma)));
    [X1,X2] = meshgrid(x1,x2);
    F = mvnpdf([X1(:) X2(:)],mu,Sigma);
    F = reshape(F,length(x2),length(x1));
    [c1,c2]=contour(ha_S,x1,x2,F,8);
    c2.LineColor = areaCols(p,:);
    c2.LineWidth=1;
end
set(ha_S,'dataaspectratio',[1 1 1]);
line(ha_S,[0 0],ylim(ha_S)); line(ha_S,xlim(ha_S),[0 0]);
xlabel(ha_S,'\Delta s_L'); ylabel(ha_S,'\Delta s_R');


CL = [linspace(1,0.1,200), linspace(0.1,0,400), zeros(1,600)];
CR = [zeros(1,600), linspace(0,0.1,400) linspace(0.1,1,200)];
%Non-laser global fit on detection trials
BL = fit.posterior.bias(:,1) ;
BR = fit.posterior.bias(:,2);
SL = fit.posterior.sens(:,1);
SR = fit.posterior.sens(:,2);
N = fit.posterior.n_exp;
NL_ph=bplot.CN(BL,SL,BR,SR,N,CL,CR);
NL_phM=mean(NL_ph,1);
NL_phQ=quantile(NL_ph,[0.025 0.975],1);

%non-laser detection trial data
idx=min([fit.data.contrastLeft fit.data.contrastRight],[],2)==0;
cDiff = fit.data.contrastRight(idx) - fit.data.contrastLeft(idx);
resp = fit.data.choice(idx);
perturb = fit.data.perturbation(idx);
sess = fit.data.sessionID(idx);
[counts,~,~,lab] = crosstab(cDiff,resp,perturb,sess);
prob = counts./sum(counts,2);%Convert to probability over choices
prob = nanmean(prob,4); %Average over sessions
cDiffUnique=cellfun(@str2num,lab(1:size(prob,1),1));

ha = tight_subplot(4,3,0.01,0.05,[0.05 0.5]);
ii=1;
for region = [1 3 2 4]
% for region = 1:3
    %Global fit
    BL = fit.posterior.bias(:,1) + fit.posterior.delta(:,1,region);
    BR = fit.posterior.bias(:,2) + fit.posterior.delta(:,2,region);
    SL = fit.posterior.sens(:,1) + fit.posterior.delta(:,3,region);
    SR = fit.posterior.sens(:,2) + fit.posterior.delta(:,4,region);
    N = fit.posterior.n_exp;
    
    ph=bplot.CN(BL,SL,BR,SR,N,CL,CR);
    phM=mean(ph,1);
    phQ=quantile(ph,[0.025 0.975],1);
    
    phDiff = ph-NL_ph;
    phDiffM=mean(phDiff,1);
    phDiffQ=quantile(phDiff,[0.025 0.975],1);
    
    jj=1;
    for r = [1 3 2]
        hold( ha(3*(ii-1) + jj), 'on');
        fx = fill(ha(3*(ii-1) + jj),[CR-CL fliplr(CR-CL)], [NL_phQ(1,:,r) fliplr( NL_phQ(2,:,r) ) ], 'k');
        fx.FaceAlpha=0.3; fx.EdgeAlpha=0;
        plot(ha(3*(ii-1) + jj), CR-CL,NL_phM(1,:,r),'k-', 'linewidth', 1);
        
        fx = fill(ha(3*(ii-1) + jj),[CR-CL fliplr(CR-CL)], [phQ(1,:,r) fliplr( phQ(2,:,r) ) ], 'r');
        fx.FaceAlpha=0.3; fx.EdgeAlpha=0; fx.FaceColor = areaCols(region,:);
        plot(ha(3*(ii-1) + jj),CR-CL,phM(1,:,r),'-','Color',areaCols(region,:), 'linewidth', 1);
        
        plot(ha(3*(ii-1) + jj),cDiffUnique,prob(:,r,1),'k.','markersize',20);
        plot(ha(3*(ii-1) + jj),cDiffUnique,prob(:,r,1+region),'.','markersize',20,'Color',areaCols(region,:));
        
        jj = jj + 1;
    end
    
    ii = ii + 1;
end
set([ha],'xcolor','none','ycolor','none');
set([ha(end-2)],'ycolor','k','xcolor','k','xticklabelmode','auto','yticklabelmode','auto')
set(ha,'ylim',[0 1],'xlim',[-1 1]);
title(ha(1),'pLeft');
title(ha(2),'pNoGo');
title(ha(3),'pRight');
set(ha,'dataaspectratio',[1 1 1]);

%% Multi-power: Calculate %correct on high contrsat conditions, averaged across sessions
pCorr = nan(max(D.sessionID),1);
for sess = 1:max(D.sessionID)
    idx = D.sessionID==sess & ((D.stimulus(:,1)==0.54 & D.stimulus(:,2)==0) | (D.stimulus(:,1)==0 & D.stimulus(:,2)==0.54));
    pCorr(sess) = mean(D.feedbackType(idx));
end
%% IIA plot?
%Get prob of choice for all possible conditions
[counts,~,~,labels] = crosstab(D.stimulus(:,1),...
    D.stimulus(:,2),...
    D.response,...
    D.laserType);

counts = counts(:,:,:,1); %only nonlaser trials
%convert to probability
prob = counts./sum(counts,3);

zl = log(prob(:,:,1)./prob(:,:,3));
zr = log(prob(:,:,2)./prob(:,:,3));
stim = cellfun(@(s) str2num(s),labels(:,1));

figure;
subplot(2,2,1);
plot(stim, zl, 'o-'); xlabel('CL'); ylabel('Log pL/pNG');
subplot(2,2,2);
plot(stim, zr', 'o-');xlabel('CR'); ylabel('Log pR/pNG');

subplot(2,2,3);
plot(stim, zl', 'o-'); xlabel('CR'); ylabel('Log pL/pNG');
hold on; plot(stim, mean(zl',2), 'ko--');
subplot(2,2,4);
plot(stim, zr, 'o-');xlabel('CL'); ylabel('Log pR/pNG');
hold on; plot(stim, mean(zr',1), 'ko--');


%% Plot old mechanistic model fit
% fit_mech = load("C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\mech_inactivation.mat");
fit_mech = load("C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\tpe9664113_3681_4f17_9c23_d4da9f50fabe.mat"); %Original mechanistic model fit
fit_mech = load("C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\tp054a2072_da39_432a_8c1b_8480a5f745a9.mat"); %Mechanistic model fit with different baselines for VIS and MOs

areaCols = [ 109 183 229;
                77 131 158;
                160 206 87;
             119 147 62]/255;
         
%Plot posterior
% bplot.plotPosterior(fit_mech.posterior)

CL = [linspace(1,0.1,200), linspace(0.1,0,400), zeros(1,600)];
CR = [zeros(1,600), linspace(0,0.1,400) linspace(0.1,1,200)];
firing_rate = nan(4, length(CL));
firing_rate(1,:) = fr_vis_contra(CR);
firing_rate(2,:) = fr_vis_contra(CL);
firing_rate(3,:) = fr_mos_contra(CR);
firing_rate(4,:) = fr_mos_contra(CL);

%Non-laser global fit on detection trials
NL_ph=bplot.MECH(fit_mech.posterior.w, firing_rate);
% NL_ph=bplot.MECH_NESTED(fit.posterior.b_G ,fit.posterior.b_LR, fit.posterior.w_G, fit.posterior.w_LR, firing_rate);

NL_phM=mean(NL_ph,1);
NL_phQ=quantile(NL_ph,[0.025 0.975],1);
%non-laser detection trial data
idx=min([fit_mech.data.contrastLeft fit_mech.data.contrastRight],[],2)==0;
cDiff = fit_mech.data.contrastRight(idx) - fit_mech.data.contrastLeft(idx);
resp = fit_mech.data.choice(idx);
perturb = fit_mech.data.perturbation(idx);
sess = fit_mech.data.sessionID(idx);
[counts,~,~,lab] = crosstab(cDiff,resp,perturb,sess);
prob = counts./sum(counts,2);%Convert to probability over choices
prob = nanmean(prob,4); %Average over sessions
cDiffUnique=cellfun(@str2num,lab(1:size(prob,1),1));

figure;
ha = tight_subplot(4,3,0.01,0.05,[0.05 0.5]);
for region = 1:4
    
    inact_firing_rate = firing_rate;
    inact_firing_rate(region,:) = 0;
    ph=bplot.MECH(fit_mech.posterior.w, inact_firing_rate);
%     ph=bplot.MECH_NESTED(fit.posterior.b_G ,fit.posterior.b_LR, fit.posterior.w_G, fit.posterior.w_LR, inact_firing_rate);

    phM=mean(ph,1);
    phQ=quantile(ph,[0.025 0.975],1);

    jj=1;
    for r = [1 3 2]
        hold( ha(3*(region-1) + jj), 'on');
        fx = fill(ha(3*(region-1) + jj),[CR-CL fliplr(CR-CL)], [NL_phQ(1,:,r) fliplr( NL_phQ(2,:,r) ) ], 'k');
        fx.FaceAlpha=0.3; fx.EdgeAlpha=0;
        plot(ha(3*(region-1) + jj), CR-CL,NL_phM(1,:,r),'k-', 'linewidth', 1);
        
        fx = fill(ha(3*(region-1) + jj),[CR-CL fliplr(CR-CL)], [phQ(1,:,r) fliplr( phQ(2,:,r) ) ], 'r');
        fx.FaceAlpha=0.3; fx.EdgeAlpha=0; fx.FaceColor = areaCols(region,:);
        plot(ha(3*(region-1) + jj),CR-CL,phM(1,:,r),'-','Color',areaCols(region,:), 'linewidth', 1);
%         
        plot(ha(3*(region-1) + jj),cDiffUnique,prob(:,r,1),'k.','markersize',20);
        plot(ha(3*(region-1) + jj),cDiffUnique,prob(:,r,1+region),'.','markersize',20,'Color',areaCols(region,:));
%         
        jj = jj + 1;
    end
    
end
set([ha],'xcolor','none','ycolor','none');
set([ha(end-2)],'ycolor','k','xcolor','k','xticklabelmode','auto','yticklabelmode','auto')
set(ha,'ylim',[0 1],'xlim',[-1 1]);
title(ha(1),'pLeft');
title(ha(2),'pNoGo');
title(ha(3),'pRight');
set(ha,'dataaspectratio',[1 1 1]);

ha = tight_subplot(1,1,0.01,[0.5 0.05],[0.55 0.05]);
hold(ha,'on');
for region = 1:4
    WL = fit_mech.posterior.w(:,2+2*(region-1)+1 );
    WR = fit_mech.posterior.w(:,2+2*(region-1)+2 );
    mu = mean([WL, WR],1);
    Sigma = cov([WL, WR]);
    x1 = (mu(1) - 100*max(diag(Sigma))):0.01:(mu(1) + 100*max(diag(Sigma)));
    x2 = (mu(2) - 100*max(diag(Sigma))):0.01:(mu(2) + 100*max(diag(Sigma)));
    [X1,X2] = meshgrid(x1,x2);
    F = mvnpdf([X1(:) X2(:)],mu,Sigma);
    F = reshape(F,length(x2),length(x1));
    [c1,c2]=contour(ha,x1,x2,F,8);
    c2.LineColor = areaCols(region,:);
    c2.LineWidth=1;
end
set(ha,'dataaspectratio',[1 1 1]);
line(ha,[0 0],ylim(ha)); line(ha,xlim(ha),[0 0]);
xlabel(ha,'w_L'); ylabel(ha,'w_R');
set(ha,'XTickLabelMode','auto','YTickLabelMode','auto');

%IIA plot
stimVals = unique(fit_mech.data.contrastLeft);
ha = tight_subplot(1,2,0.01,[0.05 0.54],[0.55 0.05]); 
hold(ha(1),'on'); hold(ha(2),'on');
%ZL(CL) with different CR
for cr = 1:length(stimVals)
    CL = [linspace(0,1,200)];
    CR = ones(1,200)*stimVals(cr);
    
    firing_rate = nan(4, length(CL));
    firing_rate(1,:) = fr_vis_contra(CR);
    firing_rate(2,:) = fr_vis_contra(CL);
    firing_rate(3,:) = fr_mos_contra(CR);
    firing_rate(4,:) = fr_mos_contra(CL);
    NL_ph=bplot.MECH(fit_mech.posterior.w, firing_rate);
    ZL = log(NL_ph(:,:,1)./NL_ph(:,:,3));
    NL_phM=mean(ZL,1);
    NL_phQ=quantile(ZL,[0.025 0.975],1);
    
    fx = fill(ha(1),[CL fliplr(CL)], [NL_phQ(1,:) fliplr( NL_phQ(2,:) ) ], 'k');
    fx.FaceAlpha=0.1; fx.EdgeAlpha=0;
    plot(ha(1),CL, NL_phM,'k-');  
    
    
    for i = 1:4
        firing_rate_inact = firing_rate;
        firing_rate_inact(i,:) = 0; %inactivate one region
        NL_ph=bplot.MECH(fit_mech.posterior.w, firing_rate_inact);
        ZL = log(NL_ph(:,:,1)./NL_ph(:,:,3));
        NL_phM=mean(ZL,1);
        NL_phQ=quantile(ZL,[0.025 0.975],1);
        plot(ha(1),CL, NL_phM,'-','Color',areaCols(i,:));  
    end
end

%ZR(CR) with different CL
for cl = 1:length(stimVals)
    CR = [linspace(0,1,200)];
    CL = ones(1,200)*stimVals(cl);
    
    firing_rate = nan(4, length(CR));
    firing_rate(1,:) = fr_vis_contra(CR);
    firing_rate(2,:) = fr_vis_contra(CL);
    firing_rate(3,:) = fr_mos_contra(CR);
    firing_rate(4,:) = fr_mos_contra(CL);
    NL_ph=bplot.MECH(fit_mech.posterior.w, firing_rate);
    ZR = log(NL_ph(:,:,2)./NL_ph(:,:,3));
    NL_phM=mean(ZR,1);
    NL_phQ=quantile(ZR,[0.025 0.975],1);
    
    fx = fill(ha(2),[CR fliplr(CR)], [NL_phQ(1,:) fliplr( NL_phQ(2,:) ) ], 'k');
    fx.FaceAlpha=0.1; fx.EdgeAlpha=0;
    plot(ha(2),CR, NL_phM,'k-');  
    
    for i = 1:4
        firing_rate_inact = firing_rate;
        firing_rate_inact(i,:) = 0; %inactivate one region
        NL_ph=bplot.MECH(fit_mech.posterior.w, firing_rate_inact);
        ZR = log(NL_ph(:,:,2)./NL_ph(:,:,3));
        NL_phM=mean(ZR,1);
        NL_phQ=quantile(ZR,[0.025 0.975],1);
        plot(ha(2),CR, NL_phM,'-','Color',areaCols(i,:));
    end
    
end

xlabel(ha(1),'CL'); ylabel(ha(1),'ZL');
xlabel(ha(2),'CR'); ylabel(ha(2),'ZR');
linkaxes(ha,'xy');

%Plot different subjects posterior
figure;
ha = tight_subplot(1,1,0.01,0.1,0.1);
hold(ha,'on');
for subj = 1:fit_mech.data.numSubjects
    for region = 1:4
        WL = fit_mech.posterior.w(:,2+2*(region-1)+1 ) + fit_mech.posterior.b_subj(:,2+2*(region-1)+1,subj);;
        WR = fit_mech.posterior.w(:,2+2*(region-1)+2 ) + fit_mech.posterior.b_subj(:,2+2*(region-1)+2,subj);;

        mu = mean([WL, WR],1);
        Sigma = cov([WL, WR]);
        x1 = (mu(1) - 100*max(diag(Sigma))):0.01:(mu(1) + 100*max(diag(Sigma)));
        x2 = (mu(2) - 100*max(diag(Sigma))):0.01:(mu(2) + 100*max(diag(Sigma)));
        [X1,X2] = meshgrid(x1,x2);
        F = mvnpdf([X1(:) X2(:)],mu,Sigma);
        F = reshape(F,length(x2),length(x1));
        [c1,c2]=contour(ha,x1,x2,F,8);
        c2.LineColor = areaCols(region,:);
        c2.LineWidth=1;
        c2.LevelList = 2;
    end
end
set(ha,'dataaspectratio',[1 1 1]);
line(ha,[0 0],ylim(ha)); line(ha,xlim(ha),[0 0]);
xlabel(ha,'w_L'); ylabel(ha,'w_R');
set(ha,'XTickLabelMode','auto','YTickLabelMode','auto');

%%%%%%%%%%%%%% full pedestal plots %%%%%%%%%%%%%%

%Transform to pedestal/diff representation
fit_mech.data.pedestal = min(fit_mech.data.contrastLeft,fit_mech.data.contrastRight);
fit_mech.data.cDiff = fit_mech.data.contrastRight - fit_mech.data.contrastLeft;

[counts,~,~,labels] = crosstab(fit_mech.data.cDiff,...
    fit_mech.data.pedestal,...
    fit_mech.data.choice,...
    fit_mech.data.perturbation,...
    fit_mech.data.sessionID);
prob = counts./sum(counts,3);%Convert to probability over choices
prob_ave = nanmean(prob,5);
cdiffs = cellfun(@str2num,labels(1:size(prob,1),1));
peds = cellfun(@str2num,labels(1:size(prob,2),2));


f=figure;
ha_background = tight_subplot( 2, 2 , [0.04 0.01], [0.01 0.04], 0.01 );
for region = 1:4

    %For each session, add new subplots
    set(ha_background(region),'units','normalized');
    bounds = get(ha_background(region),'position');
    set(ha_background(region),'xcolor','none','ycolor','none');
    title(ha_background(region),region);
    %replace with many subplots within the bounds of the original
    ha = tight_subplot(4,3,0.001, [bounds(2) 1-bounds(2)-bounds(4)], [bounds(1) 1-bounds(1)-bounds(3)]);
    set(ha(1:9),'xcolor','none');
    set(ha([2 3 5 6 8 9 11 12]),'ycolor','none');
    for i = 1:length(ha)
        hold(ha(i),'on');
    end
    
    for ped = 1:length(peds)
        CL = [linspace(1-peds(ped),0,100) zeros(1,100)] + peds(ped);
        CR = [zeros(1,100) linspace(0,1-peds(ped),100)] + peds(ped);
        firing_rate = nan(4, length(CL));
        firing_rate(1,:) = fr_vis_contra(CR);
        firing_rate(2,:) = fr_vis_contra(CL);
        firing_rate(3,:) = fr_mos_contra(CR);
        firing_rate(4,:) = fr_mos_contra(CL);
        
        ph_NL=bplot.MECH(fit_mech.posterior.w, firing_rate);
        phM_NL=mean(ph_NL,1);
        phQ_NL=quantile(ph_NL,[0.025 0.975],1);
        
        firing_rate(region,:) = 0; %inactivate one region
        ph=bplot.MECH(fit_mech.posterior.w, firing_rate);
        phM=mean(ph,1);
        phQ=quantile(ph,[0.025 0.975],1);
        
        jj=1;
        for r = [1 3 2]
            %plot datapoints
            
            plot(ha( 3*(ped-1) + jj ),cdiffs,prob_ave(:,ped,r,1),'k.','markersize',10);
            fx = fill(ha( 3*(ped-1) + jj ),[CR-CL fliplr(CR-CL)], [phQ_NL(1,:,r) fliplr( phQ_NL(2,:,r) ) ], 'k');
            fx.FaceAlpha=0.3;
            fx.EdgeAlpha=0;
            plot(ha( 3*(ped-1) + jj ),CR-CL, phM_NL(1,:,r), 'k-');
            
            
            plot(ha( 3*(ped-1) + jj ),cdiffs,prob_ave(:,ped,r,1+region),'.','markersize',10,'color',areaCols(region,:));
            fx = fill(ha( 3*(ped-1) + jj ),[CR-CL fliplr(CR-CL)], [phQ(1,:,r) fliplr( phQ(2,:,r) ) ], 'k');
            fx.FaceAlpha=0.3;
            fx.EdgeAlpha=0;
            fx.FaceColor = areaCols(region,:);
            plot(ha( 3*(ped-1) + jj ),CR-CL, phM(1,:,r), '-','color',areaCols(region,:));
            jj = jj + 1 ;
        end
        
    end
end
set(get(f,'children'),'xlim',[-1 1]*1,'ylim',[0 1],'dataaspectratio',[1 1 1]);

%Plot matrices
[counts,~,~,lab]=crosstab(fit_mech.data.contrastLeft,...
                          fit_mech.data.contrastRight,...
                          fit_mech.data.choice,...
                          fit_mech.data.perturbation,...
                          fit_mech.data.sessionID);
prob = counts ./ sum(counts,3);
prob = nanmean(prob,5);
dProb = prob(:,:,:,2:end) - prob(:,:,:,1);

%get model predictoin
cVals = [0, 0.1, 0.24 0.54];
probModel = nan(4,4,3,5);
for cl = 1:length(cVals)
    for cr = 1:length(cVals)
        firing_rate = nan(4, 1);
            firing_rate(1,:) = fr_vis_contra(cVals(cr));
            firing_rate(2,:) = fr_vis_contra(cVals(cl));
            firing_rate(3,:) = fr_mos_contra(cVals(cr));
            firing_rate(4,:) = fr_mos_contra(cVals(cl));
            ph_NL=bplot.MECH(fit_mech.posterior.w, firing_rate);
            probModel(cl,cr,:,1) = squeeze(mean(ph_NL,1));
            
        for region = 1:4
            firing_rate_inact = firing_rate;
            firing_rate_inact(region) = 0; %inactivate one region
            ph=bplot.MECH(fit_mech.posterior.w, firing_rate_inact);
            probModel(cl,cr,:,1+region) = squeeze(mean(ph,1));
        end
    end
end
dProbModel = probModel(:,:,:,2:end) - probModel(:,:,:,1);

figure;
for r=1:3
    subplot(2,3,r);
    imagesc(cVals,cVals,probModel(:,:,r,1)); xlabel('CR'); ylabel('CL')
    
    subplot(2,3,3+r);
    imagesc(cVals,cVals,prob(:,:,r,1)); ylabel('Data');
end
set(get(gcf,'children'),'ydir','normal','dataaspectratio',[1 1 1],'clim',[0 1]);

%Change in probability as you inactivate
for region = 1:4
    figure;
    for r=1:3
        subplot(2,3,r);
        imagesc(cVals,cVals,dProbModel(:,:,r,region)); xlabel('CR'); ylabel('CL')
        
        subplot(2,3,3+r);
        imagesc(cVals,cVals,dProb(:,:,r,region)); ylabel('Data');
    end

    set(get(gcf,'children'),'ydir','normal','dataaspectratio',[1 1 1],'clim',[-1 1]);
cmap = [ones(100,1) linspace(0,1,100)' linspace(0,1,100)';
    linspace(1,0,100)' linspace(1,0,100)' ones(100,1)];
    colormap(flipud(cmap));
end

%% IIA PLot

[counts,~,~,lab] = crosstab(fit_mech.data.contrastLeft,...
                            fit_mech.data.contrastRight,...
                            fit_mech.data.choice,...
                            fit_mech.data.perturbation,...
                            fit_mech.data.sessionID);
prob = counts./sum(counts,3);%Convert to probability over choices
prob = nanmean(prob,5); %Average over sessions


zl = log(prob(:,:,1,:)./prob(:,:,3,:));
zr = log(prob(:,:,2,:)./prob(:,:,3,:));
stim = cellfun(@(s) str2num(s),lab(1:4,1));

figure;
subplot(1,2,1);
plot(stim, zl(:,:,1,1), 'o-'); xlabel('CL'); ylabel('Log pL/pNG');
subplot(1,2,2);
plot(stim, zr(:,:,1,1)', 'o-');xlabel('CR'); ylabel('Log pR/pNG');

figure;
regionLabels = {'LeftVIS','RightVIS','LeftM2','RightM2'};
for region = 1:4
    subplot(4,2,2*(region-1)+1);
    plot(stim, zl(:,:,1,1+region), 'o-'); xlabel('CL'); ylabel('Log pL/pNG');
    title(regionLabels{region});
    subplot(4,2,2*(region-1)+2)
    plot(stim, zr(:,:,1,1+region)', 'o-');xlabel('CR'); ylabel('Log pR/pNG');
    title(regionLabels{region});
end
set(get(gcf,'children')','ylim',[-1 3])