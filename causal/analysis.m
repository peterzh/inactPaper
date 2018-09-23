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
    dd = loadData(expRefs{session});
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

fit = load("C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\sparse_unilateral.mat");

%% Plot map on CL=CR condition
D = getrow(D, D.stimulus(:,1) == D.stimulus(:,2) & D.stimulus(:,1)>0);

% Compute the session-averaged change in choices from inactiavtion
[counts,~,~,labels] = crosstab(D.response,...
    D.laserIdx,...
    D.sessionID);
prob = counts./sum(counts,1);%Convert to probability over choices
%Compute delta from the non-laser condition
deltaProb = prob(:,2:end,:) - prob(:,1,:);
deltaProb = nanmean(deltaProb,3);

% Compute the session-averaged change in RT from inactiavtion
% deltaRT = nan(2,52,max(D.sessionID));
% for r = 1:2
%     for laser = 1:52
%         for sess = 1:91
%             deltaRT(r,laser,sess) = median( D.RT(D.response==r & D.laserIdx==laser & D.sessionID==sess) ) - ...
%                                 median( D.RT(D.response==r & D.laserIdx==0 & D.sessionID==sess) );
%         end
%     end
% end
% deltaRT = nanmean(deltaRT,3);

%Compute null distribution based on shuffling laser identities within
%each session
numIter = 10000;
shufDeltaProb = nan(3,52,numIter);
shufDeltaRT = nan(2,52,numIter);
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
    
    
%     dRT = nan(2,52,max(D.sessionID));
%     for r = 1:2
%         for laser = 1:52
%             for sess = 1:max(D.sessionID)
%                 dRT(r,laser,sess) = median( D.RT(D.response==r & D.laserIdx(idx)==laser & D.sessionID==sess) ) - ...
%                     median( D.RT(D.response==r & D.laserIdx(idx)==0 & D.sessionID==sess) );
%             end
%         end
%     end
%     dRT = nanmean(dRT,3);
%     shufDeltaRT(:,:,i) = dRT;
end
f.delete;

%Plot change in the choices
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

%Plot change in RT
q1RT = quantile(shufDeltaRT,[0.025 0.975],3);
q2RT = quantile(shufDeltaRT,[0.005 0.995],3);
q3RT = quantile(shufDeltaRT,[0.0005 0.9995],3);
figure('color','w');
for r = 1:3
    subplot(1,3,r); hold on;
    for q = 1:numel(coords) % coords is from ctxOutlines.mat
        cx = coords(q).x/100 - 5.7;
        cy = coords(q).y/100 - 5.4;
        plot(cx,-cy, 'LineWidth', 0.5, 'Color', [1 1 1]*0.5, 'Tag', 'outline');
    end
    
    %plot NS points as small dots
    idx = q1RT(r,:,1) < deltaRT(r,:) & deltaRT(r,:) < q1RT(r,:,2);
    h=scatter(coordSet(idx,1),coordSet(idx,2),20,'k','o','filled'); axis equal;  drawnow;
    h.MarkerEdgeColor=[1 1 1]*0.75;
    h.CData = deltaRT(r,idx);
    
    %Plot *
    idx = deltaRT(r,:) < q1RT(r,:,1) |  q1RT(r,:,2) < deltaRT(r,:);
    h=scatter(coordSet(idx,1),coordSet(idx,2),80,'k','o','filled'); axis equal;  drawnow;
    h.MarkerEdgeColor=[1 1 1]*0.75;
    h.CData = deltaRT(r,idx);
    
    %Plot **
    idx = deltaRT(r,:) < q2RT(r,:,1) |  q2RT(r,:,2) < deltaRT(r,:);
    h=scatter(coordSet(idx,1),coordSet(idx,2),150,'k','o','filled'); axis equal;  drawnow;
    h.MarkerEdgeColor=[1 1 1]*0.75;
    h.CData = deltaRT(r,idx);
    
    %Plot ***
    idx = deltaRT(r,:) < q3RT(r,:,1) |  q3RT(r,:,2) < deltaRT(r,:);
    h=scatter(coordSet(idx,1),coordSet(idx,2),300,'k','o','filled'); axis equal;  drawnow;
    h.MarkerEdgeColor=[1 1 1]*0.75;
    h.CData = deltaRT(r,idx);
    
    caxis([-1 1]*0.8);
    set(gca,'xtick','','ytick','','xcolor','w','ycolor','w','xlim',[-1 1]*5.2,'ylim',[-5 4]);
end
cmap = [ linspace(0,1,100)' linspace(0,1,100)' ones(100,1);
    ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
colormap(cmap);

%%%%%%%%%%%%%%%%%%%%%REGION TESTING: shuffle
perturbationRegion = {'LeftVIS','RightVIS','LeftM2','RightM2','LeftS1','RightS1'};
D1 = getrow(D, any(D.laserRegion == perturbationRegion,2) | D.laserType==0);
D1.perturbation = zeros(size(D1.response));
for p = 1:length(perturbationRegion)
    D1.perturbation(D1.laserRegion == perturbationRegion{p}) = p;
end

numIter = 10000;
shufdContra = nan(numIter,3);
shufdIpsi = nan(numIter,3);
for i = 1:numIter
    idx = 1:length(D1.laserIdx);
    for sess = 1:max(D1.sessionID)
        sessIdx = idx(D1.sessionID==sess);
        sessIdx = sessIdx(randperm(length(sessIdx)));
        idx(D1.sessionID==sess) = sessIdx;
    end
    counts = crosstab(D1.response,...
        D1.perturbation( idx ));
    prob = counts./sum(counts,1);
    dProb = prob(:,2:end) - prob(:,1);
    shufdContra(i,:) = mean([dProb(2,1:2:end); dProb(1,2:2:end)],1);
    shufdIpsi(i,:) = mean([dProb(1,1:2:end); dProb(2,2:2:end)],1);
end
counts = crosstab(D1.response,...
                D1.perturbation);
prob = counts./sum(counts,1);
dProb = prob(:,2:end) - prob(:,1);
dContra = mean([dProb(2,1:2:end); dProb(1,2:2:end)],1); %average across hemispheres
dIpsi = mean([dProb(1,1:2:end); dProb(2,2:2:end)],1); %average across hemispheres

regions = {'VIS','MOs','SSp'};
for region = 1:length(regions)
    
    %dContra
    null = shufdContra(:,region);
    P = min( mean(dContra(region) >= null), mean(dContra(region) <= null) );
    fprintf('Inactivate %s affects contra choices: effect=%0.5g , p=%0.5g\n',regions{region},dContra(region),P);
    
    %dContra
    null = shufdIpsi(:,region);
    P = min( mean(dIpsi(region) >= null), mean(dIpsi(region) <= null) );
    fprintf('Inactivate %s affects ipsi choices: effect=%0.5g , p=%0.5g\n',regions{region},dIpsi(region),P);
    
end

%% Get pulse experiment data
clear all;

expRefs = readtable('./sessionList_pulse.csv','FileType','text','Delimiter',',');
expRefs = expRefs.expRef;
mouseName = cellfun(@(e)e{3},cellfun(@(s) strsplit(s,'_'), expRefs,'uni',0),'uni',0);
[mice,~,subjID] = unique(mouseName);

D = struct;
ii=1;
for session = 1:length(expRefs)
    dd = loadData(expRefs{session});
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
rt_ci = nan(2,1000,2);
for t = 1:length(tsteps)
%     %VIS
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
    dd = loadData(expRefs{session});
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

%% Plot multi-power model-free stuff over VIS, M2, M1
%%%%%%%%%%%%%%%%%%%%%REGION TESTING: chi2 test
perturbationRegion = {'LeftVIS','RightVIS','LeftM2','RightM2','LeftM1','RightM1'};
D1 = getrow(D, any(D.laserRegion == perturbationRegion,2) | D.laserType==0);
D1.perturbation = zeros(size(D1.response));
for p = 1:length(perturbationRegion)
    D1.perturbation(D1.laserRegion == perturbationRegion{p}) = p;
end

%%%%%%%%%%%%%%%%%%%%%REGION TESTING: shuffle
% numIter = 10000;
% shufdContra = nan(numIter,3);
% shufdIpsi = nan(numIter,3);
% for i = 1:numIter
%     idx = 1:length(D1.laserIdx);
%     for sess = 1:max(D1.sessionID)
%         sessIdx = idx(D1.sessionID==sess);
%         sessIdx = sessIdx(randperm(length(sessIdx)));
%         idx(D1.sessionID==sess) = sessIdx;
%     end
%     counts = crosstab(D1.contrastLeft>0,...
%         D1.contrastRight>0,...
%         D1.response,...
%         D1.perturbation( idx ));
%     prob = counts./sum(counts,3);
%     dProb = prob(:,2:end) - prob(:,1);
%     shufdContra(i,:) = mean([dProb(2,1:2:end); dProb(1,2:2:end)],1);
%     shufdIpsi(i,:) = mean([dProb(1,1:2:end); dProb(2,2:2:end)],1);
% end

% 1) Show that inactivation causes strong impairment in contralateral
% responses when there's only a stimulus on the contralateral side
counts = crosstab(D1.stimulus(:,1)>0,...
        D1.stimulus(:,2)>0,...
        D1.response,...
        D1.perturbation);
prob = counts./sum(counts,3);
dProb = prob(:,:,:,2:end) - prob(:,:,:,1);
    shufdContra(i,:) = mean([dProb(2,1:2:end); dProb(1,2:2:end)],1);
    shufdIpsi(i,:) = mean([dProb(1,1:2:end); dProb(2,2:2:end)],1);






perturbationRegion = {'LeftVIS','RightVIS','LeftM2','RightM2','LeftM1','RightM1'};
D1 = getrow(D, any(D.laserRegion == perturbationRegion,2) | D.laserType==0);
D1.perturbation = zeros(size(D1.response));
for p = 1:length(perturbationRegion)
    D1.perturbation(D1.laserRegion == perturbationRegion{p}) = p;
end
% 
% counts = crosstab(D1.stimulus(:,1)>0, D1.stimulus(:,2)>0, D1.response, D1.sessionID, D1.perturbation);
% prob = counts./sum(counts,3);

%1) Show that inactivation blinds the animal to the contralateral side
region = {'VIS','M2','M1'};
regionLabels = {'VIS','MOs','MOp'};
areaCols = [ 73 148 208;
             119 172 66;
             100 100 100]/255;

figure; subplot(1,2,1); hold on;
Lidx = D1.laserType==0 & D1.stimulus(:,1)==0;
Ridx = D1.laserType==0 & D1.stimulus(:,2)==0;
cho = [D1.response(Lidx)==2; D1.response(Ridx)==1];
[p, pci] = binofit(sum(cho),length(cho));
fill([0 4 4 0],pci([1 1 2 2]),'k', 'facealpha',0.3, 'edgealpha', 0);
line([0 4],[1 1]*p,'color',[0 0 0]);

for r = 1:3
    %pContra to inactivated side on Contra stim
    Lidx = (D1.laserRegion == ['Left' region{r}]) & D1.stimulus(:,1)==0 ;
    Ridx = (D1.laserRegion == ['Right' region{r}]) & D1.stimulus(:,2)==0 ;
    cho = [D1.response(Lidx)==2; D1.response(Ridx)==1];
    [p, pci] = binofit(sum(cho),length(cho));
    plot(r ,p,'k.','markersize',20,'color',areaCols(r,:));
    line([1 1]*r, [pci], 'color',areaCols(r,:));
end
set(gca,'xticklabel', regionLabels, 'xtick', 1:3);
ylabel({'Probability of choosing', 'contralateral to the inactivated side'});

%2) Show that VIS increases Ipsi choices even when there's no ipsi there
subplot(1,2,2); hold on;
Lidx = D1.laserType==0 & D1.stimulus(:,1)==0;
Ridx = D1.laserType==0 & D1.stimulus(:,2)==0;
cho = [D1.response(Lidx)==1; D1.response(Ridx)==2];
[p, pci] = binofit(sum(cho),length(cho));
fill([0 4 4 0],pci([1 1 2 2]),'k', 'facealpha',0.3, 'edgealpha', 0);
line([0 4],[1 1]*p,'color',[0 0 0]);

for r = 1:3
    %pContra to inactivated side on Contra stim
    Lidx = (D1.laserRegion == ['Left' region{r}]) & D1.stimulus(:,1)==0;
    Ridx = (D1.laserRegion == ['Right' region{r}]) & D1.stimulus(:,2)==0;
    cho = [D1.response(Lidx)==1; D1.response(Ridx)==2];
    [p, pci] = binofit(sum(cho),length(cho));
    plot(r,p,'k.','markersize',20,'color',areaCols(r,:));
    line([1 1]*r , [pci], 'color',areaCols(r,:));
end
set(gca,'xticklabel', regionLabels, 'xtick', 1:3);
ylabel({'Probability of choosing','ipsilateral to inactivated side','eventhough no ipsi stimulus there'});



%1) stats test for contra blindness
Lidx = D1.stimulus(:,1)==0;
Ridx = D1.stimulus(:,2)==0;


%stats need to average across sessions not pool over them
%Mark trials where visual stimulus is contralateral to the inactivated
%hemisphere
D1.contraStim = (contains(string(D1.laserRegion),'Left') & D1.stimulus(:,1)==0) | ...
                (contains(string(D1.laserRegion),'Right') & D1.stimulus(:,2)==0);


counts = crosstab(D1.stimulus(:,1)>0, D1.stimulus(:,2)>0, D1.response, D1.sessionID, D1.perturbation);
prob = counts./sum(counts,3);

%% Fit multi-level model, or load previous fit already
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

% fit = bfit.fitModel('Two-Level-Perturbation',dat);
fit_mech = bfit.fitModel('Two-Level-Mechanistic',dat);

% fit=load("C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\2018-06-21_1703_Two-Level-Perturbation.mat");
fit=load("C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\Two-level_removed_pre_stim_movement_trials.mat");
%% Plot multi-level model fits
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


%% Plot mechanistic model fit
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
