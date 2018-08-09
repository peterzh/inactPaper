%% Get sparse_unilateral data
clear all;
expRefs = readtable('./sessionList_sparse_unilateral.csv','FileType','text','Delimiter',',');
expRefs = expRefs.expRef;

%Load standard coordset
load('26CoordSet.mat','coordSet');
coordSet = [coordSet; -coordSet(:,1), coordSet(:,2)];

D = struct;
ii=1;
for session = 1:length(expRefs)
    dd = loadData(expRefs{session});
    dd = structfun(@(x)(x(6:(end-14),:)),dd,'uni',0); %trim first 5 trials and last 15
    dd.sessionID = ones(length(dd.response),1)*ii;
    
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

%% Plot map on CL=CR condition
D = getrow(D, D.stimulus(:,1) == D.stimulus(:,2) & D.stimulus(:,1)>0);

%Counts
[counts,~,~,labels] = crosstab(D.response,...
    D.laserIdx);
prob = counts./sum(counts,1);%Convert to probability over choices
%Compute delta from the non-laser condition
deltaProb = prob(:,2:end) - prob(:,1);

%RT
% medianRT = grpstats(D.RT,{D.laserIdx D.response}, 'median');
% deltaRT = medianRT(2:end) - medianRT(1);

figure(101);
load('ctxOutlines.mat','coords');
for r = 1:3
    subplot(1,4,r); hold on;
    for q = 1:numel(coords) % coords is from ctxOutlines.mat
        cx = coords(q).x/100 - 5.7;
        cy = coords(q).y/100 - 5.4;
        plot(cx,-cy, 'LineWidth', 0.5, 'Color', [1 1 1]*0.5, 'Tag', 'outline');
    end
    
    h=scatter(coordSet(:,1),coordSet(:,2),200,deltaProb(r,:),'o','filled'); axis equal;  drawnow;
    h.MarkerEdgeColor=[1 1 1]*0.75;
    %                         caxis([-1 1]*q);
    caxis([-1 1]*0.8);
    set(gca,'xtick','','ytick','','xcolor','w','ycolor','w','xlim',[-1 1]*5.2,'ylim',[-5 4]);
end
cmap = [ linspace(0,1,100)' linspace(0,1,100)' ones(100,1);
    ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
colormap(cmap);

% %RT plot
% subplot(1,4,4); hold on;
% for q = 1:numel(coords) % coords is from ctxOutlines.mat
%     cx = coords(q).x/100 - 5.7;
%     cy = coords(q).y/100 - 5.4;
%     plot(cx,-cy, 'LineWidth', 0.5, 'Color', [1 1 1]*0.5, 'Tag', 'outline');
% end
% h=scatter(coordSet(:,1),coordSet(:,2),200,deltaRT,'o','filled'); axis equal;  drawnow;
% h.MarkerEdgeColor=[1 1 1]*0.75;
% %                         caxis([-1 1]*q);
% caxis([-1 1]*0.1);
% set(gca,'xtick','','ytick','','xcolor','w','ycolor','w','xlim',[-1 1]*5.2,'ylim',[-5 4]);


%Now do the same thing but shuffle laser identities
shufDeltaProb = nan(3,52,5000);
for i = 1:5000
    counts = crosstab(D.response,...
        D.laserIdx(randperm(length(D.laserIdx))));
    prob = counts./sum(counts,1);
    shufDeltaProb(:,:,i) = prob(:,2:end) - prob(:,1);
end

q1 = quantile(shufDeltaProb,[0.025 0.975],3);
q2 = quantile(shufDeltaProb,[0.005 0.995],3);
q3 = quantile(shufDeltaProb,[0.0005 0.9995],3);

figure(201);
for r = 1:3
    subplot(1,3,r); hold on;
    for q = 1:numel(coords) % coords is from ctxOutlines.mat
        cx = coords(q).x/100 - 5.7;
        cy = coords(q).y/100 - 5.4;
        plot(cx,-cy, 'LineWidth', 0.5, 'Color', [1 1 1]*0.5, 'Tag', 'outline');
    end
    
    %plot NS points as small dots
    idx = q1(r,:,1) < deltaProb(r,:) & deltaProb(r,:) < q1(r,:,2);
    h=scatter(coordSet(idx,1),coordSet(idx,2),20,deltaProb(r,idx),'o','filled'); axis equal;  drawnow;
    h.MarkerEdgeColor=[1 1 1]*0.75;
    
    %Plot *
    idx = deltaProb(r,:) < q1(r,:,1) |  q1(r,:,2) < deltaProb(r,:);
    h=scatter(coordSet(idx,1),coordSet(idx,2),100,deltaProb(r,idx),'o','filled'); axis equal;  drawnow;
    h.MarkerEdgeColor=[1 1 1]*0.75;
    
    %Plot **
    idx = deltaProb(r,:) < q2(r,:,1) |  q2(r,:,2) < deltaProb(r,:);
    h=scatter(coordSet(idx,1),coordSet(idx,2),200,deltaProb(r,idx),'o','filled'); axis equal;  drawnow;
    h.MarkerEdgeColor=[1 1 1]*0.75;
    
    %Plot ***
    idx = deltaProb(r,:) < q3(r,:,1) |  q3(r,:,2) < deltaProb(r,:);
    h=scatter(coordSet(idx,1),coordSet(idx,2),300,deltaProb(r,idx),'o','filled'); axis equal;  drawnow;
    h.MarkerEdgeColor=[1 1 1]*0.75;
    
    caxis([-1 1]*0.8);
    set(gca,'xtick','','ytick','','xcolor','w','ycolor','w','xlim',[-1 1]*5.2,'ylim',[-5 4]);
end
cmap = [ linspace(0,1,100)' linspace(0,1,100)' ones(100,1);
    ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
colormap(cmap);

figure(301); set(301,'color','w')
subplot(1,2,1); hold on;
for q = 1:numel(coords) % coords is from ctxOutlines.mat
    cx = coords(q).x/100 - 5.7;
    cy = coords(q).y/100 - 5.4;
    plot(cx,-cy, 'LineWidth', 0.5, 'Color', [1 1 1]*0.5, 'Tag', 'outline');
end
r = 2;
%plot NS points as small dots
idx = q1(r,:,1) < deltaProb(r,:) & deltaProb(r,:) < q1(r,:,2);
h=scatter(coordSet(idx,1),coordSet(idx,2),20,deltaProb(r,idx),'o','filled'); axis equal;  drawnow;
h.MarkerEdgeColor=[1 1 1]*0.75;

%Plot *
idx = deltaProb(r,:) < q1(r,:,1) |  q1(r,:,2) < deltaProb(r,:);
h=scatter(coordSet(idx,1),coordSet(idx,2),100,deltaProb(r,idx),'o','filled'); axis equal;  drawnow;
h.MarkerEdgeColor=[1 1 1]*0.75;

%Plot **
idx = deltaProb(r,:) < q2(r,:,1) |  q2(r,:,2) < deltaProb(r,:);
h=scatter(coordSet(idx,1),coordSet(idx,2),200,deltaProb(r,idx),'o','filled'); axis equal;  drawnow;
h.MarkerEdgeColor=[1 1 1]*0.75;

%Plot ***
idx = deltaProb(r,:) < q3(r,:,1) |  q3(r,:,2) < deltaProb(r,:);
h=scatter(coordSet(idx,1),coordSet(idx,2),300,deltaProb(r,idx),'o','filled'); axis equal;  drawnow;
h.MarkerEdgeColor=[1 1 1]*0.75;


set(gca,'xtick','','ytick','','xcolor','w','ycolor','w','xlim',[-1 1]*5.2,'ylim',[-5 4]);
cmap = [ linspace(0,1,100)' linspace(0,1,100)' ones(100,1);
    ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
colormap(cmap);
colorbar; caxis([-1 1]*0.8);
title('delta R on CL=CR');

%% Get pulse experiment data
clear all;

expRefs = readtable('./sessionList_pulse.csv','FileType','text','Delimiter',',');
expRefs = expRefs.expRef;

D = struct;
ii=1;
for session = 1:length(expRefs)
    dd = loadData(expRefs{session});
    dd = structfun(@(x)(x(6:(end-14),:)),dd,'uni',0); %trim first 5 trials and last 15
    dd.sessionID = ones(length(dd.response),1)*ii;
    
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

%get non-laser RT and performance
NL =  getrow(D, (D.CR_gt_CL | D.CL_gt_CR) & D.laserType==0 );
VIS = getrow(D, (D.laserRegion == 'LeftVIS' & D.CR_gt_CL) | (D.laserRegion == 'RightVIS' & D.CL_gt_CR) );
M2 = getrow(D, (D.laserRegion == 'LeftM2' & D.CR_gt_CL) | (D.laserRegion == 'RightM2' & D.CL_gt_CR) );

%% Sliding window performance and RT

window_sec = 0.1;
alpha = 0.0001;

% M2.RT(M2.laserOnset>0.1)=0.4;

figure(301);
subplot(2,2,2); hold on; line([-1 1]*0.3, [1 1]*mean(NL.feedbackType),'Color',[0 0 0])
subplot(2,2,4); hold on; line([-1 1]*0.3, [1 1]*median(NL.RT(NL.feedbackType==1)),'Color',[0 0 0])


tsteps = linspace(-0.3,0.3,1000);
perf = nan(2,1000);
rt = nan(2,1000);
for t = 1:length(tsteps)
%     %VIS
    idx = WithinRanges(VIS.laserOnset, tsteps(t) + [-0.5 0.5]*window_sec);
    [counts,~,pVal] = crosstab( [NL.feedbackType; VIS.feedbackType(idx==1)], [NL.laserType; VIS.laserType(idx==1)]);
    prob = counts./sum(counts,1);
    perf(1,t) = prob(2,2);
    if pVal<alpha
        subplot(2,2,2);
        plot(tsteps(t), 0.8,'bo');
    end
    
    pVal = ranksum( NL.RT(NL.feedbackType==1), VIS.RT(idx & VIS.feedbackType==1), 'tail', 'left');
    rt(1,t) = median( VIS.RT(idx & VIS.feedbackType==1) );
    if pVal<alpha
        subplot(2,2,4);
        plot(tsteps(t), 0.31,'bo');
    end
    
    %M2
    idx = WithinRanges(M2.laserOnset, tsteps(t) + [-0.5 0.5]*window_sec);
    [counts,~,pVal] = crosstab( [NL.feedbackType; M2.feedbackType(idx==1)], [NL.laserType; M2.laserType(idx==1)]);
    prob = counts./sum(counts,1);
    perf(2,t) = prob(2,2);
    if pVal<alpha
        subplot(2,2,2);
        plot(tsteps(t), 0.7,'ro');
    end
    
    pVal = ranksum( NL.RT(NL.feedbackType==1), M2.RT(idx & M2.feedbackType==1), 'tail', 'left' );
    rt(2,t) = median( M2.RT(idx & M2.feedbackType==1) );
    if pVal<alpha
        subplot(2,2,4);
        plot(tsteps(t), 0.3,'ro');
    end
    
end

subplot(2,2,2);
plot(tsteps, perf); xlim([-1 1]*0.3); xlabel('Stim onset'); ylabel('Performance');

subplot(2,2,4);
plot(tsteps,rt); xlim([-1 1]*0.3); xlabel('Stim onset'); ylabel('Median RT');
