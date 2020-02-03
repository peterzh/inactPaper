%% 52-coordinate inactivation experiment
clear all; close all;
load('../data/inactivation/Inactivation_52Coord.mat','D');

%Load standard coordset
load('../data/inactivation/26CoordSet.mat','coordSet');
coordSet = [coordSet; -coordSet(:,1), coordSet(:,2)];

%Get only trials with equal left and right contrast
D1 = D(D.stimulus(:,1) == D.stimulus(:,2) & D.stimulus(:,1)>0,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Compute empirical change in choices across coordinates  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the session-averaged change in choices from inactiavtion
[counts,~,~,labels] = crosstab(D1.response,...
    D1.laserIdx,...
    D1.sessionID);
prob = counts./sum(counts,1);%Convert to probability over choices
deltaProb = prob(:,2:end,:) - prob(:,1,:); %Compute delta from the non-laser condition
deltaProb = nanmean(deltaProb,3); %Average over sessions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Compute sampling distribution of null: no change in choices %%%%%%%%%%%
%%%%%% Shuffling laser and non-laser trials within each session   %%%   %%%%%%
%%%%%% NOTE: not deterministic, due to random sampling %%%%%%%%%%%%%%%%%%
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
choiceLabels = {'Left Choice','Right Choice','NoGo'};
load('../data/ctxOutlines.mat','coords');
q1 = quantile(shufDeltaProb,[0.025 0.975],3);
q2 = quantile(shufDeltaProb,[0.005 0.995],3);
q3 = quantile(shufDeltaProb,[0.0005 0.9995],3);
figure('color','w','position',[308 368 1492 610]);
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
    h.CData = deltaProb(r,idx)';
    
    %Plot *
    idx = deltaProb(r,:) < q1(r,:,1) |  q1(r,:,2) < deltaProb(r,:);
    h=scatter(coordSet(idx,1),coordSet(idx,2),80,'k','o','filled'); axis equal;  drawnow;
    h.MarkerEdgeColor=[1 1 1]*0.75;
    h.CData = deltaProb(r,idx)';
    
    %Plot **
    idx = deltaProb(r,:) < q2(r,:,1) |  q2(r,:,2) < deltaProb(r,:);
    h=scatter(coordSet(idx,1),coordSet(idx,2),150,'k','o','filled'); axis equal;  drawnow;
    h.MarkerEdgeColor=[1 1 1]*0.75;
    h.CData = deltaProb(r,idx)';
    
    %Plot ***
    idx = deltaProb(r,:) < q3(r,:,1) |  q3(r,:,2) < deltaProb(r,:);
    h=scatter(coordSet(idx,1),coordSet(idx,2),300,'k','o','filled'); axis equal;  drawnow;
    h.MarkerEdgeColor=[1 1 1]*0.75;
    h.CData = deltaProb(r,idx)';
    
    caxis([-1 1]*0.8);
    set(gca,'xtick','','ytick','','xcolor','w','ycolor','w','xlim',[-1 1]*5.2,'ylim',[-5 4]);
    
    title(choiceLabels{r});
    colorbar;
end
cmap = [ linspace(0,1,100)' linspace(0,1,100)' ones(100,1);
    ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
colormap(cmap);

fig2Pdf(gcf, '../figures/Fig3_Map1.pdf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Plot map of contraversive and ipsiversive effects %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Contraversive: L hemi, change in R choices. R hemi, change in L choices
hemi = sign(coordSet(:,1)); %-1 Left, +1 Right
dContra = mean([deltaProb(2, hemi==-1); deltaProb(1, hemi==+1)],1);
dIpsi = mean([deltaProb(1, hemi==-1); deltaProb(2, hemi==+1)],1);
labels = {'contraversive','ipsiversive'};

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
    title(labels{r});
end
cmap = [ linspace(0,1,100)' linspace(0,1,100)' ones(100,1);
    ones(100,1) linspace(1,0,100)' linspace(1,0,100)'];
colormap(cmap);

fig2Pdf(gcf, '../figures/Fig3_Map2.pdf');

%% Pulse inactivation experiment
clear all; close all;
load('../data/inactivation/Inactivation_Pulse.mat','D');

%Sort by laser onset time
D.laserOnset = D.laserOnset + randn(length(D.laserOnset),1)/100000;
[~,sortIdx]=sort(D.laserOnset);
D = D(sortIdx,:);

%Laser trials are selected so that the laser is contralateral to the side
%with higher contrast.
D.CL_gt_CR = D.stimulus(:,1) > D.stimulus(:,2);
D.CR_gt_CL = D.stimulus(:,2) > D.stimulus(:,1);

%Get non-laser
NL =  D((D.CR_gt_CL | D.CL_gt_CR) & D.laserType==0, :);

VIS = D((D.laserRegion == 'LeftVIS' & D.CR_gt_CL) | (D.laserRegion == 'RightVIS' & D.CL_gt_CR),: );
M2 = D((D.laserRegion == 'LeftM2' & D.CR_gt_CL) | (D.laserRegion == 'RightM2' & D.CL_gt_CR) ,:);
M1 = D((D.laserRegion == 'LeftM1' & D.CR_gt_CL) | (D.laserRegion == 'RightM1' & D.CL_gt_CR) ,:);

window_sec = 0.1; %Box-car window width for averaging
alpha = 0.0001; %Significance threshold

figure('color','w','position',[727 610 843 340]);
subplot(1,2,1); hold on; line([-1 1]*0.3, [1 1]*mean(NL.feedbackType),'Color',[0 0 0])
subplot(1,2,2); hold on; line([-1 1]*0.3, [1 1]*median(NL.RT(NL.feedbackType==1)),'Color',[0 0 0])

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
        subplot(1,2,1);
        plot(tsteps(t), 0.8,'bo');
    end
    
    pVal = ranksum( NL.RT(NL.feedbackType==1), VIS.RT(idx & VIS.feedbackType==1), 'tail', 'left');
    rt_vals = VIS.RT(idx & VIS.feedbackType==1);
    rt(1,t) = median( rt_vals );
    rt_ci(1,t,:) = median( rt_vals ) + tinv([0.025  0.975],length(rt_vals)-1)*std(rt_vals)/sqrt(length(rt_vals));
    if pVal<alpha
        subplot(1,2,2);
        plot(tsteps(t), 0.31,'bo');
    end
    
    %M2
    idx = WithinRanges(M2.laserOnset, tsteps(t) + [-0.5 0.5]*window_sec);
    [counts,~,pVal] = crosstab( [NL.feedbackType; M2.feedbackType(idx==1)], [NL.laserType; M2.laserType(idx==1)]);
    [perf(2,t),perf_ci(2,t,:)] = binofit(counts(2,2),sum(counts(:,2)));
    if pVal<alpha
        subplot(1,2,1);
        plot(tsteps(t), 0.7,'ro');
    end

    
    pVal = ranksum( NL.RT(NL.feedbackType==1), M2.RT(idx & M2.feedbackType==1), 'tail', 'left' );
    rt_vals = M2.RT(idx & M2.feedbackType==1);
    rt(2,t) = median( rt_vals );
    rt_ci(2,t,:) = median( rt_vals ) + tinv([0.025  0.975],length(rt_vals)-1)*std(rt_vals)/sqrt(length(rt_vals));
    if pVal<alpha
        subplot(1,2,2);
        plot(tsteps(t), 0.3,'ro');
    end
end

subplot(1,2,1);
plot(tsteps, perf(1,:), 'b'); 
fill([tsteps fliplr(tsteps)],[perf_ci(1,:,1) fliplr(perf_ci(1,:,2))],'k','FaceColor','b', 'EdgeAlpha', 0, 'FaceAlpha', 0.3);
plot(tsteps, perf(2,:), 'r'); 
fill([tsteps fliplr(tsteps)],[perf_ci(2,:,1) fliplr(perf_ci(2,:,2))],'k','FaceColor','r', 'EdgeAlpha', 0, 'FaceAlpha', 0.3);
xlim([-1 1]*0.3); xlabel('Stim onset'); ylabel('Proportion correct');

subplot(1,2,2);
plot(tsteps,rt(1,:),'b');
fill([tsteps fliplr(tsteps)],[rt_ci(1,:,1) fliplr(rt_ci(1,:,2))],'k','FaceColor','b', 'EdgeAlpha', 0, 'FaceAlpha', 0.3);
plot(tsteps,rt(2,:),'r');
fill([tsteps fliplr(tsteps)],[rt_ci(2,:,1) fliplr(rt_ci(2,:,2))],'k','FaceColor','r', 'EdgeAlpha', 0, 'FaceAlpha', 0.3);
xlim([-1 1]*0.3); xlabel('Stim onset'); ylabel('Median RT');

fig2Pdf(gcf, '../figures/Fig3_Pulse.pdf');

%% Correlate inactivation with decoding
clear all; close all;