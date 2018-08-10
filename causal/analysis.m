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

%% Get multi-power inactivation data
clear all;
expRefs = readtable('./sessionList_unilateral_multiple_powers.csv','FileType','text','Delimiter',',');
expRefs = expRefs.expRef;
mouseName = cellfun(@(e)e{3},cellfun(@(s) strsplit(s,'_'), expRefs,'uni',0),'uni',0);
[names,~,subjID] = unique(mouseName);

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

%% Fit multi-level model, or load previous fit already
perturbationRegion = {'LeftVIS','RightVIS','LeftM2','RightM2','LeftM1','RightM1','LeftS1','RightS1'};
D1 = getrow(D, any(D.laserRegion == perturbationRegion,2) | D.laserType==0);
dat = struct('contrastLeft', D1.stimulus(:,1),...
              'contrastRight', D1.stimulus(:,2),...
              'choice', D1.response,...
              'sessionID', D1.sessionID,...
              'subjectID', D1.subjectID);
dat.perturbation = zeros(size(dat.choice));
for p = 1:length(perturbationRegion)
    dat.perturbation(D1.laserRegion == perturbationRegion{p}) = p;
end
fit = bfit.fitModel('Two-Level-Perturbation',dat);

% fit=load("C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\2018-06-21_1703_Two-Level-Perturbation.mat");
fit=load("C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\Two-level_removed_pre_stim_movement_trials.mat");

%% Plot multi-level 
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
set(ha,'dataaspectratio',[1.5 1 1]);