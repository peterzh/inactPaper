clear all; close all;
load('../data/Inactivation_HigherPower.mat','D');

%% Fit behavioural model to HigherPowerInact sessions
perturbationRegion = {'LeftVIS','RightVIS','LeftM2','RightM2'};
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

% fit = bfit.fitModel('Two-Level-Perturbation',dat);
% fit = load("C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\Two-level_removed_pre_stim_movement_trials.mat");
fit = load("C:\Users\Peter\Documents\MATLAB\stan2AFC\fits\tp150f8a73_eced_43f6_84c7_a5f338f2282c.mat");
%% Plot the psychometric curves for the non-laser trials only
colours = [ 0    0.4470    0.7410;
    0.8500    0.3250    0.0980; 
    0.4940    0.1840    0.5560];

%Grand average parameter estimates
BL = fit.posterior.bias(:,1) ;
BR = fit.posterior.bias(:,2) ;
SL = fit.posterior.sens(:,1) ;
SR = fit.posterior.sens(:,2) ;
N = fit.posterior.n_exp ;

%Compute empirical choice fractions for each contrast condition, choice,
%inactivation condition, and session ID
c2 = crosstab(fit.data.contrastLeft,...
    fit.data.contrastRight,...
    fit.data.choice,...
    fit.data.sessionID,...
    fit.data.perturbation);
prob = c2./sum(c2,3);%Convert to probability over choices
prob_ave = nanmean(prob,4); %average over sessions
prob_ave = prob_ave(:,:,:,:,1); %get only non-laser trials

%Plot 3D mesh psychometric curves
CLs = unique(fit.data.contrastLeft);
CRs = unique(fit.data.contrastRight);
DotSize = 32;
openDotSize = 9;

f = figure('color','w','position',[164 653 1299 255]);
ha = tight_subplot(1,4,0.1,0.1,0.1);
for i = 1:4; hold(ha(i),'on'); end;
j = 1;
for r = [1 3 2]
    CLGrid = repmat(CLs,1,4);
    CRGrid = repmat(CRs,1,4)';
    PGrid = prob_ave(:,:,r,1);
    
    plot3(ha(j),CLGrid, CRGrid, PGrid, 'k--', 'color', colours(r,:));
    plot3(ha(j),CLGrid', CRGrid', PGrid', 'k--', 'color', colours(r,:));
    xlabel(ha(j),'CL'); ylabel(ha(j),'CR'); grid(ha(j),'on');
    
    %Marker based on whether choice is correct
    if r == 1
        plot3(ha(j), CLGrid(CLGrid>CRGrid), CRGrid(CLGrid>CRGrid), PGrid(CLGrid(:)>CRGrid(:)), 'ko' ,'markersize', openDotSize,'MarkerFaceColor', [0 0 0] );
        plot3(ha(j), CLGrid(CRGrid>CLGrid), CRGrid(CRGrid>CLGrid), PGrid(CRGrid(:)>CLGrid(:)), 'ko' ,'markersize', openDotSize,'MarkerFaceColor', [1 1 1] );
        plot3(ha(j), CLGrid(CLGrid==CRGrid & CLGrid>0), CRGrid(CLGrid==CRGrid & CLGrid>0), PGrid(CLGrid(:)==CRGrid(:) & CLGrid(:)>0), 'ko' ,'markersize', openDotSize,'MarkerFaceColor', [1 1 1]*0.5 );
        plot3(ha(j), CLGrid(CLGrid==0 & CRGrid==0), CRGrid(CLGrid==0 & CRGrid==0), PGrid(CLGrid(:)==0 & CRGrid(:)==0), 'ko' ,'markersize', openDotSize,'MarkerFaceColor', [1 1 1] );

    elseif r == 2
        plot3(ha(j), CLGrid(CLGrid<CRGrid), CRGrid(CLGrid<CRGrid), PGrid(CLGrid(:)<CRGrid(:)),  'ko' ,'markersize', openDotSize,'MarkerFaceColor', [0 0 0] );
        plot3(ha(j), CLGrid(CRGrid<CLGrid), CRGrid(CRGrid<CLGrid), PGrid(CRGrid(:)<CLGrid(:)), 'ko' ,'markersize', openDotSize,'MarkerFaceColor', [1 1 1] );
        plot3(ha(j), CLGrid(CLGrid==CRGrid & CLGrid>0), CRGrid(CLGrid==CRGrid & CLGrid>0), PGrid(CLGrid(:)==CRGrid(:) & CRGrid(:)>0), 'ko' ,'markersize', openDotSize,'MarkerFaceColor', [1 1 1]*0.5 );
        plot3(ha(j), CLGrid(CLGrid==0 & CRGrid==0), CRGrid(CLGrid==0 & CRGrid==0), PGrid(CLGrid(:)==0 & CRGrid(:)==0), 'ko' ,'markersize', openDotSize,'MarkerFaceColor', [1 1 1] );

    elseif r == 3
        plot3(ha(j), CLGrid(CLGrid==0 & CRGrid==0), CRGrid(CLGrid==0 & CRGrid==0), PGrid(CLGrid(:)==0 & CRGrid(:)==0), 'ko' ,'markersize', openDotSize,'MarkerFaceColor', [0 0 0] );
        plot3(ha(j), CLGrid(CLGrid>0 | CRGrid>0), CRGrid(CLGrid>0 | CRGrid>0), PGrid(CRGrid(:)>0 | CLGrid(:)>0), 'ko' ,'markersize', openDotSize,'MarkerFaceColor', [1 1 1] );
    end

    j = j +1;
end

for cl = 1:4
    CL = ones(1,100)*CLs(cl);
    CR = linspace(0,0.54,100);
    ph=bplot.CN(BL,SL,BR,SR,N,CL,CR);
    ph = permute( mean(ph,1), [2 3 1] );
    j = 1;
    for r = [1 3 2]
        plot3(ha(j),CL,CR, ph(:,r), '-', 'color', colours(r,:));
        j = j +1;
    end
end


for cr = 1:4
    CR = ones(1,100)*CRs(cr);
    CL = linspace(0,0.54,100);
    ph=bplot.CN(BL,SL,BR,SR,N,CL,CR);
    ph = permute( mean(ph,1), [2 3 1] );
    j = 1;
    for r = [1 3 2]
        plot3(ha(j),CL,CR, ph(:,r), '-', 'color', colours(r,:));
        j = j +1;
    end
end
set(ha(1:3), 'xticklabelmode','auto',...
        'xtick', CLs,...
        'yticklabelmode','auto',...
        'ytick', CRs,...
        'xlim', [0 0.54],...
        'ylim', [0 0.54],...
        'zlim',[0 1],...
        'tickdir','out');
set(ha(1:3),'xdir','reverse'); 
set(ha(1),'view',[-160,45]);
set(ha(2),'view',[-135,45]);
set(ha(3),'view',[-110,45]);

title(ha(1),'Left');
title(ha(2),'NoGo');
title(ha(3),'Right');

%Plot pCorrect/pIncorrect/pMiss for single-sided contrast only
CL = [linspace(0.05,0.6,200), zeros(1,200)];
CR = [zeros(1,200), linspace(0.05,0.6,200)];
%Non-laser global fit on detection trials
NL_ph=bplot.CN(BL,SL,BR,SR,N,CL,CR);

pCorrect = mean(cat(4,NL_ph(:,CL>0,1),NL_ph(:,CR>0,2)),4);
q_pCorrect = quantile(pCorrect,[0.025 0.975],1);
plot(ha(4),CL(CL>0),mean(pCorrect,1),'-','color',[1 1 1]*0,'linewidth',2);
fx=fill(ha(4), [CL(CL>0) fliplr(CL(CL>0))],[q_pCorrect(1,:) fliplr(q_pCorrect(2,:))],'k'); fx.FaceAlpha=0.2; fx.EdgeAlpha=0;
fx.FaceColor=[1 1 1]*0;
plot(ha(4),CLs(CLs>0),mean(cat(2,prob_ave(2:end,1,1),prob_ave(1,2:end,2)'),2),'.','color',[1 1 1]*0,'markersize',30);

%pIncorrectOpposite
%(CL + R choice) & (CR + L choice)
pICO = mean(cat(4,NL_ph(:,CL>0,2),NL_ph(:,CR>0,1)),4);
q_pICO = quantile(pICO,[0.025 0.975],1);
plot(ha(4),CL(CL>0),mean(pICO,1),'-','color',[1 1 1]*0.8,'linewidth',2);
fx=fill(ha(4), [CL(CL>0) fliplr(CL(CL>0))],[q_pICO(1,:) fliplr(q_pICO(2,:))],'k'); fx.FaceAlpha=0.2; fx.EdgeAlpha=0;
fx.FaceColor=[1 1 1]*0.8;
plot(ha(4),CLs(CLs>0),mean(cat(2,prob_ave(2:end,1,2),prob_ave(1,2:end,1)'),2),'.','color',[1 1 1]*0.8,'markersize',30);

%pIncorrectNoGo
pING = mean(cat(4,NL_ph(:,CL>0,3),NL_ph(:,CR>0,3)),4);
q_pING = quantile(pING,[0.025 0.975],1);
plot(ha(4),CL(CL>0),mean(pING,1),'-','color',[1 1 1]*0.4,'linewidth',2);
fx=fill(ha(4), [CL(CL>0) fliplr(CL(CL>0))],[q_pING(1,:) fliplr(q_pING(2,:))],'k'); fx.FaceAlpha=0.2; fx.EdgeAlpha=0;
fx.FaceColor=[1 1 1]*0.4;
xlabel(ha(4),'Contrast (single side)');
plot(ha(4),CLs(CLs>0),mean(cat(2,prob_ave(2:end,1,3),prob_ave(1,2:end,3)'),2),'.','color',[1 1 1]*0.4,'markersize',30);

set(ha(4),'xlim',[0 0.65],'ylim',[0 1],'TickDir','out','dataaspectratio',[1 1 1],...
    'xticklabelmode','auto','xtick', CLs);