%% Plot weights of widefield-neurometric
fit_mech = load("../data/modelFits/neurometric_wf_original.mat"); %Original fit to widefield data

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
set(gca,'xlim',[-0.2 0.25],'ylim',[-0.2 0.25]);
xlabel('Weight onto Z_L'); ylabel('Weight onto Z_R'); title('Posterior dist of weights');

fig2Pdf(gcf, '../figures/Fig4_weights.pdf');

%% Plot prediction of inactivation for Correct, Incorrect & Miss rates
load('../data/inactivation/Inactivation_HigherPower.mat','D');
fit_inact = load("../data/modelFits/neurometric_inact_original.mat"); %Original fit to inactivation sessions
% fit_inact = load("C:\stanFitDump\tpaeecaac5_464c_4218_b5b9_ee96a3152748.mat");

areaCols = [ 109 183 229;
                77 131 158;
                160 206 87;
             119 147 62]/255;
         
         
         
%Get trials with inactivation in VIS and M2 only
perturbationRegion = {'LeftVIS','RightVIS','LeftM2','RightM2'};
D1 = D(any(D.laserRegion == perturbationRegion,2) | D.laserType==0,:);
D1.perturbation = zeros(size(D1.response));
for p = 1:length(perturbationRegion)
    D1.perturbation(D1.laserRegion == perturbationRegion{p}) = p;
end

%Calculate the proportion of Correct, Incorrect & Miss trials for
%single-sinded contrast
[counts,~,~,labels] = crosstab(D1.stimulus(:,1)>0,...
    D1.stimulus(:,2)>0,...
    D1.response,...
    D1.perturbation,...
    D1.sessionID);
prob = counts./sum(counts,3);%Convert to probability over choices

%pCorrect: (CL & Left choice & Right hemisphere inact) + (CR & right choice & Left hemisphere inact)
pCorrect = nanmean(cat(3,squeeze(prob(2,1,1,[1 3 5],:)),squeeze(prob(1,2,2,[1 2 4],:))),3);

%pIncorrect: (CL & Right choice & Right hemisphere inact) + (CR & Left choice & Left hemisphere inact)
pIncorrect = nanmean(cat(3,squeeze(prob(2,1,2,[1 3 5],:)),squeeze(prob(1,2,1,[1 2 4],:))),3);

%pMiss: (CL & NoGo & Right hemisphere inact) + (CR & NoGo & Left hemisphere inact)
pMiss = nanmean(cat(3,squeeze(prob(2,1,3,[1 3 5],:)),squeeze(prob(1,2,3,[1 2 4],:))),3);

figure('position',[198 259 561 545]); hold on;
errorbar(1:3, nanmean(pCorrect,2),nanstd(pCorrect,[],2)./sqrt(sum(~isnan(pCorrect),2)),'Marker','o','LineStyle','none','Markersize',8,'linewidth',2);
errorbar(4:6, nanmean(pIncorrect,2),nanstd(pIncorrect,[],2)./sqrt(sum(~isnan(pIncorrect),2)),'Marker','o','LineStyle','none','Markersize',8,'linewidth',2);
errorbar(7:9, nanmean(pMiss,2),nanstd(pMiss,[],2)./sqrt(sum(~isnan(pMiss),2)),'Marker','o','LineStyle','none','Markersize',8,'linewidth',2);
set(gca,'xtick',1:3,'xticklabel',{'Laser off','VIS','M2'},'xlim',[0 10],'ylim',[0 1]);
ylabel('Probability');

%Plot model fit for laser off trials, and also prediction for laser on trials
CL=mean(D1.stimulus(D1.stimulus(:,1)>0,1));
CR=0;

fr = load("../data/neuropixels/visMOs_allConds.mat");

VIS_window = [0.075 0.125];
MOS_window = [0.125 0.175];
%get ephys data at critical timewindow
vis = mean( mean( fr.vis(:,:,:,VIS_window(1) < fr.binsS & fr.binsS < VIS_window(2)), 4), 3);
m2 = mean( mean( fr.mos(:,:,:,MOS_window(1) < fr.binsS & fr.binsS < MOS_window(2)), 4), 3);

F=[];
F(1,:) = interp2(fr.ucl, fr.ucr, vis, CL , CR, 'linear' );
F(2,:) = interp2(fr.ucl, fr.ucr, vis', CL , CR, 'linear' );
F(3,:) = interp2(fr.ucl, fr.ucr, m2, CL , CR, 'linear' );
F(4,:) = interp2(fr.ucl, fr.ucr, m2', CL , CR, 'linear' );
NL_ph=stan_plotFcn_MECH_SYMETRICAL(fit_inact.posterior.w, F);
q_ph = quantile(NL_ph,[0.025 0.975],1);
mean_ph = mean(NL_ph,1);


%Laser off: 
    % pCorrect = (Left choice for Left contrast shown)
    % pIncorrect = (Right choice for Left contrast shown)
    % pMiss = (NoGo choice for Left contrast shown)
for r = 1:3
    line( 3*(r-1)+1 + [-1 1]*0.2, [1 1]* mean_ph(:,:,r),'linewidth',2,'color',[0 0 0]);
    fill(3*(r-1)+1 +[-0.2 0.2 0.2 -0.2],[q_ph(1,1,r) q_ph(1,1,r) q_ph(2,1,r) q_ph(2,1,r)],'k','facealpha',0.2);
end

%Laser on: Right VIS only since model is tested on left contrast.
F_inact = F; F_inact(2)=0;
ph=stan_plotFcn_MECH_SYMETRICAL(fit_inact.posterior.w, F_inact);
q_ph = quantile(ph,[0.025 0.975],1);
mean_ph = mean(ph,1);
for r = 1:3
    line( 3*(r-1)+2 + [-1 1]*0.2, [1 1]* mean_ph(:,:,r),'linewidth',2,'color',areaCols(1,:));
    fx=fill(3*(r-1)+2 +[-0.2 0.2 0.2 -0.2],[q_ph(1,1,r) q_ph(1,1,r) q_ph(2,1,r) q_ph(2,1,r)],'k','facealpha',0.2,'facecolor',areaCols(1,:));
end

%Laser on: Right MOs only since model is tested on left contrast.
F_inact = F; F_inact(4)=0;
ph=stan_plotFcn_MECH_SYMETRICAL(fit_inact.posterior.w, F_inact);
q_ph = quantile(ph,[0.025 0.975],1);
mean_ph = mean(ph,1);
for r = 1:3
    line( 3*(r-1)+3 + [-1 1]*0.2, [1 1]* mean_ph(:,:,r),'linewidth',2,'color',areaCols(3,:));
    fill(3*(r-1)+3 +[-0.2 0.2 0.2 -0.2],[q_ph(1,1,r) q_ph(1,1,r) q_ph(2,1,r) q_ph(2,1,r)],'k','facealpha',0.2,'facecolor',areaCols(3,:));
end
legend('Correct','Incorrect','Miss');
