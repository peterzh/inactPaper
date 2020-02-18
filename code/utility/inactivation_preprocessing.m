%% 52 coord dataset
clear all;
expRefs = readtable('../data/inactivation/sessionList_52Coord.csv','FileType','text','Delimiter',',');
expRefs = expRefs.expRef;
mouseName = cellfun(@(e)e{3},cellfun(@(s) strsplit(s,'_'), expRefs,'uni',0),'uni',0);
[mice,~,subjID] = unique(mouseName);

%Load standard coordset
load('../data/inactivation/26CoordSet.mat','coordSet');
coordSet = [coordSet; -coordSet(:,1), coordSet(:,2)];

D = table;
ii=1;
for session = 1:length(expRefs)
    [dd,meta] = loadBehaviouralData(expRefs{session});

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
    dd = structfun(@(f) f(find(keep),:), dd, 'uni', 0);
    
    %If there aren't any laser trials
    if ~any(dd.laserCoord(~isnan(dd.laserCoord(:,1)),1) ~= 0)
        %                         keyboard;
    end
    
    if any(dd.laserType>0)
        D = vertcat(D,struct2table(dd));
        ii=ii+1;
    end
end
D = D(D.repeatNum==1,:);

D.laserIdx = zeros(size(D.response));
%Define laser Idx
for i = 1:size(coordSet,1)
    id = D.laserCoord(:,1) == coordSet(i,1) & D.laserCoord(:,2) == coordSet(i,2);
    D.laserIdx(id) = i;
end

wheel_stimulusOn_timesteps = D.wheel_stimulusOn_timesteps(1,:);
D.wheel_stimulusOn_timesteps = [];
D.wheel_startMove_timesteps=[];
D.wheel_startMove=[];

%Save the data
save('../data/inactivation/Inactivation_52Coord.mat','D','wheel_stimulusOn_timesteps','-v7.3');

%% Pulse inactivation dataaset
clear all;
expRefs = readtable('../data/inactivation/sessionList_Pulse.csv','FileType','text','Delimiter',',');
expRefs = expRefs.expRef;
mouseName = cellfun(@(e)e{3},cellfun(@(s) strsplit(s,'_'), expRefs,'uni',0),'uni',0);
[mice,~,subjID] = unique(mouseName);

D = table;
ii=1;
for session = 1:length(expRefs)
    [dd,meta] = loadBehaviouralData(expRefs{session});

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
    
%     %add marker for whether the choice is contra or ipsi to the inactivated side is contra or ipsi to stimulus side
%     dd.ipsiversive = ( contains(string(dd.laserRegion),'Left') & (dd.response==1) ) | ...
%                ( contains(string(dd.laserRegion),'Right') & (dd.response==2) );

    if any(dd.laserType>0)
        D = vertcat(D,struct2table(dd));
        ii=ii+1;
    end
end
D = D(D.repeatNum==1,:);
% 
% %define stimulus conditions where one side is higher than the other
% D.CL_gt_CR = D.stimulus(:,1) > D.stimulus(:,2);
% D.CR_gt_CL = D.stimulus(:,2) > D.stimulus(:,1);
% 
% %prepare laser onset times for plotting later
% D.laserOnset = D.laserOnset + randn(length(D.laserOnset),1)/100000;
% [~,sortIdx]=sort(D.laserOnset);
% D = D(sortIdx,:);

%Remove trials with pre-stimulus wheel movement
keepIdx = zeros(size(D.response));
for sess = 1:max(D.sessionID)
    times = D.wheel_stimulusOn_timesteps(D.sessionID==sess,:);
    times = times(1,:);
    wheel = D.wheel_stimulusOn(D.sessionID==sess,:);
    choice = D.response(D.sessionID==sess);
    
    %Whiten
    wheel_whiten =  wheel/std(wheel(:,end));
    
    %Threshold
    wheel_whiten_window = mean( abs(wheel_whiten(:,-0.150 < times & times < 0.05 & times ~= 0)) , 2);
    keepIdx(D.sessionID==sess) = wheel_whiten_window < 0.05;
end
fprintf('Removing %0.2f%% of trials becuase of pre-stim wheel movements\n',100*mean(keepIdx==0));
D = D(keepIdx==1,:);

wheel_stimulusOn_timesteps = D.wheel_stimulusOn_timesteps(1,:);
D.wheel_stimulusOn_timesteps = [];
D.wheel_startMove_timesteps=[];
D.wheel_startMove=[];

%Save the data 
save('../data/inactivation/Inactivation_Pulse.mat','D','wheel_stimulusOn_timesteps','-v7.3');

%% Higher power inactivation dataset
clear all;
expRefs = readtable('../data/inactivation/sessionList_HigherPower.csv','FileType','text','Delimiter',',');
expRefs = expRefs.expRef;
mouseName = cellfun(@(e)e{3},cellfun(@(s) strsplit(s,'_'), expRefs,'uni',0),'uni',0);
[mice,~,subjID] = unique(mouseName);

D = table;
ii=1;
for session = 1:length(expRefs)
    [dd,meta] = loadBehaviouralData(expRefs{session});
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
        D = vertcat(D,struct2table(dd));
        ii=ii+1;
    end
end
D = D(D.repeatNum==1,:);

coordSet = unique(D.laserCoord(D.laserType==1,:),'rows');
D.laserIdx = zeros(size(D.response));
%Define laser Idx
for i = 1:size(coordSet,1)
    id = D.laserCoord(:,1) == coordSet(i,1) & D.laserCoord(:,2) == coordSet(i,2);
    D.laserIdx(id) = i;
end

%Remove trials with pre-stimulus wheel movement
keepIdx = zeros(size(D.response));
for sess = 1:max(D.sessionID)
    times = D.wheel_stimulusOn_timesteps(D.sessionID==sess,:);
    times = times(1,:);
    wheel = D.wheel_stimulusOn(D.sessionID==sess,:);
    choice = D.response(D.sessionID==sess);
    
    %Whiten
    wheel_whiten =  wheel/std(wheel(:,end));
    
    %Threshold
    wheel_whiten_window = mean( abs(wheel_whiten(:,-0.150 < times & times < 0.05 & times ~= 0)) , 2);
    keepIdx(D.sessionID==sess) = wheel_whiten_window < 0.05;
end
fprintf('Removing %0.2f%% of trials becuase of pre-stim wheel movements\n',100*mean(keepIdx==0));
D = D(keepIdx==1,:);

wheel_stimulusOn_timesteps = D.wheel_stimulusOn_timesteps(1,:);
D.wheel_stimulusOn_timesteps = [];
D.wheel_startMove_timesteps=[];
D.wheel_startMove=[];

%Save the data 
save('../data/inactivation/Inactivation_HigherPower.mat','D','wheel_stimulusOn_timesteps','-v7.3');

%% Fit neurometric model to Higher power sessions, using widefield fit as prior
clear all; close all;
load('../data/inactivation/Inactivation_HigherPower.mat','D');

fit_mech = load("../data/modelFits/neurometric_wf_original.mat"); %Original fit to widefield data
w_mu = mean( fit_mech.posterior.w, 1);
w_S = cov( fit_mech.posterior.w);

%Get trials with inactivation in VIS and M2 only
perturbationRegion = {'LeftVIS','RightVIS','LeftM2','RightM2'};
D1 = D(any(D.laserRegion == perturbationRegion,2) | D.laserType==0,:);
D1.perturbation = zeros(size(D1.response));
for p = 1:length(perturbationRegion)
    D1.perturbation(D1.laserRegion == perturbationRegion{p}) = p;
end
NL_idx = D1.laserType==0;

fr = load("../data/neuropixels/visMOs_allConds.mat");

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
dat.contrastLeft = D1.stimulus(NL_idx,1);
dat.contrastRight = D1.stimulus(NL_idx,2);
dat.choice = D1.response(NL_idx);
dat.sessionID = D1.sessionID(NL_idx);
dat.subjectID = D1.subjectID(NL_idx);
dat.firing_rate = F(NL_idx,:);
dat.PRIOR_MU = w_mu;
dat.PRIOR_S = w_S;
dat.PRIOR_S(1,1) = 1; %BROADEN PRIOR ON ALPHAS
dat.PRIOR_S(2,2) = 1; %BROADEN PRIOR ON ALPHAS
dat.PRIOR_S = dat.PRIOR_S + diag(diag(w_S)); %BROADEN PRIOR ON ALPHAS
% fit_inact = stan_fitModel('Two-Level-Mechanistic-SYMMETRICAL-PRESETW',dat);

fit_inact = load("C:\stanFitDump\tpaeecaac5_464c_4218_b5b9_ee96a3152748.mat");

%% Plot weights of inactivation-neurometric model
fit_inact = load("../data/modelFits/neurometric_inact_original.mat"); %Original fit to widefield data

areaCols = [ 109 183 229;
                77 131 158;
                160 206 87;
             119 147 62]/255;
         
%Plot posterior of Weights, copying over symmetrical values
WL = fit_inact.posterior.w(:,2+[1:4]);
WR = fit_inact.posterior.w(:,2+[2 1 4 3]);
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

%% Plot fit of widefield-neurometric model
fit_mech = load("../data/modelFits/neurometric_wf_original.mat"); %Original fit to widefield data
fr = load("..\data\neuropixels\visMOs_allConds.mat");
VIS_window = [0.075 0.125];
MOS_window = [0.125 0.175];
wf_window_delay = 0.03;
vis_ephys = mean( mean( fr.vis(:,:,:,VIS_window(1) < fr.binsS & fr.binsS < VIS_window(2)), 4), 3);
m2_ephys = mean( mean( fr.mos(:,:,:,MOS_window(1) < fr.binsS & fr.binsS < MOS_window(2)), 4), 3);

%Calculate behavioural choice proportions
uCL = unique(fit_mech.data.contrastLeft);
uCR = unique(fit_mech.data.contrastRight);
[counts,~,~,labels] = crosstab(fit_mech.data.contrastLeft,...
    fit_mech.data.contrastRight,...
    fit_mech.data.choice,...
    fit_mech.data.sessionID);
prob = counts./sum(counts,3);%Convert to probability over choices
prob = nanmean(prob,4);


figure; hold on;
CL = [linspace(0.15,1,200), zeros(1,200)];
CR = [zeros(1,200), linspace(0.15,1,200)];
F=[];
F(1,:) = interp2(fr.ucl, fr.ucr, vis_ephys, CL , CR, 'linear' );
F(2,:) = interp2(fr.ucl, fr.ucr, vis_ephys', CL , CR, 'linear' );
F(3,:) = interp2(fr.ucl, fr.ucr, m2_ephys, CL , CR, 'linear' );
F(4,:) = interp2(fr.ucl, fr.ucr, m2_ephys', CL , CR, 'linear' );
NL_ph=stan_plotFcn_MECH_SYMETRICAL(fit_mech.posterior.w, F);

pCorrect = mean(cat(4,NL_ph(:,CL>0,1),NL_ph(:,CR>0,2)),4);
q_pCorrect = quantile(pCorrect,[0.025 0.975],1);
plot(CL(CL>0),mean(pCorrect,1),'-','color',[1 1 1]*0,'linewidth',2);
fx=fill([CL(CL>0) fliplr(CL(CL>0))],[q_pCorrect(1,:) fliplr(q_pCorrect(2,:))],'k'); fx.FaceAlpha=0.2; fx.EdgeAlpha=0;
fx.FaceColor=[1 1 1]*0;
plot(uCL(uCL>0),mean(cat(2,prob(2:end,1,1),prob(1,2:end,2)'),2),'.','color',[1 1 1]*0,'markersize',30);

%pIncorrectOpposite
%(CL + R choice) & (CR + L choice)
pICO = mean(cat(4,NL_ph(:,CL>0,2),NL_ph(:,CR>0,1)),4);
q_pICO = quantile(pICO,[0.025 0.975],1);
plot(CL(CL>0),mean(pICO,1),'-','color',[1 1 1]*0.8,'linewidth',2);
fx=fill([CL(CL>0) fliplr(CL(CL>0))],[q_pICO(1,:) fliplr(q_pICO(2,:))],'k'); fx.FaceAlpha=0.2; fx.EdgeAlpha=0;
fx.FaceColor=[1 1 1]*0.8;
plot(uCL(uCL>0),mean(cat(2,prob(2:end,1,2),prob(1,2:end,1)'),2),'.','color',[1 1 1]*0.8,'markersize',30);

%pIncorrectNoGo
pING = mean(cat(4,NL_ph(:,CL>0,3),NL_ph(:,CR>0,3)),4);
q_pING = quantile(pING,[0.025 0.975],1);
plot(CL(CL>0),mean(pING,1),'-','color',[1 1 1]*0.4,'linewidth',2);
fx=fill([CL(CL>0) fliplr(CL(CL>0))],[q_pING(1,:) fliplr(q_pING(2,:))],'k'); fx.FaceAlpha=0.2; fx.EdgeAlpha=0;
fx.FaceColor=[1 1 1]*0.4;
xlabel('Contrast (single side)');
plot(uCL(uCL>0),mean(cat(2,prob(2:end,1,3),prob(1,2:end,3)'),2),'.','color',[1 1 1]*0.4,'markersize',30);

set(gca,'xlim',[0 1],'ylim',[0 1],'TickDir','out','dataaspectratio',[1 1 1],...
    'xticklabelmode','auto','xtick', uCL);
