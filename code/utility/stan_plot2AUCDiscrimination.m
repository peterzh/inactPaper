function plot2AUCDiscrimination(posterior,data)
colours = [ 0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.4940    0.1840    0.5560];

%Compute data probabilities
if ~isfield(data,'perturbation')
    data.perturbation = zeros(size(data.choice));
end

%% Plot model fit for each session

%Transform to pedestal/diff representation
data.pedestal = min(data.contrastLeft,data.contrastRight);
data.cDiff = data.contrastRight - data.contrastLeft;

[counts,~,~,labels] = crosstab(data.cDiff,...
    data.pedestal,...
    data.choice,...
    data.perturbation,...
    data.sessionID);
prob = counts./sum(counts,3);%Convert to probability over choices
prob_ave = nanmean(prob,5);


cdiffs = cellfun(@str2num,labels(1:size(prob,1),1));
peds = cellfun(@str2num,labels(1:size(prob,2),2));

for subj = 1:data.numSubjects
    
    sessionsThisSubj = find(data.subjID_session == subj);
    
    f=figure('color','w');
    ha_background = tight_subplot( ceil(sqrt( length(sessionsThisSubj) )), ceil(sqrt( length(sessionsThisSubj) )) , [0.04 0.01], [0.01 0.04], 0.01 );
    
    for sid = 1:length(sessionsThisSubj)        
        sessNum = sessionsThisSubj(sid);
        
        %Parameters for this session
        BL = posterior.bias(:,1) + posterior.b_sess(:,1,sessNum) + posterior.b_subj(:,1,subj);
        BR = posterior.bias(:,2) + posterior.b_sess(:,2,sessNum) + posterior.b_subj(:,2,subj);
        SL = posterior.sens(:,1) + posterior.b_sess(:,3,sessNum) + posterior.b_subj(:,3,subj);
        SR = posterior.sens(:,2) + posterior.b_sess(:,4,sessNum) + posterior.b_subj(:,4,subj);
        N = posterior.n_exp      + posterior.b_sess(:,5,sessNum) + posterior.b_subj(:,5,subj);
        
        %For each session, add new subplots
        set(ha_background(sid),'units','normalized');
        bounds = get(ha_background(sid),'position');
        set(ha_background(sid),'xcolor','none','ycolor','none');
        
        title(ha_background(sid),sprintf('M%d S%d',subj,sessNum));
        
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
            ph=bplot.CN(BL,SL,BR,SR,N,CL,CR);
            ph_interval = quantile(ph,[0.025 0.975]);
            ph_mean = mean(ph,1);
            
            for r = 1:3
                %plot datapoints
                plot(ha( 3*(ped-1) + r ),cdiffs,prob(:,ped,r,1,sessNum),'.','markersize',10,'color',colours(r,:));
                
                fx = fill(ha( 3*(ped-1) + r ),[CR-CL fliplr(CR-CL)], [ph_interval(1,:,r) fliplr( ph_interval(2,:,r) ) ], 'k');
                fx.FaceAlpha=0.3;
                fx.EdgeAlpha=0;
                fx.FaceColor = colours(r,:);
                plot(ha( 3*(ped-1) + r ),CR-CL, ph_mean(1,:,r), '-','color',colours(r,:));
            end
            
        end
        
    end
    
    set(get(f,'children'),'xlim',[-1 1]*1,'ylim',[0 1]);

end
% 
% 
% f=figure('color','w');
% ha_background = tight_subplot( ceil(sqrt(max(data.sessionID))), ceil(sqrt(max(data.sessionID))), [0.04 0.01], [0.01 0.04], 0.01 );
% Plot curves for each session
% for session = 1:max(data.sessionID)
%     Parameters for this session
%     BL = posterior.bias(:,1) + posterior.b_sess(:,1,session) + posterior.b_subj(:,1,data.subjID_session(session));
%     BR = posterior.bias(:,2) + posterior.b_sess(:,2,session) + posterior.b_subj(:,2,data.subjID_session(session));
%     SL = posterior.sens(:,1) + posterior.b_sess(:,3,session) + posterior.b_subj(:,3,data.subjID_session(session));
%     SR = posterior.sens(:,2) + posterior.b_sess(:,4,session) + posterior.b_subj(:,4,data.subjID_session(session));
%     N = posterior.n_exp      + posterior.b_sess(:,5,session) + posterior.b_subj(:,5,data.subjID_session(session));
%     
%     For each session, add new subplots
%     set(ha_background(session),'units','normalized');
%     bounds = get(ha_background(session),'position');
%     set(ha_background(session),'xcolor','none','ycolor','none');
%     
%     title(ha_background(session),sprintf('M%d S%d',data.subjID_session(session),session));
%     
%     replace with many subplots within the bounds of the original
%     ha = tight_subplot(4,3,0.001, [bounds(2) 1-bounds(2)-bounds(4)], [bounds(1) 1-bounds(1)-bounds(3)]);
%     set(ha(1:9),'xcolor','none');
%     set(ha([2 3 5 6 8 9 11 12]),'ycolor','none');
%     for i = 1:length(ha)
%         hold(ha(i),'on');
%     end
%     
%     for ped = 1:length(peds)
%         CL = [linspace(1-peds(ped),0,100) zeros(1,100)] + peds(ped);
%         CR = [zeros(1,100) linspace(0,1-peds(ped),100)] + peds(ped);
%         ph=bplot.CN(BL,SL,BR,SR,N,CL,CR);
%         ph_interval = quantile(ph,[0.025 0.975]);
%         ph_mean = mean(ph,1);
%         
%         for r = 1:3
%             plot datapoints
%             plot(ha( 3*(ped-1) + r ),cdiffs,prob(:,ped,r,1,session),'.','markersize',10,'color',colours(r,:));
%             
%             fx = fill(ha( 3*(ped-1) + r ),[CR-CL fliplr(CR-CL)], [ph_interval(1,:,r) fliplr( ph_interval(2,:,r) ) ], 'k');
%             fx.FaceAlpha=0.3;
%             fx.EdgeAlpha=0;
%             fx.FaceColor = colours(r,:);
%             plot(ha( 3*(ped-1) + r ),CR-CL, ph_mean(1,:,r), '-','color',colours(r,:));
%         end
%         
%     end
%     
% end

% set(get(f,'children'),'xlim',[-1 1]*1,'ylim',[0 1]);


%% Plot grand average

%Parameters of the 'average session'
BL = posterior.bias(:,1) ;
BR = posterior.bias(:,2) ;
SL = posterior.sens(:,1) ;
SR = posterior.sens(:,2) ;
N = posterior.n_exp ;

%Parameters drawn from the variation across subjects and sessions
% BL = posterior.bias(:,1) + posterior.sd_subj(:,1).*randn(2000,1) + posterior.sd_sess(:,1).*randn(2000,1);
% BR = posterior.bias(:,2) + posterior.sd_subj(:,2).*randn(2000,1) + posterior.sd_sess(:,2).*randn(2000,1);
% SL = posterior.sens(:,1) + posterior.sd_subj(:,3).*randn(2000,1) + posterior.sd_sess(:,3).*randn(2000,1);
% SR = posterior.sens(:,2) + posterior.sd_subj(:,4).*randn(2000,1) + posterior.sd_sess(:,4).*randn(2000,1);
% N = posterior.n_exp      + posterior.sd_subj(:,5).*randn(2000,1) + posterior.sd_sess(:,5).*randn(2000,1);

f=figure;
ha = tight_subplot(4,3,0.001, 0.01,0.01);
for i = 1:length(ha)
    hold(ha(i),'on');
end

numSess = size(prob,5);
sessJitter = 0.01*randn(numSess,1);
for ped = 1:length(peds)
    CL = [linspace(1-peds(ped),0,100) zeros(1,100)] + peds(ped);
    CR = [zeros(1,100) linspace(0,1-peds(ped),100)] + peds(ped);
    ph=bplot.CN(BL,SL,BR,SR,N,CL,CR);
    ph_interval = quantile(ph,[0.025 0.975]);
    ph_mean = nanmean(ph,1);
    
    for r = 1:3
        %plot datapoints
        plot(ha( 3*(ped-1) + r ),cdiffs,prob_ave(:,ped,r,1),'.','markersize',10,'color',colours(r,:));
        plot(ha( 3*(ped-1) + r ),cdiffs + sessJitter',permute( prob(:,ped,r,1,:), [1 5 2 3 4]) ,'.','markersize',1,'color',colours(r,:));
        
        
        
        fx = fill(ha( 3*(ped-1) + r ),[CR-CL fliplr(CR-CL)], [ph_interval(1,:,r) fliplr( ph_interval(2,:,r) ) ], 'k');
        fx.FaceAlpha=0.3;
        fx.EdgeAlpha=0;
        fx.FaceColor=colours(r,:);
        plot(ha( 3*(ped-1) + r ),CR-CL, ph_mean(1,:,r), '-','color',colours(r,:));
    end
    
end
set(get(f,'children'),'xlim',[-1 1]*1,'ylim',[0 1]);

%% Plot grand average mesh plots
c2 = crosstab(data.contrastLeft,...
    data.contrastRight,...
    data.choice,...
    data.perturbation,...
    data.sessionID);
prob = c2./sum(c2,3);%Convert to probability over choices
prob_ave = nanmean(prob,5); %average over sessions

CLs = unique(data.contrastLeft);
CRs = unique(data.contrastRight);

DotSize = 32;
openDotSize = 9;

f = figure('color','w');
ha = tight_subplot(1,3,0.1,0.1,0.1);
for i = 1:3; hold(ha(i),'on'); end;

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
set(ha, 'xticklabelmode','auto',...
        'xtick', CLs,...
        'yticklabelmode','auto',...
        'ytick', CRs,...
        'xlim', [0 0.54],...
        'ylim', [0 0.54],...
        'zlim',[0 1],...
        'tickdir','out');
set(ha,'xdir','reverse'); 
set(ha(1),'view',[-160,45]);
set(ha(2),'view',[-135,45]);
set(ha(3),'view',[-110,45]);

title(ha(1),'Left');
title(ha(2),'NoGo');
title(ha(3),'Right');

%% Plot CN function
c = linspace(0,1,1000);
cn = c.^posterior.n_exp;
figure; axes; hold on;
plot(c,mean(cn,1),'k-');
q = quantile(cn,[0.025 0.975],1);
fx = fill([c fliplr(c)], [q(1,:) fliplr( q(2,:) ) ], 'k');
fx.FaceAlpha=0.3;
fx.EdgeAlpha=0;
set(gca,'dataaspectratio',[1 1 1]); 
xlabel('Contrast'); ylabel('f(Contrast)');


%% Try out the other kind of curves
[counts,~,~,labels] = crosstab(data.contrastLeft,...
    data.contrastRight,...
    data.choice,...
    data.perturbation,...
    data.sessionID);
prob = counts./sum(counts,3);%Convert to probability over choices
prob = prob(:,:,:,1,:); %only laser off conditon
prob_ave = nanmean(prob,5);


cls = cellfun(@str2num,labels(1:size(prob,1),1));
crs = cellfun(@str2num,labels(1:size(prob,2),2));

colours = [ 0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.4940    0.1840    0.5560];


figure;
ha = tight_subplot(2,2,0.05,0.05,0.05);
for i = 1:length(ha)
    hold(ha(i),'on');
end
set(ha([1 3]),'xdir','reverse','yaxislocation','right');

%Varying CL while holding CR = 0
CL = linspace(0,1,100);
CR = zeros(1,100);
ph=bplot.CN(BL,SL,BR,SR,N, CL , CR );
ph_interval = quantile(ph,[0.025 0.975]);
ph_mean = mean(ph,1);
for r=1:3
    plot(ha(1),CL, ph_mean(1,:,r), '-','Color',colours(r,:));
    fx = fill(ha(1),[CL fliplr(CL)], [ph_interval(1,:,r) fliplr( ph_interval(2,:,r) ) ], 'k');
    fx.FaceAlpha=0.3;
    fx.EdgeAlpha=0;
    fx.FaceColor = colours(r,:);
    
    plot(ha(1),cls,prob_ave(:,1,r),'.','markersize',15,'color',colours(r,:));
%     plot(ha(1),cls,squeeze(prob(:,1,r,1,:)),'.','markersize',1,'color',colours(r,:));
end

%Varying CR while holding CL = 0
CR = linspace(0,1,100);
CL = zeros(1,100);
ph=bplot.CN(BL,SL,BR,SR,N, CL , CR );
ph_interval = quantile(ph,[0.025 0.975]);
ph_mean = mean(ph,1);
for r=1:3
    plot(ha(2),CR, ph_mean(1,:,r), '-','Color',colours(r,:));
    fx = fill(ha(2),[CR fliplr(CR)], [ph_interval(1,:,r) fliplr( ph_interval(2,:,r) ) ], 'k');
    fx.FaceAlpha=0.3;
    fx.EdgeAlpha=0;
    fx.FaceColor = colours(r,:);
    
    plot(ha(2),cls,prob_ave(1,:,r),'.','markersize',15,'color',colours(r,:));
%     plot(ha(2),cls,squeeze(prob(1,:,r,1,:)),'.','markersize',1,'color',colours(r,:));
end

%Varying CL while holding CR at different values, only plotting one side
for cr = 1:length(crs)
    CL = linspace(0,1,100);
    CR = ones(1,100)*crs(cr);
    ph=bplot.CN(BL,SL,BR,SR,N, CL , CR );
    ph_mean = mean(ph,1);
    
    r=1;
    plot(ha(3),CL, ph_mean(1,:,r), '-','Color',[colours(r,:) cr/length(crs)],'linewidth',1);
    
    px=scatter(ha(3), cls,prob_ave(:,cr,r),'filled','SizeData',15);
    px.MarkerFaceColor = colours(r,:);
    alpha(px,cr/length(crs))
end

%Varying CR while holding CL at different values, only plotting one side
for cl = 1:length(cls)
    CR = linspace(0,1,100);
    CL = ones(1,100)*cls(cl);
    ph=bplot.CN(BL,SL,BR,SR,N, CL , CR );
    ph_mean = mean(ph,1);
    
    r=2;
    plot(ha(4),CR, ph_mean(1,:,r), '-','Color',[colours(r,:) cl/length(cls)],'linewidth',1);
        
    px=scatter(ha(4), crs,prob_ave(cl,:,r),'filled','SizeData',15);
    px.MarkerFaceColor = colours(r,:);
    alpha(px,cl/length(cls))
end
set(ha,'xlim',[0 1],'ylim',[0 1],'dataaspectratio',[1 1 1],'xticklabelmode','auto')

xlabel(ha(3),'ContrastLeft');
xlabel(ha(4),'ContrastRight');
ylabel(ha(3),'pChoice');

title(ha(1),'CR = 0');
title(ha(2),'CL = 0');

title(ha(3),'CR = [0, 0.1, 0.24, 0.54]');
title(ha(4),'CL = [0, 0.1, 0.24, 0.54]');

set(ha,'xtick',cls);

%% Plot per-session variation in nogo (c=0) behaviour
pNG = squeeze(prob(1,1,3,1,:));
subjID = data.subjID_session;

[~,sortIdx]=sortrows([subjID pNG],[1 2 ]);
% [~,sortIdx] = sort(pNG);

subjCols = cool(max(subjID));

f=figure('color','w');
axes; hold on;
for i = 1:length(sortIdx)
    thisSessIdx = sortIdx(i);
    BL = posterior.bias(:,1) + posterior.b_sess(:,1,thisSessIdx) + posterior.b_subj(:,1,data.subjID_session(thisSessIdx));
    BR = posterior.bias(:,2) + posterior.b_sess(:,2,thisSessIdx) + posterior.b_subj(:,2,data.subjID_session(thisSessIdx));
    SL = posterior.sens(:,1) + posterior.b_sess(:,3,thisSessIdx) + posterior.b_subj(:,3,data.subjID_session(thisSessIdx));
    SR = posterior.sens(:,2) + posterior.b_sess(:,4,thisSessIdx) + posterior.b_subj(:,4,data.subjID_session(thisSessIdx));
    N = posterior.n_exp      + posterior.b_sess(:,5,thisSessIdx) + posterior.b_subj(:,5,data.subjID_session(thisSessIdx));
    
    plot(i, prob(1,1,3,1,thisSessIdx), '.', 'markersize', 20, 'color', subjCols(subjID(thisSessIdx),:));
    
    %model prediction
    ph=bplot.CN(BL,SL,BR,SR,N, 0 , 0 );
    pM = mean(ph(:,1,3),1);
    pQ = quantile(ph(:,1,3),[0.025 0.975],1);
    
    fx=fill(i + [-1 1 1 -1]*0.45, pM + [-1 -1 1 1]*0.005, 'k' ); fx.EdgeAlpha=0;
    fx.FaceColor = subjCols(subjID(thisSessIdx),:);
    
    
    fx=fill(i + [-1 1 1 -1]*0.45, [pQ(1) pQ(1) pQ(2) pQ(2)] , 'k' ); alpha(fx,0.1); fx.EdgeAlpha=0.1;
    fx.FaceColor = subjCols(subjID(thisSessIdx),:);
end

ylabel('p(NoGo) on zero contrast trials');
xlabel('Sorted session');
ylim([0 1]);

%% Plot more comprehensive inter-subj and inter-sess differences
figure('color','w');
p = [squeeze(prob(1,1,3,1,:)) squeeze(prob(4,1,1,1,:)) squeeze(prob(1,4,2,1,:)) (squeeze(prob(4,4,1,1,:)) - squeeze(prob(4,4,2,1,:)))];

ha_sess = tight_subplot(4,1,0.01,0.01,[0.01 0.3]);
for i = 1:length(ha_sess); hold(ha_sess(i),'on'); end;
for sess = 1:data.numSessions
    BL = posterior.bias(:,1) + posterior.b_sess(:,1,sess) + posterior.b_subj(:,1,data.subjID_session(sess));
    BR = posterior.bias(:,2) + posterior.b_sess(:,2,sess) + posterior.b_subj(:,2,data.subjID_session(sess));
    SL = posterior.sens(:,1) + posterior.b_sess(:,3,sess) + posterior.b_subj(:,3,data.subjID_session(sess));
    SR = posterior.sens(:,2) + posterior.b_sess(:,4,sess) + posterior.b_subj(:,4,data.subjID_session(sess));
    N = posterior.n_exp      + posterior.b_sess(:,5,sess) + posterior.b_subj(:,5,data.subjID_session(sess));
    
    %NG on C0
    ph=bplot.CN(BL,SL,BR,SR,N, [0 0.54 0 0.54] , [0 0 0.54 0.54] );
    ph = [ph(:,1,3) ph(:,2,1) ph(:,3,2) ph(:,4,1)-ph(:,4,2)];
    pM = mean(ph,1);
    pQ = quantile(ph,[0.025 0.975],1);
    
    for i = 1:4
        fx=fill(ha_sess(i),sess + [-1 1 1 -1]*0.45, pM(i) + [-1 -1 1 1]*0.005, 'k' ); fx.EdgeAlpha=0;
        fx.FaceColor = subjCols(data.subjID_session(sess),:);
        fx=fill(ha_sess(i),sess + [-1 1 1 -1]*0.45, [pQ(1,i) pQ(1,i) pQ(2,i) pQ(2,i)] , 'k' ); alpha(fx,0.1); fx.EdgeAlpha=0.1;
        fx.FaceColor = subjCols(data.subjID_session(sess),:);
        plot(ha_sess(i),sess, p(sess,i), '.', 'markersize',10,'color',subjCols(data.subjID_session(sess),:));
    end
end
set(ha_sess,'ylim',[0 1]);
set(ha_sess(end),'ylim',[-1 1]);

%% Plot curves for three sessions which span a variability axis
%plot BL and BR values for all sessions
% 
% pNOGO = nan(data.numSessions,1);
% for sess = 1:data.numSessions
%     c = data.choice(data.sessionID==sess & data.perturbation==0);
%     pNOGO(sess) = mean(c==3);
% end
% [pNOGO,sessIDs] = sort(pNOGO);
% 
% % warning('todo: choose sessions based on the overall go/nogo balance');
% % [GNGAxis, sessIdx] = sort(GNGAxis);
% % sessIDs=[ sessIdx(round([0.1 0.9]*length(sessIdx))) ];
% % sessIDs=[ 9 5 ]; %nice!
sessIDs=[ 9 1 5 ];
% 
figure;
ha = tight_subplot(length(sessIDs),2,0.05,0.05,0.05);
for i = 1:length(ha)
    hold(ha(i),'on');
end
set(ha([1:2:end]),'xdir','reverse','yaxislocation','right');
for sess = 1:length(sessIDs)

%      BL = posterior.bias(:,1) + posterior.b_subj(:,1,subjID(subj) );
%     BR = posterior.bias(:,2) + posterior.b_subj(:,2,subjID(subj) );
%     SL = posterior.sens(:,1) + posterior.b_subj(:,3,subjID(subj) );
%     SR = posterior.sens(:,2)  + posterior.b_subj(:,4,subjID(subj) );
%     N = posterior.n_exp   + posterior.b_subj(:,5,subjID(subj) );
    
    BL = posterior.bias(:,1) + posterior.b_sess(:,1,sessIDs(sess)) + posterior.b_subj(:,1,data.subjID_session(sessIDs(sess)) );
    BR = posterior.bias(:,2) + posterior.b_sess(:,2,sessIDs(sess)) + posterior.b_subj(:,2,data.subjID_session(sessIDs(sess)) );
    SL = posterior.sens(:,1) + posterior.b_sess(:,3,sessIDs(sess)) + posterior.b_subj(:,3,data.subjID_session(sessIDs(sess)) );
    SR = posterior.sens(:,2) + posterior.b_sess(:,4,sessIDs(sess)) + posterior.b_subj(:,4,data.subjID_session(sessIDs(sess)) );
    N = posterior.n_exp      + posterior.b_sess(:,5,sessIDs(sess)) + posterior.b_subj(:,5,data.subjID_session(sessIDs(sess)) );
    
%     prob_subj = mean(prob(:,:,:,1,data.subjID_session==subjID(subj)),5);
    %Varying CL while holding CR = 0
    CL = linspace(0,1,100);
    CR = zeros(1,100);
    ph=bplot.CN(BL,SL,BR,SR,N, CL , CR );
    ph_interval = quantile(ph,[0.025 0.975]);
    ph_mean = mean(ph,1);
    for r=1:3
        plot(ha( 2*(sess-1) + 1 ),CL, ph_mean(1,:,r), '-','Color',colours(r,:));
        fx = fill(ha( 2*(sess-1) + 1 ),[CL fliplr(CL)], [ph_interval(1,:,r) fliplr( ph_interval(2,:,r) ) ], 'k');
        fx.FaceAlpha=0.3;
        fx.EdgeAlpha=0;
        fx.FaceColor = colours(r,:);
        
        plot(ha( 2*(sess-1) + 1 ),cls,prob(:,1,r,1,sessIDs(sess)),'.','markersize',15,'color',colours(r,:));
        %     plot(ha(1),cls,squeeze(prob(:,1,r,1,:)),'.','markersize',1,'color',colours(r,:));
    end
    
    ylabel(ha( 2*(sess-1) + 1 ),sessIDs(sess));
   
    
    %Varying CR while holding CL = 0
    CR = linspace(0,1,100);
    CL = zeros(1,100);
    ph=bplot.CN(BL,SL,BR,SR,N, CL , CR );
    ph_interval = quantile(ph,[0.025 0.975]);
    ph_mean = mean(ph,1);
    for r=1:3
        plot(ha( 2*(sess-1) + 2 ),CR, ph_mean(1,:,r), '-','Color',colours(r,:));
        fx = fill(ha( 2*(sess-1) + 2 ),[CR fliplr(CR)], [ph_interval(1,:,r) fliplr( ph_interval(2,:,r) ) ], 'k');
        fx.FaceAlpha=0.3;
        fx.EdgeAlpha=0;
        fx.FaceColor = colours(r,:);
        
        plot(ha( 2*(sess-1) + 2 ),cls,prob(1,:,r,1,sessIDs(sess)),'.','markersize',15,'color',colours(r,:));
        %     plot(ha(2),cls,squeeze(prob(1,:,r,1,:)),'.','markersize',1,'color',colours(r,:));
    end


end
set(ha,'xlim',[0 1],'ylim',[0 1],'dataaspectratio',[1 1 1],'xticklabelmode','auto')

xlabel(ha(1),'ContrastLeft');
xlabel(ha(2),'ContrastRight');
ylabel(ha(2),'pChoice');

title(ha(1),'CR = 0');
title(ha(2),'CL = 0');

set(ha,'xtick',cls);



% 
% %% Overlay per-session curves for NoGos at different detection contrasts
% f=figure;
% ha = tight_subplot(1,1,0.001, 0.01,0.01);
% for i = 1:length(ha)
%     hold(ha(i),'on');
% end
% 
% CL = [linspace(1,0,100) zeros(1,100)];
% CR = [zeros(1,100) linspace(0,1,100)];
% 
% for sess = 1:max(data.sessionID)
%     BL = posterior.bias(:,1) + posterior.b_sess(:,1,sess) + posterior.b_subj(:,1,data.subjID_session(sess));
%     BR = posterior.bias(:,2) + posterior.b_sess(:,2,sess) + posterior.b_subj(:,2,data.subjID_session(sess));
%     SL = posterior.sens(:,1) + posterior.b_sess(:,3,sess) + posterior.b_subj(:,3,data.subjID_session(sess));
%     SR = posterior.sens(:,2) + posterior.b_sess(:,4,sess) + posterior.b_subj(:,4,data.subjID_session(sess));
%     N = posterior.n_exp      + posterior.b_sess(:,5,sess) + posterior.b_subj(:,5,data.subjID_session(sess));
%     
%     ph=bplot.CN(BL,SL,BR,SR,N, CL , CR );
%     pNG = mean(ph(:,:,3),1);
%     
%     plot( CR-CL, pNG, 'k-');
% end

% Plot bL and bR across sessions
BL=squeeze(mean(posterior.bias(:,1) + posterior.b_sess(:,1,:) + posterior.b_subj(:,1,data.subjID_session),1));
BR=squeeze(mean(posterior.bias(:,2) + posterior.b_sess(:,2,:) + posterior.b_subj(:,2,data.subjID_session),1));
figure('color','w'); axes; hold on;
cols = lines(data.numSubjects);
for subj = 1:data.numSubjects
    idx = data.subjID_session==subj;
    plot(BL(idx),BR(idx),'.', 'markersize', 20, 'color', cols(subj,:));
end
plot(BL(sessIDs),BR(sessIDs),'ko','markersize',8);
set(gca,'dataaspectratio',[ 1 1 1]);
xlabel('B_L'); ylabel('B_R');
% 
% figure('color','w'); axes; hold on;
% for sess = 1:data.numSessions
%     B=squeeze(posterior.bias + posterior.b_sess(:,1:2,sess) + posterior.b_subj(:,1:2,data.subjID_session(sess)));
%     
%     mu = mean(B,1);
%     S = cov(B);
%     error_ellipse(S,mu);
% end
end