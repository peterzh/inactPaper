%% Calculate wheel velocity and wheel stats (peak vel, etc)
regions = {'LeftVIS','RightVIS','LeftM2','RightM2','LeftM1','RightM1','LeftS1','RightS1'};
% regions = {'LeftVIS','RightVIS','LeftM2','RightM2','LeftS1','RightS1'};

%smooth wheel position
t = D.wheel_stimulusOn_timesteps(1,:);
pos = D.wheel_stimulusOn;
dt = t(2)-t(1);
smooth_t = 10/1000;  %10ms smoothing
smoothWin = myGaussWin(smooth_t, 1/dt); 
D.vel = [zeros(size(pos,1),1) conv2(diff(pos,[],2), smoothWin', 'same')]/dt;

%compute wheel stats (peak vel, and time centroid) for each trial
D.peakVel = nan(size(D.vel,1),1);
% D.centroidVel = nan(size(D.vel,1),1);
% figure;
for tr = 1:size(D.vel,1)

    p = pos(tr,:);
    v = D.vel(tr,:);
    %clip V before stimulus onset
    p(t<0) = 0;
    v(t<0) = 0;
%     figure;
%     plot(t,p,'-',t,v,'-',D.RT(tr),1,'r+');

%     %Compute peak velocity
%     if D.response(tr)==1
%         D.peakVel(tr) = max(v);
%     elseif D.response(tr)==2
%         D.peakVel(tr) = min(v);
%     end
    D.peakVel(tr) = max(abs(v));
%     
%     %Compute centroid
%     D.centroidVel(tr) = sum(t.*abs(D.vel(tr,:)))/sum(abs(D.vel(tr,:)));
%     
end
D.wheelStats = [D.RT D.peakVel];
%% 1a) Plot peak velocity and RT as a fcn of contrast for L and R choices
figure;
ha = tight_subplot(2,2,0.01,0.1,0.1); for i=1:numel(ha); hold(ha(i),'on'); end;
subj_avg = nan(4,2,2,max(D.subjectID)); %contrast x stat x response x subj
for subj = 1:max(D.subjectID)
    E = getrow(D, D.laserType==0 & D.subjectID==subj);

    %Left choices, vary CL
    CL = unique(E.stimulus(:,1));
    avg = nan(length(CL),2);
    err = nan(length(CL),2);
    for cl = 1:length(CL)
        dat=E.wheelStats(E.response==1 & E.stimulus(:,1)==CL(cl),:);
        avg(cl,:) = mean(dat,1)';
        err(cl,:) = std(dat,[],1)'/sqrt(size(dat,1));
    end
    for stat = 1:2
        errorbar(ha(2*stat - 1), CL, avg(:,stat),err(:,stat));
    end
    subj_avg(:,:,1,subj) = avg;
    
    %Right choices, vary CR
    CR = unique(E.stimulus(:,2));
    avg = nan(length(CR),2);
    err = nan(length(CR),2);
    for cr = 1:length(CR)
        dat=E.wheelStats(E.response==2 & E.stimulus(:,2)==CR(cr),:);
        avg(cr,:) = mean(dat,1)';
        err(cr,:) = std(dat,[],1)'/sqrt(size(dat,1));
    end
    for stat = 1:2
        errorbar(ha(2*stat), CR, avg(:,stat),err(:,stat));
    end
    subj_avg(:,:,2,subj) = avg;
end
%Left choices, vary CL
for stat = 1:2
    avg=mean(subj_avg(:,stat,1,:),4);
    err=std(subj_avg(:,stat,1,:),[],4)/sqrt(max(D.subjectID));
    plot(ha(2*stat - 1), CL, avg, 'k-','linewidth',3);
    fill(ha(2*stat - 1),[CL; flipud(CL)]', [avg-err; flipud(avg+err)]','k','FaceAlpha',0.2,'EdgeAlpha',0);
end

%Right choices, vary CR
for stat = 1:2
    avg=mean(subj_avg(:,stat,2,:),4);
    err=std(subj_avg(:,stat,2,:),[],4)/sqrt(max(D.subjectID));
    plot(ha(2*stat ), CR, avg, 'k-','linewidth',3);
    fill(ha(2*stat ),[CR; flipud(CR)]', [avg-err; flipud(avg+err)]','k','FaceAlpha',0.2,'EdgeAlpha',0);
end
set(ha,'xticklabelmode','auto','yticklabelmode','auto');
ylabel(ha(1),'RT');
ylabel(ha(3),'Peak velocity');
xlabel(ha(3),'Contrast');
title(ha(1),'Left choices, varying CL');
title(ha(2),'Right choices, varying CR');
linkaxes(ha,'x');
linkaxes(ha(1:2),'y');

ha_NL = ha;
%% 1b) with inactivation in each region
for loc = 1:length(regions)
    figure('name',regions{loc});
    
    ha = tight_subplot(2,2,0.01,0.1,0.1); for i=1:numel(ha); hold(ha(i),'on'); end;
    subj_avg = nan(4,2,2,max(D.subjectID)); %contrast x stat x response x subj
    for subj = 1:max(D.subjectID)
        E = getrow(D, D.laserRegion==regions{loc} & D.subjectID==subj);
        
        %Left choices, vary CL
        CL = unique(E.stimulus(:,1));
        avg = nan(length(CL),2);
        err = nan(length(CL),2);
        for cl = 1:length(CL)
            dat=E.wheelStats(E.response==1 & E.stimulus(:,1)==CL(cl),:);
            avg(cl,:) = mean(dat,1)';
            err(cl,:) = std(dat,[],1)'/sqrt(size(dat,1));
        end
        for stat = 1:2
%             errorbar(ha(2*stat - 1), CL, avg(:,stat),err(:,stat),'r');
            plot(ha(2*stat - 1),CL,avg(:,stat),'ro-');
        end
        subj_avg(:,:,1,subj) = avg;
        
        %Right choices, vary CR
        CR = unique(E.stimulus(:,2));
        avg = nan(length(CR),2);
        err = nan(length(CR),2);
        for cr = 1:length(CR)
            dat=E.wheelStats(E.response==2 & E.stimulus(:,2)==CR(cr),:);
            avg(cr,:) = mean(dat,1)';
            err(cr,:) = std(dat,[],1)'/sqrt(size(dat,1));
        end
        for stat = 1:2
            plot(ha(2*stat),CR,avg(:,stat),'ro-');
            
        end
        subj_avg(:,:,2,subj) = avg;
    end
    %Left choices, vary CL
    for stat = 1:2
        avg=nanmean(subj_avg(:,stat,1,:),4);
        err=nanstd(subj_avg(:,stat,1,:),[],4)./sqrt(sum(~isnan(subj_avg(:,stat,1,:)),4));
        plot(ha(2*stat - 1), CL, avg, 'r-','linewidth',3);
        fill(ha(2*stat - 1),[CL; flipud(CL)]', [avg-err; flipud(avg+err)]','r','FaceAlpha',0.2,'EdgeAlpha',0);
    end
    
    %Right choices, vary CR
    for stat = 1:2
        avg=nanmean(subj_avg(:,stat,2,:),4);
        err=nanstd(subj_avg(:,stat,2,:),[],4)./sqrt(sum(~isnan(subj_avg(:,stat,2,:)),4));
        plot(ha(2*stat ), CR, avg, 'r-','linewidth',3);
        fill(ha(2*stat ),[CR; flipud(CR)]', [avg-err; flipud(avg+err)]','r','FaceAlpha',0.2,'EdgeAlpha',0);
    end
    set(ha,'xticklabelmode','auto','yticklabelmode','auto');
    ylabel(ha(1),'RT');
    ylabel(ha(3),'Peak velocity');
    xlabel(ha(3),'Contrast');
    title(ha(1),'Left choices, varying CL');
    title(ha(2),'Right choices, varying CR');
    linkaxes(ha,'x');
    linkaxes(ha(1:2),'y');
    
    %copy over non-laser traces, and set those to all black
    for i=1:4
        obj = get(ha_NL(i),'children');
        copiedObj = copyobj(obj,ha(i));
        set(copiedObj(2:end),'Color',[0 0 0]);
        copiedObj(3:end).delete;
    end
end
%% 1c) Plot peak velocity and RT as a fcn of contrast for L and R choices POOLED DATA
figure;
ha = tight_subplot(2,2); for i=1:numel(ha); hold(ha(i),'on'); end;
subj_avg = nan(4,2,2); %contrast x stat x response
E = getrow(D, D.laserType==0);

%Left choices, vary CL
CL = unique(E.stimulus(:,1));
avg = nan(length(CL),2);
err = nan(length(CL),2);
for cl = 1:length(CL)
    dat=E.wheelStats(E.response==1 & E.stimulus(:,1)==CL(cl),:);
    avg(cl,:) = mean(dat,1)';
    err(cl,:) = std(dat,[],1)'/sqrt(size(dat,1));
end
for stat = 1:2
    errorbar(ha(2*stat - 1), CL, avg(:,stat),err(:,stat),'k-','linewidth',3);
end

%Right choices, vary CR
CR = unique(E.stimulus(:,2));
avg = nan(length(CR),2);
err = nan(length(CR),2);
for cr = 1:length(CR)
    dat=E.wheelStats(E.response==2 & E.stimulus(:,2)==CR(cr),:);
    avg(cr,:) = mean(dat,1)';
    err(cr,:) = std(dat,[],1)'/sqrt(size(dat,1));
end
for stat = 1:2
    errorbar(ha(2*stat), CR, avg(:,stat),err(:,stat),'k-','linewidth',3);
end

set(ha,'xticklabelmode','auto','yticklabelmode','auto');
ylabel(ha(1),'RT');
ylabel(ha(3),'Peak velocity');
xlabel(ha(3),'Contrast');
title(ha(1),'Left choices, varying CL');
title(ha(2),'Right choices, varying CR');
linkaxes(ha,'x');
linkaxes(ha(1:2),'y');

ha_NL = ha;
%% 1d) with inactivation in each region POOLED DATA
for loc = 1:length(regions)
    figure('name',regions{loc});
    
    ha = tight_subplot(2,2); for i=1:numel(ha); hold(ha(i),'on'); end;
    subj_avg = nan(4,2,2,max(D.subjectID)); %contrast x stat x response x subj
    E = getrow(D, D.laserRegion==regions{loc});
    
    %Left choices, vary CL
    CL = unique(E.stimulus(:,1));
    avg = nan(length(CL),2);
    err = nan(length(CL),2);
    for cl = 1:length(CL)
        dat=E.wheelStats(E.response==1 & E.stimulus(:,1)==CL(cl),:);
        avg(cl,:) = mean(dat,1)';
        err(cl,:) = std(dat,[],1)'/sqrt(size(dat,1));
    end
    for stat = 1:2
    	errorbar(ha(2*stat - 1), CL, avg(:,stat),err(:,stat),'r','linewidth',3);
    end
    
    %Right choices, vary CR
    CR = unique(E.stimulus(:,2));
    avg = nan(length(CR),2);
    err = nan(length(CR),2);
    for cr = 1:length(CR)
        dat=E.wheelStats(E.response==2 & E.stimulus(:,2)==CR(cr),:);
        avg(cr,:) = mean(dat,1)';
        err(cr,:) = std(dat,[],1)'/sqrt(size(dat,1));
    end
    for stat = 1:2
%         plot(ha(2*stat),CR,avg(:,stat),'ro-');
        errorbar(ha(2*stat), CR, avg(:,stat),err(:,stat),'r','linewidth',3);

    end    
   
    set(ha,'xticklabelmode','auto','yticklabelmode','auto');
    ylabel(ha(1),'RT');
    ylabel(ha(3),'Peak velocity');
    xlabel(ha(3),'Contrast');
    title(ha(1),'Left choices, varying CL');
    title(ha(2),'Right choices, varying CR');
    linkaxes(ha,'x');
    linkaxes(ha(1:2),'y');
    
    %copy over non-laser traces, and set those to all black
    for i=1:4
        obj = get(ha_NL(i),'children');
        copiedObj = copyobj(obj,ha(i));
        set(copiedObj(2:end),'Color',[0 0 0]);
        copiedObj(3:end).delete;
    end
end

%% 2a) Plot pseudocolour map of peak velocity for all contrast conditions
subj_avg = nan(4,4,2,2,max(D.subjectID),1+length(regions)); %cl x cr x stat x response x subj x inactivation location
num_trials = nan(4,4,2,2,max(D.subjectID),1+length(regions));
CL = unique(D.stimulus(:,1));
CR = unique(D.stimulus(:,2));
for r=1:2
    for cl = 1:length(CL)
        for cr = 1:length(CR)
            for subj = 1:max(D.subjectID)
                dat = D.wheelStats(D.laserType==0 & D.subjectID==subj & D.response==r & D.stimulus(:,1)==CL(cl) & D.stimulus(:,2)==CR(cr),:);
                subj_avg(cl,cr,:,r,subj,1) = mean(dat,1);
                num_trials(cl,cr,:,r,subj,1) = size(dat,1);
                
                for loc = 1:length(regions)
                    dat = D.wheelStats(D.laserRegion==regions{loc} & D.subjectID==subj & D.response==r & D.stimulus(:,1)==CL(cl) & D.stimulus(:,2)==CR(cr),:);
                    subj_avg(cl,cr,:,r,subj,1+loc) = mean(dat,1);
                    num_trials(cl,cr,:,r,subj,1+loc) = size(dat,1);
                end
            end
        end
    end
end
subj_avg = nanmean(subj_avg,5);
num_trials = round(nanmean(num_trials,5));

%Set to nan anything with fewee than 10 avg trials
subj_avg(num_trials<10)=nan;

figure;
ha = tight_subplot(2,2); for i=1:numel(ha); hold(ha(i),'on');end;
for stat=1:2
    for r=1:2
        A = subj_avg(:,:,stat,r,1);
        aplt = nan(size(A)+1);
        aplt(1:end-1,1:end-1) = A;
        pcolor(ha( 2*(stat-1)+r), aplt )
        colorbar(ha( 2*(stat-1)+r)); 
        for cl = 1:length(CL)
            for cr = 1:length(CR)
                text( ha( 2*(stat-1)+r), cr+0.5,cl+0.5, num2str(num_trials(cl,cr,stat,r,1)), 'HorizontalAlignment','center' )
            end
        end
        xlabel(ha( 2*(stat-1)+r),'CR'); ylabel(ha( 2*(stat-1)+r),'CL');
        v=subj_avg(:,:,stat,r,1);
        set(ha( 2*(stat-1)+r), 'clim',quantile(v(:),[0.2 0.8]) );
    end
end
set(ha,'dataaspectratio',[1 1 1],'ydir','normal');
set(ha,'xtick',1.5:4.5,'xticklabel',{CR},'ytick',1.5:4.5,'yticklabel',{CL});
title(ha(1),'Left choice'); title(ha(2),'Right choice');
for loc = 1:length(regions)
    figure('name',regions{loc});
    ha = tight_subplot(2,2); for i=1:numel(ha); hold(ha(i),'on');end;
    for stat=1:2
        for r=1:2
            A = subj_avg(:,:,stat,r,1+loc);
            aplt = nan(size(A)+1);
            aplt(1:end-1,1:end-1) = A;
            pcolor(ha( 2*(stat-1)+r), aplt )
            colorbar(ha( 2*(stat-1)+r));
            for cl = 1:length(CL)
                for cr = 1:length(CR)
                    text( ha( 2*(stat-1)+r), cr+0.5,cl+0.5, num2str(num_trials(cl,cr,stat,r,1+loc)), 'HorizontalAlignment','center' )
                end
            end
            xlabel(ha( 2*(stat-1)+r),'CR'); ylabel(ha( 2*(stat-1)+r),'CL');
            v=subj_avg(:,:,stat,r,1+loc); v=v(:);
            if any(~isnan(v))
                set(ha( 2*(stat-1)+r), 'clim',quantile(v,[0.2 0.8]) );
            end
        end
    end
    set(ha,'dataaspectratio',[1 1 1],'ydir','normal');
    set(ha,'xtick',1.5:4.5,'xticklabel',{CR},'ytick',1.5:4.5,'yticklabel',{CL});
    title(ha(1),'Left choice'); title(ha(2),'Right choice');
end
%% 2b) pooling data across subjects instead of averaging over subjects
% pooling data across subjects
subj_avg = nan(4,4,2,2,1+length(regions)); %cl x cr x stat x response x inactivation location
num_trials = nan(4,4,2,2,1+length(regions));
CL = unique(D.stimulus(:,1));
CR = unique(D.stimulus(:,2));
for r=1:2
    for cl = 1:length(CL)
        for cr = 1:length(CR)
                dat = D.wheelStats(D.laserType==0 & D.response==r & D.stimulus(:,1)==CL(cl) & D.stimulus(:,2)==CR(cr),:);
                subj_avg(cl,cr,:,r,1) = mean(dat,1);
                num_trials(cl,cr,:,r,1) = size(dat,1);
                
                for loc = 1:length(regions)
                    dat = D.wheelStats(D.laserRegion==regions{loc} & D.response==r & D.stimulus(:,1)==CL(cl) & D.stimulus(:,2)==CR(cr),:);
                    subj_avg(cl,cr,:,r,1+loc) = mean(dat,1);
                    num_trials(cl,cr,:,r,1+loc) = size(dat,1);
                end
            
        end
    end
end
%Set to nan anything with fewee than 10 avg trials
subj_avg(num_trials<10)=nan;

figure;
ha = tight_subplot(2,2); for i=1:numel(ha); hold(ha(i),'on');end;
for stat=1:2
    for r=1:2
        A = subj_avg(:,:,stat,r,1);
        aplt = nan(size(A)+1);
        aplt(1:end-1,1:end-1) = A;
        pcolor(ha( 2*(stat-1)+r), aplt )
        colorbar(ha( 2*(stat-1)+r)); 
        for cl = 1:length(CL)
            for cr = 1:length(CR)
                text( ha( 2*(stat-1)+r), cr+0.5,cl+0.5, num2str(num_trials(cl,cr,stat,r,1)), 'HorizontalAlignment','center' )
            end
        end
        xlabel(ha( 2*(stat-1)+r),'CR'); ylabel(ha( 2*(stat-1)+r),'CL');
        v=subj_avg(:,:,stat,r,1);
        set(ha( 2*(stat-1)+r), 'clim',quantile(v(:),[0.2 0.8]) );
    end
end
set(ha,'dataaspectratio',[1 1 1],'ydir','normal');
set(ha,'xtick',1.5:4.5,'xticklabel',{CR},'ytick',1.5:4.5,'yticklabel',{CL});
title(ha(1),'Left choice'); title(ha(2),'Right choice');
for loc = 1:length(regions)
    figure('name',regions{loc});
    ha = tight_subplot(2,2); for i=1:numel(ha); hold(ha(i),'on');end;
    for stat=1:2
        for r=1:2
            A = subj_avg(:,:,stat,r,1+loc);
            aplt = nan(size(A)+1);
            aplt(1:end-1,1:end-1) = A;
            pcolor(ha( 2*(stat-1)+r), aplt )
            colorbar(ha( 2*(stat-1)+r));
            for cl = 1:length(CL)
                for cr = 1:length(CR)
                    text( ha( 2*(stat-1)+r), cr+0.5,cl+0.5, num2str(num_trials(cl,cr,stat,r,1+loc)), 'HorizontalAlignment','center' )
                end
            end
            xlabel(ha( 2*(stat-1)+r),'CR'); ylabel(ha( 2*(stat-1)+r),'CL');
            v=subj_avg(:,:,stat,r,1+loc); v=v(:);
            if any(~isnan(v))
                set(ha( 2*(stat-1)+r), 'clim',quantile(v,[0.2 0.8]) );
            end
        end
    end
    set(ha,'dataaspectratio',[1 1 1],'ydir','normal');
    set(ha,'xtick',1.5:4.5,'xticklabel',{CR},'ytick',1.5:4.5,'yticklabel',{CL});
    title(ha(1),'Left choice'); title(ha(2),'Right choice');
end
%% 2c) pooling data across subjects: delta to parameters
% pooling data across subjects
subj_avg = nan(4,4,2,2,1+length(regions)); %cl x cr x stat x response x inactivation location
num_trials = nan(4,4,2,2,1+length(regions));
CL = unique(D.stimulus(:,1));
CR = unique(D.stimulus(:,2));
for r=1:2
    for cl = 1:length(CL)
        for cr = 1:length(CR)
                dat = D.wheelStats(D.laserType==0 & D.response==r & D.stimulus(:,1)==CL(cl) & D.stimulus(:,2)==CR(cr),:);
                subj_avg(cl,cr,:,r,1) = mean(dat,1);
                num_trials(cl,cr,:,r,1) = size(dat,1);
                
                for loc = 1:length(regions)
                    dat = D.wheelStats(D.laserRegion==regions{loc} & D.response==r & D.stimulus(:,1)==CL(cl) & D.stimulus(:,2)==CR(cr),:);
                    subj_avg(cl,cr,:,r,1+loc) = mean(dat,1);
                    num_trials(cl,cr,:,r,1+loc) = size(dat,1);
                end
            
        end
    end
end
%Set to nan anything with fewee than 10 avg trials
subj_avg = subj_avg - subj_avg(:,:,:,:,1);
subj_avg(num_trials<10)=nan;

for loc = 1:length(regions)
    figure('name',regions{loc});
    ha = tight_subplot(2,2); for i=1:numel(ha); hold(ha(i),'on');end;
    for stat=1:2
        for r=1:2
            A = subj_avg(:,:,stat,r,1+loc);
            aplt = nan(size(A)+1);
            aplt(1:end-1,1:end-1) = A;
            pcolor(ha( 2*(stat-1)+r), aplt )
            colorbar(ha( 2*(stat-1)+r));
            for cl = 1:length(CL)
                for cr = 1:length(CR)
                    text( ha( 2*(stat-1)+r), cr+0.5,cl+0.5, num2str(num_trials(cl,cr,stat,r,1+loc)), 'HorizontalAlignment','center' )
                end
            end
            xlabel(ha( 2*(stat-1)+r),'CR'); ylabel(ha( 2*(stat-1)+r),'CL');
            v=subj_avg(:,:,stat,:,1+loc); v=v(:);
            if any(~isnan(v))
                set(ha( 2*(stat-1)+r), 'clim',[-1 1]*quantile(abs(v),0.8) );
            end
        end
    end
    set(ha,'dataaspectratio',[1 1 1],'ydir','normal');
    set(ha,'xtick',1.5:4.5,'xticklabel',{CR},'ytick',1.5:4.5,'yticklabel',{CL});
    title(ha(1),'Left choice'); title(ha(2),'Right choice');
    colormap(BlueWhiteRed);
end
%% Plot wheel traces for a single  session (highest trial count)
tab=tabulate(D.sessionID);
SESS_ID = tab(tab(:,2)==max(tab(:,2)),1);
figure('color','w'); 
ha = tight_subplot(2,1); for i = 1:2; hold(ha(i),'on'); end;
idx = D.laserType==0 & D.sessionID==SESS_ID;
choiceCols = [0 0 1; 1 0 0; 0 0 0];
for r = 1:3
    plot(ha(1),t,pos(idx & D.response==r,:)','Color',[choiceCols(r,:) 0.2]);
    plot(ha(1),t,mean( pos(idx & D.response==r,:),1)','Color',choiceCols(r,:),'linewidth',3);
    
    plot(ha(2),t,D.vel(idx & D.response==r,:)','Color',[choiceCols(r,:) 0.2]);
    plot(ha(2),t,mean( D.vel(idx & D.response==r,:),1)','Color',choiceCols(r,:),'linewidth',3);
end
set(ha,'xticklabelmode','auto','yticklabelmode','auto');
% xlim([0 1]);
xlabel(ha(2),'Time from move onset');
ylabel(ha(1),'Wheel position (mm)');
ylabel(ha(2),'Wheel velocity (mm/s)');
title(ha(1),sprintf('Single session ID=%d, %d non-laser trials',SESS_ID,sum(idx))); 
%% Plot wheel traces for each subject, pooling over sessions
figure;
ha = tight_subplot(1, max(D.subjectID));
for subj = 1:max(D.subjectID)
    %non laser case,
    idx = D.laserType==0 & D.subjectID==subj;
    
    hold(ha(subj),'on');
    plot(ha(subj),t,D.vel(idx & D.response==1,:)','Color',[0 0 1 0.1])
    plot(ha(subj),t,D.vel(idx & D.response==2,:)','Color',[1 0 0 0.1])
    plot(ha(subj),t,mean( D.vel(idx & D.response==1,:),1)','Color',[0 0 1 1],'linewidth',4)
    plot(ha(subj),t,mean( D.vel(idx & D.response==2,:),1)','Color',[1 0 0 1],'linewidth',4)
end
%% Plot subject-specific avg trace, and grand average trace
figure('color','w'); hold on;
subjAvg = nan(max(D.subjectID), length(t), 2);
cols = [0 0 1; 1 0 0];
for subj = 1:max(D.subjectID)
    idx = D.laserType==0 & D.subjectID==subj;

    for r = 1:2
        dat = D.vel(idx & D.response==r,:);
        avg = mean(dat,1);
        err = 1.96*std(dat,[],1)/sqrt(size(dat,1));
        
        plot(t,avg,'Color',cols(r,:),'linewidth',2)
        fill([t fliplr(t)], [avg-err fliplr(avg+err)],cols(r,:),'FaceAlpha',0.2,'EdgeAlpha',0);
        subjAvg(subj,:,r) = mean( dat,1);
    end
end
% %add pooled avg
plot(t,mean(D.vel(D.laserType==0 & D.response==1,:),1),'Color',[0 0 0],'linewidth',3);
plot(t,mean(D.vel(D.laserType==0 & D.response==2,:),1),'Color',[0 0 0],'linewidth',3);
xlabel('Time from move onset');
ylabel('Wheel velocity (mm/sec)');
title('Session-pooled wheel velocity for each subject');
% xlim([0 1]);
%% Plot pooled grand avg for different different contrast levels
figure('color','w'); 
cols = [0.8 0.8 0.8;  0.53 0.53 0.53; 0.25 0.25 0.25; 0 0 0 ];
subplot(1,2,1); hold on; title('Varying CL (CR anything)');
CL = unique(D.stimulus(:,1)); 
for cl = 1:length(CL)
    for r = 1:2
        dat = D.vel(D.laserType==0 & D.response==r & D.stimulus(:,1)==CL(cl),:);
        avg = mean(dat,1);
        err = std(dat,[],1)/sqrt(size(dat,1));
        plot(t,avg,'Color',cols(cl,:),'linewidth',2);
        fill([t fliplr(t)], [avg-err fliplr(avg+err)],cols(cl,:),'FaceAlpha',0.2,'EdgeAlpha',0);
    end
end
subplot(1,2,2); hold on; title('Varying CR (CL anything)');
CR = unique(D.stimulus(:,2));
for cr = 1:length(CR)
    for r = 1:2
       dat = D.vel(D.laserType==0 & D.response==r & D.stimulus(:,2)==CR(cr),:);
        avg = mean(dat,1);
        err = std(dat,[],1)/sqrt(size(dat,1));
        plot(t,avg,'Color',cols(cr,:),'linewidth',2);
        fill([t fliplr(t)], [avg-err fliplr(avg+err)],cols(cr,:),'FaceAlpha',0.2,'EdgeAlpha',0); 
    end
end
xlabel('Time from move onset');
ylabel('Wheel velocity (mm/sec)');
%% Per-subject histograms of wheel stats split by contrast
figure('color','w');
ha = tight_subplot(max(D.subjectID),4,0.01,0.1,0.1); 
for i=1:numel(ha); hold(ha(i),'on');end;
for subj = 1:max(D.subjectID)
    %Left choice for different CL values
    CL = unique(D.stimulus(:,1));
    for cl = 1:length(CL)
        p = D.peakVel(D.laserType==0 & D.response==1 & D.stimulus(:,1)==CL(cl) & D.subjectID==subj);
        histogram(ha(4*(subj-1)+1), p, 'DisplayStyle','stairs','Normalization','probability');
        
        p = D.RT(D.laserType==0 & D.response==1 & D.stimulus(:,1)==CL(cl) & D.subjectID==subj);
        histogram(ha(4*(subj-1)+2), p, 'DisplayStyle','stairs','Normalization','probability');
        
    end
    %Right choice for different CR values
    CR = unique(D.stimulus(:,2));
    for cr = 1:length(CR)
        p = D.peakVel(D.laserType==0 & D.response==2 & D.stimulus(:,2)==CR(cr) & D.subjectID==subj);
        histogram(ha(4*(subj-1)+3), p, 'DisplayStyle','stairs'); 
        
        p = D.RT(D.laserType==0 & D.response==2 & D.stimulus(:,2)==CR(cr) & D.subjectID==subj);
        histogram(ha(4*(subj-1)+4), p, 'DisplayStyle','stairs'); 
    end
    
%     set(ha(4*(subj-1)+1),'yticklabelmode','auto');
    ylabel(ha(4*(subj-1)+1),sprintf('Subject %d',subj));
end
linkaxes(ha(1:2:end),'x'); 
linkaxes(ha(2:2:end),'x');
set(ha(end-3:end),'xticklabelmode','auto','yticklabelmode','auto');
xlabel(ha(end-3),'Peak velocity (LEFT)');
xlabel(ha(end-2),'RT (LEFT)');
xlabel(ha(end-1),'Peak velocity (RIGHT)');
xlabel(ha(end),'RT (RIGHT)');
%% Plot wheel traces from one subject, with and without inactivation (pooling over sessions)
figure;
tab=tabulate(D.subjectID); %get subject ID with largest nuumber of trials overall
SUBJ_ID = tab(tab(:,2)==max(tab(:,2)),1);
for loc = 1:length(regions)
    subplot(2,3,loc); hold on; title(regions{loc});
    
    for r=1:2
        %get the inactivation trials for this subject
        dat = D.vel(D.laserRegion==regions{loc} & D.subjectID==SUBJ_ID & D.response==r,:);
        avg = mean(dat,1);
        err = 1.96*std(dat,[],1)/sqrt(size(dat,1));
        if size(dat,1) > 100
        	dat = dat( randsample(size(dat,1),100), :); %get random sample, matched for number of trials
        end
        plot(t,dat,'Color',[1 0 0 0.1]); %Plot 100 trials
        plot(t,avg,'Color',[1 0 0],'linewidth',3);
        fill([t fliplr(t)], [avg-err fliplr(avg+err)],[1 0 0],'FaceAlpha',0.2,'EdgeAlpha',0); 

        
        %get the non-inactivation trials for this subject
        dat = D.vel(D.laserType==0 & D.subjectID==SUBJ_ID & D.response==r,:);
        avg = mean(dat,1);
        err = std(dat,[],1)/sqrt(size(dat,1));
        if size(dat,1) > 100
        	dat = dat( randsample(size(dat,1),100), :); %get random sample, matched for number of trials
        end

        plot(t,dat,'Color',[0 0 0 0.1]);
        plot(t,avg,'Color',[0 0 0],'linewidth',3);
        fill([t fliplr(t)], [avg-err fliplr(avg+err)],[0 0 0],'FaceAlpha',0.2,'EdgeAlpha',0); 
    end
end
linkaxes(get(gcf,'children'),'xy');
%% Plot pooled wheel trace for all inactivation conditions overlaid
figure('color','w');
areaCols = [ 109 183 229;
    77 131 158;
    160 206 87;
    119 147 62;
    50 50 50;
    100 100 100]/255;
ha = gobjects(2,length(regions));
for r = 1:2
    hold on;
    
    %non-laser pooled condition
    dat = D.vel(D.laserType==0 & D.response==r,:);
    avg = mean(dat,1);
    err = 1.96*std(dat,[],1)/sqrt(size(dat,1));
    ha(r,1)=plot(t,avg,'Color',[0 0 0],'linewidth',2);
    fill([t fliplr(t)], [avg-err fliplr(avg+err)],[0 0 0],'FaceAlpha',0.2,'EdgeAlpha',0);
    
    
    %laser pooled condition
    for loc = 1:length(regions)
        
        dat = D.vel(D.laserRegion==regions{loc} & D.response==r,:);
        avg = mean(dat,1);
        err = 1.96*std(dat,[],1)/sqrt(size(dat,1));
        ha(r,1+loc)=plot(t,avg,'Color',areaCols(loc,:),'linewidth',2);
        fill([t fliplr(t)], [avg-err fliplr(avg+err)],areaCols(loc,:),'FaceAlpha',0.2,'EdgeAlpha',0);
    end
end
legend(ha(1,:),['OFF' regions])
xlabel('Time from move onset');
ylabel('Wheel velocity (mm/sec)');
title('Pooled data across all subjects');
%% Plot histograms of the movement stats with and without inactivation
figure('color','w');
ha = tight_subplot(max(D.subjectID),4,0.01,0.1,0.1); 
for i=1:numel(ha); hold(ha(i),'on');end;
for subj = 1:max(D.subjectID)
    %Left choice with and without inactivation
	p = D.peakVel(D.laserType==0 & D.response==1 & D.subjectID==subj);
    histogram(ha(4*(subj-1)+1), p, 'DisplayStyle','stairs','EdgeColor',[0 0 0],'Normalization','probability');
    p = D.RT(D.laserType==0 & D.response==1 & D.subjectID==subj);
    histogram(ha(4*(subj-1)+2), p, 'DisplayStyle','stairs','EdgeColor',[0 0 0],'Normalization','probability');
    for loc = 1:length(regions)
        p = D.peakVel(D.laserRegion==regions{loc} & D.response==1 & D.subjectID==subj);
        histogram(ha(4*(subj-1)+1), p, 'DisplayStyle','stairs','EdgeColor',areaCols(loc,:),'Normalization','probability');
        
        p = D.RT(D.laserRegion==regions{loc} & D.response==1 & D.subjectID==subj);
        histogram(ha(4*(subj-1)+2), p, 'DisplayStyle','stairs','EdgeColor',areaCols(loc,:),'Normalization','probability');
    end
    
    %Right choice with and without inactivation
	p = D.peakVel(D.laserType==0 & D.response==2 & D.subjectID==subj);
    histogram(ha(4*(subj-1)+3), p, 'DisplayStyle','stairs','EdgeColor',[0 0 0],'Normalization','probability');
    p = D.RT(D.laserType==0 & D.response==2 & D.subjectID==subj);
    histogram(ha(4*(subj-1)+4), p, 'DisplayStyle','stairs','EdgeColor',[0 0 0],'Normalization','probability');
    for loc = 1:length(regions)
        p = D.peakVel(D.laserRegion==regions{loc} & D.response==2 & D.subjectID==subj);
        histogram(ha(4*(subj-1)+3), p, 'DisplayStyle','stairs','EdgeColor',areaCols(loc,:),'Normalization','probability');
        
        p = D.RT(D.laserRegion==regions{loc} & D.response==2 & D.subjectID==subj);
        histogram(ha(4*(subj-1)+4), p, 'DisplayStyle','stairs','EdgeColor',areaCols(loc,:),'Normalization','probability');
    end
  
    ylabel(ha(4*(subj-1)+1),sprintf('Subject %d',subj));
end
linkaxes(ha(1:2:end),'x'); 
linkaxes(ha(2:2:end),'x');
set(ha(end-3:end),'xticklabelmode','auto','yticklabelmode','auto');
xlabel(ha(end-3),'Peak velocity (LEFT)');
xlabel(ha(end-2),'RT (LEFT)');
xlabel(ha(end-1),'Peak velocity (RIGHT)');
xlabel(ha(end),'RT (RIGHT)');
%% Plot pooled wheel traces with and without inactivation
figure('color','w');
ha = tight_subplot( length(regions), 2, 0.01,0.05,0.05);
for r = 1:2
    
    %non-laser pooled
    vel_off = D.vel(D.laserType==0 & D.response==r,:);
    avg_off = mean(vel_off,1);
    err_off = 1.96*std(vel_off,[],1)/sqrt(size(vel_off,1));
    
    for loc = 1:length(regions)
        hold(ha(2*(loc-1) +r),'on');
        
        %Plot off
        
        plot( ha(2*(loc-1) +r), t, avg_off, 'k-','linewidth',4);
        fill( ha(2*(loc-1) +r), [t fliplr(t)], [avg_off-err_off fliplr(avg_off+err_off)],'k','FaceAlpha',0.2)

        %Plot on
        vel_on = D.vel(D.laserRegion==regions{loc} & D.response==r,:);

        avg_on = mean(vel_on,1);
        err_on = 1.96*std(vel_on,[],1)/sqrt(size(vel_on,1));
        plot( ha(2*(loc-1) +r), t, avg_on, 'r-','linewidth',4);
        fill( ha(2*(loc-1) +r), [t fliplr(t)], [avg_on-err_on fliplr(avg_on+err_on)],'r','FaceAlpha',0.2);
        
        if r==1
            ylabel(ha(2*(loc-1) +r), regions{loc});
        end
    end
    
end
title(ha(1),'Left choice');
title(ha(2),'Right choice');
linkaxes(ha,'xy');
%% Does inactivation affect the wheel the same way contrast does?
figure;
for loc = 1:length(regions)
    subplot(3,2,loc); hold on;
    C = unique(D.stimulus(:));
    cols_inact = [1 0.8 0.8;  1 0.53 0.53; 1 0.25 0.25; 1 0 0 ];
    for c = 1:length(C) 
        for r = 1:2
            %laser off
            if contains(regions{loc},'Left')%if Left hemi inactivation: vary CR
                dat = D.vel(D.laserType==0 & D.response==r & D.stimulus(:,2)==C(c),:);
                
            elseif contains(regions{loc},'Right') %if Right hemi inactivation: vary CL
            	dat = D.vel(D.laserType==0 & D.response==r & D.stimulus(:,1)==C(c),:);
            end
            avg = mean(dat,1);
            err = std(dat,[],1)/sqrt(size(dat,1));
            plot(t,avg,'Color',cols(c,:),'linewidth',2);
            fill([t fliplr(t)], [avg-err fliplr(avg+err)],cols(c,:),'FaceAlpha',0.2,'EdgeAlpha',0);
        end
    end
    %Plot pooled laser off across conditions
    for r=1:2
        dat = D.vel(D.laserType==0 & D.response==r,:);
        avg = mean(dat,1);
        err = std(dat,[],1)/sqrt(size(dat,1));
        plot(t,avg,'Color','k','linewidth',2,'linestyle','--');
        
        dat = D.vel(D.laserRegion==regions{loc} & D.response==r,:);
        avg = mean(dat,1);
        err = std(dat,[],1)/sqrt(size(dat,1));
        plot(t,avg,'Color','r','linewidth',2,'linestyle','--');
    end
    ylabel(regions{loc});
    title('Varying contralateral contrast, non-laser (black) and  inactivation (red)'); 
end
linkaxes(get(gcf,'children'),'xy');

% Varying two side contrast
figure;
for loc = 1:length(regions)
    subplot(3,2,loc); hold on;
    C = unique(D.stimulus(:));
    cols_inact = [1 0.8 0.8;  1 0.53 0.53; 1 0.25 0.25; 1 0 0 ];
    for c = 1:length(C) 
        for r = 1:2
            %laser off
            dat = D.vel(D.laserType==0 & D.response==r & D.stimulus(:,1)==C(c) & D.stimulus(:,2)==C(c),:);
            avg = mean(dat,1);
            err = std(dat,[],1)/sqrt(size(dat,1));
            plot(t,avg,'Color',cols(c,:),'linewidth',2);
            fill([t fliplr(t)], [avg-err fliplr(avg+err)],cols(c,:),'FaceAlpha',0.2,'EdgeAlpha',0);
        end
    end
    %Plot pooled laser off across conditions
    for r=1:2
        dat = D.vel(D.laserType==0 & D.response==r,:);
        avg = mean(dat,1);
        err = std(dat,[],1)/sqrt(size(dat,1));
        plot(t,avg,'Color','k','linewidth',2,'linestyle','--');
        
        dat = D.vel(D.laserRegion==regions{loc} & D.response==r,:);
        avg = mean(dat,1);
        err = std(dat,[],1)/sqrt(size(dat,1));
        plot(t,avg,'Color','r','linewidth',2,'linestyle','--');
    end
    ylabel(regions{loc});
    title('Varying L&R contrast, non-laser (black) and  inactivation (red)'); 
end
linkaxes(get(gcf,'children'),'xy');

%% USED IN PAPER: Testing whether laser has significant effect on wheel peak velocity
%Contrast-matched

[contrastCondSet,~,D.contrastCond] = unique(D.stimulus,'rows');
D.laserType = logical(D.laserType);

peakvelchange = nan(length(regions),2,max(D.subjectID),max(D.contrastCond));
numTrials = zeros(size(peakvelchange));

figure('color','w');
ha = tight_subplot(length(regions),max(D.subjectID));

for loc = 1:length(regions)
    for r = 1:2
        for subj = 1:max(D.subjectID)
            
            %Test LaserOn vs NoLaser for this subject
            E = getrow(D, D.response == r & D.subjectID==subj & (D.laserType==0 | D.laserRegion==regions{loc}) );
            
            for c = 1:max(D.contrastCond) %for each contrast condition
                idx = E.contrastCond==c;
                x = E.peakVel(idx & E.laserType==0);
                y = E.peakVel(idx & E.laserType==1);
                if length(x)>1 && length(y)>1
        
                    %effect size: weighted by the proportion of those
                    %trials in the data for that subject. Needs to be
                    %summed together afterwards...
                    peakvelchange(loc,r,subj,c) = mean(idx)*( mean(y)-mean(x) );
                    numTrials(loc,r,subj,c) = min(length(x),length(y));

                else
%                     fprintf('No trials of cond=%d\n',c);
                end
            end
            
            %Plot avg wheel
            hold(ha((loc-1)*max(D.subjectID) + subj),'on');
            plot(ha((loc-1)*max(D.subjectID) + subj), t, mean(E.vel(E.laserType==0,:),1), 'k-');
            plot(ha((loc-1)*max(D.subjectID) + subj), t, mean(E.vel(E.laserType==1,:),1), 'r-');
            
           if loc==1
               title(ha((loc-1)*max(D.subjectID) + subj), sprintf('SUBJ %d', subj));
           end
           if subj == 1
              ylabel(ha((loc-1)*max(D.subjectID) + subj), regions{loc});

           end
        end
 
    end
end
set(ha,'xtick','','ytick','');
linkaxes(ha,'xy');

%avg over contrast conditions
peakvelchange = nansum(peakvelchange,4); %sum because each entry is already weighted
numTrials = nansum(numTrials,4);
peakvelchange(numTrials==0)=nan;

%avg over subjects
err = nanstd(peakvelchange,[],3)./sqrt(sum(~isnan(peakvelchange),3));
avg = nanmean(peakvelchange,3);


figure;
[h1,h2]=barwitherr(err,avg); hold on;
set(h2,'CapSize',0,'Linewidth',2)
set(gca,'xtick',1:length(regions),'xticklabel',regions);
legend('Left choice','Right choice');
ylabel('contrast-matched change in peak velocity (mm/s)');

%Add t-test for each entry
for loc = 1:length(regions)
    
    for r = 1:2
        [hyp,pval] = ttest( squeeze(peakvelchange(loc,r,:)) );
        
        text(loc + h1(r).XOffset ,10, num2str(pval),'rotation',45);
        
        if pval<0.05
            text( loc + h1(r).XOffset ,9, '*', 'fontsize',20);
        end
    end
end

%% Anova on wheel velocity effect
regions = {'LeftVIS','RightVIS','LeftM2','RightM2','LeftM1','RightM1','LeftS1','RightS1'};
choices = {'Left','Right'};
[contrastCondSet,~,D.contrastCond] = unique(D.stimulus,'rows');
D.laserType = logical(D.laserType);

peakvelchange = nan(length(regions),2,max(D.subjectID),max(D.contrastCond));
numTrials = zeros(size(peakvelchange));

laser_coeff = nan(length(regions),2);
laser_coeff_pval = nan(length(regions),2);
num_tests = 0;
for loc = 1:length(regions)
    for r = 1:2
    	E = getrow(D, D.response == r & (D.laserType==0 | D.laserRegion==regions{loc}) );
        [P,T,STATS,TERMS] = anovan( E.peakVel, [E.stimulus(:,1) E.stimulus(:,2) E.sessionID E.subjectID E.laserType], 'nested',...
                            [0 0 0 0 0;0 0 0 0 0;0 0 0 1 0;0 0 0 0 0;0 0 0 0 0],'varnames',{'CL','CR','sessionID' 'subjectID', 'laser'},...
                            'display','off','continuous',[1 2]);
                               
        laser_coeff(loc,r) = STATS.coeffs(end) - STATS.coeffs(end-1);
        laser_coeff_pval(loc,r) = P(end);
        
        num_tests = num_tests+1;
    end
end

figure;
bx=bar(laser_coeff); hold on;

for loc =  1:length(regions)
    for r = 1:2
        sig = laser_coeff_pval(loc,r);
%         if sig<0.001
%             txt='***';
%         elseif sig<0.01
%            txt='**'; 
%         elseif sig<0.05
%             txt='*';
%         else
%             txt='ns';
%         end
%         tx=text(loc+(r-1.5)/3,11, laser_coeff_pval(loc,r),txt,'HorizontalAlignment','center');
%         
        %corrected for multiple comparisons
        if sig< 0.001/num_tests
            txt='***';
        elseif sig<0.01/num_tests
           txt='**'; 
        elseif sig<0.05/num_tests
            txt='*';
        else
            txt='ns';
        end
        tx=text(loc+(r-1.5)/3,10, laser_coeff_pval(loc,r),txt,'HorizontalAlignment','center');
        
    end
end
set(gca,'xtick',1:length(regions),'xticklabel',regions);
ylabel('Main effect of laser (coefficient) on wheel peak vel');
legend('Left choice','Right choice');

% Try model which compares Left and Right hemi inactivation

%Try a model which includes choices and looks at the interaction between
%choice and inactivation hemisphere
regions = {'VIS','M2','M1','S1'};
for loc = 1:length(regions)
    	E = getrow(D, D.response<3 & contains(string(D.laserRegion),regions{loc} ));
        E.laserHemisphere = contains(string(E.laserRegion),'Left') + 2*contains(string(E.laserRegion),'Right');
[P,T,STATS,TERMS] = anovan( E.peakVel, [E.stimulus(:,1) E.stimulus(:,2) E.sessionID E.subjectID E.response E.laserHemisphere], 'nested',...
                            [0 0 0 0 0 0;...
                             0 0 0 0 0 0;...
                             0 0 0 1 0 0;...
                             0 0 0 0 0 0;...
                             0 0 0 0 0 0;...
                             0 0 0 0 0 0],...
                            'varnames',{'CL','CR','sessionID' 'subjectID', 'choice','laserHemisphere'},...
                            'display','off','continuous',[1 2],'model','interaction');
      sig = P(end);
        if sig<0.001
            txt='***';
        elseif sig<0.01
           txt='**'; 
        elseif sig<0.05
            txt='*';
        else
            txt='ns';
        end                  
    fprintf('%s: interaction(choice,hemisphere) p=%0.2f %s \n',regions{loc},sig,txt);
end

%% USED IN PAPER: Get single trial examples (multi-power expt)
close all;

t = D.wheel_stimulusOn_timesteps(1,:);

f = figure; axes; hold on;
laserOff = 6375;
laserOn = 6608;

% laserOff = NaN;
% laserOn = NaN;

idx = D.laserType==0 & D.response==1 & D.feedbackType==1 & D.subjectID==3;
if isnan(laserOff)
    laserOff = randsample(find(idx),1);
end
v = D.vel(laserOff,:);
plot(t,v,'k-');
line([1 1]*(D.time_choiceMade(laserOff) - D.time_stimulusOn(laserOff)), get(gca,'ylim'),'color','k');
line([1 1]*(D.time_feedback(laserOff) - D.time_stimulusOn(laserOff)), get(gca,'ylim'),'color','k','linestyle','--');

idx = D.laserRegion=='LeftM2' & D.response==1 & D.feedbackType==1 & D.subjectID==3;
if isnan(laserOn)
    laserOn = randsample(find(idx),1);
end
v = D.vel(laserOn,:);
plot(t,v,'r-');
line([1 1]*(D.time_choiceMade(laserOn) - D.time_stimulusOn(laserOn)), get(gca,'ylim'),'color','r');
line([1 1]*(D.time_feedback(laserOn) - D.time_stimulusOn(laserOn)), get(gca,'ylim'),'color','r','linestyle','--');


xlim([0 1.5]);
fprintf('Laser Off: %d, Laser On: %d \n',laserOff,laserOn);

laserOff = NaN;
laserOn = NaN;