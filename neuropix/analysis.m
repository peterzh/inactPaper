%% Load Data, compute decoding for neurons in iso cortex
clear all; close all;
data_dir = 'E:\Steinmetz_et_al_2019_Neuropix_Data';
folders = dir('E:\Steinmetz_et_al_2019_Neuropix_Data');  folders(1:2)=[];
ctx_acronyms = load('ctx_acronyms.mat','ctx_acronyms');ctx_acronyms=ctx_acronyms.ctx_acronyms;

Stimwindow = [-0.2 0.4];
Movewindow = [-0.4 0.4];
psth_binsize = 0.001;
% psth_binsize = 0.01;

nShuf = 2000;

%define smoothing filter applied to the binned array
smoothMsec=25; %msec to smooth binned array
smoothFilt = myGaussWin(smoothMsec/1000, 1/psth_binsize); %smooth with causal filter
smoothFilt(1:round(numel(smoothFilt)/2)-1) = 0; smoothFilt = smoothFilt./sum(smoothFilt);
            
multiWaitbar( 'CloseAll' );

for sess = 1:length(folders)
    multiWaitbar('Sessions',sess/length(folders));
        
    eRef = folders(sess).name;
    fprintf('Session %d %s\n',sess,eRef);
    
    decodingFile = [ './decoding/' eRef '.mat'];
    if ~exist(decodingFile,'file')
        data_dir = fullfile('E:\Steinmetz_et_al_2019_Neuropix_Data',eRef);
        
        %Load data
        d = loadSession(data_dir);

         [cweA, cwtA] = ...
        loadSelectedTrials(eRef(1:end-11), eRef(end-9:end));
            cwtA = cwtA(cweA.inclT==1,:);
            cweA = cweA(cweA.inclT==1,:);
            
        %first construct a new "choiceInWindow" that's specifically about what "choice" the mouse makes during the analysis window. 
        %this 'choice' may be different from the choice made at the end of
        %the window
        cweA.choiceInWindow = double(cweA.hasLeftMove);  cweA.choiceInWindow(cweA.hasRightMove)=2;  cweA.choiceInWindow(cweA.hasNoMoveInWin)=3;    

        %Get clusters at cortical regions
        includedUnits = getIncludedUnits(eRef);
        channelRegions = cellstr(d.channels.brainLocation.allen_ontology);
        d.clusters.clusterRegions = channelRegions(d.clusters.peakChannel);
        d.clusters.clusterCoord = [d.channels.brainLocation.ccf_lr(d.clusters.peakChannel),...
        d.channels.brainLocation.ccf_ap(d.clusters.peakChannel)];
        goodClusterIDs = find( includedUnits & d.clusters.x_phy_annotation>1 & ismember(d.clusters.clusterRegions,ctx_acronyms) ) - 1;
        
        %Compute binned array for each cluster
        binnedArraySmoothed_stimAligned = cell(length(goodClusterIDs),1);
        binnedArraySmoothed_moveAligned = cell(length(goodClusterIDs),1);
        cluster_location = nan(length(goodClusterIDs),2);
        cluster_region = cell(length(goodClusterIDs),1);
        cluster_depth = nan(length(goodClusterIDs),1);
        fprintf('\tComputing binnedArray for each cluster\n');
        for clu = 1:length(goodClusterIDs)
            multiWaitbar('Cluster Binned Array',clu/length(goodClusterIDs));

            cluID = goodClusterIDs(clu);
            fprintf('\tCluster %d/%d [ID=%d] %s\n',clu,length(goodClusterIDs),cluID,d.clusters.clusterRegions{cluID+1});
            
            %Compute stim-aligned binned array
            [~, bins_stim, ~, ~, ~, binnedArray] = psthAndBA(d.spikes.times(d.spikes.clusters==cluID), cwtA.stimOn, Stimwindow, psth_binsize);
            binnedArraySmoothed_stimAligned{clu} = conv2(smoothFilt,1,binnedArray', 'same')'./psth_binsize;

            %Compute move-aligned binned array & smooth with same filter
            [~, bins_move, ~, ~, ~, binnedArray] = psthAndBA(d.spikes.times(d.spikes.clusters==cluID), cwtA.moveTimeAbs, Movewindow, psth_binsize);
            binnedArraySmoothed_moveAligned{clu} = conv2(smoothFilt,1,binnedArray', 'same')'./psth_binsize;
   
            cluster_location(clu,:) = d.clusters.clusterCoord(cluID+1,:);
            cluster_depth(clu) = d.clusters.depths(cluID+1);
            cluster_region{clu} = d.clusters.clusterRegions{cluID+1};
        end
        multiWaitbar('Cluster Binned Array','Close');

        
        %Create decoding variables
        labels = {'detect prob','stimL prob','stimR prob','stimR prob NOGO trials','choice prob'};    
         
        decoded_variable = {cweA.choiceInWindow<3,...   %DP
            cweA.contrastLeft>0,...             %SLP
            cweA.contrastRight>0,...            %SRP
            cweA.contrastRight(cweA.choiceInWindow==3)>0,... %SRP NOGO
            cweA.choiceInWindow(cweA.choiceInWindow<3)==1};     %CP
        
        splitting_condition = {[cweA.contrastLeft cweA.contrastRight],...
            [cweA.choiceInWindow, cweA.contrastRight],...
            [cweA.choiceInWindow, cweA.contrastLeft],...
            cweA.contrastLeft(cweA.choiceInWindow==3),...
            [cweA.contrastLeft(cweA.choiceInWindow<3) cweA.contrastRight(cweA.choiceInWindow<3)]};
        
        decoding_prob_all_stimAligned = nan(length(bins_stim),length(goodClusterIDs),length(labels));
        decoding_prob_p_all_stimAligned = nan(length(bins_stim),length(goodClusterIDs),length(labels));
        decoding_prob_all_moveAligned = nan(length(bins_move),length(goodClusterIDs),length(labels));
        decoding_prob_p_all_moveAligned = nan(length(bins_move),length(goodClusterIDs),length(labels));
        for p = 1:length(labels) %for each type of decoding
            multiWaitbar('Decoded variable',p/length(labels));

            fprintf('\tComputing %s\n',labels{p});
            
            [conds,~,trial_condition]=unique(splitting_condition{p},'rows');
            dec = decoded_variable{p};
            
            %ensure a large enough number of trials to compute the statistic
            q = crosstab(trial_condition, dec);
            if size(q,2)>1 %if dec has more than 2 values
                nComp = sum(q(:,1).*q(:,2));
            else
                nComp = 0;
            end
            if nComp>10
                
                %create shuffle labels
                shufLabels = cell(1,max(trial_condition));
                for c = 1:max(trial_condition)
                    idx = trial_condition==c;
                    chA = dec & idx;
                    nA = sum(chA);
                    chB = ~dec & idx;
                    nB = sum(chB);
                    q = arrayfun(@(x)randperm(nA+nB,nA), 1:nShuf, 'uni', false);
                    shufLabels{c} = vertcat(q{:})';
                    shufLabels{c} = [(1:nA)' shufLabels{c}];
                end
                %
                for clu = 1:length(goodClusterIDs)
                    multiWaitbar('Cluster',clu/length(goodClusterIDs));

                    %compute stim-aligned decoding
                    for t = 1:length(bins_stim) %for each timebin

                        %get the activity for this cluster
                        if p==5 %CP exclude nogo trials
                            activity = binnedArraySmoothed_stimAligned{clu}(cweA.choiceInWindow<3,t);
                        elseif p==4 %SRP NOGO excludes go trials
                            activity = binnedArraySmoothed_stimAligned{clu}(cweA.choiceInWindow==3,t);
                        else
                            activity = binnedArraySmoothed_stimAligned{clu}(:,t);
                        end
                        
                        %add tiny delta to activity to resolve ties
                        activity = activity + randn(size(activity))*10^-7;
                        
                        shuf_results =  choiceProbShuf(activity, dec, trial_condition, shufLabels);
                        decoding_prob_all_stimAligned(t,clu,p) = shuf_results(1);
                        decoding_prob_p_all_stimAligned(t,clu,p) = shuf_results(2);
                    end
                    
                    %compute move-aligned decoding
                    for t = 1:length(bins_move) %for each timebin
                        %get the activity for this cluster
                        if p==5 %CP exclude nogo trials
                            activity = binnedArraySmoothed_moveAligned{clu}(cweA.choiceInWindow<3,t);
                        elseif p==4 %SRP NOGO excludes go trials
                            activity = binnedArraySmoothed_moveAligned{clu}(cweA.choiceInWindow==3,t);
                        else
                            activity = binnedArraySmoothed_moveAligned{clu}(:,t);
                        end
                        
                        %add tiny delta to activity to resolve ties
                        activity = activity + randn(size(activity))*10^-7;
                        
                        shuf_results =  choiceProbShuf(activity, dec, trial_condition, shufLabels);
                        decoding_prob_all_moveAligned(t,clu,p) = shuf_results(1);
                        decoding_prob_p_all_moveAligned(t,clu,p) = shuf_results(2);
                    end
                    
                end
            else
                warning('too few trials to compute this decoding');
            end
            
        end
        
        %correct for p values at 0 or 1
        decoding_prob_p_all_stimAligned(decoding_prob_p_all_stimAligned==0) = (nShuf-0.5)/nShuf;
        decoding_prob_p_all_stimAligned(decoding_prob_p_all_stimAligned==1) = 1 - (nShuf-0.5)/nShuf;
        decoding_prob_p_all_moveAligned(decoding_prob_p_all_moveAligned==0) = (nShuf-0.5)/nShuf;
        decoding_prob_p_all_moveAligned(decoding_prob_p_all_moveAligned==1) = 1 - (nShuf-0.5)/nShuf;
        
        save(decodingFile,'goodClusterIDs','cluster_location','cluster_depth','cluster_region',...
                        'bins_stim','bins_move','decoding_prob_all_stimAligned','decoding_prob_p_all_stimAligned',...
                        'decoding_prob_all_moveAligned','decoding_prob_p_all_moveAligned','labels');
    end
    

end

%% Count number of custer

total = 0;
inCortex = 0;
GoodPhy = 0;
GoodResponses = 0;


for sess = 1:length(folders)
    multiWaitbar('Sessions',sess/length(folders));
    
    eRef = folders(sess).name;
    fprintf('Session %d %s\n',sess,eRef);
    
        data_dir = fullfile('E:\Steinmetz_et_al_2019_Neuropix_Data',eRef);
        
        %Load data
        d = loadSession(data_dir);
        
        total = total + length(d.clusters.x_phy_annotation);
%   
%         %Get clusters at cortical regions
        includedUnits = getIncludedUnits(eRef);
        channelRegions = cellstr(d.channels.brainLocation.allen_ontology);
        d.clusters.clusterRegions = channelRegions(d.clusters.peakChannel);
        d.clusters.clusterCoord = [d.channels.brainLocation.ccf_lr(d.clusters.peakChannel),...
            d.channels.brainLocation.ccf_ap(d.clusters.peakChannel)];
        
        inCortex = inCortex + sum(ismember(d.clusters.clusterRegions,ctx_acronyms));
        GoodPhy = GoodPhy + sum( d.clusters.x_phy_annotation>1 & ismember(d.clusters.clusterRegions,ctx_acronyms) );
        GoodResponses = GoodResponses + sum(includedUnits & d.clusters.x_phy_annotation>1 & ismember(d.clusters.clusterRegions,ctx_acronyms));  
    
end

%% Plot dorsal view of ctx, with dots showing the SP measure
clear all; close all;

%define smoothing filter applied to the decoding probability due to the tie-break introducing artificial variation
smoothMsec=25; %msec to smooth binned array
psth_binsize=0.001;
smoothFilt = myGaussWin(smoothMsec/1000, 1/psth_binsize); %smooth with causal filter
smoothFilt(1:round(numel(smoothFilt)/2)-1) = 0;  %Truncate before 0 sec, to generate causal filter
smoothFilt = smoothFilt./sum(smoothFilt);

data_dir = 'E:\Steinmetz_et_al_2019_Neuropix_Data';
folders = dir('E:\Steinmetz_et_al_2019_Neuropix_Data');  folders(1:2)=[];

load('ctxOutlines.mat','coords');
highlighted_regions = {'VISp',{'VISa','VISam','VISl','VISpm','VISrl'},'MOs','MOp','SSp'}; %regions to show traces
highlighted_regions_label = {'VISp','VIS sec','MOs','MOp','SSp'};
areaCols = [      0    0.4470    0.7410; %blue
        0.3010    0.7450    0.9330; %lightblue
    0.4660    0.6740    0.1880; %green
    0.8500    0.3250    0.0980; %orange
    0.4940    0.1840    0.5560]; %purple

alpha = 0.05; %significance threshold for showing clusters

%Collate all clusters
decoding_prob_stimAligned = cell(1,length(folders));
decoding_prob_p_stimAligned = cell(1,length(folders));
decoding_prob_moveAligned = cell(1,length(folders));
decoding_prob_p_moveAligned = cell(1,length(folders));
cluster_loc = cell(1,length(folders));
cluster_depth = cell(1,length(folders));
cluster_region = cell(1,length(folders));
cluster_session = cell(1,length(folders));
cluster_ids = cell(1,length(folders));
for sess = 1:length(folders)

    eRef = folders(sess).name;
    fprintf('Session %d %s\n',sess,eRef);
    
    decodingFile = [ './decoding/' eRef '.mat'];
    d = load(decodingFile);
    
    if ~isempty(d.decoding_prob_all_stimAligned)
        decoding_prob_stimAligned{sess} = d.decoding_prob_all_stimAligned;
        decoding_prob_p_stimAligned{sess} = d.decoding_prob_p_all_stimAligned;
        decoding_prob_moveAligned{sess} = d.decoding_prob_all_moveAligned;
        decoding_prob_p_moveAligned{sess} = d.decoding_prob_p_all_moveAligned;
        cluster_loc{sess} = d.cluster_location;
        cluster_depth{sess} = d.cluster_depth;
        cluster_region{sess} = d.cluster_region;
        cluster_session{sess} = ones(size(d.cluster_location,1),1)*sess;
        cluster_ids{sess} = d.goodClusterIDs;
    end
end
decoding_prob_stimAligned = cat(2,decoding_prob_stimAligned{:});
decoding_prob_p_stimAligned = cat(2,decoding_prob_p_stimAligned{:});
decoding_prob_moveAligned = cat(2,decoding_prob_moveAligned{:});
decoding_prob_p_moveAligned = cat(2,decoding_prob_p_moveAligned{:});
cluster_loc = cat(1,cluster_loc{:});
cluster_depth = cat(1,cluster_depth{:});
cluster_region = cat(1,cluster_region{:});
cluster_session = cat(1,cluster_session{:});
cluster_ids = cat(1,cluster_ids{:});


bins_stim = d.bins_stim;
bins_move = d.bins_move;
labels = d.labels;

%add jitter to cluster_loc
cluster_loc = cluster_loc + 100*randn(size(cluster_loc));

%Plot trace of significant decoding
f=figure('color','w','position',[130 313 1623 593]);
for p = 1:5
    
    %plot stim aligned
    subplot(2,5,p); title(labels{p}); hold on;
    
    %Plot % of significant neurons in different areas
    for loc = 1:length(highlighted_regions)
        idx = ismember(cluster_region,highlighted_regions{loc}) & ~isnan(decoding_prob_stimAligned(1,:,p)');
        sig = decoding_prob_p_stimAligned(:,idx,p)< (alpha/2) | decoding_prob_p_stimAligned(:,idx,p)>(1-alpha/2);
%         plot(bins_stim, 100*mean(sig,2),'color',areaCols(loc,:));
        plot(bins_stim,conv2(smoothFilt,1,100*mean(sig,2), 'same'),'color',areaCols(loc,:),'linewidth',2);
    end
    line([min(bins_stim) max(bins_stim)],[1 1]*alpha*100,'LineStyle','--','Color',[1 1 1]*0.5);
    set(gca,'xlim',[min(bins_stim)+2*smoothMsec/1000 max(bins_stim)],'ylim',[0 40]);
    
    
    if p == 1
        legend(gca,highlighted_regions_label,'Location','NorthWest','EdgeColor',[1 1 1]);
        xlabel('Stim onset');
        ylabel(gca,'% neurons with significant decoding');
    end
    
     %plot move aligned
    subplot(2,5,5+p); hold on;
    
    %Plot % of significant neurons in different areas
    for loc = 1:length(highlighted_regions)
        idx = ismember(cluster_region,highlighted_regions{loc}) & ~isnan(decoding_prob_moveAligned(1,:,p)');
        sig = decoding_prob_p_moveAligned(:,idx,p)< (alpha/2) | decoding_prob_p_moveAligned(:,idx,p)>(1-alpha/2);
%         plot(bins_move, 100*mean(sig,2),'color',areaCols(loc,:));
        plot(bins_move,conv2(smoothFilt,1,100*mean(sig,2), 'same'),'color',areaCols(loc,:),'linewidth',2);
    end
    line([min(bins_move) max(bins_move)],[1 1]*alpha*100,'LineStyle','--','Color',[1 1 1]*0.5);
    set(gca,'xlim',[min(bins_move)+2*smoothMsec/1000 max(bins_move)],'ylim',[0 40]);
    
    if p==1
        xlabel('Move onset');
        ylabel(gca,'% neurons with significant decoding');
    end
end
set(f,'Units','Inches','renderer','painters');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f, fullfile('C:\Users\Peter\OneDrive - University College London\Zatka-Haas et al\working\new_analyses\neuropix_decoding',['summary.pdf']),'-dpdf','-r0');

%Plot trace of SL vs SR
f=figure('color','w','position',[130 649 1623 257]);
for loc = 1:length(highlighted_regions)
    subplot(1,5,loc); title(highlighted_regions_label{loc}); hold on;
    
    idx = ismember(cluster_region,highlighted_regions{loc});
    p=2;
    sigL = decoding_prob_p_stimAligned(:,idx,p)< (alpha/2) | decoding_prob_p_stimAligned(:,idx,p)>(1-alpha/2);
    p=3;
    sigR = decoding_prob_p_stimAligned(:,idx,p)< (alpha/2) | decoding_prob_p_stimAligned(:,idx,p)>(1-alpha/2);
%     plot(bins_stim, 100*mean(sigL,2),'k-',bins_stim, 100*mean(sigR,2),'r-');

    plot(bins_stim,conv2(smoothFilt,1,100*mean(sigL,2), 'same'),'k-',bins_stim,conv2(smoothFilt,1,100*mean(sigR,2), 'same'),'r-','linewidth',2);
    
    line([min(bins_stim) max(bins_stim)],[1 1]*alpha*100,'LineStyle','--','Color',[1 1 1]*0.5);
    set(gca,'xlim',[min(bins_stim)+2*smoothMsec/1000 max(bins_stim)],'ylim',[0 25]);
    ylabel(gca,'% neurons with significant decoding');
    
    if loc == 1
        legend(gca,{'CL>0','CR>0'});
    end
end
set(f,'Units','Inches','renderer','painters');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f, fullfile('C:\Users\Peter\OneDrive - University College London\Zatka-Haas et al\working\new_analyses\neuropix_decoding',['SL_vs_SR.pdf']),'-dpdf','-r0');



%Plot maps
t_slices = [0 0.1 0.15 0.2 0.25 0.3];
alpha = 0.001;

f=figure('position',[3 81 1297 873],'color','w'); %plot outlines of ctx
ha_maps = tight_subplot(5, length(t_slices),0.02,0.05,0.05); for i = 1:numel(ha_maps); hold(ha_maps(i),'on'); end;

for p = 1:5
    
    for t = 1:length(t_slices)
        thisHa = ha_maps(length(t_slices)*(p-1) + t);
        
        for q = 1:numel(coords) % coords is from ctxOutlines.mat
            cx = coords(q).x*10;
            cy = coords(q).y*10;
            plot(thisHa,cx,cy, 'LineWidth', 0.5, 'Color', [1 1 1]*0.8, 'Tag', 'outline');
        end
        colormap(thisHa,BlueWhiteRed);
        caxis(thisHa,0.5 + [-1 1]*0.4);
        
        if p == 1
            title(thisHa,t_slices(t));
        end
        
        if t == 1
            ylabel(thisHa, labels{p});
        end
       
        %Plot decodingc
        time_idx = find(bins_stim >= t_slices(t),1,'first');
        
        sigIdx = decoding_prob_p_stimAligned(time_idx,:,p)>(1-alpha/2) | decoding_prob_p_stimAligned(time_idx,:,p)<(alpha/2);
        
        %plot non-significant clusters as black dots
        h=plot(thisHa, cluster_loc(~sigIdx,1),cluster_loc(~sigIdx,2), 'k.','color',[1 1 1]*0.7);
        h.MarkerSize=2;
        
        %plot significant clusters as coloured dots
        h=scatter(thisHa,cluster_loc(sigIdx,1),cluster_loc(sigIdx,2),10,'k','o','filled');
        h.CData = decoding_prob_stimAligned(time_idx,sigIdx,p);
        %                 h.SizeData = 1 + 10*(d.decoding_prob_p_all(time_idx,:,p)>0.99 | d.decoding_prob_p_all(time_idx,:,p)<0.01);
        h.SizeData = 5;
        h.MarkerFaceAlpha=1;
    end

end

colorbar(ha_maps(end),'Location','EastOutside');
set(ha_maps,'ydir','reverse','xtick','','ytick','','dataaspectratio',[1 1 1],'xcolor','none','ycolor','none');
set(ha_maps(1:length(t_slices):end),'ycolor','k');
set(ha_maps, 'xlim', [290.6977 1.1105e+04], 'ylim', [1.6239e+03 1.0486e+04]);

set(f,'Units','Inches','renderer','painters');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f, fullfile('C:\Users\Peter\OneDrive - University College London\Zatka-Haas et al\working\new_analyses\neuropix_decoding',['MAPS_stimAligned.pdf']),'-dpdf','-r0')


%% Plot movies

alpha = 0.001;
labels = {'Detect Probability','StimLeft Probability','StimRight Probability','StimRight Probability (NoGo trials)','Choice Probability'};
for p = 1:5
    
    %Plot stim-aligned movie
    f = figure('color','w','name',labels{p},'position',[400 306 927 625]);
    set(f,'renderer','painters');

    ha = tight_subplot(1,2);
    
    %Plot decoding traces
    title(ha(1), labels{p}); hold(ha(1),'on');
    
    %Plot % of significant neurons in different areas
    for loc = 1:length(highlighted_regions)
        idx = ismember(cluster_region,highlighted_regions{loc}) & ~isnan(decoding_prob_stimAligned(1,:,p)');
        sig = decoding_prob_p_stimAligned(:,idx,p)< (alpha/2) | decoding_prob_p_stimAligned(:,idx,p)>(1-alpha/2);
        plot(ha(1),bins_stim,conv2(smoothFilt,1,100*mean(sig,2), 'same'),'color',areaCols(loc,:),'linewidth',2);
    end
    line(ha(1),[min(bins_stim) max(bins_stim)],[1 1]*alpha*100,'LineStyle','--','Color',[1 1 1]*0.5);
    set(ha(1),'xlim',[-0.1 max(bins_stim)],'ylim',[0 25]);
    time_line =  line(ha(1),[0 0],[0 100],'LineStyle','-','Color',[1 1 1]*0.5);
    xlabel(ha(1),'Time from stimulus onset (sec)');
    ylabel(ha(1),'% of significant neurons');
    legend(ha(1),highlighted_regions_label,'box','off');
    
    %Plot map
    set(ha(2),'xlim',[253.2468 5.9805e+03],'ylim',[1000 11000],...
        'ydir','reverse','xcolor','w','ycolor','w',...
        'dataaspectratio',[1 1 1]);
%     title(ha(2), labels{p});
    hold(ha(2),'on');
    for q = 1:numel(coords) % coords is from ctxOutlines.mat
        cx = coords(q).x*10;
        cy = coords(q).y*10;
        plot(ha(2),cx,cy, 'LineWidth', 0.5, 'Color', [1 1 1]*0.8, 'Tag', 'outline');
    end
    colormap(ha(2),BlueWhiteRed);
    caxis(ha(2),0.5 + [-1 1]*0.4);
    
    %Plot dot for every neuron
    plot(ha(2), cluster_loc(:,1),cluster_loc(:,2), 'k.','color',[1 1 1]*0.7,'markersize',2);
    
    %go through time slices and plot significant neurons
    h=scatter(ha(2),cluster_loc(:,1),cluster_loc(:,2),10,'k','o','filled');    
    c=colorbar(ha(2));
    c.Label.String = labels{p};
    
    % Initialize video
    myVideo = VideoWriter(fullfile('C:\Users\Peter\OneDrive - University College London\Zatka-Haas et al\working\new_analyses\neuropix_decoding',[labels{p} '_movie_stimAligned']),...
        'MPEG-4'); %open video file
    myVideo.FrameRate = 20;  %can adjust this, 5 - 10 works well for me
    open(myVideo)

    for t = 100:3:length(bins_stim)
        time_line.XData = [1 1]*bins_stim(t);
        sigIdx = decoding_prob_p_stimAligned(t,:,p)>(1-alpha/2) | decoding_prob_p_stimAligned(t,:,p)<(alpha/2);
        sigIdx = sigIdx & (decoding_prob_p_stimAligned(t-1,:,p)>(1-alpha/2) | decoding_prob_p_stimAligned(t-1,:,p)<(alpha/2));
        sigIdx = sigIdx & (decoding_prob_p_stimAligned(t-2,:,p)>(1-alpha/2) | decoding_prob_p_stimAligned(t-2,:,p)<(alpha/2));
        
        h.CData = decoding_prob_stimAligned(t,:,p);
        
        %move non-sig clusters off the map
        h.XData = cluster_loc(:,1) + 50000*~sigIdx';
        
        pause(0.01);
        frame = getframe(f); %get frame
        writeVideo(myVideo, frame);
    end
    close(myVideo);
    close(f);
    
    
    
    
    %Plot move-aligned movie
    f = figure('color','w','name',labels{p},'position',[400 306 927 625]);
    set(f,'renderer','painters');

    ha = tight_subplot(1,2);
    
    %Plot decoding traces
    title(ha(1), labels{p}); hold(ha(1),'on');
    
    %Plot % of significant neurons in different areas
    for loc = 1:length(highlighted_regions)
        idx = ismember(cluster_region,highlighted_regions{loc}) & ~isnan(decoding_prob_stimAligned(1,:,p)');
        sig = decoding_prob_p_moveAligned(:,idx,p)< (alpha/2) | decoding_prob_p_moveAligned(:,idx,p)>(1-alpha/2);
        plot(ha(1),bins_move,conv2(smoothFilt,1,100*mean(sig,2), 'same'),'color',areaCols(loc,:),'linewidth',2);
    end
    line(ha(1),[min(bins_move) max(bins_move)],[1 1]*alpha*100,'LineStyle','--','Color',[1 1 1]*0.5);
    set(ha(1),'xlim',[-0.3 max(bins_move)],'ylim',[0 25]);
    time_line =  line(ha(1),[0 0],[0 100],'LineStyle','-','Color',[1 1 1]*0.5);
    xlabel(ha(1),'Time from movement (sec)');
    ylabel(ha(1),'% of significant neurons');
    legend(ha(1),highlighted_regions_label,'box','off');
    
    %Plot map
    set(ha(2),'xlim',[253.2468 5.9805e+03],'ylim',[1000 11000],...
        'ydir','reverse','xcolor','w','ycolor','w',...
        'dataaspectratio',[1 1 1]);
%     title(ha(2), labels{p});
    hold(ha(2),'on');
    for q = 1:numel(coords) % coords is from ctxOutlines.mat
        cx = coords(q).x*10;
        cy = coords(q).y*10;
        plot(ha(2),cx,cy, 'LineWidth', 0.5, 'Color', [1 1 1]*0.8, 'Tag', 'outline');
    end
    colormap(ha(2),BlueWhiteRed);
    caxis(ha(2),0.5 + [-1 1]*0.4);
    
    %Plot dot for every neuron
    plot(ha(2), cluster_loc(:,1),cluster_loc(:,2), 'k.','color',[1 1 1]*0.7,'markersize',2);
    
    %go through time slices and plot significant neurons
    h=scatter(ha(2),cluster_loc(:,1),cluster_loc(:,2),10,'k','o','filled');    
    c=colorbar(ha(2));
    c.Label.String = labels{p};
    
    % Initialize video
    myVideo = VideoWriter(fullfile('C:\Users\Peter\OneDrive - University College London\Zatka-Haas et al\working\new_analyses\neuropix_decoding',[labels{p} '_movie_moveAligned']),...
        'MPEG-4'); %open video file
    myVideo.FrameRate = 20;  %can adjust this, 5 - 10 works well for me
    open(myVideo)

    for t = 100:3:length(bins_move)
        time_line.XData = [1 1]*bins_move(t);
        sigIdx = decoding_prob_p_moveAligned(t,:,p)>(1-alpha/2) | decoding_prob_p_moveAligned(t,:,p)<(alpha/2);
        sigIdx = sigIdx & (decoding_prob_p_moveAligned(t-1,:,p)>(1-alpha/2) | decoding_prob_p_moveAligned(t-1,:,p)<(alpha/2));
        sigIdx = sigIdx & (decoding_prob_p_moveAligned(t-2,:,p)>(1-alpha/2) | decoding_prob_p_moveAligned(t-2,:,p)<(alpha/2));
        
        h.CData = decoding_prob_moveAligned(t,:,p);
        
        %move non-sig clusters off the map
        h.XData = cluster_loc(:,1) + 50000*~sigIdx';
        
        pause(0.01);
        frame = getframe(f); %get frame
        writeVideo(myVideo, frame);
    end
    close(myVideo);
    close(f);
end

