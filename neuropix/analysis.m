%% Load Data, compute SP for neurins in VISp and MOs
data_dir = 'E:\Steinmetz_et_al_2019_Neuropix_Data';
folders = dir('E:\Steinmetz_et_al_2019_Neuropix_Data');  folders(1:2)=[];

load('ctxOutlines.mat','coords');
window = [-0.2 0.5];
psth_binsize = 0.001;
nShuf = 2000;
smoothMsec=15; %15msec

f = figure; hold on;
colormap(BlueWhiteRed); colorbar;
caxis(0.5 + [-1 1]*0.1);
for q = 1:numel(coords) % coords is from ctxOutlines.mat
    cx = coords(q).x;
    cy = coords(q).y;
    hxx=plot(cx,-cy, 'LineWidth', 0.5, 'Color', [1 1 1]*0.8, 'Tag', 'outline');
end
set(gca,'xtick','','ytick','','dataaspectratio',[1 1 1],'xlim',[0 1200]);

for sess = 1:length(folders)
    eRef = folders(sess).name;
    fprintf('Session %d %s\n',sess,eRef);
    
    decodingFile = [ './decoding/' eRef '.mat'];
    if ~exist(decodingFile,'file')

        data_dir = fullfile('E:\Steinmetz_et_al_2019_Neuropix_Data',eRef);
        
        %Load data
        d = loadSession(data_dir);
        d.trials = getrow(d.trials, d.trials.included==1);%exclude trials (repeatNum>1, etc);
        
        %Only use good clusters from VIS and MOs
        channelRegions = cellstr(d.channels.brainLocation.allen_ontology);
        channelProbe = d.channels.probe;
        
        
        d.clusters.clusterRegions = channelRegions(d.clusters.peakChannel);
        d.clusters.clusterCoord = [d.channels.brainLocation.ccf_lr(d.clusters.peakChannel),...
            d.channels.brainLocation.ccf_ap(d.clusters.peakChannel)];
        goodClusterIdx = d.clusters.x_phy_annotation == 3 &...
            (strcmp(d.clusters.clusterRegions,'VISp') | strcmp(d.clusters.clusterRegions,'MOs'));
        goodClusters = d.clusters.originalIDs(goodClusterIdx);
        
        
        %Compute stimulus probability from spike counts
        %Split by choice
        [conds,~,trial_condition]=unique(d.trials.response_choice,'rows');
        
        %decoded variable: whether there was a stimulus on the contralateral side
        dec = d.trials.visualStim_contrastLeft>0;
        
        %ensure a large enough number of trials to compute the statistic
        q = crosstab(trial_condition, dec);
        nComp = sum(q(:,1).*q(:,2));
        if nComp<10
            error('too few trials');
        end
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
            %         shufLabels{c} = (1:nA)';
        end
        
        bins_all = cell(length(goodClusters),1);
        SP_all = cell(length(goodClusters),1);
        SP_p_all = cell(length(goodClusters),1);
        cluster_location = nan(length(goodClusters),2);
        
        for cluster = 1:length(goodClusters)
            fprintf('Cluster %d/%d\n',cluster,length(goodClusters));
            cluID = goodClusters(cluster);
            cluIdx = d.spikes.clusters==cluID;
            [psth, bins, ~, ~, ~, binnedArray] = psthAndBA(d.spikes.times(cluIdx), d.trials.visualStim_times, window, psth_binsize);
            
%             %smooth binnedArray
%             gw = gausswin(round(smoothMsec*6),3);
%             smWin = gw./sum(gw);
%             binnedArraySmoothed = conv2(smWin,1,binnedArray', 'same')'./psth_binsize;
            
            SP = nan(length(bins),1);
            SP_p = nan(size(SP));
            for t = 1:length(bins)
                activity = binnedArray(:,t);
                
                shuf_results =  choiceProbShuf(activity, dec, trial_condition, shufLabels);
                SP(t) = shuf_results(1);
                SP_p(t) = shuf_results(2);
            end
            
            
            %if significant SP, plot avg SP value over all times when it's
            %significant
            %         if any(SP_p>0.99 & SP>0.5)
            %             h=scatter(d.clusters.clusterCoord(cluster,1),d.clusters.clusterCoord(cluster,2),200,'k','o','filled');
            %             h.MarkerEdgeColor=[1 1 1]*0.75;
            %             h.SizeData=30;
            %             h.CData = max(SP);
            %             drawnow;
            %         end
            cluster_location(cluster,:) = d.clusters.clusterCoord(cluster,:);
            SP_all{cluster} = SP';
            SP_p_all{cluster} = SP_p_all';
        end
        
        bins_all = cell2num(bins_all);
        save(decodingFile,'cluster_location','bins','SP_all','SP_p_all');
    end
end

% figure;
% plot(bins, stimulus_probability)




%% Plot dorsal view of ctx, with dots showing the SP measure

%% Plot SP varying by cortical depth

