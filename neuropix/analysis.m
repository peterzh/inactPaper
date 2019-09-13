data_dir = 'E:\Steinmetz_et_al_2019_Neuropix_Data\Cori_2016-12-14';

%% Load Data, get neurons in MOs and VIS

%Behavioural data & timings (NOT PASSIVE SESSION)
b = struct;
b.contrastLeft = readNPY(fullfile(data_dir,'trials.visualStim_contrastLeft.npy'));
b.contrastRight = readNPY(fullfile(data_dir,'trials.visualStim_contrastRight.npy'));
b.choice = readNPY(fullfile(data_dir,'trials.response_choice.npy')); %-1 0 +1 coding
b.included = readNPY(fullfile(data_dir,'trials.included.npy'));
b.time_stimOnset = readNPY(fullfile(data_dir,'trials.visualStim_times.npy'));
b.time_goCue = readNPY(fullfile(data_dir,'trials.goCue_times.npy'));
b = getrow(b, b.included==1); %exclude trials (repeatNum>1, etc);

%Neural data
n = struct;
n.spikeTimes = readNPY(fullfile(data_dir,'spikes.times.npy'));
n.clusterID = readNPY(fullfile(data_dir,'spikes.clusters.npy'));
n.clusterType= readNPY(fullfile(data_dir,'clusters._phy_annotation.npy'));
n.clusterDepth = readNPY(fullfile(data_dir,'spikes.depths.npy'));

%get only good clusters
goodClusters = find(n.clusterType==3);
% n.brainRegion = dlmread(fullfile(data_dir,'channels.brainLocation.tsv'),'\t',1,0);

%% Plot rasters & PSTH aligned to stimulus onset and separated by contraStimulus

%% Compute spike counts in bins relative to stim onset

%% Compute stimulus probability from the spike counts

%% Plot dorsal view of ctx, with dots showing the SP measure

%% Plot SP varying by cortical depth

