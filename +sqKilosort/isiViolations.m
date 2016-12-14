

function isiV = isiViolations(resultsDirectory,varargin)
%isiV = isiViolations(resultsDirectory)
%  computes isi volation rate based on clustering by spike_clusters.npy 
%  (or spike_templates.npy if spike_clusters.npy does not exist)
%  returns measures for each unique cluster/template value, sorted by value
%isiV = isiViolations(resultsDirectory,spike_clusters)
%  if spike_clusters is given, this is used for spike assignments instead
%  of spike_clusters.npy or spike_templates.npy

clu=getparmC(varargin,1,0);

%% Precompute the locationsn of files to be loaded
spikeClustersPath = fullfile(resultsDirectory,'spike_clusters.npy');
spikeTemplatesPath = fullfile(resultsDirectory,'spike_templates.npy');
spikeTimesPath= fullfile(resultsDirectory,'spike_times.npy');
paramsPath= fullfile(resultsDirectory,'params.py');

%% 

refDur = 0.001;
minISI = 0.0002;

fprintf(1, 'loading data for ISI computation\n');
if exist(spikeClustersPath)
    if(isequal(clu,0))
        spike_clusters = readNPY(spikeClustersPath);
        %spike_clusters = spike_clusters + 1; % because in Python indexes start at 0
    else
        spike_clusters=clu;
    end
else
    spike_clusters = readNPY(spikeTemplatesPath);
    %spike_clusters = spike_clusters + 1; % because in Python indexes start at 0
end


spike_times = readNPY(spikeTimesPath);
params = readKSparams(paramsPath);
spike_times = double(spike_times)/params.sample_rate;

fprintf(1, 'computing ISI violations\n');

clusterIDs = unique(spike_clusters);
isiV = zeros(1,numel(clusterIDs));
for c = 1:numel(clusterIDs)
    
    [fpRate, numViolations] = ISIViolations(spike_times(spike_clusters==clusterIDs(c)), minISI, refDur);
    isiV(c) = fpRate;
    nSpikes = sum(spike_clusters==clusterIDs(c));    
    fprintf(1, 'cluster %3d: %d viol (%d spikes), %.2f estimated FP rate\n', clusterIDs(c), numViolations, nSpikes, fpRate);
    
end
