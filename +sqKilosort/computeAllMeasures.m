

function [cgs, uQ, cR, isiV] = computeAllMeasures(resultsDirectory,varargin)
%[cgs, uQ, cR, isiV] = computeAllMeasures(resultsDirectory)
%  computes measures based on clustering by spike_clusters.npy 
%  (or spike_templates.npy if spike_clusters.npy does not exist)
%  returns measures for each unique cluster/template value, sorted by value
%[cgs, uQ, cR, isiV] = computeAllMeasures(resultsDirectory,spike_clusters)
%  if spike_clusters is given, this is used for spike assignments instead
%  of spike_clusters.npy or spike_templates.npy

clu=getparmC(varargin,1,0);

clusterPath = fullfile(resultsDirectory, 'cluster_groups.tsv');
spikeClustersPath = fullfile(resultsDirectory,'spike_clusters.npy');
spikeTemplatesPath = fullfile(resultsDirectory,'spike_templates.npy');


if exist(clusterPath, 'file')
[cids, cgs] = readClusterGroupsCSV(clusterPath);
elseif exist(spikeClustersPath, 'file')
    clu2 = readNPY(spikeClustersPath);
    cgs = 3*ones(size(unique(clu2))); % all unsorted
else
    clu2 = readNPY(spikeTemplatesPath);
    cgs = 3*ones(size(unique(clu2))); % all unsorted
end
if ~isequal(clu,0)
    cids_from_clu=unique(clu);
    [~,i1]=ismember(cids_from_clu,cids);
    cgs=cgs(i1);
end
[cids2, uQ, cR] = sqKilosort.maskedClusterQuality(resultsDirectory,varargin{:});

isiV = sqKilosort.isiViolations(resultsDirectory,varargin{:});

