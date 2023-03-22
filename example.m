% This is an example of how to (1) calculated temporal clustering coefficient (also named temporal correlation coefficient) using Matlab as described in: 
% "Long Y et al. Evaluating test–retest reliability and sex‐/age‐related effects on temporal clustering coefficient of dynamic functional brain networks[J]. Human Brain Mapping, 2023." 
% and (2) how to calculate the characteristic temporal path length as described in: 
% "Long Y et al. Altered resting-state dynamic functional brain networks in major depressive disorder: Findings from the REST-meta-MDD consortium[J]. NeuroImage: Clinical, 2020, 26: 102163."
% Before calculating, you have to:
% 1. Install Matlab and the Brain Connectivity Toolbox (http://www.brain-connectivitytoolbox.net) 
% (Reference: Rubinov, M., Sporns, O., 2010. Complex network measures of brain connectivity: Uses and interpretations. Neuroimage 52, 1059–1069)
% 2. Add the functions "temporal_clustering.m" and "temporal_path_length.m" which can be found at: https://github.com/Yicheng-Long/dynamic_graph_metrics) to Matlab. 
% The function was adapted from codes in another publicly-available MATLAB toolbox.
% (Reference: Sizemore, A.E., Bassett, D.S., 2018. Dynamic graph metrics: Tutorial, toolbox, and tale.Neuroimage 180, 417–427)
% 3. Prepare the ROI-based fMRI time series files.
% If you use this code, you will be appreciated for citing the above papers.


% The following codes can be run in Matlab.


Input=['F:\HCP_timeseries\100206_REST1_AAL.mat'];   % here is the input of a (n_time_points)*(n_ROIs) fMRI time series; for example, it would be a 1200*90 matrix when there is 1200 time points and using the AAL atlas with 90 ROIs.
timeseries=importdata(Input);


% Starting sliding windows

i=139;     % Setting the window width (TRs)
j=8;   % Setting the step length (TRs)
windows=(1200-i)/j+1;   % Calculating total number of windows; the number of 1200 should be revised based on your data
windows=floor(windows); 
dynamic_network=[];
for n=1:windows;
b=timeseries(((n-1)*j+1):((n-1)*j+i),:);
dynamic_network(:,:,n)=corr(b); 
end; 
clear timeseries; 
disp(['Sliding window finished']);


%%%%%%%%%%%%  Calculating the temporal clustering coefficient and other metrics

temporal_clust=[]; 
temporal_clust_nodal=[]; 
temporal_path=[]; 
temporal_path_nodal=[]; 

for k = 0.01:0.01:0.50;     %  Setting the range of density
dynamic_network_thresholded=[]; 
disp(['Working on density #',num2str(k)]); 
for n = 1:windows; 
M=dynamic_network(:,:,n); 
a1 = threshold_proportional(M,k);  % Proportional thresholding using the Brain Connectivity Toolbox
a1(a1~=0)=1; 
dynamic_network_thresholded(:,:,n) = a1; 
end; 
%%  Calculating temporal clustering coefficient
[ C,C_vec ] = temporal_clustering(dynamic_network_thresholded,0,90);  % 90 (number of ROIS) here should be revised based on the Atlas you used
C_vec=C_vec'; 
temporal_clust_nodal=[temporal_clust_nodal;C_vec]; 
temporal_clust=[temporal_clust;C]; 
%%  Calculating characteristic temporal path length (this step may take relatively long time and you may skip it if not needed)
[ C,C_vec ] = temporal_clustering(dynamic_network_thresholded,0,90);  % 90 (number of ROIS) here should be revised based on the Atlas you used
C_vec=C_vec'; 
temporal_path_nodal=[temporal_clust_nodal;C_vec]; 
temporal_path=[temporal_clust;C]; 
end; 

% Calculation of temporal clustering coefficient finished
% The output "temporal_clust" encodes global temporal clustering coefficient at each density.
% The output "temporal_clust_nodal" encodes nodal temporal clustering coefficient of each ROI at each density.
% The output "temporal_path" encodes global characteristic temporal path length at each density.
% The output "temporal_path_nodal" encodes nodal tcharacteristic temporal path length of each ROI at each density.
