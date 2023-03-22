function [ C,C_vec ] = temporal_temporal_path( adjArray,directed,nNodes )
% this function calculates the characteristic temporal path length of a dynamic
% network as described by Long et al.
%
% Inputs:
%       adjArray = a nNodes x nNodes x nTimepoints 3-dimentional matrix encoding connections between nodes at different time points
%
%       directed = 1 if the dynamic network is directed, 0 otherwise. 
%
% Output:
%       C = characteristic temporal path length of the network
%       C_vec = nodal characteristic temporal path length of each node
%
% 
% The function was adapted from codes in https://github.com/asizemore/Dynamic-Graph-Metrics
% Reference: "Ann E. Sizemore and Danielle S. Bassett, "Dynamic Graph 
% Metrics: Tutorial, Toolbox, and Tale." Submitted. (2017)"
%




npoints = size(adjArray,3);
T = npoints;
durationShortestPaths = zeros(nNodes);
C_vec = zeros(nNodes,1);



nPathsDurationt = sum(adjArray,3);
nFastestPathsDurationt = nPathsDurationt;
nFastestPaths = nFastestPathsDurationt;
nFastestPaths(logical(eye(nNodes))) = 1;
durationShortestPaths(find(nPathsDurationt)) = 1;
durationShortestPaths(logical(eye(nNodes))) = 1;


% Note: we assume time steps are 1:nPoints 

for t = 2:npoints

    tArray = zeros(nNodes,nNodes,npoints-t+1);
    
    for p = 1:size(tArray,3)
        tArrayp = adjArray(:,:,p);
        for j = p+1:p+t-1           
            tArrayp = tArrayp*adjArray(:,:,j)+ tArrayp;
        end
        tArray(:,:,p)=tArrayp;       
    end
    
      
    % need to update shortest path matrix
    nPathsDurationt = sum(tArray,3);
	nPathsDurationt(nPathsDurationt>0) = 1;
    nFastestPathsDurationt(nFastestPaths==0) = ...
        nPathsDurationt(nFastestPaths==0);

    durationShortestPaths(nFastestPaths==0) = t*nFastestPathsDurationt(nFastestPaths==0);
    
    % update fastest paths
    nFastestPaths(nFastestPaths==0) = ...
        nFastestPathsDurationt(nFastestPaths==0);
    
       
    
    
end



% Prepare to find path members

durationShortestPaths(~durationShortestPaths) = inf;
durationShortestPaths(logical(eye(nNodes))) = 0;






L_matrix=durationShortestPaths;
L_matrix(L_matrix==0) = NaN;
L_matrix=1./L_matrix;

L_Harmonic = (nNodes*(nNodes-1))/(nansum(L_matrix(:))); 
L_nodal = (nNodes-1)./(nansum(L_matrix,2)); 
L_nodal=L_nodal';

C = L_Harmonic;
C_vec = L_nodal;



end