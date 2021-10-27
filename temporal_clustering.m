function [ C,C_vec ] = temporal_custering( adjArray,directed,nNodes )
% this function calculates the temporal clustering coefficient of a dynamic
% network as described by Long et al. 2021.
%
% Inputs:
%       adjArray = a nNodes x nNodes x nTimepoints 3-dimentional matrix encoding connections between nodes at different time points
%
%       directed = 1 if the dynamic network is directed, 0 otherwise. 
%
% Output:
%       C = temporal clustering coefficient of the network
%       C_vec = nodal temporal clustering coefficient of each node
%
% 
% The function was adapted from codes in https://github.com/asizemore/Dynamic-Graph-Metrics
% Reference: "Ann E. Sizemore and Danielle S. Bassett, "Dynamic Graph 
% Metrics: Tutorial, Toolbox, and Tale." Submitted. (2017)"
% and "Long Y et al. Evaluating test-retest reliability and sex/age-related
% effects on temporal clustering coefficient of dynamic functional brain networks (2021)"




T = size(adjArray,3);


C_vec = zeros(nNodes,1);

for node_i = 1:size(adjArray,1)

    gammaVec = zeros(1,T-1);

    for t = 1:T-1
        num = 0;

        for j = 1:nNodes
            if node_i~= j
                num = num + adjArray(node_i,j,t)*adjArray(node_i,j,t+1);
            end
        end

        den = sqrt(sum(adjArray(node_i,:,t))*sum(adjArray(node_i,:,t+1)));

        if den~=0
        gamma = num/den;
        gammaVec(t) = gamma;
        end

    end

    % sum over timepoints

    C_i = (1/(T-1))*sum(gammaVec);

    C_vec(node_i) = C_i;

end

C = sum(C_vec)/nNodes;

end

