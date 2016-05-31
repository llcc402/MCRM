function [ix, centers] = mcrm(data, alpha, maxIter)
%% init settings
if nargin < 1
    alpha = 1;
end
if nargin < 2
    maxIter = 100;
end

% all observations in one cluster
ix = ones(size(data,2)); 
centers = zeros(size(data,2)) + mean(mean(data));

% auxiliary variable
u_vecs = ones(2, maxIter + 1);
u1 = u_vecs(1,1);
u2 = u_vecs(2,1);

% learning pace for u
pace = 0.1;

%% post sampler
for iter = 1:maxIter
    