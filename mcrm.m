function [ix, centers] = mcrm(data, alpha, maxIter)
%% init settings
if nargin < 2
    alpha = 1;
end
if nargin < 3
    maxIter = 100;
end

% all observations in one cluster
ix = ones(size(data,2)); 
centers = zeros(size(data,2)) + mean(mean(data));

% auxiliary variable
u = ones(1,2);
u_vec = ones(1, maxIter);

% learning pace for u
pace = 0.1;

% the weights for each crm
w = [1 1 0; 1 0 1];

% standard deviation for likelihood
sigma1 = 0.5;
% standard deviation for base measure
sigma0 = 1;


%% post sampler
for iter = 1:maxIter
    
    % update ix(i,j) for each (i,j)
    for i = 1:2
        for j = 1:size(data,2) 
            if i == 1
                ptr1 = ix(1,:);
                ptr2 = ix(2,:);
                
                % delete (i,j)
                ptr1(j) = [];
            else
                ptr1 = ix(1,:);
                ptr2 = ix(2,:);
                
                % delete (i,j)
                ptr1(j) = [];
            end
                
            % the number of cluters
            K = max(max(ptr1), max(ptr2));

            % the frequencies for each group
            q1 = histcounts(ptr1, 1:(K+1));
            q2 = histcounts(ptr2, 1:(K+1));
            q = [q1; q2];
            
            % prob
            prob = zeros(1, K+1);
            
            % prob for existing clusters
            for k = 1:K
                q_p = q(:,k);
                q_p(i) = q_p(i) + 1;
                if sum(q_p) > 0
                    prob(k) = log(sum(q_p)) ...
                        - 1/2/sigma1^2 * (data(i,j) - centers(k))^2 ...
                        + log(get_tau(w, u, q_p)) ...
                        - log(get_tau(w, u, q(:,k)));
                end
            end
            % prob for new cluster
            tau1 = 0;
            for r = 1:size(w,2)
                tau1 = tau1 + w(i,r) / (u * w(:,r) + 1);
            end
            prob(K+1) = log(alpha) + log(sigma1) ...
                - log(sigma1^2 + sigma0^2)/2 ...
                - data(i,j)^2 / 2 / (sigma1^2 + sigma0^2) ...
                + log(tau1);
            
            % normalize prob
            prob = prob - max(prob);
            prob = exp(prob);
            prob = prob / sum(prob);
            
            [~,~,ix(i,j)] = histcounts(rand(1), [0, cumsum(prob)]);
        end
    end
    
    % update centers
    d_vec = data(:);
    ix_vec = ix(:);
    B = accumarray(ix_vec, 1:length(ix_vec), [], @(x){x});
    for i = 1:length(B)
        if ~isempty(B{i})
            if length(B{i}) == 1
                centers(i) = d_vec(B{i});
            else
                centers(i) = mean(d_vec(B{i}));
            end
        end
    end
    
    % update u1 and u2
    for i = 1:2
        term1 = 1 / u(i);
        
        term2 = 0;
        for r = 1:3
            term2 = term2 + w(i,r)/(1 + u * w(:,r));
        end
        term2 = alpha * term2;
        
        term3 = 0;
        K = max(ix(:));
        q1 = histcounts(ix(1,:), 1:(K+1));
        q2 = histcounts(ix(2,:), 1:(K+1));
        q = [q1; q2];
        for k = 1:K
            numerator = 0;
            for r = 1:3
                numerator = numerator + sum(q(:,k)) * w(i,r) *...
                    prod(w(:,r).^q(:,k)) / (u * w(:,r) + 1) ^ (sum(q(:,k)) + 1);
            end
            
            denominator = 0;
            for r = 1:3
                denominator = denominator + prod(w(:,r).^q(:,k)) ...
                    / (u * w(:,r) + 1) ^ sum(q(:,k));
            end
            term3 = term3 + numeritor / denominator;
        end
        
        grad = term1 - term2 - term3;
        u(i) = u(i) + pace * grad;
    end
    fprintf(['iter ', num2str(iter), ' done'])
    
    % monitor u(1)
    u_vec(iter) = u(1);
    plot(u_vec)
end


            
    
    