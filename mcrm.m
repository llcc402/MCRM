function [ix, centers, K_vec, u_vec, w, w_1_1_vec, w_1_2_vec,w_1_3_vec,w_2_1_vec,w_2_2_vec,w_2_3_vec] = mcrm(data, alpha, maxIter)
%% init settings
if nargin < 2
    alpha = 1;
end
if nargin < 3
    maxIter = 100;
end

% all observations in one cluster
ix = ones(size(data)); 
centers = zeros(1, size(data,2)) + mean(mean(data));

% auxiliary variable
u = ones(1,2) * 100;
u_vec = ones(1, maxIter);

% learning pace for u
pace = 1;

% learning pace for w
pace_w = 0.003;

% the weights for each crm
w = ones(2,3);

% standard deviation for likelihood
sigma1 = 3;
% standard deviation for base measure
sigma0 = 6;

% the number of clusters for each iteration
K_vec = zeros(1, maxIter);

% monitor w
w_1_1_vec = zeros(1,maxIter);
w_1_2_vec = zeros(1,maxIter);
w_1_3_vec = zeros(1,maxIter);
w_2_1_vec = zeros(1,maxIter);
w_2_2_vec = zeros(1,maxIter);
w_2_3_vec = zeros(1,maxIter);

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
                ptr2(j) = [];
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
                if sum(q(:,k)) > 0
                    prob(k) = log(sum(q(:,k))) ...
                        - 1/2/sigma1^2 * (data(i,j) - centers(k))^2 ...
                        + log(get_tau_frac(w, u, q(:,k), i));
%                     fprintf(['tau_k = ', num2str(get_tau_frac(w, u, q(:,k), i)), '\n'])
                end
            end
            % prob for new cluster
            tau1 = 0;
            for r = 1:size(w,2)
                tau1 = tau1 + w(i,r) / (u * w(:,r) + 1);
            end
%             fprintf(['ta1 = ', num2str(tau1), '\n'])
            prob(K+1) = log(alpha) + log(sigma1) ...
                - log(sigma1^2 + sigma0^2)/2 ...
                - data(i,j)^2 / 2 / (sigma1^2 + sigma0^2) ...
                + log(tau1);
            
            % normalize prob
            prob = prob - max(prob);
            prob = exp(prob);
            prob = prob / sum(prob);
            
            [~,~,ix(i,j)] = histcounts(rand(1), [0, cumsum(prob)]);
            
            if ix(i,j) == K+1
                centers(K+1) = data(i,j);
            end
        end
    end
    
    % update centers
    d_vec = data(:);
    ix_vec = ix(:);
    K = 0;
    B = accumarray(ix_vec, 1:length(ix_vec), [], @(x){x});
    for i = 1:length(B)
        if ~isempty(B{i})
            % count number of clusters
            K = K + 1;
            
            % update centers
            if length(B{i}) == 1
                centers(i) = d_vec(B{i});
            else
                centers(i) = mean(d_vec(B{i}));
            end
        end
    end
    K_vec(iter) = K;
    
    % the cluster frequencies
    K = max(ix(:));
    q1 = histcounts(ix(1,:), 1:(K+1));
    q2 = histcounts(ix(2,:), 1:(K+1));
    q = [q1; q2];
    
    % update u1 and u2
    for i = 1:2
        term1 = (size(data,2) - 1) / u(i);
        
        term2 = 0;
        for r = 1:size(w, 2)
            term2 = term2 + w(i,r)/(1 + u * w(:,r));
        end
        term2 = alpha * term2;
        
        term3 = 0;
        
        for k = 1:K
            term3 = term3 + sum(q(:,k)) * get_tau_frac(w, u, q(:,k), i);
        end
        
        grad = term1 - term2 - term3;
        u(i) = u(i) + pace * grad;
    end
    u_vec(iter) = u(1);

    % update w
    for i = 1:size(data, 1)
        for r = 1:size(w,2)
            w_new = w(i,r) + pace_w * update_w_i_r(alpha, w, q, u, i, r);
            if w_new > 0
                w(i,r)=  w_new;
            end
        end
    end
    w_1_1_vec(iter) = w(1,1);
    w_1_2_vec(iter) = w(1,2);
    w_1_3_vec(iter) = w(1,3);
    w_2_1_vec(iter) = w(2,1);
    w_2_2_vec(iter) = w(2,2);
    w_2_3_vec(iter) = w(2,3);
    
    fprintf(['iter ', num2str(iter), ' done\n'])
    
end


% monitor u(1)
figure(1)
plot(u_vec)

% monitor w(1)
figure(2)
plot(w_1_1_vec)

end


            
    
    