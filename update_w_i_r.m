% Input:
%     w       A matrix of order d * R. 'd' is the number of mixed random
%             measures and 'R' is the number of crms 
%     u       A row vector of length 'd'.
%     q       A matrix of order d * K, where K is the number of clusters.
%     i       A scalar. 
%     alpha   A scalar.
function grad = update_w_i_r(alpha, w, q, u, i, r)
uw = u * w + 1;
term1 = - alpha * u(i) / uw(i);

term2 = 0;
for k = 1:size(q,2)
    numeritor = q(i,k) / w(i,r) - sum(q(:,k)) * u(i) / uw(i);
    denominator = 0;
    for l = 1:size(w,2)
        summard = q(:,k)' * log(w(:,l) ./ w(:,r)) ...
            + sum(q(:,k)) * log(uw(l) / uw(r));
        summard = exp(summard);
        denominator = denominator + summard;
    end
    term2 = term2 + numeritor / denominator;
end

grad = term1 + term2;

end
        
