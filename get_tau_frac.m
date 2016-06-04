% Input:
%     w       A matrix of order d * R. 'd' is the number of mixed random
%             measures and 'R' is the number of crms 
%     u       A row vector of length 'd'.
%     q       A column vector of length 'd'.
%     i       A scalar. 
function frac = get_tau_frac(w, u, q, i)

% the frequencies in atom k over all group
t = sum(q);

% wu(r) =  sum_{i=1}^d w_{i,r} * u_i + 1
wu = u * w + 1;

% init frac
frac = 0;

for r = 1:size(w,2)    
    % compute the denominator
    s = 0;
    for l = 1:size(w,2)
        s = s + exp(t * (log(wu(r)) - log(wu(l))) + q' * log(w(:,l) ./ w(:,r)));
    end
    
    frac = frac + w(i, r) / wu(r) / s;
end

end