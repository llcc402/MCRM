% Input:
%     w       a matrix. The element w(i, r) is the weight for the i-th
%             mixed random measure in the r-th crm.
%     u       a row vector. The length of u is equal to the number of rows 
%             of w.
%     q       a column vector. The length of q is equal to the number of 
%             rows of w.
function tau = get_tau(w, u, q)
tau = 0;
t = sum(q);
for r = 1:size(w, 2)
    numeritor = prod(w(:,r) .^ q);
    denominator = (u * w(:,r) + 1) .^ t;
    tau = tau + numeritor / denominator;
end
end