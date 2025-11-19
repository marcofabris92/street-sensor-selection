function z = vech(Z)

[m,n] = size(Z);
if m ~= n
    error('The argument must be a symmetric square matrix.');
end

N = round(n*(n-1)/2);
if norm(Z-Z',"fro") > 2e-15*N
    error('The argument must be a symmetric square matrix.');
end

z = zeros(N+n,1);
k = 0;
for j = 1:n
    for i = j:n
        k = k + 1;
        z(k) = Z(i,j);
    end
end

end