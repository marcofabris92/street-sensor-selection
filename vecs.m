function z = vecs(Z)

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
    k = k + 1;
    z(k) = Z(j,j);
    for i = j+1:n
        k = k + 1;
        z(k) = 2*Z(i,j);
    end
end

end