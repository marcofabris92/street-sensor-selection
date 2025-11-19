function M = force_symmetry(M)

[m,n] = size(M);
if m ~= n
    error('The argument must be a symmetric square matrix.');
end


if norm(M-M',"fro") > 2e-15*round(n*(n-1)/2)
    error('The argument must be a symmetric square matrix.');
end

for i = 1:n-1
    for j = i+1:n
        m = (M(i,j)+M(j,i))/2;
        M(i,j) = m;
        M(j,i) = m;
    end
end

end