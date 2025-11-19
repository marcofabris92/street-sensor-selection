function Z = vech_1(z)

[m,n] = size(z);
if m > 1 && n > 1
    error('The argument must be a vector.')
end

N = length(z);
n = round((-1+sqrt(1+8*N))/2);
if 2*N ~= n*(n+1)
    error('The argument must be the half-vectorized form of a square matrix.')
end

Z = zeros(n);
k = 0;
for j = 1:n
    k = k + 1;
    Z(j,j) = z(k);
    for i = j+1:n
        k = k + 1;
        Z(i,j) = z(k);
        Z(j,i) = z(k);
    end
end

end