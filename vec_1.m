function Z = vec_1(z)

[m,n] = size(z);
if m > 1 && n > 1
    error('The argument must be a vector.')
end

n = round(sqrt(length(z)));
if n^2 ~= length(z)
    error('The argument must be the vectorized form of a square matrix.')
end
Z = reshape(z,n,n);

end