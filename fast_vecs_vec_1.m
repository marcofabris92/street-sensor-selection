function z = fast_vecs_vec_1(z)

[m,n] = size(z);
if m > 1 && n > 1
    error('The argument must be a vector.')
end

n = round(sqrt(length(z)));
if n^2 ~= length(z)
    error('The argument must be the vectorized form of a square matrix.')
end


j = 0;
reset = 0;
pivot = 1;
while pivot <= length(z)
    i = pivot;
    j = j + 1;
    z(j) = z(i);
    reset = reset + 1;
    while i < n*reset
        i = i + 1;
        j = j + 1;
        z(j) = 2*z(i);
    end
    pivot = pivot + n+1;
end

z = z(1:round(n*(n+1)/2));

end