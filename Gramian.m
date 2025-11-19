function W = Gramian(A,C,order)

if order == Inf
    % A'*W_inf*A - W_inf + C'*C = 0
    % X = dlyap(A,Q) solves A*X*A' - X + Q = 0
    W = dlyap(A',C'*C);
else
    W = obs(A,C,order,0,0,0,0);
    W = W'*W;
end
W = force_symmetry(W);

end