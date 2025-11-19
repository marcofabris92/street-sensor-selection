function value = f(metric,sys)

n = size(sys.C,2);

switch metric
    case 1      % rank
        value = my_rank(Gramian(sys.A,sys.C,n));
    case 2      % trace/n
        value = trace(Gramian(sys.A,sys.C,n))/n;
    case 3      % inverse condition number
        value = 1/cond(Gramian(sys.A,sys.C,n));
    case 4      % smallest eigenvalue
        value = min(svd(Gramian(sys.A,sys.C,n)));
    case 5      % determinant^1/n
        value = nthroot(abs(det(Gramian(sys.A,sys.C,n))),n);
    case 6      % H2 norm 
        value = trace(sys.B'*Gramian(sys.A,sys.C,Inf)*sys.B);
    case 7      % log det W_inf
        value = log(det(sys.B'*Gramian(sys.A,sys.C,Inf)*sys.B));
    case 8      % discounted H2 norm 
        value = trace(sys.B'*Gramian(sys.eta*sys.A,sys.C,Inf)*sys.B);
    otherwise   % rank by default
        value = my_rank(Gramian(sys.A,sys.C,n));
end

end