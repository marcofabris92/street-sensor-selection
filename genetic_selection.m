function [Q_star,detectable,observable] =...
    genetic_selection(sys,metric,p_star)


C = sys.C;
[p,n] = size(C);


% x = ga(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,intcon,options)
par.metric = metric;
par.sys = sys;
par.p_star = p_star;
fun = @(beta)cost(par,beta);
% nvars == p
A = [ones(1,p); -ones(1,p)];
b = [p_star; -p_star];
% Aeq == ones(1,nvars);
% beq == p_star;
lb = zeros(p,1);
ub = ones(p,1);
nonlcon = @(beta)constraint(par,beta);
intcon = (1:p)';
options = optimoptions('ga');
% options.CreationFcn = 'gacreationuniform';
%options.FunctionTolerance = 1e-6;
%options.ConstraintTolerance = 1e-6;
beta_star = ga(fun,p,A,b,[],[],lb,ub,nonlcon,intcon,options);

Q_star = var2index(beta_star);
[~,observable,detectable] = obs(sys.A,C(Q_star,:),n,0,1,0,1);

end

function value = cost(par,beta)

par.sys.C = par.sys.C(var2index(beta),:);
value = -f(par.metric,par.sys);

end

function [c,ceq] = constraint(par,beta)

n = size(par.sys.C,2);
%[O,isObs,isDet,r,lag] = obs(A,C,order,getObs,getDet,getLag,check)
[~,~,isDet,r] = obs(par.sys.A,par.sys.C(var2index(beta),:),n,0,1,0,1);
c = 0;
if ~isDet
    c = n-r;
end
ceq = []; 
    
end

function Q = var2index(beta)

Q = zeros(sum(beta),1);
k = 0;
for i = 1:length(beta)
    if beta(i)
        k = k + 1;
        Q(k) = i;
    end
end

end
