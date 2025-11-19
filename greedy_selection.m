function [Q_star,detectable,observable] =...
    greedy_selection(sys,f,metric,p_star)


C = sys.C;
[p,n] = size(C);

Q_star = zeros(p_star,1);
S = 1:p;

for k = 1:p_star
    max_fQi = -Inf;
    max_j = 0;
    for j = 1:p
        i = S(j);
        sys.C = C([Q_star(1:k-1); i],:);
        fQi = f(metric,sys);
        if fQi > max_fQi
            max_fQi = fQi;
            max_j = j;
            Q_star(k) = i;
        end
    end
    S(max_j) = S(p);
    p = p - 1;
end

Q_star = sort(Q_star,'ascend');
[~,observable,detectable] = obs(sys.A,C(Q_star,:),n,0,1,0,1);

end