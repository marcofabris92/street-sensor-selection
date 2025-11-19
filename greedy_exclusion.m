function [Q_star,detectable,observable,partial] =...
    greedy_exclusion(sys,f,metric,p_star)


C = sys.C;
[p,n] = size(C);


Q_star = [];
partial = 0;
[~,observable,detectable] = obs(sys.A,C,n,0,1,0,1);
if ~detectable
    return
end

Q = (1:p)';
q = zeros(p,1);
k = 1;
while k <= p-p_star && partial == 0
    min_fQi = +Inf;
    for j = 1:length(Q)
        i = Q(j);
        sys.C = C(setdiff(Q,i),:);
        [~,isObs,isDet] = obs(sys.A,sys.C,n,0,1,0,1);
        if isDet
            fQi = f(metric,sys);
            if fQi < min_fQi
                min_fQi = fQi;
                q(p-k+1) = i;
            end
            if fQi == min_fQi && isObs
                min_fQi = fQi;
                q(p-k+1) = i;
            end
        end
    end
    if q(p-k+1) > 0
        Q = setdiff(Q,q(p-k+1));
        k = k + 1;
    else
        partial = k;
    end
end

Q_star = sort(Q,'ascend');
if partial > 0
    detectable = 0;
    observable = 0;
else
    [~,observable] = obs(sys.A,C(Q_star,:),n,1,0,0,1);
end

end