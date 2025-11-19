function [Q_star,detectable,observable] =...
    random_selection(sys,f,metric,p_star,alpha)


C = sys.C;
[p,n] = size(C);

Q = zeros(p_star,1);
Q_star = [];
score = -Inf;
detectable = 0;
observable = 0;
Q0 = [];
score0 = -Inf;

for j = 1:ceil(alpha*nchoosek(p,p_star))
    idx = 1:p;
    for k = 1:p_star
        i = unidrnd(p-k+1);
        Q(k) = idx(i);
        idx(i) = idx(p-k+1);
    end
    sys.C = C(Q,:);
    fQ = f(metric,sys);
    % obs(A,C,order,getObs,getDet,getLag,check)
    [~,isObs,isDet] = obs(sys.A,sys.C,n,0,1,0,1);
    if fQ > score && isDet
        Q_star = Q;
        score = fQ;
        detectable = 1;
    end
    if fQ == score && ~observable && isObs
        Q_star = Q;
        score = fQ;
        observable = 1;
        detectable = 1;
    end
    if fQ > score0
        Q0 = Q;
        score0 = fQ;
    end
end

if ~detectable 
    Q_star = Q0;
end

Q_star = sort(Q_star,'ascend');

end