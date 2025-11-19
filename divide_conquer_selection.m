function [Q_star,detectable,observable] =...
    divide_conquer_selection(sys,f,metric,p_star)


C = sys.C;
[p,n] = size(C);

ff = zeros(p,1);
for j = 1:p
    sys.C = C(j,:);
    ff(j) = f(metric,sys);
end
[~,idx] = sort(ff,'descend');
Q = zeros(p_star,1);
Q_star = idx(1:p_star);
detectable = 0;
observable = 0;

ii = p_star;
depth = 1;
cw = [ones(1,ii) zeros(1,p-ii)];
done = 0;
while ~detectable && ~done
    
    k = 1;
    for i = 1:p
        if cw(i) == 1
            Q(k) = idx(i);
            k = k + 1;
        end
    end
    % obs(A,C,order,getObs,getDet,getLag,check)
    [~,isObs,isDet] = obs(sys.A,C(Q,:),n,0,1,0,1);
    if isDet
        Q_star = Q;
        detectable = 1;
    end
    if ~observable && isObs
        Q_star = Q;
        observable = 1;
        detectable = 1;
    end
        
    % handling the while loop to search all sensor configurations
    if ii < p
        cw(ii) = 0;
        ii = ii + 1;
        cw(ii) = 1;
    else
        jj = p;
        kk = 1;
        stack = zeros(depth,1);
        stack(kk) = p;
        found = 0;
        while ~found && jj > 1
            jj = jj - 1;
            if cw(jj) == 1
                if cw(jj+1) == 0
                    found = 1;
                    ii = jj+1+kk;
                    depth = p-jj;
                    cw(jj) = 0;
                    cw(jj+1:ii) = ones(1,kk+1);
                    for hh = 1:kk
                        jj = stack(hh);
                        if jj > ii
                            cw(jj) = 0;
                        end
                    end
                else
                    kk = kk + 1;
                    stack(kk) = jj;
                end
            end
        end
        if ~found
            done = 1;
        end
    end
end

Q_star = sort(Q_star,'ascend');

end