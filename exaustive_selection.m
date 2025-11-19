function [Q_star,detectable,observable] =...
    exaustive_selection(sys,f,metric,p_star)


C = sys.C;
[p,n] = size(C);

Q = zeros(p_star,1);
Q_star = [];
score = -Inf;
detectable = 0;
observable = 0;
Q0 = [];
score0 = -Inf;

ii = p_star;
depth = 1;
cw = [ones(1,ii) zeros(1,p-ii)];
done = 0;
count = 0;
total_iterations = nchoosek(p,p_star);
completion = -1;
while ~done
    count = count + 1;
    new_completion = floor(100*count/total_iterations);
    if new_completion > completion
        %clc
        completion = new_completion;
        %fprintf(num2str(completion));
        %fprintf('\n')
    end
    k = 1;
    for i = 1:p
        if cw(i) == 1
            Q(k) = i;
            k = k + 1;
        end
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

if ~detectable 
    Q_star = Q0;
end

%clc

end