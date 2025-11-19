function [O,isObs,isDet,r,lag] = ...
    obs(A,C,order,getObs,getDet,getLag,check)

% O: observability matrix of the given order
% isObs: flag == 1, if (A,C) happens to be observable 
% isDet: flag == 1, if (A,C) happens to be detectable
% r: rank of the observability matrix of order size(C,2)
% lag: an integer in {1,...,min(n,order)} equal to the minimum number of
% steps required to obtain a full column rank observability matrix

% To construct the observability matrix, set: check == 0
% To speed up the construction of the observability matrix O_order, set:
% getObs == 0, getDet == 0, getLag == 0
% Important: if check == 1 then order == size(C,2) is forced.

% To compute the lag, set:
% order == size(C,2), getObs == 0, getDet == 0, getLag == 1, check == 1

% To speed up observability check, set: check == 1.
% Important: if check == 1 then getObs = 1 is forced.

% To check detectability, set: check == 1
% To speed up detectability check, set: order == size(C,2)
% Important: if getDet == 1 then getObs = 1 is forced.

% Note: for certain couples (A,C), computing the lag may or may not help 
% determining whether observability or detectability hold...


p = size(C,1);
n = size(C,2);
O = zeros(p*order,n);
O(1:p,:) = C;

if check
    getObs = 1;
    order = n;
end

isDet = NaN;
if getDet
    getObs = 1;
    isDet = 0;
end

isObs = NaN;
if getObs
    isObs = 0;
end

lag = NaN;
r = NaN;
if getLag && my_rank(C) == n 
    lag = 1;
    isObs = 1;
    r = n;
    if check && ~getDet
        return
    end
end

for i = 1:min(order,n)-1
    O(i*p+1:(i+1)*p,:) = O((i-1)*p+1:i*p,:)*A;
    if getLag && isnan(lag) && my_rank(O(1:(i+1)*p,:)) == n 
        lag = i+1;
        isObs = 1;
        r = n;
        isDet = 1;
        if check && ~getDet
            return
        end
    end
end

Or = zeros(0,n);
if order < n && getObs && ~isObs
    Or = zeros(p*(n-order),n);
end

if order >= n
    if getDet || ~getLag && getObs
        r = my_rank(O(1:p*n,:));
    end
    if ~getLag && getObs && r == n 
        isObs = 1;
        isDet = 1;
    end
    if check && ~getDet
        return
    end
    for i = n:order-1
        O(i*p+1:(i+1)*p,:) = O((i-1)*p+1:i*p,:)*A;
    end
elseif getObs && ~isObs
    Or(1:p,:) = O((order-1)*p+1:order*p,:)*A;
    for i = 1:n-order-1
        Or(i*p+1:(i+1)*p,:) = Or((i-1)*p+1:i*p,:)*A;
    end
    r = my_rank([O; Or]);
    if r == n 
        isObs = 1;
        isDet = 1;
    end
end

if getDet && ~isDet
    [VT,~] = qr([O(1:p*n,:); Or]',0);
    V = VT';
    A = VT*A*V;
    C = C*V;
    while r >= 1 && my_rank(obsv(A(1:r,1:r),C(:,1:r))) < r
        r = r-1;
    end
    lambda = abs(eig(A(r+1:end,r+1:end)));
    isDet = 1;
    k = 0;
    while isDet && k < length(lambda)
        k = k + 1;
        if lambda(k) >= 1
            isDet = 0;
        end
    end
end

end

