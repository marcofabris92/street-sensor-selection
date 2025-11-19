function sys = system_setup(sys)

rng(3)

A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;


%% check dimensions
[p,n] = size(C);
[p_prime,m_prime] = size(D);
[n_prime,m] = size(B);
if size(A,1) ~= n || size(A,2) ~= n ||...
        n ~= n_prime || m ~= m_prime || p ~= p_prime
    error('\nDimensions of system matrices are not consistent.')
end

%% default initial condition
sys.x0 = zeros(n,1);

%% find Q_tilde
p_tilde = 0; % <--- can be changed manually (test starts from p_tilde+1)
Q_tilde = [];
avg_rank = n; % <--- huristics to skip some values of p_tilde
while p_tilde < p-1 && isempty(Q_tilde)
    p_tilde = p_tilde + 1;
    %p_tilde = min(p-1, p_tilde + max(1,round((n-avg_rank)/2)));
    fprintf(['testing p_tilde = ' num2str(p_tilde) '\n'])
    [Q_tilde,avg_rank] = get_Q_tilde(sys,p_tilde);
end
if p_tilde == p-1 && isempty(Q_tilde)
    p_tilde = p;
    Q_tilde = get_Q_tilde(sys,p_tilde);
    if ~isempty(Q_tilde)
        fprintf('Warning: p == p_tilde!\n')
    end
end
Q_tilde
p_tilde

if isempty(Q_tilde)
    error('Asm. 1 is violated.\n')
end

%% lag computation w.r.t. (A,C_Q_tilde)
%               obs(A,C,order,getObs,getDet,getLag,check)
[~,~,~,~,lag] = obs(A,C(Q_tilde,:),n,0,0,1,1);
if isnan(lag)
    error('The current sensor configuration makes the system unobservable.')
end
sys.lag = lag;
lag

%% -----------------------------------------------------------------------
% SLIGHT VARIATION ON FILIPPOS FOTIADIS CODE TO COMPUTE pinv(Phi')

C_tilde = C(Q_tilde,:);
N = 1*lag;
L = N*(m+p_tilde);
t0 = N+1;
tK = round(L*(L+1));

% Controller to stabilize the system output
R_weight = 1e0*eye(m);
Q_weight = C_tilde'*C_tilde;
[~,~,K] = dare(A,B,Q_weight,R_weight);

% Data gathering
x = zeros(n,tK);
u = randn(m,tK);
for t = 1:tK-1
    x(:,t+1) = A*x(:,t)+B*u(:,t);
end
y_tilde_1 = C_tilde*x; % <-- used to estimate eta

noisegain = 1e2; % the higher the more persistency of excitation 
x(:,1) = 1e-3*randn(n,1);
u = noisegain*randn(m,tK);
for t = 1:tK-1
    x(:,t+1) = A*x(:,t)+B*u(:,t);
    u(:,t+1) = u(:,t+1) - K*x(:,t+1);
end
y_tilde = C_tilde*x;   % <-- used to compute H2 norm in a data-driven way
    
%%% ********************************************************
% select eta to get a stable A_tilde = eta*A 
eta = 1;
rho_A = max(abs(eig(A)));
if rho_A >= 1
    eta = (1-eps)/rho_A; % <--- set the numerator in (0,1)
end
eta % <--- true eta

% estimate eta from data 
yy = zeros(tK,1);
tEND = 0;
while tEND < tK && norm(y_tilde_1(:,tEND+1)) < 1e50
    tEND = tEND + 1;
    yy(tEND) = norm(y_tilde_1(:,tEND))^2;
end
tSS = ceil(0.01*tEND);
expfit = fit((tSS:tEND-1)',yy(tSS+1:tEND),'exp1');

stbl_thr = 0.95; % in (0,1]: the lowest the more stable & inaccurate the
                 % H2 norm will be
eta_est = 1;
if expfit.a > 0 && expfit.b >= 0
    eta_est = stbl_thr/exp(expfit.b/2); % <--- set the numerator in (0,1)
end
eta_est % <--- estimated eta

eta_assigned = eta_est; % <--- assigned eta
eta_assigned
sys.eta = eta_assigned;      
sys.eta2 = sys.eta^2;
%fprintf('Press any key to continue\n')
%pause
%%% ********************************************************

% Data matrix formation
Eu = [eye(m); zeros(L-m, m)];
k = 1;
tsamples = t0:tK; 
Nsamples = length(tsamples);
Phi_rows = round(L*(L+1)/2);
Phi = zeros(Phi_rows,Nsamples);

t = tsamples(1);
tmpy = y_tilde(:,t-1:-1:t-N);
tmpu = u(:,t-1:-1:t-N);
z = [tmpu(:); tmpy(:)];

for t = tsamples
    tmpy = y_tilde(:,t:-1:t-N+1);
    tmpu = u(:,t:-1:t-N+1);
    z1 = [tmpu(:); tmpy(:)];
    
    %Compute phi(t) and add to Phi(k)
    z1Eu = z1-Eu*u(:,t);
    Phi(:,k) = fast_vecs_vec_1(sys.eta2*kron(z1Eu,z1Eu)-kron(z,z));
    
    z = z1;
    k = k + 1;
end

% SMALLER SINGULAR VALUE == LARGER ERRORS, CHECK PE IN DATA
sing_Phi = svd(Phi);
smin = sing_Phi(Phi_rows);
fprintf('Min. sing. val.: %d. \n', smin);
% if smin < 1e-12  % CAREFUL!!! This check should be enabled
%     error('RETRY: Phi is not full row rank.')
% end

sys.tsamples = tsamples;
sys.Nsamples = Nsamples;
sys.y = C*x;
sys.Eu = Eu;
sys.PhiTinv = pinv(Phi',1e-15);

% -----------------------------------------------------------------------

end



function [Q_tilde,avg_rank] = get_Q_tilde(sys,p_tilde)


C = sys.C;
[p,n] = size(C);

Q = zeros(p_tilde,1);
Q_tilde = [];

ii = p_tilde;
depth = 1;
cw = [ones(1,ii) zeros(1,p-ii)];
done = 0;
avg_rank = 0;
count = 0;
while ~done
    k = 1;
    for i = 1:p
        if cw(i) == 1
            Q(k) = i;
            k = k + 1;
        end
    end
    % obs(A,C,order,getObs,getDet,getLag,check)
    [~,observable,~,r] = obs(sys.A,C(Q,:),n,1,0,0,1);
    avg_rank = avg_rank + r;
    count = count + 1;
    if observable
        Q_tilde = Q;
        return
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

avg_rank = round(avg_rank / count);

end
