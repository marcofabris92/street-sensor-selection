clearvars
close all
clc
RHO_MAX = 2;   % model order (or better, maximum ARX lag)

%% --- Load data ---
path = '.\SANCARLOcE\';
%path = '/home/marco/Desktop/OneDrive_2025-10-28/Street_sensor_selection/Code/v05/SANCARLOcE/';
A = table2array(struct2table(load([path '_A_sancarlo_c.mat'])));
B = table2array(struct2table(load([path '_B_sancarlo_c.mat'])));
load([path 'flux_data.mat']) % loads vars u and y

y00 = y(1,:)';

T = size(y,1) - 1;
n = size(y,2);

% Find K such that each node can be reached in at most K steps
K = 0;
AA = eye(n);
while K < n
    K = K + 1;
    AA = A*AA;
    if isempty(find(AA == 0,1))
        break
    end
end


AA = zeros(n,n,K-1);
AA(:,:,1) = (A ~= 0);
for k = 2:K-1
    AA(:,:,k) = A*AA(:,:,k-1);
    AA(:,:,k) = (AA(:,:,k) ~= 0);
end
AA(:,:,K) = ones(n);


NA = zeros(1,K);
for k = 1:K-1
    NA(k) = sum(sum(AA(:,:,k)));
end
NA(K) = n^2;

BBB = 1*(B~=0);
CCC = ones(n,1)-BBB;
NB = sum(BBB);
T_MIN = 1+K+NA(K)+NB;
if T < T_MIN
    error(['Not enough data contained in y: T should be increased to ' num2str(T_MIN) ' at the very least.'])
end
MAX_ORDER = ceil((T+K*n^2-NA(K)-NB)/(1+n^2))-1;
if RHO_MAX > MAX_ORDER
    error(['RHO_MAX must be an integer selected within 1 and ' num2str(MAX_ORDER)])
end

AICc = zeros(1,RHO_MAX);
A_est_ = cell(RHO_MAX,1);
B_est_ = cell(RHO_MAX,1);
for rho = 1:RHO_MAX
    %% --- Prepare adjacency and free variables ---
    Adj = (A ~= 0);   % free A entries
    
    % Build regressor matrix Phi for multiple lags
    Phi = zeros(n*rho, T - rho + 1);
    for k = 1:rho
        Phi((k-1)*n+1:k*n, :) = y(rho - k + 1:T - k + 1, :)';
    end
    Y = y(rho+1:T+1,:)';   % target outputs
    U = u(rho+1:T+1)';     % inputs
    T_eff = T - rho + 1;   % effective number of samples
    
    %% --- Map free variables (extended to rho lags) ---
    A_map = zeros(n,n,rho);
    col = 0;
    for i = 1:n
        for j = 1:n
            if Adj(i,j)
                for k = 1:rho
                    col = col + 1;
                    A_map(i,j,k) = col;
                end
            end
        end
    end
    
    B_map = zeros(n,1);
    for i = 1:n
        if B(i) ~= 0
            col = col + 1;
            B_map(i) = col;
        end
    end
    n_vars = col;
    
    %% --- Build H matrix ---
    H = zeros(n*T_eff, n_vars);
    y_vec = Y(:);
    
    for i = 1:n
        for j = 1:n
            for k = 1:rho
                if A_map(i,j,k) > 0
                    idx = A_map(i,j,k);
                    for t = 1:T_eff
                        row = (i-1)*T_eff + t;
                        H(row,idx) = Phi((k-1)*n + j, t);
                    end
                end
            end
        end
        if B_map(i) > 0
            idx = B_map(i);
            for t = 1:T_eff
                row = (i-1)*T_eff + t;
                H(row,idx) = U(t);
            end
        end
    end
    
    %% --- Equality constraints (only first lag + B) ---
    Aeq = zeros(n, n_vars);
    beq = ones(n,1);
    for j = 1:n
        free_rows = find(Adj(:,j));
        for k = 1:numel(free_rows)
            Aeq(j, A_map(free_rows(k), j, 1)) = 1; % only first lag
        end
        if B_map(j) > 0
            Aeq(j, B_map(j)) = 1;
        end
    end
    
    %% --- Display info ---
    fprintf('Model order (rho): %d\n', rho);
    fprintf('Number of variables: %d\n', n_vars);
    fprintf('Equality constraints: %d\n', size(Aeq,1));
    fprintf('H size: %d x %d\n', size(H,1), size(H,2));
    fprintf('\n')



    %% --- Solve regularized QP ---
    lambda_vec = [0 logspace(-5,5,100)];
    best_mse = +Inf; 
    best_theta = []; 
    best_lambda = NaN;
    lb = zeros(n_vars,1);      % lower bounds
    ub = ones(n_vars,1);       % upper bounds
    
    for lam = lambda_vec
        Q = H'*H + lam*eye(n_vars); 
        f = -H'*y_vec;
    
        theta = quadprog(2*Q,2*f,[],[],Aeq,beq,lb,ub,[],...
                optimoptions('quadprog','Display','off'));
    
        % Recover A_est and B_est
        A_est = zeros(rho*n);
        for r = 1:rho
            for i = 1:n
                for j = 1:n
                    if A_map(i,j,r) > 0
                        A_est(i,j+(r-1)*n) = theta(A_map(i,j,r));
                    end
                end
            end
        end
        B_est = zeros(rho*n,1);
        for i = 1:n
            if B_map(i) > 0
                B_est(i) = theta(B_map(i));
            end
        end
    
        % Compute MSE and AICc
        for i = n+1:rho*n
            A_est(i,i-n) = 1;
        end
        Y_pred = A_est*Phi + B_est*U;
        Y_pred = Y_pred(1:n,:);
        mse = mean((Y(:)-Y_pred(:)).^2);
    
        if mse < best_mse
            best_mse = mse;
            best_theta = theta;
            A_est_{rho} = A_est;
            B_est_{rho} = B_est;
            best_lambda = lam;
            AICc(rho) = ( T_eff*log(det((reshape(y_vec - H*best_theta, n, T_eff) * reshape(y_vec - H*best_theta, n, T_eff)')/T_eff)) + 2*n_vars + T_eff*( n*(log(2*pi)+1) ) ) + (2*n_vars*(n_vars+1))/(T_eff - n_vars - 1);
        end
    end
    
    fprintf('Selected lambda = %.4g, val MSE = %.4f\n', best_lambda, best_mse);
    fprintf('\n')

end


AICc
[~,rho] = min(AICc);
A_est = A_est_{rho};
B_est = B_est_{rho};

%% --- Simulation y_hat ---
y_hat = zeros(rho*n,T+1);
y_hat(1:n,1) = y00; 
for t = 1:T
    y_hat(:,t+1) = A_est*y_hat(:,t) + B_est*u(t);
end
y_hat = y_hat(1:n,:);

%% --- Plot ---
ftsz = 20;
intr = 'latex';
lw = 1.5;

figure('position',[100 100 1000 600])
tt = (0:T)';
y_plot = y';
grid on
hold on
set(gca, 'ColorOrder', lines(n), 'NextPlot', 'add');
colors = get(gca,'ColorOrder');
for i = 1:n
    plot(tt,60*y_plot(i,:),'Color',colors(i,:),'LineStyle','-','LineWidth',lw)
    plot(tt,60*y_hat(i,:),'Color',colors(i,:),'LineStyle','--','LineWidth',lw)
end
xlabel('$t$ [s]', 'Interpreter',intr)
ylabel('Trajectories [veh/min]','Interpreter',intr)
xaxisproperties = get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = intr;
xaxisproperties.FontSize = ftsz;
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = intr;   
yaxisproperties.FontSize = ftsz;
xlim([0 T])

figure('position',[100 100 1000 600])
grid on
hold on
y_tilde = 60*(y_plot - y_hat);
for i = 1:n
    plot(tt,y_tilde(i,:),'Color',colors(i,:),'LineWidth',lw)
end
xlabel('$t$ [s]', 'Interpreter',intr)
ylabel('Estimation errors [veh/min]','Interpreter',intr)
xaxisproperties = get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = intr;
xaxisproperties.FontSize = ftsz;
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = intr;   
yaxisproperties.FontSize = ftsz;
xlim([0 T])

%% --- Save results ---
save([path '_A_est_sancarlo_c.mat'],'A_est')
save([path '_B_est_sancarlo_c.mat'],'B_est')
