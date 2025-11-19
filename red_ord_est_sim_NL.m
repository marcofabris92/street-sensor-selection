function [x,y,x_hat,ys,y_hat,e_x,norm_e_x,os_ref] = ...
    red_ord_est_sim_NL(sys,est,S,u,x_hat_0)

% Useful dimensions
Tsim = size(u,2)-1;
n = length(sys.x0);
p = size(sys.C,1);
p_star = size(S,1);
p_bar = n-p_star;

% Initializations
Ds = S*sys.D;
x = zeros(2*n,Tsim+1);
y = zeros(p,Tsim+1);
ys = zeros(p_star,Tsim+1);
xi_hat = zeros(p_bar,Tsim+1);
x_hat = zeros(n,Tsim+1);
y_hat = zeros(p,Tsim+1);
e_x = zeros(n,Tsim+1);
norm_e_x = zeros(1,Tsim+1);
os_ref = NaN;

%% Bounded METANET parameters

pp.Ts = 1;
pp.rho_crit = 30 /1000; % > 10 and < 150
% pp.kappa = 10 /1000; % >0 small
% pp.anti = 0; % 50 % >= 0
% pp.conv = 0; % 1 % either 0 or 1
pp.vfree = 50 /3.6; % speed limits
pp.tau = 5; % > 0
pp.a = 2; % between 1 and 3
pp.eta_tilde = 0.5;
pp.kappa_tilde = 150;

street_lengths = [
    158;
    114;
    264;
    153;
    464;
    43;
    467;
    694;
    900;
    90;
    125;
    174;
    171;
    120;
    212;
    331;
    149;
    167;
    213;
    119;
    193;
    409;
    74;
    40;
    68];

pp.DL = inv(diag(street_lengths));

pp.n = size(sys.A,1);
pp.nx = 0;
pp.PU = sys.A-eye(pp.n); 
pp.PD = sys.A'-eye(pp.n);
pp.B = sys.B;

%% Dynamics of the reduced order estimator

% 1) Start [2b)]
t = 1;
x(:,t) = sqrt([sys.x0; sys.x0]);
x_bar_0 = est.V*sys.x0;
x_bar_hat_0 = est.V*x_hat_0;
xi_hat(:,t) = x_bar_hat_0(1:p_bar) + est.L*x_bar_0(p_bar+1:end);

% 2) Updates
for t = 1:Tsim
    % 2a) Outputs, estimates and errors
    [xx,yy] = METANET(pp,x(:,t),u(:,t:t+1));
    y(:,t) = yy(1,1:n)';
    ys(:,t) = S*y(:,t);
    ys_hat = ys(:,t) - Ds*u(:,t);
    x_hat(:,t) = est.Vinv*[xi_hat(:,t)-est.L*ys_hat; ys_hat];
    y_hat(:,t) = sys.C*x_hat(:,t) + sys.D*u(:,t);
    e_x(:,t) = x_hat(:,t) - y(:,t);
    norm_e_x(t) = norm(e_x(:,t));
    if isnan(os_ref) && norm_e_x(t) > 0
        os_ref = t;
    end
    % 2b) Next system state, next estimator state
    x(:,t+1) = xx(2,1:2*n)'; 
    xi_hat(:,t+1) = est.E*xi_hat(:,t) + est.N*u(:,t) + est.M*ys_hat;
end

% 3) Conclusion [2a)]
t = Tsim+1;
[~,yy] = METANET(pp,x(:,t),u(:,t));
y(:,t) = yy(1,1:n)';
ys(:,t) = S*y(:,t);
ys_hat = ys(:,t) - Ds*u(:,t);
x_hat(:,t) = est.Vinv*[xi_hat(:,t)-est.L*ys_hat; ys_hat];
y_hat(:,t) = sys.C*x_hat(:,t) + sys.D*u(:,t);
e_x(:,t) = x_hat(:,t) - y(:,t);
norm_e_x(t) = norm(e_x(:,t));
if isnan(os_ref) && norm_e_x(t) > 0
    os_ref = t;
end


end