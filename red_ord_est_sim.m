function [x,y,x_hat,ys,y_hat,e_x,norm_e_x,os_ref] = ...
    red_ord_est_sim(sys,est,S,u,x_hat_0)

% Useful dimensions
Tsim = size(u,2)-1;
n = length(sys.x0);
p = size(sys.C,1);
p_star = size(S,1);
p_bar = n-p_star;

% Initializations
Ds = S*sys.D;
x = zeros(n,Tsim+1);
y = zeros(p,Tsim+1);
ys = zeros(p_star,Tsim+1);
xi_hat = zeros(p_bar,Tsim+1);
x_hat = zeros(n,Tsim+1);
y_hat = zeros(p,Tsim+1);
e_x = zeros(n,Tsim+1);
norm_e_x = zeros(1,Tsim+1);
os_ref = NaN;

%% Dynamics of the reduced order estimator

% 1) Start [2b)]
t = 1;
x(:,t) = sys.x0;
x_bar_0 = est.V*sys.x0;
x_bar_hat_0 = est.V*x_hat_0;
xi_hat(:,t) = x_bar_hat_0(1:p_bar) + est.L*x_bar_0(p_bar+1:end);

% 2) Updates
for t = 1:Tsim
    % 2a) Outputs, estimates and errors
    y(:,t) = sys.C*x(:,t) + sys.D*u(:,t);
    ys(:,t) = S*y(:,t);
    ys_hat = ys(:,t) - Ds*u(:,t);
    x_hat(:,t) = est.Vinv*[xi_hat(:,t)-est.L*ys_hat; ys_hat];
    y_hat(:,t) = sys.C*x_hat(:,t) + sys.D*u(:,t);
    e_x(:,t) = x_hat(:,t) - x(:,t);
    norm_e_x(t) = norm(e_x(:,t));
    if isnan(os_ref) && norm_e_x(t) > 0
        os_ref = t;
    end
    % 2b) Next system state, next estimator state
    x(:,t+1) = sys.A*x(:,t) + sys.B*u(:,t); 
    xi_hat(:,t+1) = est.E*xi_hat(:,t) + est.N*u(:,t) + est.M*ys_hat;
end

% 3) Conclusion [2a)]
t = Tsim+1;
y(:,t) = sys.C*x(:,t) + sys.D*u(:,t);
ys(:,t) = S*y(:,t);
ys_hat = ys(:,t) - Ds*u(:,t);
x_hat(:,t) = est.Vinv*[xi_hat(:,t)-est.L*ys_hat; ys_hat];
y_hat(:,t) = sys.C*x_hat(:,t) + sys.D*u(:,t);
e_x(:,t) = x_hat(:,t) - x(:,t);
norm_e_x(t) = norm(e_x(:,t));
if isnan(os_ref) && norm_e_x(t) > 0
    os_ref = t;
end


end