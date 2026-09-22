function [est,best_rho,best_os,best_t_os] = red_ord_est_design(sys,S)

% Useful dimensions
p_star = size(S,1);
n = size(sys.A,1);
p_bar = n-p_star;

% Change of basis
Cs = S*sys.C;
[Csbar,~] = qr(Cs');   
Csbar = Csbar(:,p_star+1:end)';

est.V = [Csbar; Cs];
est.Vinv = pinv(est.V);

% Handle the case where all sensors have been selected
if p_bar == 0
    best_rho = 0;
    best_os = 0;
    best_t_os = 0;
    est.L = zeros(0,n); 
    est.E = [];
    est.M = zeros(0,n);
    est.N = zeros(0,size(sys.B,2));
    return
end

F = est.V*sys.A*est.Vinv;
G = est.V*sys.B;
% Also: H = [zeros(p_star,p_bar) eye(p_star)];

F11 = F(1:p_bar,1:p_bar);
F12 = F(1:p_bar,p_bar+1:end);
F21 = F(p_bar+1:end,1:p_bar);
F22 = F(p_bar+1:end,p_bar+1:end);

G1 = G(1:p_bar,:);
G2 = G(p_bar+1:end,:);

%% Dealing with the dual reachable subset to find an estimator gain
Phi = F11';
Gamma = F21';
Reach = ctrb(Phi,Gamma);
r = rank(Reach);
[T,~] = qr(Reach,0);
Tinv = T';
Phi = Tinv*Phi*T;
Gamma = Tinv*Gamma;
Phi_11 = Phi(1:r,1:r);
Phi_12 = Phi(1:r,r+1:end);
Gamma_1 = Gamma(1:r,:);


while r > 1 && my_rank(ctrb(Phi_11,Gamma_1)) < r
    r = r - 1;
    Phi_11 = Phi(1:r,1:r);
    Phi_12 = Phi(1:r,r+1:end);
    Gamma_1 = Gamma(1:r,:);
end

%% Distinguish between stable and unstable poles of Phi_11:
% < k means that the values in poles_11 are asym. stable poles
% or poles on the unit circle
poles_11 = sort(eig(Phi_11));
k = 1;
while k <= length(poles_11)
    if abs(poles_11(k)) > 1
        break
    end
    k = k + 1;
end

%% Pole placement method: mirroring & alpha-attenuation

% Mirroring
for h = k:length(poles_11)
    poles_11(h) = 1/conj(poles_11(h));
end
rho_11 = max(abs(poles_11));

% Attenuating - preliminary phase
warning('off','all')
u = ones(1,n+1); % input
x_hat_0 = zeros(n,1);
best_K_bar_1 = zeros(size(Gamma_1,2),size(Gamma_1,1));
K_bar_2 = zeros(size(Gamma_1,2),size(Phi_12,2));
one = 1; % careful: this 'one' is what enforces estimation stability!
est.L = -([best_K_bar_1 K_bar_2]*Tinv)'; 
est.E = F11+est.L*F21;
best_rho = max(abs(eig(est.E)));
search_interval = -1:1e-2:1;
X = zeros(length(search_interval)+1,4);
X(1,:) = [NaN +Inf +Inf +Inf];
if best_rho < one
    est.M = F12+est.L*F22-F11*est.L-est.L*F21*est.L;
    est.N = G1+est.L*G2;
    [~,~,~,~,~,~,norm_e_x,os_ref] = red_ord_est_sim(sys,est,S,u,x_hat_0);
    [max_norm_e_x,t_os]= max(norm_e_x);
    num_os = max_norm_e_x-norm_e_x(os_ref);
    if num_os <= 0
        X(1,:) = [NaN best_rho -Inf t_os];
    else
        X(1,:) = [NaN best_rho log10(num_os/norm_e_x(os_ref)) t_os];
    end
end

% Attenuating - optimization phase (search)
best_K_bar_2 = pinv(Gamma_1)*Phi_12;
if r > 0
    i = 1;
    for alpha = search_interval
        i = i + 1;
        try
            K_bar_1 = place(Phi_11,Gamma_1,alpha/rho_11*poles_11);
            
            est.L = -([K_bar_1 best_K_bar_2]*Tinv)'; 
            est.E = F11+est.L*F21;
            rho = max(abs(eig(est.E)));
            if rho < one
                est.M = F12+est.L*F22-F11*est.L-est.L*F21*est.L;
                est.N = G1+est.L*G2;
                [~,~,~,~,~,~,norm_e_x,os_ref] =...
                    red_ord_est_sim(sys,est,S,u,x_hat_0);

                % fig1 = figure;
                % grid on
                % hold on
                % tt = 0:Tsim;
                % plot(tt,norm_e_x,'r')
                % set(gca, 'YScale', 'log')
                % pause
                % close(fig1)

                [max_norm_e_x,t_os] = max(norm_e_x);
                num_os = max_norm_e_x-norm_e_x(os_ref);
                if num_os <= 0
                    X(i,:) = [alpha rho -Inf t_os];
                else
                    X(i,:) = [alpha rho log10(num_os/norm_e_x(os_ref)) t_os];
                end
            else
                X(i,:) = [alpha +Inf +Inf +Inf];
            end
        catch ME
            X(i,:) = [alpha +Inf +Inf +Inf];
        end
    end
elseif p_star == n
    best_rho = 0;
end

% Attenuating - optimization phase (comparison)
best_alpha = NaN;
best_t_os = +Inf;
best_os = min(X(:,3));
min_os_indexes = find(X(:,3) == best_os);
for i_ = 1:length(min_os_indexes)
    i = min_os_indexes(i_);
    if isnan(X(i,1)) || X(i,2) < min(best_rho,one) ||...
            X(i,2) <= best_rho && X(i,2) < one && X(i,4) < best_t_os
        best_alpha = X(i,1);
        best_rho = X(i,2);
        best_t_os = X(i,4);
    end
end
if ~isnan(best_alpha) && p_star < n
    best_K_bar_1 = place(Phi_11,Gamma_1,best_alpha/rho_11*poles_11);
else
    best_K_bar_2 = K_bar_2;
end

% best_K_bar_1 
warning('on','all')

%% Gains
est.L = -([best_K_bar_1 best_K_bar_2]*Tinv)'; 
est.E = F11+est.L*F21;
est.M = F12+est.L*F22-F11*est.L-est.L*F21*est.L;
est.N = G1+est.L*G2;


end