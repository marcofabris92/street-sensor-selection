% MAIN SANCARLOc
% Location: Italy, Padova, Arcella, San Carlo church, east area 
% (lat,lon) = (45.430373615189424, 11.888818838634402)°

clearvars
close all
clc

M = 0;
U = 10 /60;
switch M
    case 0
        path = '.\SANCARLOc';
        % path = '/home/marco/Desktop/OneDrive_2025-10-28/Street_sensor_selection/Code/v05/SANCARLOc'; % Lenovo
        % path = '/home/marcof/MEGA/v05/SANCARLOc'; % remote server
    case 1
        path = '.\SANCARLOcE';
        % path = '/home/marco/Desktop/OneDrive_2025-10-28/Street_sensor_selection/Code/v05/SANCARLOcE'; % Lenovo
        % path = '/home/marcof/MEGA/v05/SANCARLOcE'; % remote server
end
if U > 1/6
    path = [path '_cng'];
end
path = [path '\'];
% path = [path '/']; % Lenovo & remote server


%% Model generation/ selection
switch M
    case 0
        sys.A = table2array(struct2table(load([path '_A_sancarlo_c.mat'])));
    case 1
        sys.A = table2array(struct2table(load([path '_A_est_sancarlo_c.mat'])));
end
n = size(sys.A,1);
sys.C = eye(n);

[~,~,isDet] = obs(sys.A,sys.C,n,1,1,0,1);
if ~isDet
    error('The given system model is not detectable')
end

switch M
    case 0
        sys.B = table2array(struct2table(load([path '_B_sancarlo_c.mat'])));
    case 1
        sys.B = table2array(struct2table(load([path '_B_est_sancarlo_c.mat'])));
end
sys.D = zeros(n,size(sys.B,2));

sys = system_setup(sys); % sets further parameters, such as the lag

% change the default initial condition (=0), if you wish
% don't use 0.5 because it returns zero error... must be caused by 
% some very nice subspace
sys.x0 = 0.1/60*ones(n,1); % 0.1/60 % in August, 2025 the value was set to 1


eta = sys.eta;
save([path 'eta.mat'],'eta')
%error('stop')

%% Plot maker - method 6: "data-driven" is excluded AT THE MOMENT
% 1) exaustive_selection(sys,@f,metric,p_star);
% 2) random_selection(sys,f,metric,p_star,alpha)
% 3) greedy_selection(sys,@f,metric,p_star);
% 4) greedy_exclusion(sys,@f,metric,p_star);
% 5) genetic_selection(sys,metric,p_star);
% 6) divide_conquer_selection(sys,@f,metric,p_star);
% 7) data_driven_selection(sys,p_star);

Tsim = 3*60*60; % simulation time
fprintf('START\n')
for method = 3:7 %1:7
    compute_data(path,U,sys,Tsim,method);
end
fprintf('END')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [] = compute_data(path,U,sys,Tsim,method)

%% selection of the optimization method (opt. method is fixed)
switch method
    case 1 % exhaustive
        metrics = 1:7;
    case 2 % random
        metrics = 1:7;
    case 3 % greedy selection
        metrics = [1 2 5 6 7];
    case 4 % greedy exclusion
        metrics = [1 2 5 6 7];
    case 5 % genetic
        metrics = 1:7;
    case 6 % divide & conquer
        metrics = [2 6];
    case 7 % data-driven
        metrics = 8;
end


%% Initialization
n = length(sys.x0);
p = size(sys.C,1);
N_metrics = length(metrics);
rate = zeros(N_metrics,p-1);  
obsdet = zeros(N_metrics,p-1);
error = zeros(N_metrics,p-1,Tsim+1);
errorNL = zeros(N_metrics,p-1,Tsim+1);
conf = cell(N_metrics,p-1);

%% Main cycle: sensor placement
for metric = metrics
    parfor p_star = 1:p

        fprintf(['mthd, mtrc, p_str: ' num2str(method) ', ' num2str(metric) ', ' num2str(p_star) '\n'])

        % Sensor selection is carried out depending on the method and metric
        switch method
           case 1
            [Q_star,detectable,observable] = exaustive_selection(sys,@f,metric,p_star);
           case 2
            [Q_star,detectable,observable] = random_selection(sys,@f,metric,p_star,1);
           case 3
            [Q_star,detectable,observable] = greedy_selection(sys,@f,metric,p_star);
           case 4
            [Q_star,detectable,observable] = greedy_exclusion(sys,@f,metric,p_star);
           case 5
            [Q_star,detectable,observable] = genetic_selection(sys,metric,p_star);
           case 6
            [Q_star,detectable,observable] = divide_conquer_selection(sys,@f,metric,p_star);
           case 7
            [Q_star,detectable,observable] = data_driven_selection(sys,p_star);
        end
   
        conf{metric,p_star} = Q_star;
     
        % in case of detectability, find the performance
        if detectable
    
            obsdet(metric,p_star) = 1; % = 1: just detectable
            if observable
                obsdet(metric,p_star) = 2; % = 2: observable
            end

            %% Selection matrix computation
            S = zeros(p_star,p);
            for k = 1:p_star
                S(k,Q_star(k)) = 1;
            end

            % Luenberger gain + other estimator matrices
            [est,best_rho,~,~] = red_ord_est_design(sys,S);
            rate(metric,p_star) = best_rho;
   
            % System and estimator trajectories 
            % Measurements from linear model
            u = U*ones(1,Tsim+1); % input
            x_hat_0 = zeros(n,1); % initial guess for x_hat(0) 
            [~,~,~,~,~,~,norm_e_x] = red_ord_est_sim(sys,est,S,u,x_hat_0);
            error(metric,p_star,:) = norm_e_x;

            % System and estimator trajectories 
            % Measurements from nonlinear model
            [~,~,~,~,~,~,norm_e_x] = red_ord_est_sim_NL(sys,est,S,u,x_hat_0);
            errorNL(metric,p_star,:) = norm_e_x;
        
        end

    end
end

%% Saving data
rate

save([path 'SANCARLOc_errorNL_' num2str(method) '.mat'],'errorNL')
save([path 'SANCARLOc_error___' num2str(method) '.mat'],'error')
save([path 'SANCARLOc_obsdet__' num2str(method) '.mat'],'obsdet')
save([path 'SANCARLOc_rate____' num2str(method) '.mat'],'rate')
save([path 'SANCARLOc_conf____' num2str(method) '.mat'],'conf')

end

%%
% close all
% load([path 'SANCARLOc_errorNL_3.mat'])
% e = 60*squeeze(errorNL(5,9,:))/25;
% 
% figure
% grid on
% hold on
% plot(0:length(e)-1,e)