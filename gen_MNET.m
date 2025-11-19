%% data generation - METANET MODEL

clearvars
close all
clc

path = '.\SANCARLOcE\';
% path = '/home/marcof/MEGA/v03/SANCARLOcE/'; % remote server
% path = '/home/marco/Desktop/OneDrive_2025-10-28/Street_sensor_selection/Code/v04/SANCARLOcE/';

A = table2array(struct2table(load([path '_A_sancarlo_c.mat'])));
B = table2array(struct2table(load([path '_B_sancarlo_c.mat'])));

n = size(A,1);
p.n = n;
p.B = B;

p.Ts = 1;
p.rho_crit = 30 /1000; % > 10 and < 150
% p.kappa = 10 /1000; % >0 small
% p.anti = 0; % 50 % >= 0
% p.conv = 0; % 1 % either 0 or 1
p.vfree = 50 /3.6; % speed limits
p.tau = 5; % > 0
p.a = 2; % between 1 and 3
p.eta_tilde = 0.5; % this is \nu in the monograph
p.kappa_tilde = 150;

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

p.DL = inv(diag(street_lengths));
p.PU = A-eye(n); 
p.PD = A'-eye(n); 
nx = 0;
%x0 = sqrt(0.1/60)*ones(2*n,1);
x0 = zeros(2*n,1);

HOURS = 1;
Tsim = round(HOURS*3600/p.Ts);
tspan = 0:p.Ts:Tsim;

u = (10 /60)*ones(length(tspan),1);
%u = u+(10 /60)*0.1*randn(size(u));
u(1:nx) = zeros(nx,1);

[x,y] = METANET(p,x0,u);


% Final time:
% rho in veh/m | v in km/h | flow in veh/min
[x(end,1:25)' 3.6*x(end,26:end)' 60*y(end,1:25)']


figure
grid on
hold on
plot(tspan,60*y)
plot(tspan,60*sum(y,2),'k')

save([path 'flux_data.mat'],'tspan','x','y','u','p','x0')

