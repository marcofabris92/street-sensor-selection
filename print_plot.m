%% function print_plot_SANCARLOc()

clearvars
close all
clc

%% Plot maker. Methods:
% 1) exhaustive_selection(sys,@f,metric,p_star);
% 2) random_selection(sys,f,metric,p_star,alpha)
% 3) greedy_selection(sys,@f,metric,p_star);
% 4) greedy_exclusion(sys,@f,metric,p_star);
% 5) genetic_selection(sys,metric,p_star);
% 6) divide_conquer_selection(sys,@f,metric,p_star);
% 7) data_driven_selection(sys,p_star);

folder = 'SANCARLOcE';
path = ['.\' folder '\'];
% path = ['./' folder '/'];
% path = [];
A = table2array(struct2table(load([path '_A_sancarlo_c.mat'])));
if strcmp(folder,'SANCARLOcE')
    A = table2array(struct2table(load([path '_A_est_sancarlo_c.mat'])));
end
n = size(A,1);
n_2 = floor((n-1)/2);

Tsim_prime = 11; % must be shorter (or equal) than Tsim used in simulations
tt = 0:Tsim_prime-1;
ttt = 0:3*3600;

methods = 1:7; 
for method = methods
    load([path 'SANCARLOc_error___' num2str(method) '.mat'],'error')
    load([path 'SANCARLOc_errorNL_' num2str(method) '.mat'],'errorNL')
    load([path 'SANCARLOc_obsdet__' num2str(method) '.mat'],'obsdet')
    load([path 'SANCARLOc_rate____' num2str(method) '.mat'],'rate')

    % if method == 1 || method == 2 || method == 3 || method == 4 || method == 5
    %     metric = 5; % det(W)^(1/n)
    %     errorFig(1:n_2,n,n_2,method,metric,error,tt,path,0);
    %     errorFig(n_2+1:n-1,n,n_2,method,metric,error,tt,path,0);
    %     errorFig(1:n_2,n,n_2,method,metric,errorNL,ttt,path,1);
    %     errorFig(n_2+1:n-1,n,n_2,method,metric,errorNL,ttt,path,1);
    % end
    
    for metric = getMetricsFromMethod(method)
        errorFig(1:n_2,n,n_2,method,metric,error,tt,path,0);
        errorFig(n_2+1:n-1,n,n_2,method,metric,error,tt,path,0);
        errorFig(1:n_2,n,n_2,method,metric,errorNL,ttt,path,1);
        errorFig(n_2+1:n-1,n,n_2,method,metric,errorNL,ttt,path,1);
    end

    rateFig(method,(1:n)',obsdet,rate,path);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = errorFig(p_stars,n,n_2,method,metric,error,tt,path,NL)

ftsz = 20;
intr = 'latex';

fig = figure();
grid on
hold on
% seed = 123;
c = 0;
h = gobjects(length(p_stars),1);
l = cell(length(p_stars),1);
for p_star = p_stars
    c = c + 1;
    scale = 60/(n-p_star);
    h(c) = plot(tt,scale*squeeze(error(metric,p_star,1+tt)),...
        'Color',getColor(),'LineWidth',1.5);
    l{c} = ['$p^{\star} = ' num2str(p_star) '$'];
end
xlabel('$t$ [s]','Interpreter',intr,'FontSize',ftsz)
ylabel('$\left\| e^{s}_{x}(t)  \right\| \quad $ [veh/min]','Interpreter',intr,'FontSize',ftsz)
set(gca, 'YScale', 'log')
xaxisproperties = get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = intr;
xaxisproperties.FontSize = ftsz;
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = intr;   
yaxisproperties.FontSize = ftsz;
ylim([1e-18 1e3])
location = 'northeast';
if NL
    location = 'southwest';
    xlim([0 8000])
    ylim([1e-7 1e3])
end
My_LGD = legend(h,l,'Location',location,'Interpreter',intr,'FontSize',ftsz);

My_LGD.NumColumns = ceil(n_2/2);


NLstr = [];
if NL
    NLstr = ' (NL)';
end
title(['Error dynamics' NLstr '. Approach: ' getMethod(method) '. Metric: '...
    getMetric(metric) '.'],'FontSize',ftsz,'Interpreter',intr)

if ~isempty(path)
    currentFolder = pwd;
    cd([currentFolder path(2:end)])
    filename = ['.\figs\fig_' num2str(metric) '_error_'];
    % filename = ['./figs/fig_' num2str(metric) '_error_'];
    if NL
        filename = [filename 'NL_'];
    end
    if p_stars(end) <= n_2
        saveas(fig,[filename 'mthd_' num2str(method) '_m1-' num2str(n_2)],'fig')
    else
        saveas(fig,[filename 'mthd_' num2str(method) '_m' num2str(n_2+1) '-' num2str(n-1)],'fig')
    end
    cd(currentFolder)
end

function color = getColor()
    % % Index-based pseudo-random generator
    % % Returns ONLY the c-th random number in [0,1]
    % 
    % % classical LCG parameters (Numerical Recipes)
    % % a = 1664525;
    % % q = 1013904223;
    % % m = 2^32;
    % 
    % % compute x_c directly with modular arithmetic
    % color = seed;
    % for k = 1:c
    %     color = mod(1664525*color + 1013904223, 4294967296);
    % end
    % color = color / 4294967296;   % normalized to [0,1]
    % 
    % color = hsv2rgb([color 0.7 0.9]);        % converts to RGB

    normalization = 255;
    crimson = [220 20 60]/normalization;
    % red = [255 0 0]/normalization;
    salmon = [250 128 114]/normalization;
    % orangered = [255 69 0]/normalization;
    % darkorange = [255 165 0]/normalization;
    gold = [255 215 0]/normalization;
    green = [0 128 0]/normalization;
    % olive = [128 128 0]/normalization;
    lime = [0 255 0]/normalization;
    % teal = [0 128 128]/normalization;
    % cyan = [0 255 255]/normalization;
    dodgerblue = [20 144 255]/normalization;
    % navy = [0 0 128]/normalization;
    blue = [0 0 255]/normalization;
    darkviolet = [148 0 211]/normalization;
    % purple = [128 0 128]/normalization;
    magenta = [255 0 255]/normalization;
    maroon = [128 0 0]/normalization;
    % saddlebrown = [139 69 19]/normalization;
    slategray = [112 128 144]/normalization;
    black = [0 0 0]/normalization;
    switch c
        case 1
            color = crimson;
        case 2
            color = green;
        case 3
            color = gold;
        case 4
            color = darkviolet;
        case 5
            color = dodgerblue;
        case 6
            color = lime;
        case 7
            color = salmon;
        case 8
            color = maroon;
        case 9
            color = slategray;
        case 10
            color = blue;
        case 11
            color = magenta;
        otherwise
            color = black;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = rateFig(method,p_stars,obsdet,rate,path)

ftsz = 20;
intr = 'latex';

fig = figure();
grid on
hold on

metrics = getMetricsFromMethod(method);
h = gobjects(length(metrics),1);
l = cell(length(metrics),1);
c = 0;
rank_flag = 0;
if metrics(1) == 1
    rank_flag = 1;
    metrics = [metrics(2:end) 1];
    c = c + 1;
end
for metric = metrics
    if rank_flag && metric == 1
        c = 1;
    else
        c = c + 1;
    end
    [mrk,mrks,mrkec,mrkfc] = getMarker(metric,obsdet);
    h(c) = plot(-100,rate(metric,1),'color',[1 1 1],'Marker',mrk,...
        'MarkerSize',ftsz,...
        'MarkerEdgeColor',mrkec,'MarkerFaceColor',mrkfc);
    scatter(p_stars,rate(metric,:),...
        mrks,'MarkerFaceAlpha',0.45,...
        'Marker',mrk,'MarkerEdgeColor',mrkec,'MarkerFaceColor',mrkfc);
    l{c} = getMetric(metric);
end

xticks(p_stars)
xlim([1 p_stars(end)])
ylim([0 1])
xlabel('$p^{\star}$','Interpreter',intr,'FontSize',ftsz)
ylabel('$\rho_{\xi}$','Interpreter',intr,'FontSize',ftsz)
xaxisproperties = get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = intr;
xaxisproperties.FontSize = ftsz;
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = intr;   
yaxisproperties.FontSize = ftsz;
legend(h,l,'Interpreter',intr,'FontSize',ftsz);
%title(['Rates. Approach: ' getMethod(method) '.'],'FontSize',ftsz,'Interpreter',intr)

if ~isempty(path)
    currentFolder = pwd;
    cd([currentFolder path(2:end)])
    saveas(fig,['.\figs\fig_rate_' num2str(method)],'fig')
    % saveas(fig,['./figs/fig_rate_' num2str(method)],'fig')
    cd(currentFolder)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mrk,mrks,mrkec,mrkfc] = getMarker(metric,obsdet)

normalization = 255;

% crimson = [220 20 60]/normalization;
red = [255 0 0]/normalization;
% salmon = [250 128 114]/normalization;
% orangered = [255 69 0]/normalization;
% darkorange = [255 165 0]/normalization;
gold = [255 215 0]/normalization;
green = [0 128 0]/normalization;
% olive = [128 128 0]/normalization;
% lime = [0 255 0]/normalization;
% teal = [0 128 128]/normalization;
cyan = [0 255 255]/normalization;
% dodgerblue = [20 144 255]/normalization;
% navy = [0 0 128]/normalization;
blue = [0 0 255]/normalization;
% darkviolet = [148 0 211]/normalization;
% purple = [128 0 128]/normalization;
magenta = [255 0 255]/normalization;
% maroon = [128 0 0]/normalization;
% saddlebrown = [139 69 19]/normalization;
% slategray = [112 128 144]/normalization;
black = [0 0 0]/normalization;

circle = 'o';
% plus = '+';
% asterisk = '*';
% point = '.';
cross = 'x';
% horizontal_line = '_';
% vertical_line = '|';
square = 's';
diamond = 'd';
upward_triangle = '^';
downward_triangle = 'v';
% rightward_triangle = '>';
% leftwrd_triangle = '<';
pentagram = 'p'; % star with 5 spikes
hexagram = 'h'; % star with 6 spikes
none = 'none'; % not applicable

mrk = none;
mrks = 30*4.^(obsdet(metric,:));
mrkec = black;
mrkfc = black;

switch metric
    case 1
        mrk = cross;
    case 2
        mrk = circle;
        mrkfc = magenta;
    case 3
        mrk = upward_triangle;
        mrkfc = blue;
    case 4
        mrk = downward_triangle;
        mrkfc = red;
    case 5
        mrk = square;
        mrkfc = cyan;
    case 6
        mrk = pentagram;
        mrkfc = gold;
    case 7
        mrk = diamond;
        mrkfc = green;
    case 8
        mrk = hexagram;
        mrkfc = gold;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = getMetric(metric)

switch metric
    case 1
        str = '$\mathrm{rank}[\mathcal{W}_n]$';
    case 2
        str = '$\mathrm{tr}_n[\mathcal{W}_n]$';
    case 3
        str = '$\mathcal{K}^{-1}[\mathcal{W}_n]$';
    case 4
        str = '$\lambda_{min}[\mathcal{W}_n]$';
    case 5
        str = '$\mathrm{det}_n[\mathcal{W}_n]$';
    case 6
        str = '$H_2[\mathcal{W}_{\infty}]$';
    case 7
        str = '$\ell d[\mathcal{W}_{\infty}]$';
    case 8
        str = '$H_2[\tilde{\mathcal{W}}_{\infty}]$'; 
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = getMethod(method)

switch method
    case 1
        str = 'exhaustive search';
    case 2
        str = 'random selection';
    case 3
        str = 'greedy selection';
    case 4
        str = 'greedy exclusion';
    case 5
        str = 'genetic algorithm';
    case 6
        str = 'divide \& conquer';
    case 7
        str = 'data-driven selection';
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function metrics = getMetricsFromMethod(method)

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

end