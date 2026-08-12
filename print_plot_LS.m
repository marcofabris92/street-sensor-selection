%% function print_plot_ARCELLA()

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

folder = 'ARCELLA';
path = ['.\' folder '\'];
% path = ['./' folder '/'];
% path = [];
sys.A = table2array(struct2table(load([path '_A_arcella.mat'])));
n = size(sys.A,1);
n_2 = floor((n-1)/2);
% sys.B = table2array(struct2table(load([path '_B_arcella.mat'])));
% sys.C = eye(n);
% sys.D = zeros(n,size(sys.B,2));
% sys.x0 = zeros(n,1);
% sys = heuristic_system_setup(sys);

Tsim_prime = 11; % must be shorter (or equal) than Tsim used in simulations
tt = 0:Tsim_prime-1;
ttt = 0:500;

p_star = 23;
load([path 'ARCELLA_error___0.mat'],'error')
error0 = squeeze(error(8,p_star,:));


methods = [2 3 4 5 6 7]; 
for method = methods
    load([path 'ARCELLA_conf____' num2str(method) '.mat'],'conf')
    load([path 'ARCELLA_error___' num2str(method) '.mat'],'error')
    load([path 'ARCELLA_obsdet__' num2str(method) '.mat'],'obsdet')
    load([path 'ARCELLA_rate____' num2str(method) '.mat'],'rate')

    fprintf(['Method: ' getMethod(method) '\n'])
    conf = cell2mat(conf);
    reshape(conf,p_star,[])

    eval(['error' num2str(method) ' = error;']);

    % if method == 1 || method == 2 || method == 3 || method == 4 || method == 5
    %     metric = 5; % det(W)^(1/n)
    %     errorFig(1:n_2,n,n_2,method,metric,error,tt,path,0);
    %     errorFig(n_2+1:n-1,n,n_2,method,metric,error,tt,path,0);
    %     errorFig(1:n_2,n,n_2,method,metric,errorNL,ttt,path,1);
    %     errorFig(n_2+1:n-1,n,n_2,method,metric,errorNL,ttt,path,1);
    % end
    
    % for metric = getMetricsFromMethod(method)
    %     errorFig(1:n_2,n,n_2,method,metric,error,tt,path,0);
    %     errorFig(n_2+1:n-1,n,n_2,method,metric,error,tt,path,0);
    %     errorFig(1:n_2,n,n_2,method,metric,errorNL,ttt,path,1);
    %     errorFig(n_2+1:n-1,n,n_2,method,metric,errorNL,ttt,path,1);
    % end

    %rateFig(method,(1:n)',obsdet,rate,path);
end

ftsz = 20;
intr = 'latex';

fig = figure();
scale = 60/(n-p_star);
grid on
hold on


K = 20;

hh = [];
for metric = 1:7
    h = plot([0 0],[0 0],getMarker(metric),'Color','k','LineWidth',2);
    hh = [hh h];
end

h = plot(ttt,scale*error0(1+ttt),'--','Color',getColor(0),'LineWidth',1.5);
hh = [hh h];

% RANDOM approach
h = plot(ttt,scale*squeeze(error2(1,p_star,1+ttt)),...
        'Color',getColor(2),'LineWidth',1.5);
hh = [hh h];
for metric = 1:7 % red
    plot(ttt,scale*squeeze(error2(metric,p_star,1+ttt)),...
        'Color',getColor(2),'LineWidth',1.5)
    ttt_ = downsample(ttt,50+(metric-1)*K);
    plot(ttt_,scale*squeeze(error2(metric,p_star,1+ttt_)),...
        getMarker(metric),'Color',getColor(2),'LineWidth',2)
end
% GREEDY SELECTION approach
h = plot(ttt,scale*squeeze(error3(1,p_star,1+ttt)),...
        'Color',getColor(3),'LineWidth',1.5);
hh = [hh h];
for metric = [1 2 5 6 7] % cyan
    plot(ttt,scale*squeeze(error3(metric,p_star,1+ttt)),...
        'Color',getColor(3),'LineWidth',1.5)
    ttt_ = downsample(ttt,50+(metric-1)*K+30);
    plot(ttt_,scale*squeeze(error3(metric,p_star,1+ttt_)),...
        getMarker(metric),'Color',getColor(3),'LineWidth',2)
end
% GREEDY EXCLUSION approach 
h = plot(ttt,scale*squeeze(error4(1,p_star,1+ttt)),...
        'Color',getColor(4),'LineWidth',1.5);
hh = [hh h];
for metric = [1 2 5 6 7] % magenta
    plot(ttt,scale*squeeze(error4(metric,p_star,1+ttt)),...
        'Color',getColor(4),'LineWidth',1.5)
    ttt_ = downsample(ttt,50+(metric-1)*K+60);
    plot(ttt_,scale*squeeze(error4(metric,p_star,1+ttt_)),...
        getMarker(metric),'Color',getColor(4),'LineWidth',2)
end
% GENETIC approach
h = plot(ttt,scale*squeeze(error5(1,p_star,1+ttt)),...
        'Color',getColor(5),'LineWidth',1.5);
hh = [hh h];
for metric = 1:7 % blue
    plot(ttt,scale*squeeze(error5(metric,p_star,1+ttt)),...
        'Color',getColor(5),'LineWidth',1.5)
    ttt_ = downsample(ttt,50+(metric-1)*K+90);
    plot(ttt_,scale*squeeze(error5(metric,p_star,1+ttt_)),...
        getMarker(metric),'Color',getColor(5),'LineWidth',2)
end
% DIVIDE AND CONQUER approach
h = plot(ttt,scale*squeeze(error6(2,p_star,1+ttt)),...
        'Color',getColor(6),'LineWidth',1.5);
hh = [hh h];
for metric = [2 6] % green
    plot(ttt,scale*squeeze(error6(metric,p_star,1+ttt)),...
        'Color',getColor(6),'LineWidth',1.5)
    ttt_ = downsample(ttt,50+(metric-1)*K+120);
    plot(ttt_,scale*squeeze(error6(metric,p_star,1+ttt_)),...
        getMarker(metric),'Color',getColor(6),'LineWidth',2)
end
% DATA-DRIVEN approach
h = plot(ttt,scale*squeeze(error7(1,p_star,1+ttt)),...
        'Color',getColor(7),'LineWidth',1.5);
hh = [hh h];
for metric = 8 % gold
    plot(ttt,scale*squeeze(error7(metric,p_star,1+ttt)),...
        'Color',getColor(7),'LineWidth',1.5)
    ttt_ = downsample(ttt,50+(metric-1)*K+150);
    plot(ttt_,scale*squeeze(error7(metric,p_star,1+ttt_)),...
        getMarker(metric),'Color',getColor(7),'LineWidth',2)
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
ylim([1e-16 1e-2])
location = 'northeast';
l = {getMetric(1),...
    getMetric(2),...
    getMetric(3),...
    getMetric(4),...
    getMetric(5),...
    getMetric(6),...
    getMetric(7),...
    'Manual',...
    'Random','Greedy sel.','Greedy excl.','Genetic','D. \& C.', 'Data-driven'};
My_LGD = legend(hh,l,'Location',location,'Interpreter',intr,'FontSize',ftsz);

function color = getColor(c)
  

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
    cyan = [0 255 255]/normalization;
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
            color = '';
        case 2
            color = crimson; 
        case 3
            color = cyan;
        case 4
            color = magenta;
        case 5
            color = blue;
        case 6
            color = green;
        case 7
            color = gold; 
        otherwise
            color = black;
    end
end

function mrkr = getMarker(c)

switch c
    case 1
        mrkr = 'x';
    case 2
        mrkr = 'o';
    case 3
        mrkr = '^';
    case 4
        mrkr = 'v';
    case 5
        mrkr = 'square';
    case 6
        mrkr = 'pentagram';
    case 7
        mrkr = 'hexagram';
    otherwise
        mrkr = 'pentagram';
end

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        str = 'divide & conquer';
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