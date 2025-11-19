clear all
close all
clc

path1 = '.\SANCARLOc\figs';
%path1 = '/home/marco/Desktop/OneDrive_2025-10-28/Street_sensor_selection/Code/v04/SANCARLOc/figs';
path2 = '.\SANCARLOcE\figs';
%path2 = '/home/marco/Desktop/OneDrive_2025-10-28/Street_sensor_selection/Code/v04/SANCARLOcE/figs';;

folder_path = path2; %uigetdir(pwd, path);
if folder_path == 0
    error('No folder selected.');
end

fig_files = dir(fullfile(folder_path, '*.fig'));

if isempty(fig_files)
    error('No file .fig found in the selected folder.');
end

fig_width = 1000;  
fig_height = 600;  


for k = 1:length(fig_files)

    
  
    fig_path = fullfile(folder_path, fig_files(k).name);
    
    fig = openfig(fig_path, 'invisible');  % don't show on the screen
   
    set(fig, 'Position', [100, 100, fig_width, fig_height]);
    
    set(fig, 'Color', 'white');
    set(fig, 'InvertHardcopy', 'off');
    

    lgd = findobj(fig, 'Type', 'Legend');
    if ~isempty(lgd)
        set(lgd, 'FontSize', 12);
    end

    
    if length(fig_files(k).name) >= 13 &&...
            strcmp(fig_files(k).name(13), 'm')
        ax = findobj(fig, 'Type', 'Axes');
        xticks(ax,0:1:10);
    end
    
    
    if length(fig_files(k).name) == 14
        set(fig, 'Position', [100, 100, fig_width+100, fig_height]);
        ax = findobj(fig, 'Type', 'Axes');
        yticks(ax,0:0.1:1);
    end

    if strcmp(folder_path,path2) && ...
            strcmp(fig_files(k).name,'fig_2_error_NL_mthd_5_m13-24.fig')
        ylim([10^-7 10^6])
    end

    [~, name, ~] = fileparts(fig_files(k).name);
    eps_name = fullfile(folder_path, [name, '.eps']);
    
    % Save as EPS
    print(fig, eps_name, '-depsc', '-r300');  % 300 dpi EPS format
    
    close(fig);
    
    fprintf(['k = ' num2str(k) ', '])
    fprintf('Saved: %s\n', eps_name);

    
   
end

disp('All files have been saved as .eps.');
