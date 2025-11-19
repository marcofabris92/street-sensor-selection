function [A,B,C,D] = model_generation()

% Instructions:
% To generate real model matrices:
% -- disable [A,B,C] = get_system_random(); in MAIN.m
% -- enable setenv('TCL_LIBRARY' ...
% -- enable setenv('TK_LIBRARY' ...
% -- enable pyrunfile("main.py")
% -- launch
% -- close the window with the map image
% -- done after 3 seconds

setenv('TCL_LIBRARY', 'C:\Users\Marco Fabris\AppData\Local\Programs\Python\Python310\tcl\tcl8.6')
setenv('TK_LIBRARY', 'C:\Users\Marco Fabris\AppData\Local\Programs\Python\Python310\tcl\tk8.6')

fprintf('Generation of the model\n')
pyrunfile("main.py")
A = table2array(struct2table(load("A.mat")));
B = table2array(struct2table(load("B.mat")));
C = table2array(struct2table(load("C.mat")));
D = zeros(size(C,1),size(B,2));
fprintf('Model has been generated.\n') 


end