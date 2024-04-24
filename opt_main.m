%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ---- Optimize nonlinear function f(x) using MADS_SMF_CON -------- %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Workflow:
%   Step 1) Initialize problem to be solved using 'problem_function'
%           This will contain solver options as well!
%
%   Step 2) Solve using MADS-CON
%
%   Step 3) Postprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---- Initialize Problem ------------------------------------------------
clear all;
close all;
addpath('MADS_SMF_CON','nurbs')
global figcount history
figcount = 1;

%%%% Rosenbrock function
problem_rosenbrock;

%%%% Curvematching problem
% problem_curvematch;

tic;
%% ---- Solve -------------------------------------------------------------
[x_sol, f_min, history, options] = madsqn(f, x0,options);
time_to_solve = toc;
f_count = sum(history.fCount); % total objective function count

disp('MADS-QN Solution:')
disp(x_sol)
disp('Total number of objective function calls:')
disp(f_count)


% -- Postprocessing
%%%% Rosenbrock function
plot_rosenbrock;

%%%% Curvematching problem
% plot_curvematch;

