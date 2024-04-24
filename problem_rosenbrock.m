%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ---- Extended Rosenbrock, Problem Set-Up ------------------------ %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the objective function, initial guess, derivative information,
% and any desired constraint functions.
%
% Then, define the options structure to be input into madsqn.
%
% Goal to minimize the Extended Rosenbrock function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---- Objective Function and Initial Guess ------------------------------
% -- Specify the problem dimension, i.e. the number of control parameters
n = 2; % dimension, must be an even number for rosenbrock.

% -- Specify function to minimize
f = @(x) rosenbrock(n,x);

% -- Specify surrogate function
f_sur = @(x) 1;

% -- Specify initial guess
% (this is a test x0, see "Dennis and Schnabel. Numerical methods 
% for unconstrained optimization and nonlinear equations. 1983")
if mod(n,2) == 0
    temp = repmat([-1 2 1],1,ceil(n/3));
    x0 = temp(1:n);
else
    error('n must be an even positive integer.')
end
x0 = x0'; % make sure correct vector dimensions

% -- Can use other initial guesses, if desired
% x0 = 100*x0';
% x0 = [1 0]';

% -- Define true minimizer
x_min = ones(n,1);


% -- Specify gradient and hessian of f (if known)
if n == 2
    grad_f = @(x) [-400*x(1).*(x(2)-x(1).^2)-2*(1-x(1)); 200*(x(2)-x(1).^2)];
    hess_f = @(x) [802*x(1) - 400*(x(2)-x(1).^2) , -400*x(1);
                   -400*x(1) , 200*x(2)];
end

% -- Define constraint functions (if desired)
con1 = @(x) x(1) - 1.1;
con2 = @(x) 0.25 - x(2);


%% ---- Options -----------------------------------------------------------
% -- Build options (input problem dimension)
options = build_options(n);

% -- Specify optimization options
% Unspecified options will remain their default values
options.maxIter = 1000; % maximum iterations
options.f_grad = grad_f; % gradient of f
options.f_hess = hess_f; % hessian of f

options.displayIterInfo = true; % display iteration info during optimization

% -- Specify which search steps to use, with any desired options
% options.searchStep.samplingSearch = true; % sampling search
% options.numSurrogateSamples = 10;
options.searchStep.quasiNewton = true; % QN search
options.derivativeMethod = 'CS_1'; % type of FD scheme to use (if needed)
options.searchStep.quasiNewtonStartTol = 1e-4; % when QN search kicks in

% -- Specify constraints in options structure
% Can specify unrelaxable or relaxable constraints.
% options.constraints.nlineqConstraint{1} = con1;
% options.constraints.nlineqConstraint{2} = con2;

% options.constraints.nlineqConstraint_relax{1} = con2



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---- EXTENDED ROSENBROCK FUNCTION --------------------------------------
% A common test problem for minimization solvers
function out = rosenbrock(n,x)
% INPUTS:
% n -- dimension of space
% x -- point to evaluate on
%
% OUTPUTS:
% out -- value of rosenbrock function at each x

% Define rosenbrock function:
F = zeros(n,1);
for i = 1:n/2 
    F(2*i-1) = 10*(x(2*i) - x(i).^2);
    F(2*i) = 1 - x(2*i-1);
end

out = sum(F.^2);

end



