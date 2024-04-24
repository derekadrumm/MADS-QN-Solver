%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ---- OPTIONS STRUCTURE INITIALIZATION --------------------------- %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define options structure and set default options values
%
% Should be no need to mess with this, unless want to modify solver.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function options = build_options(n)
% INPUTS:
% n -- dimension of parameter search space
%
% OUTPUTS:
% options -- options structure

%% ---- Define default optimization parameters ----------------------------
options = struct();

% display iteration information
if ~isfield(options, 'displayIterInfo')
    options.displayIterInfo = false;
end
% stopping function tolerance
if ~isfield(options, 'functionTolerance')
    options.functionTolerance = 1e-6;
end
% stopping mesh tolerance
if ~isfield(options, 'meshTolerance')
    options.meshTolerance = 1e-6;
end
% stopping frame tolerance
if ~isfield(options, 'frameTolerance')
    options.frameTolerance = 1e-6;
end
% stopping relative mesh-gradient tolerance
if ~isfield(options, 'relativeMeshGradientTolerance')
    options.relativeMeshGradientTolerance = 1e-18;
end
% maximum number of iterations
if ~isfield(options, 'maxIter')
    options.maxIter = 100;
end
% maximum step size
if ~isfield(options, 'maxStepSize')
    options.maxStepSize = Inf;
end
% frame size parameter (in (0,Inf))
if ~isfield(options, 'frameSize')
    options.frameSize = 1;
end
% mesh size adjustment parameter (in (0,1))
if ~isfield(options, 'meshSizeAdjust')
    options.meshSizeAdjust = 0.5;
end
% sampling method to use during MADS
if ~isfield(options, 'sampleType')
    options.sampleType = 'SES';
end
% number of surrogate function search sample points on each step
if ~isfield(options, 'numSurrogateSamples')
    options.numSurrogateSamples = 0;
end
% maximum times the surrogate search step can fail
if ~isfield(options, 'maxSurrogateFailure')
    options.maxSurrogateFailure = 100;
end
% define search direction matrix (collumns are search directions)
if ~isfield(options, 'searchDirections')
    options.searchDirections = [eye(n) -eye(n)];
end


%% ---- Define function info ----------------------------------------------
% function to optimize (may not be used from options, maybe nice to have)
if ~isfield(options, 'f')
    options.f = @(x) 0;
end
% function gradient (if empty, must specify numerical method)
if ~isfield(options, 'f_grad')
    options.f_grad = @(x) zeros(n,1);
end
% function hessian (if empty, must specify numerical method)
if ~isfield(options, 'f_hess')
    options.f_hess = @(x) zeros(n);
end


% surrogate function of objective function
if ~isfield(options, 'f_sur')
    options.f_sur = @(x) 0;
end
% surrogate function gradient (if empty, must specify numerical method)
if ~isfield(options, 'f_sur_grad')
    options.f_sur_grad = @(x) zeros(n,1);
end
% surrogate function hessian (if empty, must specify numerical method)
if ~isfield(options, 'f_sur_hess')
    options.f_sur_hess = @(x) zeros(n);
end


% derivative calculation method (can maybe implement multiple methods)
if ~isfield(options, 'derivativeMethod')
    options.derivativeMethod = 'FS_1'; % default 1st order fwd diff
end


%% ---- Define search step methods ----------------------------------------
if ~isfield(options, 'searchStep')
    options.searchStep = struct();
end
% line search
if ~isfield(options.searchStep, 'lineSearch')
    options.searchStep.lineSearch = false;
end
% line search counter
if ~isfield(options.searchStep, 'lineSearchCount')
    options.searchStep.lineSearchCount = 0;
end
% trust region (hook)
if ~isfield(options.searchStep, 'hook')
    options.searchStep.hook = false;
end
% trust region (hook) counter
if ~isfield(options.searchStep, 'hookCount')
    options.searchStep.hookCount = 0;
end
% trust region (double dog leg)
if ~isfield(options.searchStep, 'doubleDogLeg')
    options.searchStep.doubleDogLeg = false;
end
% trust region (double dog leg) counter
if ~isfield(options.searchStep, 'doubleDogLegCount')
    options.searchStep.doubleDogLegCount = 0;
end
% sampling search
if ~isfield(options.searchStep, 'samplingSearch')
    options.searchStep.samplingSearch = false;
end
% sampling search counter
if ~isfield(options.searchStep, 'samplingSearchCount')
    options.searchStep.samplingSearchCount = 0;
end


% QuasiNewton Search and parameters
if ~isfield(options.searchStep, 'quasiNewton')
    options.searchStep.quasiNewton = false;
end
% QuasiNewton search size tolerance (need mesh size < than tol to us QN)
if ~isfield(options.searchStep, 'quasiNewtonStartTol')
    options.searchStep.quasiNewtonStartTol = Inf;
end
% QuasiNewton search alpha condition
if ~isfield(options.searchStep, 'quasiNewtonAlpha')
    options.searchStep.quasiNewtonAlpha = 1e-4;
end
% QuasiNewton search counter
if ~isfield(options.searchStep, 'quasiNewtonCount')
    options.searchStep.quasiNewtonCount = 0;
end


%% ---- Define optimization constraints: ----------------------------------
if ~isfield(options, 'constraints')
    options.constraints = struct();
end
% lower bounds of x, default unbounded
if ~isfield(options.constraints, 'lb')
    options.constraints.lb = -Inf*ones(n,1);
end
% upper bounds of x, default unbounded
if ~isfield(options.constraints, 'ub')
    options.constraints.ub = Inf*ones(n,1);
end
% enforced linear inequality constraints: Cx <= 0
if ~isfield(options.constraints, 'ineqConstraint')
    options.constraints.ineqConstraint = [];
end
% relaxable linear inequality constraints: Cx <= 0
if ~isfield(options.constraints, 'ineqConstraint_relax')
    options.constraints.ineqConstraint_relax = [];
end
% enforced nonlinear inequality constraints: c(x) <= 0
% if multiple, order in a cell array by rows
if ~isfield(options.constraints, 'nlineqConstraint')
    options.constraints.nlineqConstraint = [];
end
% relaxable nonlinear inequality constraints: c(x) <= 0
% if multiple, order in a cell array by rows
if ~isfield(options.constraints, 'nlineqConstraint_relax')
    options.constraints.nlineqConstraint_relax = [];
end


%% ---- Define variable options -------------------------------------------
if ~isfield(options, 'variable')
    options.variable = struct();
end
% typical values of variables
if ~isfield(options.variable, 'xTypical')
    options.variable.xTypical = ones(n,1);
end
% typical value of function
if ~isfield(options.variable, 'fTypical')
    options.variable.fTypical = 1;
end
% periodic variables on [lb,ub]
if ~isfield(options.variable, 'periodic')
    options.variable.periodic = false(n,1);
end
% discrete variables
if ~isfield(options.variable, 'discrete')
    options.variable.discrete = false(n,1);
end


end 