%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- Mesh Adaptive Direct Search (MADS) with Quasi Newton (QN) Search -- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Search for the minimizers of a nonlinear function f(x) using MADS, with a
% Sampling search and a QN search.
%
% Allow for linear and nonlinear constraints. Constraints can be relaxable
% or unrelaxable.
%
% Allows for implementation of periodic variables.
% 
%
% For more details about MADS, see 
% [1] "Audet. Derivative-Free and Blackbox
% Optimization. 2017." , Chapter 8 "Mesh Adaptive Direct Search"
%
% For more details about QN, see 
% [2] "Dennis and Schnabel. Numerical methods 
% for unconstrained optimization and nonlinear equations. 1983" Chapter 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_sol, varargout] = madsqn(f, x_0, varargin)
% INPUTS:
% f -- function to optimize
% x0 -- initial point
% options -- options structure with various properties
% history -- history structure (if running madscon multiple times)

% OUTPUTS:
% x_min -- minimum location
% f_min -- minimum value
% history -- history structure with opt information
% options -- return options structure


%% ---- Load/Build structures ---------------------------------------------
% Load/build the solver structures
% -- Build initial structure for solver.
global fCount history
n = size(x_0,1);

% -- Define or Build options structure
if size(varargin) < 1
    options = build_options(n);
else
    options = varargin{1};
end

% -- Define or Build history structure
if size(varargin) < 2
    history = build_history();
else
    history = varargin{2};
end
% retrieve previous iteration count (if running multiple simulations)
if isempty(history.iterCount)
    lastIter = 0;
else
    lastIter = history.iterCount(end);
end

% -- Define or Build sampleParameters structure
if size(varargin) < 3
    sampleParameters = build_sample(n);
else
    sampleParameters = varargin{3};
end


%% ---- Initialize parameter names ----------------------------------------
% Renaming some structure parameters.
% -- Algorithm initial variables
del = options.frameSize; % frame size
tau = options.meshSizeAdjust; % size adjustment value

% -- Tolerance and error variables
fTol = options.functionTolerance;
meshTol = options.meshTolerance;
frameTol = options.frameTolerance;
relmeshgradTol = options.relativeMeshGradientTolerance;
typx = options.variable.xTypical; % typical x
typf = options.variable.fTypical; % typical f

% -- Misc. parameters/variables
D = options.searchDirections;
p = size(D,2); % number of search directions
lower = options.constraints.lb; % lower bound of space
upper = options.constraints.ub; % upper bound of space
maxstep = options.maxStepSize; % maximum step size

% -- Sampling parameters
% Only if using sampling search step.
if options.searchStep.samplingSearch
    num_samples = options.numSurrogateSamples; % number of LHS samples
    lhs_delta = min(upper-lower)/2; % LHS search width
    % in case no ub,lb
    if isnan(lhs_delta)
        lhs_delta = del;
    end

    sampleType = options.sampleType;
    sampleParameters.numSamples = num_samples;
    sampleParameters.maxSamples = max(num_samples,100);
    sampleParameters.lowerBound = lower; % lower bound of sample space
    sampleParameters.upperBound = upper; % upper bound of sample space
    sampleParameters.lowerRadius = min(norm(lower),1); % inner radius
    sampleParameters.upperRadius = min(norm(upper),8); % outer radius
    sampleParameters.x_center = x_0; % center of sample space
    sampleParameters.searchDirections = D;
end


%% ---- Preliminary Calculations ------------------------------------------
% -- Determine period of periodic variables
periodic = options.variable.periodic;
period = zeros(n,1);
period(periodic) = upper(periodic) - lower(periodic);


% -- Check Constraints to make sure sizes match up/ are valid
constraint_check(x_0,options.constraints);


%% ---- Define Functions ------------------------------------------
% -- Define constraint, objective, and surrogate functions
mu = 1; % constraint violation paramter (fixed, can be tweaked)
h = @(x) constraint_violation(x,options.constraints);
f = @(x) f(x) + mu*h(x);
f_sur = @(x) options.f_sur(x) + mu*h(x);


% -- Define gradient and hessian functions
% (Precision parameters fixed, can be tweaked).
grad_f = @(x,f_c) numerical_gradient_infnan(x,f,options.derivativeMethod,f_c,ones(n,1),1e-10);
hess_f = @(x,f_c) numerical_hessian_infnan(x,f,f_c,ones(n,1),1e-10);


%% ---- Main Algorithm ----------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Where the magic happens :)
% Iteration Codes:
% 1 - Refine Mesh
% 2 - Poll Step Success
% 3 - QuasiNewton Success
% 4 - QuasiNewton-LineSearch Success
% 5 - (reserved)
% 6 - Sampling Success

% -- Solver variable names
x_c = x_0; % current point
x_k = x_0; % point at iteration k
x_p = x_0; % newest point
x_prev = x_0; % previous best
f_c = f(x_c); % f at current point
f_k = f_c; % f at iteration k
f_p = f_c; % f at newest point
f_prev = f_c; % previous best function evaluation
xsur_c = x_0; % current surrogate location
fsur_c = f_sur(x_0); % current surrogate value

del_p = del; % frame size
del_k = del; % frame size at iteration k

% -- Quality of life variable names
iter = 0 + lastIter; % iteration count
fCount = 1; % function evaluation count
prev_fCount = 0; % function count on previous iteration
f_surCount = 1; % surrogate function evaluation counter
prev_f_surCount = 0; % surrogate function count on previous iteration
surfailCount = 0; % counts how many smf steps are taken
x_error = Inf; % initialize error
f_error = Inf; % initialize error
FOUND_IMPROVEMENT = false; % determine if solution is updated
iterType = ''; % what step solver is on

% -- Initialize errors and logicals
grad_error = Inf;
relgrad_error = Inf;
relmeshgrad_error = Inf;
QN_SEARCH = false;
FOUND_IMPROVEMENT_PREV = false; % if solution updated on previous iteration
FOUND_IMPROVEMENT_QN = true;
QN_COMPUTATIONS = false; % determine if need to calculate grad and hess
QN_counter = 0; % count how many iterations QN is active


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while iter <= options.maxIter + lastIter
    
% -------------------------------------------------------------------------
% ---- UPDATE STEP 1 ------------------------------------------------------
% -------------------------------------------------------------------------

    % -- Update counters
    iter = iter+1;
    prev_fCount = fCount;
    prev_f_surCount = f_surCount;
    
    % -- Store previous solution, update current
    x_c = x_k;
    f_c = f_k;
    % -- Mesh Size and Parameter update
    del_c = del_p;
    delta_c = min(del_c,del_c^2);
    
    % -- Reset logicals
    FOUND_IMPROVEMENT = false;
    FOUND_IMPROVEMENT_QN = false;

    
% -------------------------------------------------------------------------
% ---- PRELIMINARY STEPS --------------------------------------------------
% -------------------------------------------------------------------------
    % This step checks if using QN or Sampling searches.
    % If they are used, update and calculate relevant values.

    % -- Update sampling parameters
    if options.searchStep.samplingSearch
        sampleParameters.x_center = x_c;
        sampleParameters.meshSize = delta_c;
        sampleParameters.lowerBound = xsur_c - lhs_delta;
        sampleParameters.upperBound = xsur_c + lhs_delta;
    end
    
    % -- Determine if using QN search
    % turn on if mesh is fine enough, and x_c is feasible
    % once on, only turn off if mesh gets too course
    if options.searchStep.quasiNewton && (delta_c < options.searchStep.quasiNewtonStartTol) && (~isinf(f_c))
        QN_SEARCH = true;
        QN_counter = QN_counter + 1;
        
    % turn off, reset counter
    elseif (delta_c >= 1)
        QN_SEARCH = false;
        QN_counter = 0;
    end
    
    % -- Determine if need to calculate grad/hess
    % If first time QN is active, OR
    % If found an improvement, need to compute grad/hess again
    if (QN_counter == 1) || (QN_SEARCH && FOUND_IMPROVEMENT_PREV)
        QN_COMPUTATIONS = true;
        QN_counter = QN_counter + 1;
    % If no improvement (refined mesh) can keep old grad/hess
    else     
        QN_COMPUTATIONS = false;
        % Check previous grad/hess for NaN and Inf
        % Need to do this here to make sure QN_SEARCH wont run using bad vals
        if QN_counter > 1
            gradcheck = (any(isnan(grad_f_c),'all') || any(isinf(grad_f_c),'all'));
            hesscheck = any(isnan(H_c),'all') || any(isinf(H_c),'all');
            if gradcheck || hesscheck
                QN_SEARCH = false;
            end
        end
    end
    
    % -- If using QN search, compute the grad/hess
    % Also check that grad/hess are feasible, if not, skip QN step
    if options.searchStep.quasiNewton && QN_SEARCH && QN_COMPUTATIONS
        grad_f_c = grad_f(x_c,f_c); 
        % check that grad is finite
        if any(isnan(grad_f_c),'all') || any(isinf(grad_f_c),'all')
            QN_SEARCH = false;
        else       
            H_c = hess_f(x_c,f_c);        
            % check that hess is finite
            if any(isnan(H_c),'all') || any(isinf(H_c),'all')
                QN_SEARCH = false;
            else
                % check hessian is spd (add small amount to diagonal)
                mu_k = 0;
                if ~all(eig(H_c) > 0)
                    mu_k = abs(min(eig(H_c))) + 1e-4;
                end
                H_c = H_c + mu_k*eye(n);
            end
        end
    end
    
    
% -------------------------------------------------------------------------
% ---- SEARCH STEPS -------------------------------------------------------
% -------------------------------------------------------------------------

    % ---- Search Step (Quasi Newton) -------------------------------------
    % ---------------------------------------------------------------------
    % Find the QN iterate of x_c,
    % then project x_c onto the mesh,
    % then search set is the 1st nearest neighbors to projected point.
    if options.searchStep.quasiNewton && QN_SEARCH

        % -- Only do QN until relgrad is within tolerance
        if relgrad_error > 1e-6
            % -- Calculate Newton point
            s_N = H_c\-grad_f_c; % calculate Newton step
            max_step_size = 1;
            if norm(s_N) > max_step_size
                s_N = max_step_size*s_N/norm(s_N);
            end

            x_N = x_c + s_N; % travel along Newton direction
            
            % check feasibilty of Newton point
            feasible = ~isinf(h(x_N)); % if h(t) = Inf, then t breaks enforced conditions
            if feasible
                f_N = f(x_N);
                fCount = fCount + 1;
            else
                f_N = Inf;
            end

            % -- Check alpha condition to determine if x_N is a good step
            alpha = options.searchStep.quasiNewtonAlpha;
            ALPHA_CONDITION = (f_N <= f_c + alpha*grad_f_c'*s_N);

            % -- If alpha_condition satisfied, do Newton Step
            if ALPHA_CONDITION
                % project Newton point onto mesh
                x_proj = mesh_projection(x_c,x_N,delta_c,D,1);
                for i = 1:size(x_proj,2)
                    f_proj = f(x_proj(:,i));
                    fCount = fCount + 1;

                    % check feasibilty of trial point
                    feasible = ~isinf(h(x_proj(:,i))); % if h(t) = Inf, then t breaks enforced conditions
                    if feasible
                        f_test = f_proj;
                    else
                        f_test = Inf;
                    end

                    % check if trial point is an improvement
                    if f_test < f_c
                        x_k = x_proj(:,i);
                        f_k = f_test;
                        FOUND_IMPROVEMENT = true;
                        FOUND_IMPROVEMENT_QN = true;
                        % display number of mesh points searched through
                        % S-0 --> took proj(x_N) as new point
                        iter_type = ['Search Step (QN:S-', num2str(i-1),':LS-0)'];
                        history.iterTypeCode(iter) = 3;
                        options.searchStep.quasiNewtonCount = options.searchStep.quasiNewtonCount + 1;
                        break
                    end
                end
                
            % -- If alpha_condition not satisfied, do Backtracking Line Search (LS):
            % only if computed grad/hess on this iteration
            else
                % parameters (fixed values, see [2])
                maxcount = 10; % max number of backtracking steps (can be tweaked)
                alpha = 1e-4; % alpha in (0,1/2) (should be < 1/4)
                u_lam = 1/2; % lam < u_lam
                l_lam = 1/10; % l_lam < lam

                descent_direction = s_N;

                init_slope = grad_f_c'*descent_direction; % search in descent direction
                lam = 1;
                f_p = f_N;

                LS_count = 0;
                LS_SUCCESS = false;
                % -- Begin LS algorithm
                while ~LS_SUCCESS && (LS_count < maxcount)
                    LS_count = LS_count + 1;
                    % on first backtrack, use a quadratic model
                    if lam == 1
                        lam_temp = -init_slope/(2*(f_p - f_c - init_slope));

                    % otherwise, use a cubic model    
                    else
                        A = [1/(lam^2) -1/(lam_prev^2); -lam_prev/(lam^2) lam/(lam_prev^2)];
                        b = [f_p - f_c - lam*init_slope; f_p_prev - f_c - lam_prev*init_slope];
                        ab = (1/(lam - lam_prev))*A*b;
                        a = ab(1);
                        b = ab(2);
                        disc = b^2 - 3*a*init_slope;

                        % is actually a quadratic
                        if a == 0
                           lam_temp = -init_slope/(2*b);

                        % is actually a cubic    
                        else
                            lam_temp = (-b + sqrt(disc))/(3*a);
                        end
                    end

                    % update lam
                    lam_prev = lam;
                    f_p_prev = f_p;

                    % check that lam in [lower*lam,upper*lam]
                    % lam too small
                    if lam_temp < l_lam*lam
                        lam = l_lam*lam;
                    % lam too big or Nan/Inf
                    elseif (lam_temp > u_lam*lam) || (isnan(lam_temp))
                        lam = u_lam*lam;
                    % lam just right
                    else
                        lam = lam_temp;
                    end

                    % check alpha condition and TR condition
                    s_c = lam*descent_direction; % current step
                    x_p = x_c + s_c;
                    f_p = f(x_p);
                    fCount = fCount + 1;
                    ALPHA_CONDITION = (f_p <= f_c + alpha*lam*init_slope);
                    STEPSIZE_CONDITION = (norm(s_c) <= maxstep);

                    LS_SUCCESS = (ALPHA_CONDITION && STEPSIZE_CONDITION);
                end

                % -- If LS succeeded, update
                if LS_SUCCESS
                    % project newton point onto mesh
                    x_proj = mesh_projection(x_c,x_p,delta_c,D,1);
                    for i = 1:size(x_proj,2)

                        % check feasibilty of trial point
                        feasible = ~isinf(h(x_proj(:,i))); % if h(t) = Inf, then t breaks enforced conditions
                        if feasible
                            f_proj = f(x_proj(:,i));
                            fCount = fCount + 1;
                            f_test = f_proj;
                        else
                            f_test = Inf;
                        end

                        % check if trial point is an improvement
                        if f_test < f_c
                            x_k = x_proj(:,i);
                            f_k = f_proj;
                            FOUND_IMPROVEMENT = true;
                            FOUND_IMPROVEMENT_QN = true;
                            % display number of mesh points searched through
                            % S-0 --> took proj(x_N) as new point
                            iter_type = ['Search Step (QN:S-', num2str(i-1),':LS-', num2str(LS_count), ')']; 
                            history.iterTypeCode(iter) = 4;
                            options.searchStep.quasiNewtonCount = options.searchStep.quasiNewtonCount + 1;
                            break
                        end
                    end
                end
            end
            
        end
    
    end
        
    
    
    % ---- Search Step (using sampling method) ----------------------------
    % ---------------------------------------------------------------------
    % Use desired sampling method to heuristically search.
    % (can modify to order sampled points based on the surrogate function)
    if options.searchStep.samplingSearch && ~FOUND_IMPROVEMENT && ~QN_SEARCH
        sample = sampling_method(sampleType,sampleParameters,h);
        p_out = size(sample,2);

        % -- Check each point in sample set
        for i = 1:p_out
            t = mod(sample(:,i),period);

            % check feasibilty of trial point
            feasible = ~isinf(h(t)); % feasible determined previously
            if feasible
                f_test = f(t);
                fCount = fCount + 1;
            else
                f_test = Inf;
            end

            % check if trial point is an improvement
            if f_test < f_c
                x_k = t;
                f_k = f_test;
                FOUND_IMPROVEMENT = true;
                iter_type = 'Search Step (Sampling)';
                history.iterTypeCode(iter) = 6;
                options.searchStep.samplingSearchCount = options.searchStep.samplingSearchCount + 1;
                break
            end
        end
    end
    
    
    
% -------------------------------------------------------------------------
% ---- POLL STEP ----------------------------------------------------------
% -------------------------------------------------------------------------
    % If search steps unsuccessful, search using mesh via polling.
    % Evaluate surrogate function first, reorder points,
    % then evaluate true function in that order.
    if ~FOUND_IMPROVEMENT
        % -- Generate Polling Directions
        DD = poll_direction(n,del_c,delta_c);

        % -- Check poll points on surrogate function
        % if no surrogate, f_sur(x) := 0, so this has no effect
        test_pts = [];
        fsur_test_pts = [];
        for i = 1:p
            t = mod(x_c + delta_c*DD(:,i),period);

            % check feasibilty of trial point
            feasible = ~isinf(h(t));
            if feasible
                fsur_t = f_sur(t);
                f_surCount = f_surCount + 1;
                test_pts = [test_pts t];
                fsur_test_pts = [fsur_test_pts fsur_t];
            else
                fsur_t = Inf;
                test_pts = [test_pts t];
                fsur_test_pts = [fsur_test_pts fsur_t];
            end
        end

        % -- Sort surrogate function points
        [fsur_test_pts,I] = sort(fsur_test_pts);
        test_pts = test_pts(:,I);

        % -- Check each point in poll set
        for i = 1:p
            t = mod(test_pts(:,i),period);

            % check feasibilty of trial point
            feasible = ~isinf(fsur_test_pts(i)); % feasible determined previously
            if feasible
                f_test = f(t);
                fCount = fCount + 1;
            else
                f_test = Inf;
            end

            % check if trial point is an improvement
            if f_test < f_c
                x_k = t;
                f_k = f_test;
                FOUND_IMPROVEMENT = true;
                iter_type = 'Poll Step';
                history.iterTypeCode(iter) = 2;
                break

            end
        end
    end
    

% -------------------------------------------------------------------------
% ---- UPDATE STEP 2 ------------------------------------------------------
% -------------------------------------------------------------------------

    % ---- Frame Size Update ----------------------------------------------
    % ---------------------------------------------------------------------
    % -- Store mesh size and counts
    del_k = del_c;
    delta_k = delta_c;
    fCount_k = fCount - prev_fCount;
    f_surCount_k = f_surCount - prev_f_surCount;
    
    % -- Update frame size and update errors
    % if improvement found, coarsen frame
    % otherwise, refine frame.
    if FOUND_IMPROVEMENT
        FOUND_IMPROVEMENT_PREV = true;
        x_prev = x_c;
        f_prev = f_c;
        del_max = 8; % maximum possible frame size
        
        del_p = min(del_max,del_k/tau);
        
        x_error = norm(x_k - x_prev);
        f_error = norm(f_k - f_prev);
        
    else      
        FOUND_IMPROVEMENT_PREV = false;
        del_p = del_k*tau;
        iter_type = 'Refine Mesh';
        history.iterTypeCode(iter) = 1;
    end
     
    
    % ---- Error Calculation ----------------------------------------------
    % ---------------------------------------------------------------------
    % -- QN Errors
    delta_pow = log(delta_k)/log(1/tau);
    if options.searchStep.quasiNewton && QN_SEARCH
        % gradient error
        grad_f_k = grad_f(x_k,f_k);
        grad_error = norm(grad_f_k);
        
        % relative gradient error
        relgrad = (grad_f_k.*max(abs(x_k),typx))/max(abs(f_k),typf);
        relgrad_error = norm(relgrad,Inf);

        % relative mesh-gradient error
        relmeshgrad_error = relgrad_error*(10^delta_pow);
        
    end
    
    % relative mesh size
    relmesh_error = delta_k*norm(typx);

    
    % -- History Saving
    history.iterCount(iter) = iter;
    history.xVal(:,iter) = x_k;
    history.fVal(iter) = f_k;
    history.fCount(iter) = fCount_k;
    history.f_surCount(iter) = f_surCount_k;
    history.iterType{iter} = iter_type;
    history.meshSize(iter) = delta_k;
    history.frameSize(iter) = del_k;
    
    history.error.xError(iter) = x_error;
    history.error.fError(iter) = f_error;
    history.error.fGradError(iter) = grad_error;
    history.error.relGradError(iter) = relgrad_error;
    history.error.relMeshGradError(iter) = relmeshgrad_error;
    history.error.relMeshError(iter) = relmesh_error;
    
    % -- Display Progress
    if options.displayIterInfo
        message = [num2str(iter) ' '...
                   iter_type ' '...
                   num2str(f_k) ' '...
                   num2str(delta_pow) ' '...
                   num2str(fCount_k)];
        disp(message)
    end
    
% -------------------------------------------------------------------------
% ---- TERMINATION STEP ---------------------------------------------------
% -------------------------------------------------------------------------
    % terminate if within mesh tolerance
    if delta_k < meshTol
        if options.displayIterInfo
            disp('Mesh size within tolerance')
        end
        break

    % terminate if within mesh tolerance (power)
    elseif delta_pow < -Inf
        if options.displayIterInfo
            disp('Mesh size power within tolerance')
        end
        break

    % terminate if within frame tolerance
    elseif del_k < frameTol
        if options.displayIterInfo
            disp('Frame size within tolerance')
        end
        break

    % terminate if within function tolerance    
    elseif f_error < -Inf
        if options.displayIterInfo
            disp('Function value within tolerance')
        end        
        break

    % terminate if within relgrad-mesh tolerance   
    elseif relmeshgrad_error < relmeshgradTol
        if options.displayIterInfo
            disp('Relative mesh-gradient within tolerance')
        end
        break

    % terminate if maximum iterations
    elseif iter >= options.maxIter
        if options.displayIterInfo
            disp('Maximum number of iterations reached')
        end
        break
    end  
    
    
end


%% ---- Outputs -----------------------------------------------------------
x_sol = x_k; % approximate solution
varargout{1} = f_k; % function value at x_sol
varargout{2} = history; % history info
varargout{3} = options; % options info
disp(' ')
end



%% ------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SUBROUTINES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---- Check Constraints -------------------------------------------------
% Perform sanity check on constraints to make sure sizes, etc, match up
% with initial conditions
function constraint_check(x,constraints)

lb = constraints.lb;
ub = constraints.ub;
ineq_constraint = constraints.ineqConstraint;
ineq_constraint_relax = constraints.ineqConstraint_relax;

num_ineq = size(ineq_constraint,1);
num_ineq_relax = size(ineq_constraint_relax,1);


% Check Bounds
if ~all(size(lb) == size(x))
    msg = 'Size Mismatch: Lowerbound and x0.';
    error(msg)
end
if ~all(size(ub) == size(x))
    msg = 'Size Mismatch: Upperbound and x0.';
    error(msg)
end

% Check linear constraints
for i = 1:num_ineq
    if ~all(size(ineq_constraint{i},2) == size(x,1))
        msg = ['Size Mismatch: Cols of Inequality Constraint ',...
            num2str(i), ' must equal rows of x0.'];
        error(msg)
    end
end
for i = 1:num_ineq_relax
    if ~all(size(ineq_constraint_relax{i},2) == size(x,1))
        msg = ['Size Mismatch: Cols of Relaxale Inequality Constraint ',...
            num2str(i), ' must equal rows of x0.'];
        error(msg)
    end
end
end



%% ---- Sampling Method ---------------------------------------------------
% choose which sampling method to use based on 'options'
% and output the sample points
%
% INPUT:
% sampleType - which sampling method to use
% sampleParameters - the parameters necessary for the sampling mehtod
% cvf - constraint violation function for checking if sample points lie in
%       domain space
function sample = sampling_method(sampleType,sampleParameters,cvf)

switch sampleType
    % Latin Hypercube Sampling
    case 'LHS'
        sample = latin_hypercube_sampling(...
            sampleParameters.dimension,...
            sampleParameters.numSamples,...
            sampleParameters.lowerBound,...
            sampleParameters.upperBound);
        
    % Spherical Exterior Sampling    
    case 'SES'
    [sample,~] = spherical_exterior_sampling(...
        sampleParameters.x_center,...
        sampleParameters.numSamples,...
        sampleParameters.maxSamples,...
        sampleParameters.lowerRadius,...
        sampleParameters.upperRadius,...
        sampleParameters.meshSize,...
        sampleParameters.searchDirections,...
        cvf);
end
end



%% ---- Mesh Projection ---------------------------------------------------
% Project a given point onto the mesh centered at x0
% Mesh points are of the form x0 + delta*D*y
% where y in N^p (y is a vector of naturals of size p=size(D,2))
% 
% NOTE: Only works if D is plus-minus cartesian directions
%       i.e. D = [eye(n) -eye(n)]
function out = mesh_projection(x0,x,delta,D,step_num)
% INPUTS
% x0 -- point mesh is centered at
% x -- point to project onto mesh
% delta -- mesh size parameter
% D -- direction matrix
% step_num -- number of steps to take away from x0
% 
% OUTPUTS:
% out -- point proj(x) on mesh

p = size(D,2);
n = p/2;
% determine which mesh elements xN is between
y = (1/delta)*(x - x0);
ytemp = zeros(p,1);
for i = 1:p/2
    if sign(y(i)) == 1
        ytemp(i) = y(i);
        ytemp(i+n) = 0;
    else
        ytemp(i) = 0;
        ytemp(i+n) = -y(i);
    end
end
yround = round(ytemp);

% Generate vectors step_num away from x0
% Generate all vectors with combinations of 0,...,step_num
% which sum to step_num.
yshift = [];
for i = 1:step_num
    K = p+i; % The required sum
    c = nchoosek(2:K,p-1);
    m = size(c,1);
    A = zeros(p,m);
    for ix = 1:m
        A(:,ix) = diff([1,c(ix,:),K+1]) - 1;   
    end
    yshift = [yshift yround + A];
end

out = x0 + delta*D*([yround yshift]);
end



%% ---- Poll Direction Generation -----------------------------------------
% An intermediate function to generate polling directions in MADS
%
% Requires construction of a Householder Matrix from a normal row vector v
% (see [1] for details)
function D = poll_direction(n,del,delta)
v = randi([-10000 10000],n,1);
v = v/norm(v);

% Householder matrix
H = eye(length(v)) - 2*(v*v');

% generate Poll directions
B = zeros(size(H));
for i = 1:length(v)
    B(:,i) = (del/delta)*(H(:,i)/max(abs(H(:,i))));
end

D = [B -B];
end


%% ---- Numerical Gradient (PRECISE) (CHECK INFS/NANS) --------------------
%(assuming f cannot take vector inputs)
% Want to generically compute numerical gradient of f:R^n -> R
% in n variables, at the point x, using various finite difference methods.
%
% Actively checking for Inf/NaN values in grad
% If Inf/NaN values appear, set val to be 100f(x).
function out = numerical_gradient_infnan(x,f,dtype,f_c,varargin)
% INPUTS
% x -- point to calculate at (only one point at a time!) (row vector)
% f -- function to evaluate on (needs to be able to take vector inputs!)
% dtype -- what type of differentiation to perform
% f_c -- value of f(x) (if already calculated)
% typx -- typical value of x
% eta -- parameter of accuracy (10^-DIGITS)
% 
% OUTPUTS:
% out -- gradient of f at x
global fCount

n = size(x,1);
% check f_c exists
if isempty(f_c)
    f_c = f(x);
end

% check varargin parameters
if size(varargin) < 1
    typx = zeros(n,1);
    eta = 1e-6;
elseif size(varargin) < 2
    typx = varargin{1};
    eta = 1e-6;
else
    typx = varargin{1};
    eta = varargin{2};
end

cubeta = eta^(1/3);
I = eye(n);

% Compute step sizes
stepsize = cubeta*max(abs(x),typx).*sign(x);
stepsize(abs(stepsize) < eta) = cubeta; % in case x = 0

% Compute i-th row of gradient
grad_f = zeros(n,1);
% grad_inf = Inf*ones(n,1); % Inf vector for outputting if any Inf/NaN vals

switch dtype
    % Forward Space, 1st order
    case 'FS_1'        
        for j = 1:n
            tempj = x(j) - stepsize(j);
            stepsize(j) = tempj - x(j);
            
            grad_f(j) = (f(x + I(:,j)*stepsize(j)) - f_c)./stepsize(j);
            
            fCount = fCount + 1;
            
            % check for Inf/NaN
            check = any(isnan(grad_f(j)) | isinf(grad_f(j)));
            if check
                grad_f(j) = 100*f_c;
            end
        end
        
        out = grad_f;
        
    % Central Space, 1st order    
    case 'CS_1'
        for j = 1:n
            tempj = x(j) - stepsize(j);
            stepsize(j) = tempj - x(j);
            
            fp = f(x + I(:,j)*stepsize(j));
            fm = f(x - I(:,j)*stepsize(j));
            
            grad_f(j) = (fp - fm)./(2*stepsize(j));
            
            fCount = fCount + 2;
            
            % check for Inf/NaN
            check = any(isnan(grad_f(j)) | isinf(grad_f(j)));
            if check
                grad_f(j) = 100*f_c;
            end
        end
        out = grad_f;
        
end
end



%% ---- Numerical Hessian Matrix (PRECISE) (CHECK INFS/NANS) --------------
% Calculate Hessian as in [2] Algorithm 5.6.2 (pg 321)
% Modified to more closely follow the Finite Difference formula.
%
% Actively checking for Inf/NaN values in H.
% If Inf/NaN values appear, stop calculation and output Inf/NaN matrix.
function H = numerical_hessian_infnan(x,f,f_c,varargin)
% INPUTS
% x -- point to calculate at (only one point at a time!) (row vector)
% f -- function to evaluate on (needs to be able to take vector inputs!)
% f_c -- function value at x
% typx -- typical value of x
% eta -- parameter of accuracy (10^-DIGITS)
% 
% OUTPUTS:
% H -- hessian matrix of f at x
global fCount

n = size(x,1);
% check f_c exists
if isempty(f_c)
    f_c = f(x);
    fCount = fCount + 1;
end

% check varargin parameters
if size(varargin) < 1
    typx = ones(n,1);
    eta = 1e-6;
elseif size(varargin) < 2
    typx = varargin{1};
    eta = 1e-6;
else
    typx = varargin{1};
    eta = varargin{2};
end


cubeta = eta^(1/3);
I = eye(n);
% Hinf = Inf*ones(n); % inf matrix to return if Inf/NaN appear

% Compute step sizes
stepsize = cubeta*max(abs(x),typx).*sign(x);
stepsize(abs(stepsize) < eta) = cubeta; % in case x = 0

% Compute f_i vals, store in n-vector
F_N = zeros(n,1);
for i = 1:n
    F_N(i) = f(x + stepsize(i)*I(:,i));
    fCount = fCount + 1;
    
    % check for Inf/NaN
    check = any(isnan(F_N(i)) | isinf(F_N(i)));
    if check
        F_N(i) = 100*f_c;
    end
end

% Compute i-th row of H
H = zeros(n);
for i = 1:n
    f_i = F_N(i);
    f_ii = f(x + 2*stepsize(i)*I(:,i));
    
    fCount = fCount + 1;
    
    % check for Inf/NaN
    check = any(isnan(f_ii) | isinf(f_ii));
    if check
        f_ii = 100*f_c;
    end
    
    H(i,i) = (f_ii - 2*f_i + f_c)/(stepsize(i)*stepsize(i));
    
    % Compute ij-th entry of H
    for j = i+1:n
        f_ij = f(x + stepsize(i)*I(:,i) + stepsize(j)*I(:,j));
        f_j = F_N(j);
        
        fCount = fCount + 1;
        
        % check for Inf/NaN
        check = any(isnan(f_ij) | isinf(f_ij));
        if check
            f_ij = 100*f_c;
        end
        
        H(i,j) = ((f_ij - f_i) - (f_j - f_c))/(stepsize(i)*stepsize(j));
    end
end

% Build lower triangle of H (symmetric)
D = diag(diag(H));
UT = H - D;
LT = UT';
H = H + LT;
end

