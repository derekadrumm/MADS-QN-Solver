%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ---- Constraint Violation Function (CVF) ------------------------ %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the cvf of the constraints given by the constraints structure.
% Constraints structure is defined in the options struct.
%
% See [1] "Audet. Derivative-Free and Blackbox
% Optimization. 2017." , Chapter 8 "Mesh Adaptive Direct Search"
function h = constraint_violation(x,constraints)
% INPUTS:
% x - point to evaluate cvf at
% constraints - the structure of linear and nonlinear constraints
%
% OUTPUTS:
% h -- the value of the cvf

lb = constraints.lb;
ub = constraints.ub;
ineq_constraint = constraints.ineqConstraint;
ineq_constraint_relax = constraints.ineqConstraint_relax;
nlineq_constraint = constraints.nlineqConstraint;
nlineq_constraint_relax = constraints.nlineqConstraint_relax;


num_ineq = size(ineq_constraint,1);
num_ineq_relax = size(ineq_constraint_relax,1);
num_nlineq = size(nlineq_constraint,2);
num_nlineq_relax = size(nlineq_constraint_relax,2);

h = 0;
% if domain constraints are violated, return Inf
if ~all((lb <= x) & (x <= ub))
    h = Inf;
    return
end


%If nonrelaxable constraints are violated, return Inf
% enforced linear inequality constraints: Cx <= 0
for i = 1:num_ineq
    if ~all(ineq_constraint{i}*x > 0)
        h = Inf;
        return
    end
end

% enforced nonlinear inequality constraints: Cx <= 0
for i = 1:num_nlineq
    if (nlineq_constraint{i}(x) > 0)
        h = Inf;
        return
    end
end


% If relaxable constraints are not violated, compute CVF
% relaxable linear inequality constraints
for i = 1:num_ineq_relax
    h = h + max(ineq_constraint_relax{i}*x,0).^2;
end

% relaxable nonlinear inequality constraints
for i = 1:num_nlineq_relax
    h = h + max(nlineq_constraint_relax{i}(x),0).^2;
end

end
