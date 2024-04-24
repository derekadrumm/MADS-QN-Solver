%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ---- Curve Matching, Problem Set-Up ----------------------------- %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the objective function, initial guess, derivative information,
% and any desired constraint functions.
%
% Then, define the options structure to be input into madsqn.
%
% Goal to determine the control point (cpt) locations and weights which define a
% periodic nurbs (pnurbs) curve that closely matches a given input curve.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---- Objective Function and Initial Guess ------------------------------

% -- Build curve to match
% try to match an elipse
r = 2;
theta = linspace(0,2*pi,200);
x = 1.5*r*cos(theta) + 2*r*sin(theta);
y = 0.8*r*sin(theta) - 0.4*r*cos(theta);
curve = [x; y]';


% -- Define control point parameters and initial nurbs curve.
% Initial condition
lb = -10; % minimimum cpt location (x,y)
ub = 10; % maximum cpt location (x,y)

N = 4; % number of cpts

% Initialize N cpts along a circle
R = 2;
Rx = 0;
Ry = 0;
theta = linspace(0,2*pi,N+1);
theta = theta(1:end-1);
cpts_x_init = R*cos(theta) + Rx;
cpts_y_init = R*sin(theta) + Ry;

cpts_init = [cpts_x_init; cpts_y_init];
weights_init = ones(N,1);

% initial parameter to input into MADS (column vector)
x0 = [cpts_x_init cpts_y_init weights_init']';
n =  length(x0); % dimension of problem (number of params to solve for)

% -- Objective function
% Calculates the sum of square distances between points along the curves.
p = 5; % degree splines should be during solve
f = @(x) curve_matching(curve,x,p); % NURBS

% Can fix search to only Bsplines if fixed weights to be 1.
% f = @(x) curve_matching(curve,[x; ones(N,1)],p); % Bspline

% -- Create and plot spline
[pnurbs_init,t_pnurbs,cpts_new] = pnurbspline(cpts_init,weights_init,p,100,"none",[0,1,250]);

% plot
figure(figcount)
figcount = figcount + 1;
plot(pnurbs_init(1,1),pnurbs_init(2,1),'ko',pnurbs_init(1,end),pnurbs_init(2,end),'ko')
hold on
plot(pnurbs_init(1,:),pnurbs_init(2,:),'k','Linewidth',1)
axis([lb ub lb ub])
plot([cpts_new(1,:) cpts_new(1,1)],[cpts_new(2,:) cpts_new(2,1)],'r*--')
text(cpts_new(1,:)+0.03, cpts_new(2,:), string(1:numel(cpts_new(1,:))))
plot(curve(:,1),curve(:,2),'b','Linewidth',1)
drawnow
hold off


%% ---- Options -----------------------------------------------------------
% Build options
options = build_options(n);

% -- Specify optimization options
options.maxIter = 300;
options.relativeMeshGradientTolerance = 1e-12;

% Specify options parameters
options.displayIterInfo = true;

% -- Specify which search steps to use, with any desired options
% options.searchStep.samplingSearch = true;
options.searchStep.quasiNewton = true;
options.derivativeMethod = 'CS_1';
options.searchStep.quasiNewtonStartTol = 1e-3;

% variable options
options.variable.xTypical = x0;
options.variable.fTypical = 1;


% -- Specify constraints in options structure
% bounds
options.constraints.lb = [lb*ones(N,1); lb*ones(N,1); 0.5*ones(N,1)];
options.constraints.ub = [ub*ones(N,1); ub*ones(N,1); 5*ones(N,1)];

% for use in closed curves, checks spline for self intersection
options.constraints.nlineqConstraint{1} = @(x) poly_selfintersection(x(1:2*N),2,N);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---- CURVE MATCHING FUNCTION -------------------------------------------
% Want to determine spline which closely matches to curve

% Objective function to minimize
function out = curve_matching(curve,cpts_in,p)
% INPUT:
% curve:   curve to match
% cpts:    control points of spline as a nx1 collumn vector [x | | | y]
% p:       degree of spline
%
% OUTPUT:
% out:      objective function value to minimize

N = length(cpts_in)/3;
cpts_x = cpts_in(1:N)';
cpts_y = cpts_in(N+1:2*N)';
weights = cpts_in(2*N+1:end);
cpts = [cpts_x; cpts_y];
[pnurbs,~] = pnurbspline(cpts,weights,p,100,"none",[0,1,size(curve,1)]);

% Minimize sum of square distances
out = sum(sum((curve-pnurbs').^2));

end

