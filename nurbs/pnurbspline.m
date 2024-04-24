%% Periodic NonUniform Rational B-Spline (PNURBS)
% Create a periodic cubic nonuniform rational B-Spline (PNURBS) from control points
% INPUT:
% cpts:     control points as row vector (2xn) [x;y]
% weights:  weight of each cpt (nx1)
% d:        degree of spline (not sure if higher degrees will work)
% N:        desired number of points along spline
%
% sorting = {alg,Cx,Cy}:  indicate if should sort control points
% alg:      "poly", "geo", "center", "none" or []
% [Cx,Cy]:  center to sort around [x;y]
%
% [a,b,M]:  indicate if should interpolate points so curve is equispaced
%           along parameter t in [a,b], M total points
%
% OUTPUT:
% PBS:      periodic bspline as row vector [PBSx PBSy]
% cpts:     control points (order may change from input)
% u:        parameterization vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PNURBS,varargout] = pnurbspline(cpts_in,weights_in,p,N,varargin)

%% ---- Initialization ----------------------------------------------------
% sort the control points in the specified manner
sorting = varargin{1};
if sorting{1} == "none" || isempty(sorting)
    cpts = cpts_in;
else
    if sorting{1} == "poly"
        % center of polygonal center of cpts
        polyintersect = polyselfintersect(cpts_in(1,:),cpts_in(2,:));
        center = mean(polyintersect,2);     
    elseif sorting{1} == "geo"
        % center on geometric center of cpts
        center = mean(cpts_in,2);
    elseif sorting{1} == "center"
        % Sort control points around specified center
        center = [sorting{2}; sorting{3}];
    else
        msg = ['Did not specify sorting algorithm.'...
            '\nMust specify "poly", "geo", "center", or "none"'];
        error(msg);
    end
    vec = cpts_in - center ; % vectors connecting the centroid and cpts
    theta = atan2(vec(2,:),vec(1,:)); % angles of cpts
    [~, idx] = sort(theta);   % sorting the angles 
    cpts = cpts_in(:,idx); % sorting the given points
end


%% ---- Cpt Wrapping ------------------------------------------------------
n = size(cpts,2)-1; % n+1 cpts

% Wrapping p cpts to make periodic (use a loop in case p > n+1)
cpts_new = [cpts zeros(2,p)];
weights_new = [weights_in; zeros(p,1)];
for i = 1:p
    cpts_new(:,n+1+i) = cpts_new(:,i);
    weights_new(n+1+i) = weights_new(i);
    
%     cpts_new = [cpts_new cpts_new(:,i)];
%     weights_new = [weights_new weights_new(:,i)];
end
nwrap = n + p + 1; % number of total cpts after wrapping


%% ---- Define Knots ------------------------------------------------------
% define the knots in a way so that curve is parameterized on [0,1]
m = nwrap + p;
% knots = -p/(p+1):1/(p+1):1+(p/p+1);
knots = -p/(m-2*p):1/(m-2*p):1+(p/(m-2*p));
u = linspace(0,1,N);


%% ---- Rational Basis Functions ------------------------------------------
% calculate rational basis function R_i^p
% first find denominator
denom = 0;
for j = 0:nwrap-1
    jp1 = j+1;
    denom = denom + weights_new(jp1)*BSpline_Basis(u,knots,j,p); 
end

% now find PNURBS using RBF
PNURBS = zeros(2,length(u));
for i = 0:nwrap-1
    
    ip1 = i+1;
    RBF = (weights_new(ip1)*BSpline_Basis(u,knots,i,p))./denom;
    PNURBS = PNURBS + cpts_new(:,ip1)*RBF;
    
end


%% ---- Parameter Interpolation -------------------------------------------
if length(varargin) >= 2
    interval = varargin{2};
    [PNURBS,u] = equispaced(PNURBS,interval(1),interval(2),interval(3));
end


% output parameterization u
varargout{1} = u;

% output control points
varargout{2} = cpts;


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ---- Rational Basis Functions ------------------------------------------
% generate the i,n-th rational basis function
% INPUTS:
% x : domain space of knot parameter (min(knots),max(knots))
% knots: a knot vector of length m
% i: index of B-Spline Basis function (up to number of control points)
%    0 <= i <= numel(knots) - n
% n: degree of B-Spline Basis function
function RBF = Rational_Basis(u,weights,knots,i,n)

denom = 0;
for j = 0:n
    denom = denom + weights(j)*BSpline_Basis(u,knots,j,p); 
end

RBF = weights(i)*denom;



end

%% ---- Recursive BSpline function ----------------------------------------
% generate the i,n-th BSpline basis function
%
% NOTE: knot vector must be of size i+n+1 !!!
function y = BSpline_Basis(u,knots,i,n)
% INPUTS:
% x : domain space of knot parameter (min(knots),max(knots))
% knots: a knot vector of length m
% i: index of B-Spline Basis function (up to number of control points)
%    0 <= i <= numel(knots) - n
% n: degree of B-Spline Basis function
if (n==0)
    ip1 = i+1;
    check = (knots(ip1)<=u) & (u<knots(ip1+1));
    y(check) = 1;
    y(~check) = 0;
    
else
    omega_i = BSpline_Coef(u,knots,i,n);
    omega_ip1 = BSpline_Coef(u,knots,i+1,n);
    
    y = omega_i.*BSpline_Basis(u,knots,i,n-1) + (1-omega_ip1).*BSpline_Basis(u,knots,i+1,n-1);
    
end
end


% calculate the i,r-th BSpline Coefficients
function omega = BSpline_Coef(u,knots,i,r)
ip1 = i+1;
if knots(ip1) == knots(ip1+r)
    omega = 0;
else
    omega = (u - knots(ip1))./(knots(ip1+r)-knots(ip1));
end
end


%% ---- Equispaced Interpolation ------------------------------------------
% Interpolate curve along t in [a,b] with M equally spaced points
function [C,t] = equispaced(C_in,a,b,M)
% find the parameter s
dxdy_C = [diff(C_in') ; C_in(1,end)-C_in(1,1) C_in(2,end)-C_in(2,1)]';

ds_C = sqrt(dxdy_C(1,:).^2 + dxdy_C(2,:).^2);

s = zeros(1,length(ds_C));
for i = 1:length(ds_C)-1
    s(i+1) = s(i) + ds_C(i);
end
s = (s./s(end))*(b - a) + a;


% create new parameter t and interpolate points from s to t
t = linspace(a,b,M);

C_x = interp1(s,C_in(1,:),t);
C_y = interp1(s,C_in(2,:),t);

C = [C_x ; C_y];

end