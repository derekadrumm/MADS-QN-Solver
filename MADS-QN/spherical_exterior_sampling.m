%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ---- Spherical Exterior Sampling -------------------------------- %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For use during the MADS search step.
%
% Create a sample search set using SES Sampling, to sample points on the
% mesh centered at x, outside a sphere of radius rho.
%
% NOTE: may not give exactly 'p' sample points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sample_out,p_out] = spherical_exterior_sampling(x,p,pmax,rho_l,rho_u,delta,D,h)
% INPUTS:
% x - center point
% p - preferred number of sample points
% pmax - maximum number of possible sample points
% rho_l - minimum distance of sample points from x
% rho_u - distance from x to boundary of space (max distance)
% delta - mesh coarseness (must be rational < 1)
% D - direction matrix
% h - constraint violation function
%
% OUTPUTS:
% sample - sampled points outside sphere

% only run if delta <= 1
if delta <= 1
    n = length(x); % dimension of sample space 
    
    Nu = floor(rho_u/(delta*min(vecnorm(D))));
    
    
    y = randi([0 Nu],size(D,2),pmax); % random perm of naturals in [0, Nu]
    
    % check length of y (to ensure outside sphere)
    check = ((vecnorm(delta*D*y) > rho_l));
    y = y(:,check);
    
    % sample points on the mesh
    sample_temp = unique((x + delta*D*y)','rows')';
    
    % choose only points satisfying constraints
    check = zeros(1,size(sample_temp,2));
    for i = 1:size(sample_temp,2)
       
        check(i) = (h(sample_temp(:,i)) < Inf);
        
    end
    check = logical(check);
    
    index = 1:size(sample_temp,2);
    sample = sample_temp(:,index(check));
    

    % only output p points from sample (or entire sample if < p)
    p_out = min(p,size(sample,2));
    if p_out <= 0
        sample_out = [];
        p_out = 0;
    else
        sample_out = sample(:,randi([1 size(sample,2)],1,p_out));
    end
    
    
% if delta > 1, return empty sample
else
    sample_out = [];
    p_out = 0;
end
end