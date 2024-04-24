%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ---- Latin Hypercube Sampling ----------------------------------- %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For use during the MADS search step.
%
% Create a sample search set using Latin Hypercube Sampling on the
% rectangle with bounds [lower, upper].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sample = latin_hypercube_sampling(n,p,lower,upper)
% INPUTS:
% n - dimension of space
% p - number of sample points
% lower - lower bound (row vector)
% upper - upper bound (row vector)
%
% OUTPUTS:
% sample -- sampled points

sample = zeros(n,p);
for i = 1:n
    sample(i,:) = lower(i) + ((randperm(p) - rand(1,p))/p)*(upper(i) - lower(i));
end

end