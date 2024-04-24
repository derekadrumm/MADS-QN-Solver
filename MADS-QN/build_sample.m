%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ---- SAMPLING STRUCTURE INITIALIZATION -------------------------- %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define sampling structure.
%
% Used when want to sample points in solution space.
%
% Should be no need to mess with this, unless want to modify solver.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sampleParameters = build_sample(n)
% INPUTS:
% n -- dimension of parameter search space
%
% OUTPUTS:
% sampleParameters -- updated sampleParameters structure

%% ---- Initialize history structure --------------------------------------
sampleParameters = struct();


%% ---- Sample Parameters -------------------------------------------------
% spacial dimension
if ~isfield(sampleParameters, 'dimension')
    sampleParameters.dimension = n;
end
% number of samples
if ~isfield(sampleParameters, 'numSamples')
    sampleParameters.numSamples = 0;
end
% maximum number of samples
if ~isfield(sampleParameters, 'maxSamples')
    sampleParameters.maxSamples = 0;
end
% sampling mesh size
if ~isfield(sampleParameters, 'meshSize')
    sampleParameters.meshSize = 0;
end
% solver search directions
if ~isfield(sampleParameters, 'searchDirections')
    sampleParameters.searchDirections = [eye(n) -eye(n)];
end
% average number of steps to take in D directions
if ~isfield(sampleParameters, 'averageSteps')
    sampleParameters.averageSteps = 10;
end


% sampling lower bound
if ~isfield(sampleParameters, 'lowerBound')
    sampleParameters.lowerBound = -Inf;
end
% sampling upper bound
if ~isfield(sampleParameters, 'upperBound')
    sampleParameters.upperBound = Inf;
end
% sampling center
if ~isfield(sampleParameters, 'x_center')
    sampleParameters.x_center = 0;
end
% sampling lower radius
if ~isfield(sampleParameters, 'lowerRadius')
    sampleParameters.lowerRadius = 1;
end
% sampling upper radius
if ~isfield(sampleParameters, 'upperRadius')
    sampleParameters.upperRadius = 20;
end



end 