%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ---- HISTORY STRUCTURE INITIALIZATION/UPDATE -------------------- %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define history structure.
% Contains all of the solution procedure information we may want.
%
% Should be no need to mess with this, unless want to add solver outputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function history = build_history()
% INPUTS:
% n -- dimension of parameter search space
%
% OUTPUTS:
% history -- updated history structure

%% ---- Initialize history structure --------------------------------------
history = struct();


%% ---- Solver iteration info ---------------------------------------------
% iteration count
if ~isfield(history, 'iterCount')
    history.iterCount = [];
end
% iteration type
if ~isfield(history, 'iterType')
    history.iterType = {};
end
% iteration type code
if ~isfield(history, 'iterTypeCode')
    history.iterTypeCode = [];
end


%% ---- Mesh/frame info ---------------------------------------------------
% mesh size
if ~isfield(history, 'meshSize')
    history.meshSize = [];
end
% frame size
if ~isfield(history, 'frameSize')
    history.frameSize = [];
end


%% ---- Function info -----------------------------------------------------
% function evaluation count at each iteration
if ~isfield(history, 'fCount')
    history.fCount = [];
end
% surrogate function evaluation count at each iteration
if ~isfield(history, 'f_surCount')
    history.f_surCount = [];
end
% solution position history
if ~isfield(history, 'xVal')
    history.xVal = [];
end
% solution function history
if ~isfield(history, 'fVal')
    history.fVal = [];
end
% solution surrogate function history
if ~isfield(history, 'f_surVal')
    history.f_surVal = [];
end


%% ---- Error measure info ------------------------------------------------
if ~isfield(history, 'error')
    history.error = struct();
end
% solution error
if ~isfield(history.error, 'xError')
    history.error.xError = [];
end
% function error
if ~isfield(history.error, 'fError')
    history.error.fError = [];
end
% function gradient error
if ~isfield(history.error, 'fGradError')
    history.error.fGradError = [];
end
% relative gradient error
if ~isfield(history.error, 'relGradError')
    history.error.relGradError = [];
end
% relative meshing error
if ~isfield(history.error, 'relMeshError')
    history.error.relMeshError = [];
end


end 