function R = run_taskSpacePV_interference_sym(ratNames, varargin)
% Thin wrapper: run your existing pipeline but force symmetric "across" (cross-halves).
R = run_taskSpacePV_interference(ratNames, 'AcrossMode','cross-halves', varargin{:});
end
