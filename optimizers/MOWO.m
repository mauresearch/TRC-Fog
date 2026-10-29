function [ParetoSet, ParetoFront] = MOWO(N, T, lb, ub, dim, fobj, numObjectives, opts)
% Multi-objective walrus optimizer (wrapper of the published code).
if nargin < 8 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'integer'), opts.integer = false; end
addpath(fullfile(fileparts(mfilename('fullpath')), 'MOWO_published'));
lb = lb(:)'; ub = ub(:)';
if numel(lb) == 1, lb = repmat(lb, 1, dim); end
if numel(ub) == 1, ub = repmat(ub, 1, dim); end
if opts.integer
    fun = @(x) safe_eval(fobj, round(x), numObjectives);
else
    fun = @(x) safe_eval(fobj, x, numObjectives);
end
MultiObj = struct('name', 'problem', 'numOfObj', numObjectives, 'fun', fun, ...
                  'nVar', dim, 'var_min', lb, 'var_max', ub);
params = struct('Np', N, 'Nr', N, 'maxgen', T, 'showPlot', false, 'verbose', false);
[ParetoSet, ParetoFront] = MOWO_published_core(params, MultiObj);
if opts.integer, ParetoSet = round(ParetoSet); end
end

function f = safe_eval(fobj, x, M)
try
    f = real(fobj(x)); f = f(:)';
    if numel(f) ~= M || any(~isfinite(f)), f = 1e6 * ones(1, M); end
catch
    f = 1e6 * ones(1, M);
end
end
