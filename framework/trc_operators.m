function ops = trc_operators(ctx, model)
% TRC-specific operators (thermal repair, FTRS, DRFA).
P = ctx.P;
ops = {};
if model.opThermal, ops{end+1} = struct('fn', @(x, pr) op_thermal(x, pr, ctx, model), 'prob', P.pThermal); end
if model.opFTRS && model.useRep, ops{end+1} = struct('fn', @(x, pr) op_ftrs(x, pr, ctx, model), 'prob', P.pFTRS); end
if model.opDRFA, ops{end+1} = struct('fn', @(x, pr) op_drfa(x, pr, ctx, model), 'prob', P.pDRFA); end
end

function [a, r] = split(x, ctx, model)
n = ctx.nT; x = round(x(:)');
a = x(1:n);
if model.useRep && numel(x) >= 2 * n, r = x(n+1:2*n); else, r = zeros(1, n); end
end
function x = join(x, a, r, ctx, model)
n = ctx.nT; x(1:n) = a;
if model.useRep && numel(x) >= 2 * n, x(n+1:2*n) = r; end
end

function x = op_thermal(x, progress, ctx, model)
[a, r] = split(x, ctx, model);
S = trc_proxy_state(a, r, ctx, model);
hot = find(S.active & S.margin < 0);
if isempty(hot), return; end
psiN = ctx.psi / max(ctx.psi);
tasks = find(ismember(a, hot));
tasks = tasks(randperm(numel(tasks)));
for i = tasks
    if S.margin(a(i)) >= 0, continue; end
    pm = 0.1 + 0.3 * (1 - progress) + 0.2 * psiN(i);
    if rand() < pm
        cand = S.margin; cand(a(i)) = -inf; cand(S.margin < 0) = -inf;
        [best, j] = max(cand);
        if isfinite(best)
            a(i) = j;
            if r(i) == j, r(i) = 0; end
            S = trc_proxy_state(a, r, ctx, model);
        end
    end
end
x = join(x, a, r, ctx, model);
end

function x = op_ftrs(x, ~, ctx, model)
[a, r] = split(x, ctx, model);
S = trc_proxy_state(a, r, ctx, model);
risk = ctx.psi .* S.p(a);
if max(risk) <= 0, return; end
thr = (0.3 + 0.4 * rand()) * max(risk);
score = S.margin / 10 + (1 - S.share) - 5 * S.p;
for i = 1:ctx.nT
    if risk(i) >= thr
        s = score; s(a(i)) = -inf;
        [~, ord] = sort(s, 'descend');
        r(i) = ord(randi(min(2, numel(ord))));
    elseif risk(i) < 0.1 * max(risk)
        r(i) = 0;
    end
end
x = join(x, a, r, ctx, model);
end

function x = op_drfa(x, ~, ctx, model)
P = ctx.P;
[a, r] = split(x, ctx, model);
S = trc_proxy_state(a, r, ctx, model);
low = find(S.W > 0 & S.share < P.tauLow);
low = low(randperm(numel(low)));
for j = low
    if any(r == j), continue; end
    targets = find(S.active & (1:numel(S.active)) ~= j);
    if isempty(targets), continue; end
    [~, o] = sort(S.share(targets), 'descend');
    for k = targets(o)
        a2 = a; r2 = r;
        mv = (a == j);
        a2(mv) = k;
        r2(mv & r2 == k) = 0;
        S2 = trc_proxy_state(a2, r2, ctx, model);
        if consolidation_ok(S, S2, k, ctx, model)
            a = a2; r = r2; S = S2;
            break;
        end
    end
end
x = join(x, a, r, ctx, model);
end

function ok = consolidation_ok(S, S2, k, ctx, model)
noThrottle = S2.T(k) <= ctx.fogEnv.Tmax(k);
dE = S.E - S2.E;
dL = S2.L - S.L;
gain = dE > 0 && (dL <= 0 || dE / max(S.E, eps) >= dL / max(S.L, eps));
switch model.drfaRule
    case 'thermal'
        ok = noThrottle;
    case 'failure'
        ok = gain;
    otherwise
        ok = gain && (noThrottle || ~model.throttle);
end
end
