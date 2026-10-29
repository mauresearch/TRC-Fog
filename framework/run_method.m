function out = run_method(name, ctx, PopSize, MaxIter)
% Plans a schedule with one method.
model = trc_models(name);
n = ctx.nT; m = ctx.fogEnv.numNodes;
t0 = tic;

switch model.special
    case 'EFT'
        [~, bg] = trc_plan_params(ctx, model);
        assign = eft_schedule(ctx.wf, ctx.fogEnv, bg);
        rep = zeros(1, n);
        F = trc_fobj(assign, ctx, model);

    case 'PSB'
        [~, bg] = trc_plan_params(ctx, model);
        [assign, rep] = psb_schedule(ctx.wf, ctx.fogEnv, bg);
        F = trc_fobj([assign, rep], ctx, model);

    case 'TDRM'
        [~, bg] = trc_plan_params(ctx, model);
        assign = tdrm_schedule(ctx.wf, ctx.fogEnv, bg, ctx.P, ctx.deadline);
        rep = zeros(1, n);
        F = trc_fobj(assign, ctx, model);

    case 'EqE'
        [X, F] = optimise(model, ctx, n, m, PopSize, MaxIter);
        [~, bg] = trc_plan_params(ctx, model);
        [aP, rP] = psb_schedule(ctx.wf, ctx.fogEnv, bg);
        fP = trc_fobj([aP, rP], ctx, model);
        ok = F(:, 2) <= fP(2);
        if any(ok)
            idx = find(ok); [~, q] = min(F(idx, 4)); q = idx(q);
        else
            [~, q] = min(F(:, 2));
        end
        xb = round(X(q, :));
        assign = xb(1:n); rep = xb(n+1:2*n);

    case 'Sequential'
        [X, F] = optimise(model, ctx, n, m, PopSize, MaxIter);
        xb = X(compromise_point(F, eft_reference(model, ctx, n)), :);
        assign = round(xb(1:n));
        [~, bg] = trc_plan_params(ctx, model);
        Tpred = ctx.fogEnv.Tambient + ctx.P.ambientShift + ctx.fogEnv.eta .* bg;
        score = (1 - bg) + 0.02 * (ctx.fogEnv.Tmax - Tpred);
        rep = zeros(1, n);
        for i = find(ctx.crit == 3)
            s = score; s(assign(i)) = -inf;
            [~, rep(i)] = max(s);
        end
        mdl = model; mdl.useRep = true; mdl.opThermal = false; mdl.opFTRS = false;
        mdl.opDRFA = true; mdl.drfaRule = 'thermal';
        ops = trc_operators(ctx, mdl);
        x = ops{1}.fn([assign, rep], 1);
        assign = x(1:n); rep = x(n+1:2*n);
        model.useRep = true;
        F = trc_fobj([assign, rep], ctx, model);

    otherwise
        [X, F] = optimise(model, ctx, n, m, PopSize, MaxIter);
        xb = round(X(compromise_point(F, eft_reference(model, ctx, n)), :));
        assign = xb(1:n);
        if model.useRep, rep = xb(n+1:2*n); else, rep = zeros(1, n); end
end

rep(rep == assign) = 0;
fsel = trc_fobj([assign, rep], ctx, setfield(model, 'useRep', true));
out.assign = assign;
out.rep = rep;
out.promised = fsel(1);
out.F = F;
out.time = toc(t0);
end

function [X, F] = optimise(model, ctx, n, m, PopSize, MaxIter)
if model.useRep
    dim = 2 * n; lb = [ones(1, n), zeros(1, n)]; ub = m * ones(1, 2 * n);
else
    dim = n; lb = ones(1, n); ub = m * ones(1, n);
end
fob = @(x) trc_fobj(x, ctx, model);
if strcmp(model.optimizer, 'IMOWO')
    opts = model.imowoOpts;
    opts.integer = true;
    opts.operators = trc_operators(ctx, model);
    [X, F] = IMOWO(PopSize, MaxIter, lb, ub, dim, fob, 4, opts);
else
    [X, F] = feval(model.optimizer, PopSize, MaxIter, lb, ub, dim, fob, 4, struct('integer', true));
end
ok = all(isfinite(F), 2) & all(F < 1e7, 2);
X = X(ok, :); F = F(ok, :);
if isempty(F), error('run_method: %s returned no feasible solution', model.name); end
end

function fref = eft_reference(model, ctx, n)
ref = trc_models('TRC');
[~, bg] = trc_plan_params(ctx, ref);
a = eft_schedule(ctx.wf, ctx.fogEnv, bg);
fref = trc_fobj([a, zeros(1, n)], ctx, ref);
if strcmp(model.objective, 'Thermal')
    [~, bgm] = trc_plan_params(ctx, model);
    am = eft_schedule(ctx.wf, ctx.fogEnv, bgm);
    if model.useRep, xm = [am, zeros(1, n)]; else, xm = am; end
    fown = trc_fobj(xm, ctx, model);
    fref(4) = fown(4);
end
end
