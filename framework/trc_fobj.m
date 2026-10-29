function [f, D] = trc_fobj(x, ctx, model)
% Objective vector of a chromosome.
n = ctx.nT;
x = round(x(:)');
assign = x(1:n);
if model.useRep && numel(x) >= 2 * n
    rep = x(n+1:2*n);
else
    rep = zeros(1, n);
end
[Pp, bg] = trc_plan_params(ctx, model);
[f4, D] = trc_evaluate_schedule(assign, rep, ctx.wf, ctx.fogEnv, bg, ctx.psi, Pp, model.mode);
if strcmp(model.objective, 'Thermal')
    over = D.Tnode - ctx.fogEnv.Tmax;
    over = over(D.active & over > 0);
    Trisk = sum(exp(over ./ ctx.fogEnv.theta0));
    f = [f4(1:3), Trisk];
else
    f = f4;
end
f(~isfinite(f)) = 1e8;
end
