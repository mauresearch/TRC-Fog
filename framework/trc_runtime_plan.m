function rt = trc_runtime_plan(out, ctx, model)
% Plan information for the predictive SLA controller (PSC).
rt = struct('mode', model.runtime);
if strcmp(rt.mode, 'none'), return; end
rt.thermal = numel(rt.mode) > 2 && strcmp(rt.mode(end-1:end), '@T');
if rt.thermal, rt.mode = rt.mode(1:end-2); end
rt.reliable = strcmp(rt.mode, 'PSC-R');
P = ctx.P; wf = ctx.wf; n = ctx.nT;
[Pp, bg] = trc_plan_params(ctx, model);
[~, D] = trc_evaluate_schedule(out.assign, out.rep, wf, ctx.fogEnv, bg, ctx.psi, Pp, model.mode);
rt.plannedFinish = D.finishTime;
rt.plannedExec = D.execTime;
rt.pPlan = D.pFail;
tail = zeros(1, n);
for q = numel(wf.topoOrder):-1:1
    i = wf.topoOrder(q);
    c = wf.children{i};
    if ~isempty(c), tail(i) = max(D.execTime(c) + tail(c)); end
end
rt.tail = tail;
psi = ctx.psi(:)';
psiN = (psi - min(psi)) / max(max(psi) - min(psi), eps);
rt.psiN = psiN;
switch rt.mode
    case {'PSC', 'PSC-R'}
        rt.eps = P.epsMax - (P.epsMax - P.epsMin) .* psiN;  rt.propagate = true;
    case 'PSC-U'
        rt.eps = P.epsFixed * ones(1, n);                    rt.propagate = true;
    case 'Local'
        rt.eps = P.epsFixed * ones(1, n);                    rt.propagate = false;
    otherwise
        error('trc_runtime_plan: unknown runtime mode %s', rt.mode);
end
rt.bound = (1 + P.epsSLA) * out.promised;
end
