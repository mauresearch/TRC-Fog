function S = trc_proxy_state(assign, rep, ctx, model)
% Closed-form estimate of the node state.
[Pp, bg] = trc_plan_params(ctx, model);
fog = ctx.fogEnv; m = fog.numNodes; n = ctx.nT;
assign = round(min(max(assign(:)', 1), m));
rep = round(min(max(rep(:)', 0), m)); rep(rep == assign) = 0;

W = accumarray(assign(:), ctx.wf.taskFLOPs(:), [m 1])';
capEff = fog.nodeCapacity .* max(0.15, 1 - bg);
active = W > 0;
active(rep(rep > 0)) = true;
speed = ones(1, m);
wf = ctx.wf;
for it = 1:2
    t = W ./ (capEff .* speed);
    tau = wf.taskFLOPs ./ (capEff(assign) .* speed(assign));
    fin = zeros(1, n);
    for q = 1:numel(wf.topoOrder)
        i = wf.topoOrder(q);
        r0 = 0;
        for pp = wf.parents{i}(:)'
            d = fin(pp);
            if assign(pp) ~= assign(i)
                d = d + (wf.edgeDataSize(pp, i) / fog.nodeBandwidth(assign(i))) * 8;
            end
            r0 = max(r0, d);
        end
        fin(i) = r0 + tau(i);
    end
    M = max([max(t), max(fin), eps]);
    u = t / M;
    T = fog.Tambient + Pp.ambientShift + fog.eta .* (bg + (1 - bg) .* u);
    T(~active) = fog.Tambient(~active) + Pp.ambientShift;
    newSpeed = ones(1, m);
    if model.throttle
        over = T - fog.Tmax; hot = over > 0;
        newSpeed(hot) = max(Pp.minSpeed, 1 - Pp.throttleSlope .* over(hot));
    end
    if max(abs(newSpeed - speed)) < 1e-3, break; end
    speed = newSpeed;
end
if strcmp(model.mode, 'coupled'), AF = arrhenius_af(T, Pp); else, AF = ones(1, m); end
p = 1 - (1 - Pp.p0) .^ (AF .* (M / Pp.Href));
p(~active) = 0;
Pn = fog.Pidle + (fog.Ppeak - fog.Pidle) .* min(1, bg + (1 - bg) .* u);
Pn(~active) = Pp.sleepFrac .* fog.Pidle(~active);
pt = p(assign);
has = rep > 0;
pt(has) = pt(has) .* p(rep(has));

S.W = W; S.T = T; S.share = u; S.active = active; S.p = p; S.M = M;
S.E = sum(Pn .* M);
S.L = sum(ctx.psi .* pt) / max(sum(ctx.psi), eps);
S.margin = fog.Tmax - T;
S.n = n;
end
