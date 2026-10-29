function R = simulate_TRC_truth2(assign, rep, wf, fogEnv, bgTrue, psi, P, deadline, seed, promised)
% Second ground truth (transient thermal model, step-wise DTM, Weibull failures).

dflt = struct('numMC2', 100, 'dt2', 1, 'tau2', [120 90 60], 'hyst2', 3, ...
              'levels2', [1 0.85 0.70 0.55 0.40], 'weibullK2', 1.5, ...
              'bgSigma2', 0.08, 'bgTau2', 60);
fn = fieldnames(dflt);
for q = 1:numel(fn), if ~isfield(P, fn{q}), P.(fn{q}) = dflt.(fn{q}); end, end

n = wf.numTasks; m = fogEnv.numNodes;
assign = round(min(max(assign, 1), m));
rep = round(min(max(rep, 0), m)); rep(rep == assign) = 0;
S = P.numMC2; dt = P.dt2;
tau = P.tau2(fogEnv.nodeTier);
aT = 1 - exp(-dt ./ tau);
rhoB = exp(-dt / P.bgTau2);
kW = P.weibullK2;
thetaW = P.Href / (-log(1 - P.p0))^(1 / kW);
levels = P.levels2; nL = numel(levels);
pos = zeros(1, n); pos(wf.topoOrder) = 1:n;
tMax = 50 * max(deadline, promised);

mk = zeros(S,1); en = zeros(S,1); ru = zeros(S,1); mg = zeros(S,1);
ur = zeros(S,1); cl = zeros(S,1); nf = zeros(S,1); mt = zeros(S,1);
th = zeros(S,1); an = zeros(S,1); ag = zeros(S,1); maf = zeros(S,1);

s0 = rng(); rng(seed, 'twister');
for s = 1:S
    bgBase = min(0.95, max(0, bgTrue + 0.05 * randn(1, m)));
    etaT = fogEnv.eta .* max(0.5, 1 + P.etaError * randn(1, m));
    PT = P; PT.Ea = P.Ea * max(0.3, 1 + P.EaError * randn());
    E = -log(rand(1, m));
    x = P.bgSigma2 * randn(1, m);

    awake = false(1, m); awake(assign) = true; awake(rep(rep > 0)) = true;
    initAwake = awake;
    alive = true(1, m);
    Tamb = fogEnv.Tambient + P.ambientShift;
    bg = min(0.95, max(0, bgBase + x));
    T = Tamb + etaT .* bg; T(~awake) = Tamb(~awake);
    L = ones(1, m); everThr = false(1, m);
    D = zeros(1, m);

    host = assign; remain = wf.taskFLOPs; state = zeros(1, n);
    earliest = zeros(1, n); finish = nan(1, n); recovered = false(1, n); lost = false(1, n);
    running = zeros(1, m);
    queue = cell(1, m);
    for j = 1:m, qj = find(assign == j); [~, o] = sort(pos(qj)); queue{j} = qj(o); end

    t = 0; Etot = 0; peak = -inf; agingAcc = 0; afAcc = 0; afCnt = 0;
    nRep = 0; nMig = 0; nUnr = 0; nFail = 0; loss = 0;
    while any(state < 2) && t < tMax
        for j = find(alive & running == 0)
            if isempty(queue{j}), continue; end
            i = queue{j}(1);
            if t + 1e-9 < earliest(i) || ~inputs_ready(i, j), continue; end
            queue{j}(1) = []; running(j) = i; state(i) = 1;
            if ~awake(j), awake(j) = true; end
        end
        busy = running > 0; runSnap = running;
        u = bg + (1 - bg) .* busy;
        capEff = fogEnv.nodeCapacity .* max(0.15, 1 - bg) .* levels(L);
        for j = find(busy)
            i = running(j);
            remain(i) = remain(i) - capEff(j) * dt;
            if remain(i) <= 0
                finish(i) = t + dt + remain(i) / capEff(j);
                state(i) = 2; running(j) = 0;
            end
        end
        Pw = fogEnv.Pidle + (fogEnv.Ppeak - fogEnv.Pidle) .* u;
        Pw(~awake | ~alive) = P.sleepFrac * fogEnv.Pidle(~awake | ~alive);
        Etot = Etot + sum(Pw) * dt;
        for j = find(busy)
            i = runSnap(j);
            if rep(i) > 0 && ~recovered(i) && alive(rep(i))
                Etot = Etot + P.standbyFrac * (fogEnv.Ppeak(rep(i)) - fogEnv.Pidle(rep(i))) * dt;
            end
        end
        Tss = Tamb + etaT .* u;
        on = awake & alive;
        T(on) = T(on) + (Tss(on) - T(on)) .* aT(on);
        T(~on) = T(~on) + (Tamb(~on) - T(~on)) .* aT(~on);
        if any(on), peak = max(peak, max(T(on))); end
        hot = on & T > fogEnv.Tmax; cool = T < fogEnv.Tmax - P.hyst2;
        L(hot) = min(L(hot) + 1, nL); L(cool) = max(L(cool) - 1, 1);
        everThr = everThr | (on & L > 1);
        AF = arrhenius_af(T, PT);
        D(on) = D(on) + AF(on) * dt;
        agingAcc = agingAcc + sum(AF(on)) * dt;
        afAcc = afAcc + sum(AF(on)); afCnt = afCnt + sum(on);
        x = rhoB * x + P.bgSigma2 * sqrt(1 - rhoB^2) * randn(1, m);
        bg = min(0.95, max(0, bgBase + x));
        t = t + dt;
        newFail = find(on & (D / thetaW).^kW >= E);
        for j = newFail
            alive(j) = false; nFail = nFail + 1;
            affected = queue{j}; queue{j} = [];
            wasRun = running(j);
            if wasRun > 0, affected = [wasRun, affected]; running(j) = 0; end
            [~, o] = sort(pos(affected)); affected = affected(o);
            for i = affected
                reroute(i, t, i == wasRun);
            end
        end
    end
    if any(state < 2)
        k = find(state < 2); nUnr = nUnr + numel(k);
        for i = k, if ~lost(i), loss = loss + psi(i); lost(i) = true; end, end
        finish(k) = tMax;
    end

    Mact = max(max(finish), eps);
    mk(s) = Mact; en(s) = Etot / 3600; ru(s) = nRep; mg(s) = nMig; ur(s) = nUnr;
    cl(s) = loss / max(sum(psi), eps); nf(s) = nFail; mt(s) = peak;
    th(s) = sum(everThr); an(s) = sum(awake);
    ag(s) = agingAcc / 3600; maf(s) = afAcc / max(afCnt, 1);
end
rng(s0);

R.makespan = mean(mk); R.makespanP95 = prctile_simple(mk, 95);
R.deadlineMiss = mean(mk > deadline);
R.slaViolation = mean(mk > (1 + P.epsSLA) * promised);
R.slowdown = mean(mk) / max(promised, eps);
R.energy = mean(en); R.replicaUsed = mean(ru); R.migrated = mean(mg);
R.unrecovered = mean(ur); R.critLoss = mean(cl); R.nodeFailures = mean(nf);
R.maxTemp = mean(mt); R.throttledNodes = mean(th); R.activeNodes = mean(an);
R.alarms = 0; R.moves = 0; R.aging = mean(ag); R.meanAF = mean(maf);

    function ok = inputs_ready(i, j)
        ok = true;
        for p = wf.parents{i}(:)'
            if state(p) < 2, ok = false; return; end
            arr = finish(p);
            if host(p) ~= j, arr = arr + wf.edgeDataSize(p, i) * 8 / fogEnv.nodeBandwidth(j); end
            if arr > t + 1e-9, ok = false; return; end
        end
    end

    function reroute(i, tf, wasRunning)
        remain(i) = wf.taskFLOPs(i); state(i) = 0;
        b = rep(i);
        if ~recovered(i) && b > 0 && alive(b) && b ~= host(i)
            host(i) = b; earliest(i) = tf + P.detectDelay + P.replicaDelay;
            recovered(i) = true; nRep = nRep + 1;
            enqueue(i, b);
            return;
        end
        if recovered(i) && wasRunning
            nUnr = nUnr + 1;
            if ~lost(i), loss = loss + psi(i); lost(i) = true; end
            state(i) = 2; finish(i) = tf + 3 * wf.taskFLOPs(i) / fogEnv.nodeCapacity(host(i));
            return;
        end
        cand = find(alive);
        if isempty(cand)
            nUnr = nUnr + 1;
            if ~lost(i), loss = loss + psi(i); lost(i) = true; end
            state(i) = 2; finish(i) = tf + 3 * wf.taskFLOPs(i) / fogEnv.nodeCapacity(host(i));
            return;
        end
        backlog = zeros(1, numel(cand));
        for q = 1:numel(cand)
            j2 = cand(q);
            w = sum(remain(queue{j2}));
            if running(j2) > 0, w = w + remain(running(j2)); end
            backlog(q) = w / fogEnv.nodeCapacity(j2) + P.wakeDelay * ~awake(j2);
        end
        [~, q] = min(backlog); j2 = cand(q);
        extra = P.wakeDelay * ~awake(j2);
        if ~recovered(i), nMig = nMig + 1; end
        host(i) = j2; recovered(i) = true;
        if ~lost(i), loss = loss + psi(i); lost(i) = true; end
        earliest(i) = tf + P.detectDelay + extra + wf.taskDataSize(i) * 8 / fogEnv.nodeBandwidth(j2);
        enqueue(i, j2);
    end

    function enqueue(i, j)
        qj = [queue{j}, i]; [~, o] = sort(pos(qj)); queue{j} = qj(o);
    end
end

function v = prctile_simple(x, p)
x = sort(x(:)); k = max(1, min(numel(x), ceil(p / 100 * numel(x)))); v = x(k);
end
