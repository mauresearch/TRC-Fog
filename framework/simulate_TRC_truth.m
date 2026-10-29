function R = simulate_TRC_truth(assign, rep, wf, fogEnv, bgForecast, psi, P, deadline, seed, promised, rt)
% Monte-Carlo ground-truth simulator (steady-state thermal model).

if nargin < 11 || isempty(rt), rt = struct('mode', 'none'); end
rtOn = ~strcmp(rt.mode, 'none');
n = wf.numTasks;
m = fogEnv.numNodes;
assign = round(min(max(assign, 1), m));
rep    = round(min(max(rep, 0), m));
rep(rep == assign) = 0;

S = P.numMC;
mk = zeros(S,1); en = zeros(S,1); ru = zeros(S,1); mg = zeros(S,1);
ur = zeros(S,1); cl = zeros(S,1); nf = zeros(S,1); mt = zeros(S,1);
th = zeros(S,1); an = zeros(S,1); al = zeros(S,1); mv = zeros(S,1);
ag = zeros(S,1); maf = zeros(S,1);

s0 = rng();
rng(seed, 'twister');

for s = 1:S
    bg = min(0.95, max(0, bgForecast + P.bgNoise * randn(1, m)));

    envT = fogEnv;
    envT.eta = fogEnv.eta .* max(0.5, 1 + P.etaError * randn(1, m));
    PT = P;
    PT.Ea = P.Ea * max(0.3, 1 + P.EaError * randn());

    [~, D] = trc_evaluate_schedule(assign, rep, wf, envT, bg, psi, PT, 'coupled');
    cap = D.cap;
    Mplan = D.makespan;

    failTime = inf(1, m);
    Ufail = rand(1, m);
    fails = (Ufail < D.pFail) & D.active;
    failTime(fails) = rand(1, sum(fails)) * Mplan;
    nf(s) = sum(fails);

    initActive = D.active;
    busyP = D.busy;
    baseCap = fogEnv.nodeCapacity .* max(0.15, 1 - bg);
    woken = false(1, m);
    nodeReady = zeros(1, m);
    busyAct = zeros(1, m);
    finish = zeros(1, n);
    hostOf = zeros(1, n);
    nRep = 0; nMig = 0; nUnr = 0; loss = 0; nAlarm = 0; nMove = 0;

    for k = 1:numel(wf.topoOrder)
        i = wf.topoOrder(k);
        a = assign(i); b = rep(i);

        readyA = ready_on(i, a);
        extra = 0;
        if failTime(a) > readyA
            host = a;  rdy = readyA;
            if rtOn, [host, rdy, extra] = psc_control(i, host, rdy, extra, b); end
        elseif b > 0 && failTime(b) > ready_on(i, b)
            host = b;  rdy = ready_on(i, b);  extra = P.replicaDelay;  nRep = nRep + 1;
        else
            [host, extra] = pick_alive(max(readyA, failTime(a)) + P.detectDelay);
            if host == 0
                finish(i) = readyA + 3 * wf.taskFLOPs(i) / cap(a);  hostOf(i) = a;
                nUnr = nUnr + 1;  loss = loss + psi(i);
                continue;
            end
            rdy = max(ready_on(i, host), failTime(a) + P.detectDelay);
            nMig = nMig + 1;  loss = loss + psi(i);
        end

        st = max(rdy + extra, nodeReady(host));
        ex = wf.taskFLOPs(i) / cap(host);
        fin = st + ex;

        if failTime(host) < fin
            busyAct(host) = busyAct(host) + max(0, failTime(host) - st);
            t0 = failTime(host) + P.detectDelay;
            if (host == a || (rtOn && host ~= b)) && b > 0 && failTime(b) > t0
                host2 = b;  st2 = max(t0 + P.replicaDelay, nodeReady(b));
                nRep = nRep + 1;
            else
                [host2, extra2] = pick_alive(t0);
                if host2 == 0
                    finish(i) = fin + 3 * ex;  hostOf(i) = host;
                    nUnr = nUnr + 1;  loss = loss + psi(i);
                    continue;
                end
                st2 = max(t0 + extra2 + wf.taskDataSize(i) / fogEnv.nodeBandwidth(host2) * 8, nodeReady(host2));
                nMig = nMig + 1;  loss = loss + psi(i);
            end
            ex = wf.taskFLOPs(i) / cap(host2);
            fin = st2 + ex;
            host = host2;
            if failTime(host2) < fin
                nUnr = nUnr + 1;  loss = loss + psi(i);
                fin = fin + 3 * ex;
            end
        end

        nodeReady(host) = fin;
        busyAct(host) = busyAct(host) + ex;
        finish(i) = fin;
        hostOf(i) = host;
    end

    Mact = max(max(finish), eps);
    act = initActive | woken;
    uTask = min(busyAct ./ Mact, 1);
    uNode = min(1, bg + (1 - bg) .* uTask);
    Pw = fogEnv.Pidle + (fogEnv.Ppeak - fogEnv.Pidle) .* uNode;
    Pw(~act) = P.sleepFrac .* fogEnv.Pidle(~act);
    E = sum(Pw .* Mact) / 3600;
    hasRep = rep > 0;
    if any(hasRep)
        idx = find(hasRep);
        E = E + sum(P.standbyFrac .* (fogEnv.Ppeak(rep(idx)) - fogEnv.Pidle(rep(idx))) .* D.execTime(idx)) / 3600;
    end

    mk(s) = Mact; en(s) = E; ru(s) = nRep; mg(s) = nMig; ur(s) = nUnr;
    cl(s) = loss / max(sum(psi), eps);
    mt(s) = max(D.Tnode(D.active));
    th(s) = sum(D.speed < 1 & D.active);
    an(s) = sum(act);
    al(s) = nAlarm; mv(s) = nMove;
    AFt = arrhenius_af(D.Tnode, PT);
    ag(s) = sum(AFt(act)) * Mact / 3600;
    maf(s) = mean(AFt(act));
end
rng(s0);

R.makespan       = mean(mk);
R.makespanP95    = prctile_simple(mk, 95);
R.deadlineMiss   = mean(mk > deadline);
R.slaViolation   = mean(mk > (1 + P.epsSLA) * promised);
R.slowdown       = mean(mk) / max(promised, eps);
R.energy         = mean(en);
R.replicaUsed    = mean(ru);
R.migrated       = mean(mg);
R.unrecovered    = mean(ur);
R.critLoss       = mean(cl);
R.nodeFailures   = mean(nf);
R.maxTemp        = mean(mt);
R.throttledNodes = mean(th);
R.activeNodes    = mean(an);
R.alarms         = mean(al);
R.moves          = mean(mv);
R.aging          = mean(ag);
R.meanAF         = mean(maf);

    function r = ready_on(i, h)
        r = 0;
        pr = wf.parents{i};
        for p = pr(:)'
            d = finish(p);
            if hostOf(p) ~= h
                d = d + (wf.edgeDataSize(p, i) / fogEnv.nodeBandwidth(h)) * 8;
            end
            r = max(r, d);
        end
    end

    function [h, r, ex] = psc_control(i, h, r, ex, bkp)
        finA = max(r + ex, nodeReady(h)) + wf.taskFLOPs(i) / cap(h);
        alarm = false;
        for p = wf.parents{i}(:)'
            dev = (finish(p) - rt.plannedFinish(p)) / max(rt.plannedExec(p), eps);
            if dev > rt.eps(p), alarm = true; end
        end
        if rt.propagate
            Mhat = max(max(finish), finA + rt.tail(i));
            trig = alarm || Mhat > rt.bound;
        else
            trig = alarm;
        end
        if ~trig, return; end
        nAlarm = nAlarm + 1;
        best = finA; bh = 0; bEx = 0; bR = 0;
        for hh = find(failTime > r & (initActive | woken))
            if hh == h || hh == bkp, continue; end
            if isfield(rt, 'reliable') && rt.reliable && ...
                    rt.pPlan(hh) > rt.pPlan(h) * (1 + (1 - rt.psiN(i)) * P.pscRiskTol)
                continue;
            end
            rr = ready_on(i, hh);
            f = max(rr, nodeReady(hh)) + wf.taskFLOPs(i) / cap(hh);
            if f < best, best = f; bh = hh; bR = rr; end
        end
        if bh > 0 && best < finA - P.pscGain * (finA - r)
            h0 = h;
            h = bh; r = bR; ex = bEx; nMove = nMove + 1;
            if isfield(rt, 'thermal') && rt.thermal
                busyP(h0) = max(0, busyP(h0) - wf.taskFLOPs(i) / cap(h0));
                busyP(h) = busyP(h) + wf.taskFLOPs(i) / cap(h);
                update_node(h0, r); update_node(h, r);
            end
        end
    end

    function update_node(j, tnow)
        u = min(1, bg(j) + (1 - bg(j)) * min(busyP(j) / Mplan, 1));
        T = fogEnv.Tambient(j) + PT.ambientShift + envT.eta(j) * u;
        sp = 1;
        if T > fogEnv.Tmax(j), sp = max(PT.minSpeed, 1 - PT.throttleSlope * (T - fogEnv.Tmax(j))); end
        cap(j) = baseCap(j) * sp;
        pj = 1 - (1 - PT.p0) ^ (arrhenius_af(T, PT) * Mplan / PT.Href);
        if failTime(j) > tnow
            if Ufail(j) < pj
                if isinf(failTime(j))
                    failTime(j) = tnow + (Ufail(j) / pj) * max(Mplan - tnow, 0);
                    nf(s) = nf(s) + 1;
                end
            elseif ~isinf(failTime(j))
                failTime(j) = inf;
                nf(s) = nf(s) - 1;
            end
        end
    end

    function [h, extraDelay] = pick_alive(t)
        alive = find(failTime > t);
        h = 0; extraDelay = 0;
        if isempty(alive), return; end
        wakePen = zeros(1, numel(alive));
        sleeping = ~(initActive(alive) | woken(alive));
        wakePen(sleeping) = P.wakeDelay;
        [~, q] = min(max(nodeReady(alive), t) + wakePen);
        h = alive(q);
        extraDelay = wakePen(q);
        if sleeping(q), woken(h) = true; end
    end
end

function v = prctile_simple(x, p)
x = sort(x(:));
k = max(1, min(numel(x), ceil(p / 100 * numel(x))));
v = x(k);
end
