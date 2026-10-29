function assign = tdrm_schedule(wf, fogEnv, bgUtil, P, deadline)
% Temperature- and deadline-aware resource management (TDRM), Niknafs et al., J. Syst. Archit. 167 (2025) 103506.

n = wf.numTasks;
m = fogEnv.numNodes;
cap = fogEnv.nodeCapacity .* max(0.15, 1 - bgUtil);
W = wf.taskFLOPs(:) ./ cap;
E = W .* fogEnv.Ppeak;
K = deadline;
Kres = K * ones(1, m);
busy = zeros(1, m);
nextPt = zeros(1, m);
finish = zeros(1, n);
assign = zeros(1, n);
mapped = false(1, n);
Tamb = fogEnv.Tambient + P.ambientShift;

for step = 1:n
    ready = find(~mapped & arrayfun(@(i) all(mapped(wf.parents{i})), 1:n));
    ready = wf.topoOrder(ismember(wf.topoOrder, ready));
    dbest = inf; jstar = ready(1); Fstar = [];
    for j = ready
        Fj = find(Kres >= W(j, :));
        if isempty(Fj), d = inf; elseif numel(Fj) == 1, d = 0; else, d = deadline; end
        if d < dbest, dbest = d; jstar = j; Fstar = Fj; end
    end
    j = jstar;

    cand = Fstar;
    placed = false;
    while ~isempty(cand)
        [~, q] = min(E(j, cand));
        i = cand(q);
        [ok, st] = schedulable(j, i);
        if ok
            commit(j, i, st); placed = true; break;
        end
        cand(q) = [];
    end
    if ~placed
        Tnew = zeros(1, m);
        for i = 1:m
            u = min(1, (busy(i) + W(j, i)) / K);
            Tnew(i) = Tamb(i) + fogEnv.eta(i) * (bgUtil(i) + (1 - bgUtil(i)) * u);
        end
        [~, i] = min(Tnew - fogEnv.Tmax);
        commit(j, i, start_on(j, i));
    end
end

    function [ok, st] = schedulable(j, i)
        st = start_on(j, i);
        fin = st + W(j, i);
        u = min(1, (busy(i) + W(j, i)) / K);
        T = Tamb(i) + fogEnv.eta(i) * (bgUtil(i) + (1 - bgUtil(i)) * u);
        ok = fin <= deadline && T <= fogEnv.Tmax(i);
    end

    function st = start_on(j, i)
        r = 0;
        for p = wf.parents{j}(:)'
            d = finish(p);
            if assign(p) ~= i, d = d + wf.edgeDataSize(p, j) * 8 / fogEnv.nodeBandwidth(i); end
            r = max(r, d);
        end
        st = max(r, nextPt(i));
    end

    function commit(j, i, st)
        assign(j) = i; mapped(j) = true;
        finish(j) = st + W(j, i);
        nextPt(i) = finish(j);
        busy(i) = busy(i) + W(j, i);
        Kres(i) = Kres(i) - W(j, i);
    end
end
