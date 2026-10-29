function [assign, rep] = psb_schedule(wf, fogEnv, bgUtil)
% Primary-secondary backup scheduling (PSB), She et al., Future Gener. Comput. Syst. 168 (2025) 107738.

n = wf.numTasks;
m = fogEnv.numNodes;
cap = fogEnv.nodeCapacity .* max(0.15, 1 - bgUtil);
W = wf.taskFLOPs(:) ./ cap;
wbar = mean(W, 2)';
bbar = mean(fogEnv.nodeBandwidth);

urank = zeros(1, n);
for q = numel(wf.topoOrder):-1:1
    i = wf.topoOrder(q);
    best = 0;
    for j = wf.children{i}(:)'
        cij = wf.edgeDataSize(i, j) * 8 / bbar;
        best = max(best, cij + urank(j));
    end
    urank(i) = wbar(i) + best;
end
[~, order] = sort(urank, 'descend');
order = topo_consistent(order, wf);

avail = zeros(1, m);
AFT = zeros(1, n);
assign = zeros(1, n);
for i = order
    best = inf; bj = 1;
    for k = 1:m
        eat = 0;
        for p = wf.parents{i}(:)'
            d = AFT(p);
            if assign(p) ~= k, d = d + wf.edgeDataSize(p, i) * 8 / fogEnv.nodeBandwidth(k); end
            eat = max(eat, d);
        end
        est = max(avail(k), eat);
        eft = est + W(i, k);
        if eft < best, best = eft; bj = k; end
    end
    assign(i) = bj; AFT(i) = best; avail(bj) = best;
end

rep = zeros(1, n);
if m < 2, return; end
for i = order
    best = inf; bj = 0;
    for k = 1:m
        if k == assign(i), continue; end
        eat = 0;
        for p = wf.parents{i}(:)'
            d = AFT(p);
            if assign(p) ~= k, d = d + wf.edgeDataSize(p, i) * 8 / fogEnv.nodeBandwidth(k); end
            eat = max(eat, d);
        end
        eft = max(avail(k), eat) + W(i, k);
        if eft < best, best = eft; bj = k; end
    end
    rep(i) = bj; avail(bj) = best;
end
end

function order = topo_consistent(order, wf)
n = numel(order); placed = false(1, n); out = zeros(1, n); c = 0;
while c < n
    for q = 1:n
        i = order(q);
        if ~placed(i) && all(placed(wf.parents{i}))
            c = c + 1; out(c) = i; placed(i) = true; break;
        end
    end
end
order = out;
end
