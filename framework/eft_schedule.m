function assign = eft_schedule(wf, fogEnv, bgUtil)
% Earliest-finish-time list scheduling.
n = wf.numTasks;
m = fogEnv.numNodes;
cap = fogEnv.nodeCapacity .* max(0.15, 1 - bgUtil);
assign = zeros(1, n);
finishTime = zeros(1, n);
nodeReady = zeros(1, m);
for k = 1:numel(wf.topoOrder)
    i = wf.topoOrder(k);
    bestFin = inf; bestJ = 1;
    for j = 1:m
        r = 0;
        for p = wf.parents{i}(:)'
            d = finishTime(p);
            if assign(p) ~= j
                d = d + (wf.edgeDataSize(p, i) / fogEnv.nodeBandwidth(j)) * 8;
            end
            r = max(r, d);
        end
        fin = max(r, nodeReady(j)) + wf.taskFLOPs(i) / cap(j);
        if fin < bestFin, bestFin = fin; bestJ = j; end
    end
    assign(i) = bestJ;
    finishTime(i) = bestFin;
    nodeReady(bestJ) = bestFin;
end
end
