function [criticality, psi, cp] = CARS_forecast(wf, fogEnv, bgF, P)
% Criticality-aware ranking of the workflow tasks (CARS).

wD = P.carsW(1); wCP = P.carsW(2); wS = P.carsW(3); wPhi = P.carsW(4);

n = wf.numTasks;
capF = mean(fogEnv.nodeCapacity .* max(0.15, 1 - bgF));
bwM  = mean(fogEnv.nodeBandwidth);
w = wf.taskFLOPs / capF;

order = wf.topoOrder(:)';
up = zeros(1, n);
for i = fliplr(order)
    best = 0;
    for c = wf.children{i}(:)'
        best = max(best, wf.edgeDataSize(i, c) / bwM * 8 + up(c));
    end
    up(i) = w(i) + best;
end
down = zeros(1, n);
for i = order
    best = 0;
    for p = wf.parents{i}(:)'
        best = max(best, down(p) + w(p) + wf.edgeDataSize(p, i) / bwM * 8);
    end
    down(i) = best;
end
cp = (up + down) / max(up + down);

Di = wf.taskDeadline; Si = wf.taskDataSize;
reach = false(n);
for i = fliplr(order)
    for c = wf.children{i}(:)'
        reach(i, c) = true;
        reach(i, :) = reach(i, :) | reach(c, :);
    end
end
phi = sum(reach, 2)' / max(n - 1, 1);
Dmax = max(max(Di), eps); Smax = max(max(Si), eps);
psi = wD .* (Dmax - Di) ./ Dmax + wCP .* cp + wS .* (Si ./ Smax) + wPhi .* phi;

if isfield(P, 'carsTop'), top = P.carsTop; else, top = [0.5 0.2]; end
[~, ord] = sort(psi, 'descend');
rk = zeros(1, n); rk(ord) = 1:n;
criticality = ones(1, n);
criticality(rk <= ceil(top(1) * n)) = 2;
criticality(rk <= ceil(top(2) * n)) = 3;
end
