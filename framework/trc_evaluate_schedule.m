function [f, D] = trc_evaluate_schedule(assign, rep, wf, fogEnv, bgUtil, psi, P, mode)
% Thermal-reliability coupled planning model.

n = wf.numTasks;
m = fogEnv.numNodes;

assign = round(min(max(assign, 1), m));
rep    = round(min(max(rep, 0), m));
rep(rep == assign) = 0;

headroom = max(0.15, 1 - bgUtil);
baseCap  = fogEnv.nodeCapacity .* headroom;

active = false(1, m);
active(assign) = true;
active(rep(rep > 0)) = true;

speed = ones(1, m);
fpIters = 0; fpConv = false; fpRes = inf;
fpTol = 1e-3; if isfield(P, 'fpTol'), fpTol = P.fpTol; end
for it = 1:P.thermalIters
    fpIters = it;
    cap = baseCap .* speed;
    [finishTime, busy, execTime] = list_schedule(assign, wf, fogEnv, cap);
    M = max(max(finishTime), eps);
    uTask = min(busy ./ M, 1);
    uHeat = min(1, bgUtil + (1 - bgUtil) .* uTask);
    Tnode = fogEnv.Tambient + P.ambientShift + fogEnv.eta .* uHeat;
    Tnode(~active) = fogEnv.Tambient(~active) + P.ambientShift;
    newSpeed = ones(1, m);
    over = Tnode - fogEnv.Tmax;
    hot = over > 0;
    if isfield(P, 'modelThrottle') && ~P.modelThrottle
        hot(:) = false;
    end
    newSpeed(hot) = max(P.minSpeed, 1 - P.throttleSlope .* over(hot));
    fpRes = max(abs(newSpeed - speed));
    if fpRes < fpTol
        speed = newSpeed; fpConv = true;
        break;
    end
    speed = newSpeed;
end
cap = baseCap .* speed;
[finishTime, busy, execTime] = list_schedule(assign, wf, fogEnv, cap);
M = max(max(finishTime), eps);

if strcmpi(mode, 'coupled')
    AF = arrhenius_af(Tnode, P);
else
    AF = ones(1, m);
end
exposure = M / P.Href;
pFail = 1 - (1 - P.p0) .^ (AF .* exposure);
pFail(~active) = 0;

pTask = pFail(assign);
hasRep = rep > 0;
pTask(hasRep) = pTask(hasRep) .* pFail(rep(hasRep));
w = psi(:)' / max(sum(psi), eps);
Loss = sum(w .* pTask);

uNode = min(1, bgUtil + (1 - bgUtil) .* min(busy ./ M, 1));
Pnode = fogEnv.Pidle + (fogEnv.Ppeak - fogEnv.Pidle) .* uNode;
Pnode(~active) = P.sleepFrac .* fogEnv.Pidle(~active);
E = sum(Pnode .* M) / 3600;
if any(hasRep)
    idx = find(hasRep);
    standbyJ = P.standbyFrac .* (fogEnv.Ppeak(rep(idx)) - fogEnv.Pidle(rep(idx))) .* execTime(idx);
    E = E + sum(standbyJ) / 3600;
end

Cst = sum(fogEnv.nodeCostRate(assign) .* execTime);
if any(hasRep)
    idx = find(hasRep);
    Cst = Cst + P.replicaCostFrac * sum(fogEnv.nodeCostRate(rep(idx)) .* execTime(idx));
end

f = [M, E, Cst, Loss];
f(~isfinite(f)) = 1e8;

if nargout > 1
    D.finishTime = finishTime;
    D.busy       = busy;
    D.execTime   = execTime;
    D.makespan   = M;
    D.Tnode      = Tnode;
    D.speed      = speed;
    D.cap        = cap;
    D.pFail      = pFail;
    D.active     = active;
    D.assign     = assign;
    D.rep        = rep;
    D.fpIters    = fpIters;
    D.fpConv     = fpConv;
    D.fpRes      = fpRes;
end
end

function [finishTime, busy, execTime] = list_schedule(assign, wf, fogEnv, cap)
n = wf.numTasks;
m = fogEnv.numNodes;
finishTime = zeros(1, n);
execTime   = zeros(1, n);
nodeReady  = zeros(1, m);
busy       = zeros(1, m);
for idx = 1:numel(wf.topoOrder)
    i = wf.topoOrder(idx);
    j = assign(i);
    preds = wf.parents{i};
    ready = 0;
    for p = preds(:)'
        d = finishTime(p);
        if assign(p) ~= j
            d = d + (wf.edgeDataSize(p, i) / fogEnv.nodeBandwidth(j)) * 8;
        end
        ready = max(ready, d);
    end
    st = max(ready, nodeReady(j));
    execTime(i) = wf.taskFLOPs(i) / cap(j);
    finishTime(i) = st + execTime(i);
    nodeReady(j) = finishTime(i);
    busy(j) = busy(j) + execTime(i);
end
end
