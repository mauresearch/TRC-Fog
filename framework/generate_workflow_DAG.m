function wf = generate_workflow_DAG(workflowName, numTasks, seed)
% Synthetic scientific workflows (Montage, Epigenomics, CyberShake, Inspiral, SIPHT).

if nargin < 3 || isempty(seed)
    seed = 42;
end

s = rng();
rng(seed, 'twister');

adjMatrix = zeros(numTasks, numTasks);

switch workflowName
    case 'Montage'
        nRoot = 1;
        nFinal = 1;
        nMiddle = numTasks - nRoot - nFinal;
        if nMiddle < 2
            nMiddle = max(2, numTasks - 2);
            nFinal = numTasks - nRoot - nMiddle;
        end
        nLayer1 = max(1, round(0.6 * nMiddle));
        nLayer2 = nMiddle - nLayer1;
        if nLayer2 < 1
            nLayer2 = 1;
            nLayer1 = nMiddle - nLayer2;
        end

        rootIdx   = 1;
        layer1Idx = (rootIdx+1) : (rootIdx+nLayer1);
        layer2Idx = (layer1Idx(end)+1) : (layer1Idx(end)+nLayer2);
        finalIdx  = numTasks;

        for j = layer1Idx
            adjMatrix(rootIdx, j) = 1;
        end
        for j = layer2Idx
            srcCount = min(numel(layer1Idx), randi([2, max(2, min(4, numel(layer1Idx)))]));
            srcs = layer1Idx(randperm(numel(layer1Idx), srcCount));
            for i = srcs
                adjMatrix(i, j) = 1;
            end
        end
        for i = layer2Idx
            adjMatrix(i, finalIdx) = 1;
        end
        if isempty(layer2Idx)
            for i = layer1Idx
                adjMatrix(i, finalIdx) = 1;
            end
        end

        taskFLOPs    = 1e10 * (0.5 + 1.5 * rand(1, numTasks));
        taskDataSize = 20 + 80 * rand(1, numTasks);
        taskDeadlineBase = 8;

    case 'CyberShake'
        numStages = 5;
        tasksPerStage = floor(numTasks / numStages) * ones(1, numStages);
        remainder = numTasks - sum(tasksPerStage);
        tasksPerStage(1:remainder) = tasksPerStage(1:remainder) + 1;
        tasksPerStage(tasksPerStage < 1) = 1;
        while sum(tasksPerStage) > numTasks
            [~, idxMax] = max(tasksPerStage);
            tasksPerStage(idxMax) = tasksPerStage(idxMax) - 1;
        end

        stageStart = cumsum([1, tasksPerStage(1:end-1)]);
        stageIdx = cell(1, numStages);
        for st = 1:numStages
            stageIdx{st} = stageStart(st) : (stageStart(st) + tasksPerStage(st) - 1);
        end

        for st = 1:numStages-1
            srcNodes = stageIdx{st};
            dstNodes = stageIdx{st+1};
            for i = srcNodes
                nLinks = min(numel(dstNodes), randi([1, max(1, min(2, numel(dstNodes)))]));
                dsts = dstNodes(randperm(numel(dstNodes), nLinks));
                for j = dsts
                    adjMatrix(i, j) = 1;
                end
            end
            for j = dstNodes
                if ~any(adjMatrix(srcNodes, j))
                    i = srcNodes(randi(numel(srcNodes)));
                    adjMatrix(i, j) = 1;
                end
            end
        end

        taskFLOPs    = 1e10 * (2 + 8 * rand(1, numTasks));
        taskDataSize = 5 + 30 * rand(1, numTasks);
        taskDeadlineBase = 5;

    case 'Epigenomics'
        layerRatios = [0.32, 0.26, 0.20, 0.12, 0.06, 0.04];
        nLayersE = numel(layerRatios);
        layerSizes = max(1, round(layerRatios * numTasks));
        diffN = numTasks - sum(layerSizes);
        layerSizes(1) = layerSizes(1) + diffN;
        layerSizes(layerSizes < 1) = 1;
        while sum(layerSizes) > numTasks
            [~, idxMax] = max(layerSizes);
            layerSizes(idxMax) = layerSizes(idxMax) - 1;
        end
        while sum(layerSizes) < numTasks
            layerSizes(1) = layerSizes(1) + 1;
        end

        layerIdxE = cell(1, nLayersE);
        cursor = 1;
        for L = 1:nLayersE
            layerIdxE{L} = cursor : (cursor + layerSizes(L) - 1);
            cursor = cursor + layerSizes(L);
        end

        for L = 1:nLayersE-1
            srcNodes = layerIdxE{L};
            dstNodes = layerIdxE{L+1};
            for j = dstNodes
                srcCount = min(numel(srcNodes), randi([2, max(2, min(5, numel(srcNodes)))]));
                srcs = srcNodes(randperm(numel(srcNodes), srcCount));
                for i = srcs
                    adjMatrix(i, j) = 1;
                end
            end
            for i = srcNodes
                if ~any(adjMatrix(i, dstNodes))
                    j = dstNodes(randi(numel(dstNodes)));
                    adjMatrix(i, j) = 1;
                end
            end
        end

        taskFLOPs    = 1e10 * (3 + 9 * rand(1, numTasks));
        taskDataSize = 30 + 90 * rand(1, numTasks);
        taskDeadlineBase = 9;

    case 'Inspiral'
        nHalf1 = floor(numTasks / 2);
        nHalf2 = numTasks - nHalf1;
        halves = {1:nHalf1, (nHalf1+1):numTasks};

        for h = 1:2
            members = halves{h};
            nM = numel(members);
            if nM < 3
                continue;
            end
            nTop = max(1, round(0.45 * nM));
            nMid = max(1, round(0.40 * nM));
            nBot = nM - nTop - nMid;
            if nBot < 1
                nBot = 1;
                nMid = nM - nTop - nBot;
                if nMid < 1
                    nMid = 1;
                    nTop = nM - nMid - nBot;
                end
            end
            topIdx = members(1:nTop);
            midIdx = members(nTop+1 : nTop+nMid);
            botIdx = members(nTop+nMid+1 : end);

            for j = midIdx
                srcCount = min(numel(topIdx), randi([2, max(2, min(4, numel(topIdx)))]));
                srcs = topIdx(randperm(numel(topIdx), srcCount));
                for i = srcs
                    adjMatrix(i, j) = 1;
                end
            end
            for i = midIdx
                for j = botIdx
                    adjMatrix(i, j) = 1;
                end
            end
        end

        taskFLOPs    = 1e10 * (1 + 4 * rand(1, numTasks));
        taskDataSize = 4 + 15 * rand(1, numTasks);
        taskDeadlineBase = 6;

    case 'SIPHT'
        nRootsS = max(2, min(4, max(1, round(0.06 * numTasks))));
        nFinalChain = min(4, max(1, floor(numTasks * 0.02)));
        nHub = 1;
        remaining = numTasks - nRootsS - nHub - nFinalChain;
        if remaining < 4
            remaining = max(2, numTasks - nRootsS - nHub - 1);
            nFinalChain = numTasks - nRootsS - nHub - remaining;
            if nFinalChain < 1
                nFinalChain = 1;
            end
        end
        nGrid = remaining;

        rootIdxS = 1:nRootsS;
        hubIdx   = nRootsS + 1;
        gridIdx  = (hubIdx+1) : (hubIdx + nGrid);
        chainIdx = (hubIdx + nGrid + 1) : numTasks;

        for i = rootIdxS
            adjMatrix(i, hubIdx) = 1;
        end
        for j = gridIdx
            adjMatrix(hubIdx, j) = 1;
        end
        if ~isempty(chainIdx)
            for i = gridIdx
                adjMatrix(i, chainIdx(1)) = 1;
            end
            for k = 1:numel(chainIdx)-1
                adjMatrix(chainIdx(k), chainIdx(k+1)) = 1;
            end
        end

        taskFLOPs    = 1e10 * (1.5 + 6 * rand(1, numTasks));
        taskDataSize = 10 + 40 * rand(1, numTasks);
        taskDeadlineBase = 7;

    otherwise
        rng(s);
        error('Workflow not recognized: %s (use one of: Montage, Epigenomics, CyberShake, Inspiral, SIPHT)', workflowName);
end

edgeDataSize = zeros(numTasks, numTasks);
[srcList, dstList] = find(adjMatrix);
for e = 1:numel(srcList)
    i = srcList(e); j = dstList(e);
    edgeDataSize(i,j) = taskDataSize(i) * (0.3 + 0.4*rand());
end

depth = compute_topo_depth(adjMatrix, numTasks);
maxDepth = max(depth);
if maxDepth == 0, maxDepth = 1; end
taskDeadline = taskDeadlineBase * (1 + (maxDepth - depth) ./ maxDepth) .* (0.8 + 0.4*rand(1, numTasks));
taskFailImpact = rand(1, numTasks);
outDeg = sum(adjMatrix, 2)';
inDeg  = sum(adjMatrix, 1);
taskFailImpact = min(1, taskFailImpact + 0.1*(outDeg > 2) + 0.1*(inDeg > 2));

parents = cell(1, numTasks);
children = cell(1, numTasks);
for i = 1:numTasks
    parents{i} = find(adjMatrix(:,i))';
    children{i} = find(adjMatrix(i,:));
end
topoOrder = topological_sort(adjMatrix, numTasks);

rng(s);

wf.name          = workflowName;
wf.numTasks      = numTasks;
wf.adjMatrix     = adjMatrix;
wf.parents       = parents;
wf.children      = children;
wf.taskFLOPs     = taskFLOPs;
wf.taskDataSize  = taskDataSize;
wf.edgeDataSize  = edgeDataSize;
wf.taskDeadline  = taskDeadline;
wf.taskFailImpact = taskFailImpact;
wf.topoOrder     = topoOrder;

end
