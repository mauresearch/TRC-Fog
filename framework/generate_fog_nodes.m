function fogEnv = generate_fog_nodes(numNodes, seed)

if nargin < 2 || isempty(seed)
    seed = 2024;
end

s = rng();
rng(seed, 'twister');

nTier1 = round(0.4 * numNodes);
nTier2 = round(0.4 * numNodes);
nTier3 = numNodes - nTier1 - nTier2;
if nTier3 < 0
    nTier3 = 0;
    nTier2 = numNodes - nTier1;
end
nodeTier = [ones(1,nTier1), 2*ones(1,nTier2), 3*ones(1,nTier3)];
nodeTier = nodeTier(randperm(numNodes));

tierCapacityRange     = [3e9, 5e9; 1.5e9, 3e9; 0.5e9, 1.5e9];
tierCostRange         = [0.06, 0.10; 0.03, 0.06; 0.01, 0.03];
tierPpeakRange        = [200, 280; 130, 200; 80, 130];
tierPidleRange        = [60, 85; 40, 60; 25, 40];
tierThermalGainRange  = [20, 30; 32, 44; 42, 55];

nodeCapacity = zeros(1, numNodes);
nodeCostRate = zeros(1, numNodes);
Ppeak        = zeros(1, numNodes);
Pidle        = zeros(1, numNodes);
eta          = zeros(1, numNodes);
for j = 1:numNodes
    tr = nodeTier(j);
    nodeCapacity(j) = tierCapacityRange(tr,1) + diff(tierCapacityRange(tr,:)) * rand();
    nodeCostRate(j) = tierCostRange(tr,1) + diff(tierCostRange(tr,:)) * rand();
    Ppeak(j)        = tierPpeakRange(tr,1) + diff(tierPpeakRange(tr,:)) * rand();
    Pidle(j)        = tierPidleRange(tr,1) + diff(tierPidleRange(tr,:)) * rand();
    eta(j)          = tierThermalGainRange(tr,1) + diff(tierThermalGainRange(tr,:)) * rand();
end

nodeBandwidth = 10 + 90 * rand(1, numNodes);
Tambient      = 28 + 12 * rand(1, numNodes);
Tmax          = 50 + 25 * rand(1, numNodes);

fogEnv.numNodes      = numNodes;
fogEnv.nodeTier      = nodeTier;
fogEnv.nodeCapacity  = nodeCapacity;
fogEnv.nodeBandwidth = nodeBandwidth;
fogEnv.nodeCostRate  = nodeCostRate;
fogEnv.Pidle         = Pidle;
fogEnv.Ppeak         = Ppeak;
fogEnv.Tambient      = Tambient;
fogEnv.Tmax          = Tmax;
fogEnv.eta           = eta;
fogEnv.theta0        = 10;

rng(s);
end
