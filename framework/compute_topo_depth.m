function depth = compute_topo_depth(adjMatrix, numTasks)
    indeg = sum(adjMatrix, 1);
    depth = zeros(1, numTasks);
    order = topological_sort(adjMatrix, numTasks);
    for idx = 1:numel(order)
        i = order(idx);
        preds = find(adjMatrix(:,i))';
        if isempty(preds)
            depth(i) = 0;
        else
            depth(i) = max(depth(preds)) + 1;
        end
    end
end
