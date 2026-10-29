function order = topological_sort(adjMatrix, numTasks)
    indeg = sum(adjMatrix, 1);
    order = zeros(1, numTasks);
    pos = 1;
    queue = find(indeg == 0);
    indegRemaining = indeg;
    visited = false(1, numTasks);

    while ~isempty(queue)
        n = queue(1);
        queue(1) = [];
        if visited(n), continue; end
        visited(n) = true;
        order(pos) = n; pos = pos + 1;

        succ = find(adjMatrix(n,:));
        for s = succ
            indegRemaining(s) = indegRemaining(s) - 1;
            if indegRemaining(s) == 0 && ~visited(s)
                queue(end+1) = s;
            end
        end
    end

    if pos <= numTasks
        remaining = find(~visited);
        order(pos:pos+numel(remaining)-1) = remaining;
    end
end
