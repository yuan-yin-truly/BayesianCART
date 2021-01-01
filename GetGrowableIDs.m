function nodeIDs = GetGrowableIDs(tree, n)
% Return ID of nodes who have more than n samples

allIDs = cell2mat(tree.lchild_ids.keys);
leafIDs = allIDs(isnan(cell2mat(tree.lchild_ids.values)));

nodeIDs = [];

for i = 1:length(leafIDs)
    leafID = leafIDs(i);
    if length(tree.nodes(leafID).sample_ids) >= n
        nodeIDs = [nodeIDs, leafID];
    end
end

end