function emptyLeaf = HasEmptyLeaf(tree, n)
% return true if the tree has leaf with <n sample_ids

allIDs = cell2mat(tree.lchild_ids.keys);
leafIDs = allIDs(isnan(cell2mat(tree.lchild_ids.values)));

emptyLeaf = false;

for i = 1:length(leafIDs)
    leafID = leafIDs(i);
    if length(tree.nodes(leafID).sample_ids) < n
        emptyLeaf = true;
        break
    end
end

end