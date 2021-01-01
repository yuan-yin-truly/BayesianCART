function depth = NodeDepth(tree, nodeID)
% Return depth of nodeID (count number of parents)
% Depth of root is 0

allID = cell2mat(tree.parent_ids.keys);
root_id = allID(isnan(cell2mat(tree.parent_ids.values)));

if nodeID == root_id
    depth = 0;
    return
end

depth = 1+ NodeDepth(tree, tree.parent_ids(nodeID));

end