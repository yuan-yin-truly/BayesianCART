function nodeIDs = GetPrunableIDs(tree)
% Return the IDs of nodes whose left and right childs are both leaves

allIDs = cell2mat(tree.lchild_ids.keys);
leafIDs = allIDs(isnan(cell2mat(tree.lchild_ids.values)));
leafParentIDs = values(tree.parent_ids,num2cell(leafIDs));

% left children of these leafParentID
l = cell2mat(values(tree.lchild_ids, leafParentIDs));

% right children of these leafParentID
r = cell2mat(values(tree.rchild_ids, leafParentIDs));

% leafParentIDs whose l and r are both leafIDs

l_leaf = ismember(l,leafIDs);
r_leaf = ismember(r, leafIDs);
nodeIDs = cell2mat(leafParentIDs);
nodeIDs = unique(nodeIDs(l_leaf & r_leaf));


end