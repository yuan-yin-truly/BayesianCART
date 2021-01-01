function prob = TreeProb(tree, alpha_var, beta_var)
% return the log probability of a tree

% find root id
allIDs = cell2mat(tree.lchild_ids.keys);
rootID = allIDs(isnan(cell2mat(tree.parent_ids.values)));

prob = 0;

for i = 1:length(allIDs)
    nodeID = allIDs(i);
    if isnan(tree.lchild_ids(nodeID)) % no children
        prob = prob + log(1-alpha_var*(1+NodeDepth(tree, nodeID))^-beta_var);
    else
        prob = prob + log(alpha_var*(1+NodeDepth(tree, nodeID))^-beta_var);
    end


end