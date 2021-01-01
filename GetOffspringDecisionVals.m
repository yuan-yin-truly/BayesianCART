function used_decision_vals = GetOffspringDecisionVals(tree, nodeID, decision_feature)
% Return vector of decision_vals in all offspring nodes of nodeID that has
% used decision_feature

used_decision_vals = [];
l = tree.lchild_ids(nodeID);

if isnan(l) 
    % nodeID is root
    return
end

r = tree.rchild_ids(nodeID);
queue = [l,r];


while ~isempty(queue)
    
    thisID = queue(1);
    queue = queue(2:end);
    
    l = tree.lchild_ids(thisID);
    r = tree.rchild_ids(thisID);
    
    % if thisID is not leaf:
    % add left and right child to queue,
    % add decision_vals
    if ~isnan(l) 
        queue = [queue, l, r];
    
        if tree.nodes(thisID).decision_feature == decision_feature
            used_decision_vals = [used_decision_vals, tree.nodes(thisID).decision_val];
        end
        
    end
    
end