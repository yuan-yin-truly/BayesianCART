function used_decision_vals = GetAncestralDecisionVals(tree, nodeID, decision_feature)
% Return a vector of decision_values that are used by nodeID's ancestral
% nodes that used decision_feature

used_decision_vals = [];

while true
    if isnan(tree.parent_ids(nodeID)) 
        % nodeID is root
        return
    end
    
    parentNode = tree.nodes(tree.parent_ids(nodeID)); 
    
    if parentNode.decision_feature == decision_feature
        used_decision_vals = [used_decision_vals, parentNode.decision_val];
    end
    
    % update nodeID to parent of input node
    nodeID = tree.parent_ids(nodeID);

end