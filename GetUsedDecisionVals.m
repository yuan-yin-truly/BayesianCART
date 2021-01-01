function used_decision_vals = GetUsedDecisionVals(tree, nodeID, decision_feature)
% Return a vector of decision_values that are used by nodeID's ancestral
% and offspring nodes that used decision_feature

used_val_ances = GetAncestralDecisionVals(tree, nodeID, decision_feature);
used_val_offsp = GetOffspringDecisionVals(tree, nodeID, decision_feature);

used_decision_vals = [used_val_ances, used_val_offsp];
end
        