function nodeSpace = NodeSpace(T, current_nodeID, nodeSpace, TF_ids)

% Return the space of each node as [TF1min, TF1max; TF2min, TF2max]
%   If decision node: the space where node.decision_val line to be drawn
%   If leaf node: leaf space
% (Only works for 2 TFs(2D))
%
% Input
%   nodeSpace: container.Map with key=decision_nodeID,
%       value=[TF1min,TF1max;TF2min,TF2max];
%       At initialization nodeSpace should already contain space of 
%       current node.
%   TF_ids: vector of TF ids. eg. in T the only 2 TFs involved may have id
%       of 4, and 5, then TF_ids should be [4, 5]
%
% Output
%   updated nodeSpace



% Find space for lchild and rchild of current_node
lchildID = T.lchild_ids(current_nodeID);
rchildID = T.rchild_ids(current_nodeID);

% Stop iteration at leaf node 
if ~isnan(lchildID)
    if T.nodes(current_nodeID).decision_feature == TF_ids(1) % TF1 is decision feature
        lspace = nodeSpace(current_nodeID);
        lspace(1,2) = T.nodes(current_nodeID).decision_val;
    elseif T.nodes(current_nodeID).decision_feature == TF_ids(2)
        lspace = nodeSpace(current_nodeID);
        lspace(2,2) = T.nodes(current_nodeID).decision_val;
    end
    
    nodeSpace(lchildID) = lspace;
    
    nodeSpace = NodeSpace(T, lchildID, nodeSpace, TF_ids);

end

% rchild's space should be min's modified
if ~isnan(rchildID)
    if T.nodes(current_nodeID).decision_feature == TF_ids(1) % TF1 is decision feature
        rspace = nodeSpace(current_nodeID);
        rspace(1,1) = T.nodes(current_nodeID).decision_val;
    elseif T.nodes(current_nodeID).decision_feature == TF_ids(2)
        rspace = nodeSpace(current_nodeID);
        rspace(2,1) = T.nodes(current_nodeID).decision_val;
    end
    
    nodeSpace(rchildID) = rspace;
    
    nodeSpace = NodeSpace(T, rchildID, nodeSpace, TF_ids);
end
        


end