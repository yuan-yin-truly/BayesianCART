function col_ranges = PlotTreeStructure_helper(tree)

% To plot the tree, each row is a depth in the tree, identify col_ranges of each
% node
%
% Rule:
% If a node is lchild, then range is 1st half of parent_range
% If a node is rchild, then range is 2nd half of parent_range
% If the node is root, then range is all (1:nrow)
%
% Output:
% col_ranges: container map with input nodeID output vector of cols


allIDs = cell2mat(tree.parent_ids.keys);
leafIDs = allIDs(isnan(cell2mat(tree.lchild_ids.values)));

% Find depth of tree as number of rows
node_depths = [];
for i = 1:length(leafIDs)
    node_depths = [node_depths, NodeDepth(tree, leafIDs(i))];
end

ncol = 2^max(node_depths);



for i = 1:length(allIDs)
    node_id = allIDs(i);
    
    if isnan(tree.parent_ids(node_id)) % root
        col_ranges = containers.Map(node_id,1:ncol);

    elseif ismember(node_id, cell2mat(tree.lchild_ids.values)) % this node lchild
        parent_range = col_ranges(tree.parent_ids(node_id));
        
        col_ranges(node_id) = parent_range(1):...
            (parent_range(1)+0.5*(parent_range(end)-parent_range(1)-1));
        
    elseif ismember(node_id, cell2mat(tree.rchild_ids.values)) % this node rchild
        parent_range = col_ranges(tree.parent_ids(node_id));
        
        col_ranges(node_id) = ...
            (parent_range(1)+0.5*(parent_range(end)-parent_range(1)+1)):...
            parent_range(end);
    end
end


end