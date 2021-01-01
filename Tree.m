classdef Tree < matlab.mixin.Copyable
    % When initialize a tree, we initialize a binary tree consists of nodes
    
    properties
        % At leaf node i, lchild(i) and rchild(i) = NaN
        lchild_ids = containers.Map('KeyType','double', 'ValueType','double');
        rchild_ids = containers.Map('KeyType','double', 'ValueType','double');
        % Root node parent = NaN
        parent_ids = containers.Map('KeyType','double', 'ValueType','double');
        
        nodes = containers.Map('KeyType','double', 'ValueType','any');
        
    end
    
    methods (Access = protected)
        function objCopy = copyElement(obj)
        % copy this tree to make a new object (not referenced to originals)
        
            objCopy = copyElement@matlab.mixin.Copyable(obj);
            objCopy.lchild_ids = containers.Map(obj.lchild_ids.keys, obj.lchild_ids.values);
            objCopy.rchild_ids = containers.Map(obj.rchild_ids.keys, obj.rchild_ids.values);
            objCopy.parent_ids = containers.Map(obj.parent_ids.keys, obj.parent_ids.values);

            objCopy.nodes = containers.Map('KeyType','double', 'ValueType','any');
            for i = obj.nodes.keys
                objCopy.nodes(i{1}) = copy(obj.nodes(i{1}));
            end
        end
    end
    
    
    methods
        
        
        %---- Object Constructor ----%
        
        function obj = Tree(lchild_ids, rchild_ids, decision_features, decision_vals) % Inputs: container maps
           
            obj.lchild_ids = lchild_ids;
            obj.rchild_ids = rchild_ids;
            
            
            % Iterate all nodes ID and initialize nodes
            all_nodes_id = keys(obj.lchild_ids);
            
            for i = 1:length(all_nodes_id)
                thisNode = all_nodes_id{i};
                
                % Initialize nodes
                obj.nodes(thisNode) = Node(thisNode, decision_features(thisNode), decision_vals(thisNode));
                
                if isnan(obj.lchild_ids(thisNode))
                    continue
                end
                
                % Recover parent
                obj.parent_ids(obj.lchild_ids(thisNode)) = thisNode;
                obj.parent_ids(obj.rchild_ids(thisNode)) = thisNode;
                
                
            end
            
            % add root node's parent NaN
            withChildrenButNoParent = ~ismember(cell2mat(all_nodes_id), cell2mat(keys(obj.parent_ids)));
            root_id = all_nodes_id{withChildrenButNoParent};
            
            obj.parent_ids(root_id) = NaN;
            
        end
        
        
        
        function assign_data(obj, training_set)
        % training_set = nsample x nfeature
        % This function assign sample id to each nodes of the tree
            
            % Find root id
            all_ids = keys(obj.parent_ids);
            root_id = all_ids{isnan(cell2mat(values(obj.parent_ids)))};
            
            root_node = obj.nodes(root_id);
            root_node.sample_ids = 1:size(training_set,1);
            
            
            obj.assign_data_helper(root_id, training_set, root_node.sample_ids);
        
        end
        
        function assign_data_helper(obj, current_id, training_set, current_sample_ids)
        % Recursively assign data index to nodes
            if isnan(obj.lchild_ids(current_id))
                return
            else
                
                all_ids = 1:size(training_set,1);
                
                feature = obj.nodes(current_id).decision_feature;
                val = obj.nodes(current_id).decision_val;
                % assign left child
                lid = obj.lchild_ids(current_id);
                
                lnode = obj.nodes(lid);
                lnode.sample_ids = intersect(...
                    all_ids(training_set(:,feature) < val),...
                    current_sample_ids);
                
                
                % assign right child
                rid = obj.rchild_ids(current_id);
                
                rnode = obj.nodes(rid);
                rnode.sample_ids = intersect(...
                    all_ids(training_set(:,feature) >= val),...
                    current_sample_ids);
                
                
                % recursive call
                assign_data_helper(obj, lid, training_set, lnode.sample_ids);
                assign_data_helper(obj, rid, training_set, rnode.sample_ids);
                
            end
        end
                
        
        %---- Prediction ----%
        
        function update_leaf_value(obj, y)
        % Convert the sample_id at leaves to mean and variance and assign
        % to leaf nodes' node_mean and node_var.
        
        % Input
        % y: target genes expression: nsamples x ngenes
        % Hence nrow(y) should = nrow(training_set)
            
            % If y is a vector, ensure y is a col vector
            if size(y,1) == 1
                y = y(:);
            end
            
            allIDs = cell2mat(obj.lchild_ids.keys);
            leafIDs = allIDs(isnan(cell2mat(obj.lchild_ids.values)));
            
            for i = 1:length(leafIDs)
                leafID = leafIDs(i);
                thisNode = obj.nodes(leafID);
                thisSet = y(thisNode.sample_ids,:);
                
                % convert thisSet to vector for easy computation of
                % sufficient stats
                thisSet = thisSet(:);
                
                thisNode.node_mean = mean(thisSet);
                thisNode.node_var = var(thisSet);
            end
        
        end
        
        function update_leaf_value_linear(obj, training_set, y)
        % Find linear regression coefficients at leaf nodes
        % Regress gene level on level of TFs of ancestors of leaf
        
        % Input
        % training_set: all TF_level training set: nsamples x ngenes
        % y: target genes expression: nsamples x ngenes

            
            allIDs = cell2mat(obj.lchild_ids.keys);
            leafIDs = allIDs(isnan(cell2mat(obj.lchild_ids.values)));
            
            for i = 1:length(leafIDs)
                leafID = leafIDs(i);
                thisNode = obj.nodes(leafID);
                
                % TFs that are ancestors to this leaf
                features = obj.get_ancestral_decision_features(leafID);
                if isempty(features) % root is the only node of tree
                    features = 1:size(training_set,2);
                end
                features = sort(unique(features));
                thisNode.ancestral_decision_features = features;
                
                thisX = training_set(thisNode.sample_ids, features);
                thisy = y(thisNode.sample_ids,:);
                
                % stack up all col
                thisy = thisy(:);
                thisX = [ones(length(thisy),1),repmat(thisX,size(y,2),1)];
                
                thisNode.linear_reg_coef = thisX\thisy;
            end
        end
        
        
        function features = get_ancestral_decision_features(obj, nodeID)
        % return decision_features of all ancestral nodes of nodeID
        
            features = [];
            while true
                if isnan(obj.parent_ids(nodeID)) 
                    % nodeID is root
                    return
                end

                parentNode = obj.nodes(obj.parent_ids(nodeID)); 
                features = [features, parentNode.decision_feature];
                
                % update nodeID to parent of input node
                nodeID = obj.parent_ids(nodeID);
            end
        end
        
        
        function predicted_vals = predict(obj,test_set)
        % return node_mean of leaf where samples in test_set lands
            
            % Find root id
            allIDs = keys(obj.parent_ids);
            rootID = allIDs{isnan(cell2mat(values(obj.parent_ids)))};
            
            predicted_vals = [];
            
            for i = 1:size(test_set,1)
                test_set_entry = test_set(i,:);
                leafID = obj.predict_val_helper(test_set_entry, rootID);
                predicted_vals = [predicted_vals; obj.nodes(leafID).node_mean];
            end
            
        end
        
        function predicted_vals = predict_linear(obj,test_set)
        % return linear regression at leaf where samples in test_set lands
            
            % Find root id
            allIDs = keys(obj.parent_ids);
            rootID = allIDs{isnan(cell2mat(values(obj.parent_ids)))};
            
            predicted_vals = [];
            
            for i = 1:size(test_set,1)
                test_set_entry = test_set(i,:);
                leafID = obj.predict_val_helper(test_set_entry, rootID);
                predicted_vals = [predicted_vals;...
                    dot([1,test_set_entry(obj.nodes(leafID).ancestral_decision_features)],...
                    obj.nodes(leafID).linear_reg_coef)];
            end
            
        end
        
        function nodeID = predict_val_helper(obj, test_set_entry, startingID)
        % recursively find the leaf for test_set_entry
        
            if isnan(obj.lchild_ids(startingID))
                nodeID = startingID;
                return
            end
            
            dfture = obj.nodes(startingID).decision_feature;
            dval = obj.nodes(startingID).decision_val;
        
            if test_set_entry(dfture) >= dval
                startingID = obj.rchild_ids(startingID);
            else
                startingID = obj.lchild_ids(startingID);
            end
            
            nodeID = predict_val_helper(obj,test_set_entry, startingID);
        
        end
        
        
        %---- tree manipulations ----%
        
        function grow(obj, node_id, decision_feature, decision_val, training_set)
        % split at node_id (must be a leaf)
            
            max_id = max(cell2mat(keys(obj.lchild_ids)));
            
            obj.lchild_ids(node_id) = max_id + 1;
            obj.rchild_ids(node_id) = max_id + 2;
            
            obj.lchild_ids(max_id + 1) = NaN;
            obj.lchild_ids(max_id + 2) = NaN;
            
            obj.rchild_ids(max_id + 1) = NaN;
            obj.rchild_ids(max_id + 2) = NaN;
            
            obj.parent_ids(max_id + 1) = node_id;
            obj.parent_ids(max_id + 2) = node_id;
            
            % Update decision at node_id
            thisNode = obj.nodes(node_id);
            thisNode.decision_feature = decision_feature;
            thisNode.decision_val = decision_val;
            
            % Create nodes and ...
            % update sample_ids at 2 new nodes
            
            obj.nodes(max_id + 1) = Node(max_id + 1);
            obj.nodes(max_id + 2) = Node(max_id + 2);
            
            obj.assign_data_helper(node_id, training_set, thisNode.sample_ids)

        end
        
        
        function prune(obj, node_id)
        % Remove two leaves under node_id
            
            % Leaf id to remove:
            rm_id = {obj.lchild_ids(node_id),obj.rchild_ids(node_id)};
            
            remove(obj.lchild_ids, rm_id);
            remove(obj.rchild_ids, rm_id);
            remove(obj.parent_ids, rm_id);
            
            % Remove nodes
            remove(obj.nodes, rm_id);
            
            % Remove link to the nodes
            obj.lchild_ids(node_id) = NaN;
            obj.rchild_ids(node_id) = NaN;
            
            % Remove decision feature and val at this node
            thisNode = obj.nodes(node_id);
            thisNode.decision_feature = NaN;
            thisNode.decision_val = NaN;
        end
        
        
        function change(obj, node_id, decision_feature, decision_val, training_set)
        % reassign splitting rule of node_id (must be internal)
            
            thisNode = obj.nodes(node_id);
            thisNode.decision_feature = decision_feature;
            thisNode.decision_val = decision_val;
            
            % Update sample_ids at downstream nodes
            obj.assign_data_helper(node_id, training_set, thisNode.sample_ids)
        
        end
        
        
        function swap(obj,node_id, training_set)
        % swap the decision rule between node_id and a child (both must be internal)
        % if both children have same decision rule, swap both with parent
        
            % both children are leaves
            if isnan(obj.lchild_ids(obj.lchild_ids(node_id))) &&...
                    isnan(obj.rchild_ids(obj.rchild_ids(node_id)))
                return
            
                % lchild is leaf
            elseif isnan(obj.lchild_ids(obj.lchild_ids(node_id))) &&...
                    ~isnan(obj.rchild_ids(obj.rchild_ids(node_id)))
                
                obj.swap_helper(node_id, obj.rchild_ids(node_id));
                
                                
                %rchild is leaf
            elseif ~isnan(obj.lchild_ids(obj.lchild_ids(node_id))) &&...
                    isnan(obj.rchild_ids(obj.rchild_ids(node_id)))
                
                obj.swap_helper(node_id, obj.lchild_ids(node_id));
                
                
            else
                
                % if both children have same decision rule
                if (obj.nodes(obj.lchild_ids(node_id)).decision_feature == ...
                        obj.nodes(obj.rchild_ids(node_id)).decision_feature) &&...
                        (obj.nodes(obj.lchild_ids(node_id)).decision_val == ...
                            obj.nodes(obj.rchild_ids(node_id)).decision_val)
                    
                    % swap parent and right child
                    obj.swap_helper(node_id, obj.rchild_ids(node_id));
                    
                    
                    % assign right child's decision to left child
                    lnode = obj.nodes(obj.lchild_ids(node_id));
                    lnode.decision_feature = ...
                        obj.nodes(obj.rchild_ids(node_id)).decision_feature;
                    lnode.decision_val = ...
                        obj.nodes(obj.rchild_ids(node_id)).decision_val;
                
                % if different dicision rule, randomly choose one to swap
                % with parent
                else
                    if rand(1) > 0.5
                        obj.swap_helper(node_id, obj.rchild_ids(node_id));
                        
                    else
                        obj.swap_helper(node_id, obj.lchild_ids(node_id));
                        
                    end
                end
            end
            
            % Update sample_ids from this node downwards
            obj.assign_data_helper(node_id, training_set, obj.nodes(node_id).sample_ids);
            
        end
        
        function swap_helper(obj, node1_id, node2_id)
        % swap the decision_feature and decision_val of node1 and node2
            
            parent_feature = obj.nodes(node1_id).decision_feature;
            parent_val = obj.nodes(node1_id).decision_val;

            node1 = obj.nodes(node1_id);
            node1.decision_feature = obj.nodes(node2_id).decision_feature;
            node1.decision_val = obj.nodes(node2_id).decision_val;

            node2 = obj.nodes(node2_id);
            node2.decision_feature = parent_feature;
            node2.decision_val = parent_val;
        
        end
        
        
    end
    
    
end