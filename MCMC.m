function [chain, lnprob, acceptance_rate] = MCMC(T, mu, a, nu, lambda, alpha_var, beta_var, n, training_set, Y, iteration, BLR, weights)

% Parameters
% T: initial tree. Recommend to start with just one root node
% mu: normal prior's mean of leaf means' 
% a: normal prior's variance's scaling factor
% nu: inverse-gamma prior's hyperparam of leaf variance
% lambda: inverse-gamma prior's hyperparam of leaf variance
% alpha: tree prior's hyperparam
% beta: tree prior's hyperparam
% n: min num of samples in each leaf
% training_set: expression level of TFs
% Y: expression level of target gene
% iteration: number of iterations for MCMC
% BLR: if true, use bayesian linear regression to evaluate marginal L
% weights: probability of sampling each TF. Use 1 for uniform sampling

if BLR
    mu_b = mu; % mu_b is prior's mean(vector). If scalar, scalar*one_vector
    V_b = a; % V_b is prior's covariance up to sigma^2. If scalar, scalar*I
    a = nu; % inverse-gamma prior's first param
    b = lambda; % inverse-gamma prior's second param
end

if weights == 1 % Uniform sampling
    weights = 1/size(training_set,2)*ones(1,size(training_set,2));
end


chain = [];
lnprob = [];
acceptance_rate = [];

naccept = 0;
for i = 1:iteration
    % propose an operation on old_tree
    % 1: grow
    % 2: prune
    % 3: change

    operation = randperm(3,1);
    
    if length(T.lchild_ids) == 1 % if only a root, only allow grow
        operation = 1;
    end

    if operation == 1
        % grow
        
        % find growable nodes:
        % nodes with training_set size >= n
        growableIDs = GetGrowableIDs(T, n);
        growID = datasample(growableIDs, 1);
        
        
        % sample decision_feature with weights
        decision_feature = datasample(1:size(training_set,2),1,'Weights',weights);
        
        % sample decision_val: should not have appeared in ancestor nodes
        used_decision_vals = GetAncestralDecisionVals(T, growID, decision_feature);
        
%         % also should not be in top/lower 1% under that
%         % deicion_feature
%         used_decision_vals_2 = GetTopBotVals(training_set(:,decision_feature), 1);

        unique_decision_vals = unique(...
            setdiff(training_set(T.nodes(growID).sample_ids,decision_feature),used_decision_vals));
        
        % further screen down to those that do not produce empty leaf tree
        usable_decision_vals = [];
        for j = 1:length(unique_decision_vals)
            
            % to test grow, we have to make new copy of temp tree each time
            T_temp = T.copy();
            
            unique_decision_val = unique_decision_vals(j);
            T_temp.grow(growID, decision_feature, unique_decision_val, training_set);
            
            % if both l and r child growID HAVE >= n samples, then growID is
            % legal
            if length(T_temp.nodes(T_temp.lchild_ids(growID)).sample_ids)>=n &&...
                    length(T_temp.nodes(T_temp.rchild_ids(growID)).sample_ids)>=n
                
                usable_decision_vals = [usable_decision_vals, unique_decision_val];
                
            end
        end
        
        if isempty(usable_decision_vals)
            continue
        end
        
        decision_val = datasample(usable_decision_vals, 1);
        
        
        % propose new tree
        T_star = T.copy();
        T_star.grow(growID, decision_feature, decision_val, training_set);
        
        
        % tree prior ratio
        growID_depth = NodeDepth(T, growID);
        tree_ratio = ...
            -log(1-alpha_var*(1+growID_depth)^-beta_var)...
            +log(alpha_var*(1+growID_depth)^-beta_var)...
            +log(1-alpha_var*(2+growID_depth)^-beta_var)...
            +log(1-alpha_var*(2+growID_depth)^-beta_var);
        
        
        % likelihood ratio
        if BLR
            L = LogLikelihood_linear(T, mu_b, V_b, a, b, training_set, Y);
            L_star = LogLikelihood_linear(T_star, mu_b, V_b, a, b, training_set, Y);
        else
            L = LogLikelihood(T, mu, a, nu, lambda, Y);
            L_star = LogLikelihood(T_star, mu, a, nu, lambda, Y);
        end
        likelihood_ratio = L_star - L;
        
        
%         % gamma function ratio
%         % after growing, there are one more node in T_star
%         gamma_star = sum(gammaln(nsample_star/2+nu/2))-length(nsample_star)*gammaln(nu/2);
%         gamma_old = sum(gammaln(nsample/2+nu/2))-length(nsample)*gammaln(nu/2);
%         gamma_ratio = gamma_star - gamma_old;
%         
%         
%         % in case of NaN (Inf-Inf): assume the gamma ratio is 1
%         if isnan(gamma_ratio)
%             gamma_ratio = 0;
%         end
        
        
        % transition ratio
        transition_ratio = ...
            -log(length(GetPrunableIDs(T_star)))... % grown tree to old tree
            ...
            +log(length(growableIDs))+log(sum(weights)/weights(decision_feature))...
            +log(length(usable_decision_vals)); % old tree to grown tree
        
        
        % posterior ratio
        posterior_ratio = tree_ratio + likelihood_ratio + transition_ratio;
        

    elseif operation == 2
        % prune
        
        % if T now is only a root, continue to next iteration
        if length(T.lchild_ids.keys) == 1
            continue
        end
        
        
        % find prunable nodes: both left and right child are leaves
        prunableIDs = GetPrunableIDs(T);
        pruneID = datasample(prunableIDs, 1);

        
        % propose new tree
        T_star = T.copy();
        T_star.prune(pruneID);
        
        
        % tree prior ratio
        pruneID_depth = NodeDepth(T, pruneID);
        tree_ratio = ...
            log(1-alpha_var*(1+pruneID_depth)^-beta_var)...
            -log(alpha_var*(1+pruneID_depth)^-beta_var)...
            -log(1-alpha_var*(2+pruneID_depth)^-beta_var)...
            -log(1-alpha_var*(2+pruneID_depth)^-beta_var);

        
        % likelihood ratio
        if BLR
            L = LogLikelihood_linear(T, mu_b, V_b, a, b, training_set, Y);
            L_star = LogLikelihood_linear(T_star, mu_b, V_b, a, b, training_set, Y);
        else
            L = LogLikelihood(T, mu, a, nu, lambda, Y);
            L_star = LogLikelihood(T_star, mu, a, nu, lambda, Y);
        end
        
        likelihood_ratio = L_star - L;
        
        
        % transition ratio
        
        % if to grow back, what is the probability
        growableIDs = GetGrowableIDs(T_star, n);
        % to get the original tree's decision_feature
        decision_feature = T.nodes(pruneID).decision_feature;
        
        used_decision_vals = GetAncestralDecisionVals(T_star, pruneID, decision_feature);
                
        unique_decision_vals = unique(...
            setdiff(training_set(T_star.nodes(pruneID).sample_ids,decision_feature),used_decision_vals));
        
        
        % further screen down to those that do not produce empty leaf tree
        usable_decision_vals = [];
        for j = 1:length(unique_decision_vals)
            
            % to test grow, we have to make new copy of temp tree each time
            T_temp = T.copy();
            
            unique_decision_val = unique_decision_vals(j);
            T_temp.grow(pruneID, decision_feature, unique_decision_val, training_set);
            
            % if both l and r child growID HAVE >=n samples, then growID is
            % legal
            if length(T_temp.nodes(T_temp.lchild_ids(pruneID)).sample_ids)>=n &&...
                    length(T_temp.nodes(T_temp.rchild_ids(pruneID)).sample_ids)>=n
                usable_decision_vals = [usable_decision_vals, unique_decision_val];
            end
        end
        
        
        transition_ratio = ...
            -log(length(growableIDs))-log(sum(weights)/weights(decision_feature))...
            -log(length(usable_decision_vals))... % pruned tree to old tree
            ...
            +log(length(prunableIDs)); % old tree to pruned tree
        
        
        % posterior ratio
        posterior_ratio = tree_ratio + likelihood_ratio + transition_ratio;


    elseif operation == 3
        % change
        
        
        % find changable nodes: must be internal
        allID = cell2mat(T.lchild_ids.keys);
        internalIDs = allID(~isnan(cell2mat(T.lchild_ids.values)));
        
        % ID of node to change
        if isempty(internalIDs) && length(T.lchild_ids.keys) == 1
            changeID = cell2mat(T.lchild_ids.keys);
        else
            changeID = datasample(internalIDs, 1);
        end
        
        
        % sample decision_feature with weights
        decision_feature = datasample(1:size(training_set,2),1,'Weights',weights);
        
        
        % decision_val must have not been used in ancestral or children nodes
        used_decision_vals = GetUsedDecisionVals(T, changeID, decision_feature);
        
        
        unique_decision_vals = unique(...
            setdiff(training_set(T.nodes(changeID).sample_ids,decision_feature),used_decision_vals));
        
        % further screen down to those that do not produce >=n leaf sample tree
        usable_decision_vals = [];
        T_temp = T.copy();
        for j = 1:length(unique_decision_vals)
            unique_decision_val = unique_decision_vals(j);
            T_temp.change(changeID, decision_feature, unique_decision_val, training_set);
            if ~HasEmptyLeaf(T_temp,n)
                usable_decision_vals = [usable_decision_vals, unique_decision_val];
            end
        end
        
        if isempty(usable_decision_vals)
            continue
        end
        
        decision_val = datasample(usable_decision_vals, 1);
        
        
        % propose new tree
        T_star = T.copy();
        T_star.change(changeID, decision_feature, decision_val, training_set);
        
        
        % likelihood ratio
        if BLR
            L = LogLikelihood_linear(T, mu_b, V_b, a, b, training_set, Y);
            L_star = LogLikelihood_linear(T_star, mu_b, V_b, a, b, training_set, Y);
        else
            L = LogLikelihood(T, mu, a, nu, lambda, Y);
            L_star = LogLikelihood(T_star, mu, a, nu, lambda, Y);
        end
        
        likelihood_ratio = L_star - L;
        
        
        % transition ratio
        % consider new tree to old tree        
        decision_feature_old = T.nodes(changeID).decision_feature;
        
        used_decision_vals_star = GetUsedDecisionVals(T_star, changeID, decision_feature_old);
        
        unique_decision_vals_star = unique(...
            setdiff(training_set(T_star.nodes(changeID).sample_ids,decision_feature_old),used_decision_vals_star));
        
        % find total num of unique vals that will NOT give nsample < n leaves
        usable_decision_vals_star = [];
        T_temp = T_star.copy();
        for k = 1:length(unique_decision_vals_star)
            unique_decision_val_star = unique_decision_vals_star(k);
            T_temp.change(changeID, decision_feature_old, unique_decision_val_star, training_set);
            
            if ~HasEmptyLeaf(T_temp, n)
                usable_decision_vals_star = [usable_decision_vals_star, unique_decision_val_star];
            end
        end
        
        transition_ratio = ...
            -log(length(internalIDs))-log(sum(weights)/weights(decision_feature_old))...
            -log(length(usable_decision_vals_star))... % changed tree to old tree
            ...
            +log(length(internalIDs))+log(sum(weights)/weights(decision_feature))...
            +log(length(usable_decision_vals)); % old tree to changed tree

        


        % posterior ratio
        posterior_ratio = likelihood_ratio + transition_ratio;

    end
    
    % MCMC
    randnum = rand;
    accept = false;
    if randnum < min(1,exp(posterior_ratio))
        T = T_star.copy();
        accept = true;
        naccept = naccept + accept;
    end
    
    % update data
    chain = [chain, T];
    disp(strcat(num2str(i),'Current tree has',num2str(sum(isnan(cell2mat(T.lchild_ids.values)))),'leaves.'));
    
    acceptance_rate = [acceptance_rate, naccept/i];
    if accept
        lnprob = [lnprob,L_star+TreeProb(T,alpha_var,beta_var)];
        disp(strcat('Current tree log-posterior',num2str(L_star+TreeProb(T,alpha_var,beta_var))));
    
    else
        lnprob = [lnprob,L+TreeProb(T,alpha_var,beta_var)];
        disp(strcat('Current tree log-posterior',num2str(L+TreeProb(T,alpha_var,beta_var))));
        
    end

end


end