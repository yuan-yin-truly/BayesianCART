function PlotTreeScatter(T,TF_level, gene_level, training_idx, TF_name, gene_name, Error, alpha_var, beta_var, Linear, bootstrp)

% Plot data as scatter; tree as boundary lines & regression lines/planes
% Work for up to 2 TFs


% Input
% T: ensure leaf nodes have node_mean and node_var filled after MCMC search. 
%   Call Tree.update_leaf_val with input genes training set.

% TF_level: all TF expression vector(matrix): nsamples x 1(2)
% gene_level: all gene expression matrix: nsamples x ngenes
% training_idx: index to find TF and gene's training set. Only used when
%   bootstrap is required.

% TF_name: cell array length 1(2). To be displayed on plot title.
% gene_name: cell array length >= 1. To be displayed on plot title.

% Error: Optional. CV error for display only.
% alpha_var: Optional. alpha used to find the tree. For display only.
% beta_var: Optional. beta used to find the tree. For display only.

% Linear: bool. Use false if unsure. If true, instead of  plotting node_mean
%   at each leaf, plot the linear regression line/plane at the leaf.

% bootstrp: int. 0 if unsure. 100 recommended. If 0, no bootstrap/CI; 
%   If 100, do 100 resampling and plot CI.

    

    % number of TFs in the tree
    allID = cell2mat(T.lchild_ids.keys);
    decisionIDs = allID(~isnan(cell2mat(T.lchild_ids.values)));
    
    TF_ids = [];
    for i = 1:length(decisionIDs)
        decisionID = decisionIDs(i);
        TF_ids = [TF_ids, T.nodes(decisionID).decision_feature];
    end
    
    TF_ids = unique(TF_ids);


    figure;
 
    % scatter plot each gene vs. TFs
    % only 1 TF
    if length(TF_ids) == 1
        for i = 1:size(gene_level,2)
            % s as handle to scatter plots only (exclude line plots below)
            % for legends of genes only
            s(i) = scatter(TF_level(:,TF_ids), gene_level(:,i),'.');
            hold on
        end
        
        xlabel(strjoin(TF_name(TF_ids)))
        ylabel('genes')
        
    % 2 TFs
    elseif length(TF_ids) == 2
        for i = 1:size(gene_level,2)
            s(i) = scatter3(TF_level(:,TF_ids(1)), TF_level(:,TF_ids(2)), gene_level(:,i),'.');
            hold on
        end
        
        xlabel(TF_name{TF_ids(1)})
        ylabel(TF_name{TF_ids(2)})
        zlabel('genes')
        
    end
    title(strcat(['genes vs. ',strjoin(TF_name(TF_ids),', '),'  Error_{CV}=',num2str(Error),' \alpha=',num2str(alpha_var),' \beta=',num2str(beta_var)]));
    

    hold on

    
    
    
    % Plot decision_vals and regression line
    
    
    
    % 1 TF
    if length(TF_ids) == 1

        decision_vals = [];

        for i = 1:length(decisionIDs)
            decisionID = decisionIDs(i);
            decision_vals=[decision_vals,T.nodes(decisionID).decision_val];
        end

        decision_vals = sort(decision_vals);
        decision_vals = [min(TF_level),decision_vals]; % add in min of TF_level_log


        if Linear  % plot linear regression

            for i = 2:length(decision_vals)
                % vertical line (decision_val)
                decision_val = decision_vals(i);
                line([decision_val,decision_val],[min(gene_level(:)),max(gene_level(:))],...
                    'Color','red','LineWidth',3);
                hold on

                % regression line
                line([decision_vals(i-1),decision_val],...
                    [T.predict_linear(decision_vals(i-1)+0.01),T.predict_linear(decision_val-0.01)],...
                    'Color','green','LineWidth',2);
                hold on
            end

            % regression line on right of right-most decision_val
            line([decision_vals(end),max(TF_level)],...
                [T.predict_linear(decision_vals(end)+0.01), T.predict_linear(max(TF_level))],...
                    'Color','green','LineWidth',2);

        else  % not linear regression: plot node_mean

            for i = 2:length(decision_vals)
                % vertical line
                decision_val = decision_vals(i);
                line([decision_val,decision_val],[min(gene_level(:)),max(gene_level(:))],...
                    'Color','red','LineWidth',3);
                hold on

                % horizontal line
                node_mean = T.predict(mean([decision_vals(i-1),decision_val]));
                line([decision_vals(i-1),decision_val],[node_mean, node_mean],...
                    'Color','green','LineWidth',2);
                hold on
            end

            % horizonal line on right of right-most decision val(vertical line)
            node_mean = T.predict(max(TF_level));
            line([decision_vals(end),max(TF_level)],[node_mean, node_mean],...
                    'Color','green','LineWidth',2);
        end
        
        
    % 2TFs
    elseif length(TF_ids) == 2
        
        % Plot decision_vals as lines on bottom plane
        
        rootID = allID(isnan(cell2mat(T.parent_ids.values)));
        % Get X, Y coordinates where each node govern/can plot
        nodeSpaces = containers.Map(rootID,[get(gca,'XLim');get(gca,'YLim')]);
        nodeSpaces = NodeSpace(T, rootID, nodeSpaces, TF_ids);
        
        for i = 1:length(decisionIDs)
            decisionID = decisionIDs(i);
            nodeSpace = nodeSpaces(decisionID);
            
            if T.nodes(decisionID).decision_feature == TF_ids(1)
                xval = T.nodes(decisionID).decision_val;
                zvals = get(gca,'ZLim');
                plot3([xval,xval], nodeSpace(2,:),[zvals(1),zvals(1)],'k');
                hold on
            elseif T.nodes(decisionID).decision_feature == TF_ids(2)
                yval = T.nodes(decisionID).decision_val;
                zvals = get(gca,'ZLim');
                plot3(nodeSpace(1,:),[yval,yval],[zvals(1),zvals(1)],'k');
                hold on
            end
        end
        
        
        % Plot regression plane in each leaf
        if Linear
            
            leafIDs = allID(isnan(cell2mat(T.lchild_ids.values)));
            for i = 1:length(leafIDs)
                
                leafID = leafIDs(i);
                nodeSpace = nodeSpaces(leafID);

                [X,Y] = meshgrid(...
                    linspace(nodeSpace(1,1),nodeSpace(1,2),20),...
                    linspace(nodeSpace(2,1),nodeSpace(2,2),20));
                
                if T.nodes(leafID).ancestral_decision_features==TF_ids(1)
                    Z = T.nodes(leafID).linear_reg_coef(1) +...
                        T.nodes(leafID).linear_reg_coef(2)*X;
                elseif T.nodes(leafID).ancestral_decision_features==TF_ids(2)
                    Z = T.nodes(leafID).linear_reg_coef(1) +...
                        T.nodes(leafID).linear_reg_coef(2)*Y;
                else
                    Z = T.nodes(leafID).linear_reg_coef(1) +...
                        T.nodes(leafID).linear_reg_coef(2)*X +...
                        T.nodes(leafID).linear_reg_coef(3)*Y;
                end
                
                if ~bootstrp % Don't plot CI
                    surf(X,Y,Z,'FaceAlpha',0.5,'EdgeColor','none');
                    hold on
            
                else % Plot CI from bootstraping
                    
                    % bootstrap resample is from training_set
                    TF_training = TF_level(training_idx,:);
                    gene_training = gene_level(training_idx,:);
                     
                    % Resample at this leaf
                    bstrp_param = [];
                    for j = 1:bootstrp
                        
                        resample_id = datasample(T.nodes(leafID).sample_ids,...
                            length(T.nodes(leafID).sample_ids));
                        x = TF_training(resample_id,...
                            T.nodes(leafID).ancestral_decision_features);
                        y = gene_training(resample_id,:);
                        
                        thisy = y(:);
                        thisx = [ones(length(thisy),1),repmat(x,size(y,2),1)];
                        
                        % regression param from bootstrap sample
                        bstrp_param = [bstrp_param,thisx\thisy];
                    end
                    
                    
                    % calc and find 95% CI in plotting space
                    X_vec = X(:);
                    Y_vec = Y(:);
                    Z_upper = [];
                    Z_lower = [];
                    for k = 1:length(X_vec) % at each coordinates
                        if T.nodes(leafID).ancestral_decision_features==TF_ids(1)
                            Z_range = 1*bstrp_param(1,:)+...
                                X_vec(k)*bstrp_param(2,:);
                        elseif T.nodes(leafID).ancestral_decision_features==TF_ids(2)
                            Z_range = 1*bstrp_param(1,:)+...
                                Y_vec(k)*bstrp_param(2,:);
                        else
                            Z_range = 1*bstrp_param(1,:)+...
                                X_vec(k)*bstrp_param(2,:)+...
                                Y_vec(k)*bstrp_param(3,:);
                        end
                        Z_upper = [Z_upper;prctile(Z_range,97.5)];
                        Z_lower = [Z_lower;prctile(Z_range,2.5)];
                    end
                    
                    Z_upper = reshape(Z_upper,size(X));
                    Z_lower = reshape(Z_lower,size(X));
                    surf(X,Y,Z_upper,'FaceAlpha',0.5,'EdgeColor','none');
                    hold on
                    surf(X,Y,Z_lower,'FaceAlpha',0.5,'EdgeColor','none');
                    hold on
                    
                    % Plot vertical planes to join CI planes
                    patch([X(1,:),flip(X(1,:))],...
                        [Y(1,:),flip(Y(1,:))],...
                        [Z_upper(1,:),flip(Z_lower(1,:))],...
                        [Z_upper(1,:),flip(Z_lower(1,:))],... %colormap val
                        'FaceAlpha',0.5);
                    hold on
                    patch([X(:,1)',flip(X(:,1)')],...
                        [Y(:,1)',flip(Y(:,1)')],...
                        [Z_upper(:,1)',flip(Z_lower(:,1)')],...
                        [Z_upper(:,1)',flip(Z_lower(:,1)')],...
                        'FaceAlpha',0.5);
                    hold on
                    patch([X(end,:),flip(X(end,:))],...
                        [Y(end,:),flip(Y(end,:))],...
                        [Z_upper(end,:),flip(Z_lower(end,:))],...
                        [Z_upper(end,:),flip(Z_lower(end,:))],...
                        'FaceAlpha',0.5);
                    hold on
                    patch([X(:,end)',flip(X(:,end)')],...
                        [Y(:,end)',flip(Y(:,end)')],...
                        [Z_upper(:,end)',flip(Z_lower(:,end)')],...
                        [Z_upper(:,end)',flip(Z_lower(:,end)')],...
                        'FaceAlpha',0.5);
                    hold on
                    
                    
                    % Plot intersection of regression and vertical planes
                    plot3(X(1,:),Y(1,:),Z(1,:),'LineWidth',2,'Color','k');
                    plot3(X(end,:),Y(end,:),Z(end,:),'LineWidth',2,'Color','k');
                    plot3(X(:,end),Y(:,end),Z(:,end),'LineWidth',2,'Color','k');
                    plot3(X(:,1),Y(:,1),Z(:,1),'LineWidth',2,'Color','k');
                    
                    
                end
            end
        
        else % % not linear regression: plot node_mean
            
            leafIDs = allID(isnan(cell2mat(T.lchild_ids.values)));
            for i = 1:length(leafIDs)
                
                leafID = leafIDs(i);
                nodeSpace = nodeSpaces(leafID);
                
                node_mean = T.nodes(leafID).node_mean;

                if ~bootstrp % Don't plot CI
                    
                    patch([nodeSpace(1,:),flip(nodeSpace(1,:))],...
                        [repmat(nodeSpace(2,1),1,2),repmat(nodeSpace(2,2),1,2)],...
                        repmat(node_mean,1,4),...
                        repmat(node_mean,1,4),... %colormap val
                        'FaceAlpha',0.5);
                    hold on
            
                else % Plot CI from bootstraping
                    
                    % bootstrap resample is from training_set
                    TF_training = TF_level(training_idx,:);
                    gene_training = gene_level(training_idx,:);
                     
                    % Resample at this leaf
                    bstrp_param = [];
                    for j = 1:bootstrp
                        
                        resample_id = datasample(T.nodes(leafID).sample_ids,...
                            length(T.nodes(leafID).sample_ids));
                        y = gene_training(resample_id,:);
                        thisy = y(:);
                        
                        % param from bootstrap sample
                        bstrp_param = [bstrp_param,mean(thisy)];
                    end
                    
                    
                    % calc and find 95% CI in plotting space
                    patch([nodeSpace(1,:),flip(nodeSpace(1,:))],...
                        [repmat(nodeSpace(2,1),1,2),repmat(nodeSpace(2,2),1,2)],...
                        repmat(prctile(bstrp_param,97.5),1,4),...
                        repmat(prctile(bstrp_param,97.5),1,4),... %colormap val
                        'FaceAlpha',0.5);
                    hold on
                    patch([nodeSpace(1,:),flip(nodeSpace(1,:))],...
                        [repmat(nodeSpace(2,1),1,2),repmat(nodeSpace(2,2),1,2)],...
                        repmat(prctile(bstrp_param,2.5),1,4),...
                        repmat(prctile(bstrp_param,2.5),1,4),... %colormap val
                        'FaceAlpha',0.5);
                    hold on
                    
                    % Plot vertical planes to join CI planes
                    patch([nodeSpace(1,:),flip(nodeSpace(1,:))],...
                        repmat(nodeSpace(2,1),1,4),...
                        [repmat(prctile(bstrp_param,2.5),1,2),repmat(prctile(bstrp_param,97.5),1,2)],...
                        [repmat(prctile(bstrp_param,2.5),1,2),repmat(prctile(bstrp_param,97.5),1,2)],... %colormap val
                        'FaceAlpha',0.5);
                    hold on
                    patch([nodeSpace(1,:),flip(nodeSpace(1,:))],...
                        repmat(nodeSpace(2,2),1,4),...
                        [repmat(prctile(bstrp_param,2.5),1,2),repmat(prctile(bstrp_param,97.5),1,2)],...
                        [repmat(prctile(bstrp_param,2.5),1,2),repmat(prctile(bstrp_param,97.5),1,2)],... %colormap val
                        'FaceAlpha',0.5);
                    hold on
                    patch(repmat(nodeSpace(1,1),1,4),...
                        [nodeSpace(2,:),flip(nodeSpace(2,:))],...
                        [repmat(prctile(bstrp_param,2.5),1,2),repmat(prctile(bstrp_param,97.5),1,2)],...
                        [repmat(prctile(bstrp_param,2.5),1,2),repmat(prctile(bstrp_param,97.5),1,2)],... %colormap val
                        'FaceAlpha',0.5);
                    hold on
                    patch(repmat(nodeSpace(1,2),1,4),...
                        [nodeSpace(2,:),flip(nodeSpace(2,:))],...
                        [repmat(prctile(bstrp_param,2.5),1,2),repmat(prctile(bstrp_param,97.5),1,2)],...
                        [repmat(prctile(bstrp_param,2.5),1,2),repmat(prctile(bstrp_param,97.5),1,2)],... %colormap val
                        'FaceAlpha',0.5);
                    hold on
                    
                    
                    % Plot intersection of regression and vertical planes
                    plot3(nodeSpace(1,:),repmat(nodeSpace(2,1),1,2),repmat(node_mean,1,2),'LineWidth',2,'Color','k');
                    plot3(nodeSpace(1,:),repmat(nodeSpace(2,2),1,2),repmat(node_mean,1,2),'LineWidth',2,'Color','k');
                    plot3(repmat(nodeSpace(1,1),1,2),nodeSpace(2,:),repmat(node_mean,1,2),'LineWidth',2,'Color','k');
                    plot3(repmat(nodeSpace(1,2),1,2),nodeSpace(2,:),repmat(node_mean,1,2),'LineWidth',2,'Color','k');
                    
                end
            end
            
        end
    end
        
    % add legend only to the scatter plot
    legend(s(1:size(gene_level,2)),gene_name);
    hold off
    
end