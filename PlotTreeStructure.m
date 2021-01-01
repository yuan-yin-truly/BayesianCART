function PlotTreeStructure(tree, TF_level, gene_level, training_idx, TF_name, gene_name, median_factor, plotSlope)

% Plot decision_nodes of tree where at each nodes is a violin plot of the gene.

% Input 
%   tree: Tree object
%   TF_level: all TF expression matrix: nsamples x nTFs
%   TF_name: cell array. TF names
%   gene_level: all gene expression matrix: nsamples x ngene
%   gene_name: cell array. gene names
%   median_factor: ratio of first gene's median over later gene's median in
%       a TU. Applicable to median mediated learning (main_cv_median.m) output.
%   plotSlope: boolean. If true, plot linear regression at each leaf (heatmap);
%                       if false, plot violin plot at each leaf



    if nargin == 6
        median_factor = NaN;
        plotSlope = false;
    end
    

    fig=figure;
    
    
    % Dim check: row as genes, cols as samples
    if size(TF_level,1) ~= length(TF_name)
        TF_level = TF_level';
    end
    if size(gene_level,1) ~= length(gene_name)
        gene_level= gene_level';
    end


    TF_training = TF_level(:,training_idx);
    gene_training = gene_level(:,training_idx);
    
    
    % Find node IDs
    allIDs = cell2mat(tree.parent_ids.keys);
    
    % Depth of each nodes
    node_depths = [];
    for i = 1:length(allIDs)
        node_depths = [node_depths, NodeDepth(tree, allIDs(i))];
    end
    
    % nrow of subplot: max depths (+1 to make 1 based)
    figure_nrow = max(node_depths)+1;
    % ncol of subplot 
    figure_ncol = 2^max(node_depths);
    
    
    
    % Get col_ranges of subplots
    col_ranges = PlotTreeStructure_helper(tree);
    
    
    % Plot decision nodes
    decision_ids = allIDs(~isnan(cell2mat(tree.lchild_ids.values)));
    for i = 1:length(decision_ids)
        
        decision_id = decision_ids(i);
        
        all_gene_level = gene_training(:); % enveloped by all gene's data
        
        nbins = floor(length(all_gene_level)/20)+5;
        [cnts,ctrs] = hist(all_gene_level, nbins); % bin counts and bin centers
        % cnts = cnts./max(cnts);
        cnts = smooth(cnts);
        cnts = cnts(:)';
        
        
        figure_id = NodeDepth(tree,decision_id)*figure_ncol+...
            col_ranges(decision_id);
        
        subplot(figure_nrow, figure_ncol, figure_id);
        patch([ctrs,ctrs(end:-1:1)],...
            [cnts,-cnts(end:-1:1)],...
            [ctrs,ctrs(end:-1:1)],'FaceAlpha',0);
        
        set(gca,'ytick',[]);
        set(gca,'yticklabel',[]);
        
        hold on
        
        thisGene_level = gene_training(:,tree.nodes(decision_id).sample_ids);
        
        [cnts,ctrs] = hist(thisGene_level(:), nbins); % bin counts and bin centers
        cnts = smooth(cnts);
        cnts = cnts(:)';
        
        patch([ctrs,ctrs(end:-1:1)],...
            [cnts,-cnts(end:-1:1)],...
            [ctrs,ctrs(end:-1:1)]);
        
        
        % Add mean as vertical line 
        x_coor = mean(thisGene_level,'all');
        xline(x_coor,'--', strcat('mean=',num2str(x_coor,3)),'LineWidth',2,'LabelOrientation','horizontal');
            
              
        TFDecision = strcat(TF_name{tree.nodes(decision_id).decision_feature},'\geq',...
            num2str(tree.nodes(decision_id).decision_val, 3));
        
        if i == 1
            xlabel(sprintf('Gene Expression Level\nNode %d: %s',decision_id, TFDecision));
            ylabel('Distribution');
        else
            xlabel(sprintf('Node %d: %s',decision_id, TFDecision));
        end
        
        % Save subplot x position
        pos = get(gca,'Position');
        subplot_lx_pos(decision_id) = pos(1); % from subplot left edge to figure left boundary
        subplot_rx_pos(decision_id) = pos(1)+pos(3); % from right edge to left boundary
        subplot_by_pos(decision_id) = pos(2); % from bottom edge to bottom boundary
        subplot_ty_pos(decision_id) = pos(2)+pos(4); % from top edge to bottom boundary
    end
    
    
    
    % Gather data at leaves
    leaf_ids = allIDs(isnan(cell2mat(tree.lchild_ids.values)));
    
    if plotSlope
        coeff_store = NaN(length(leaf_ids),length(TF_name));
        r2_store = NaN(length(leaf_ids),length(TF_name));
        pval_store = NaN(length(leaf_ids),length(TF_name));

        SST = 0; % total sum of squares across all leaves 
        SSR = 0; % total residual sum of squares

        for i = 1:length(leaf_ids)
            leaf_id = leaf_ids(i);

            % Find R^2 and F-test pval at this leaf
            y = gene_training(:,tree.nodes(leaf_id).sample_ids);
            y = y';
            thisy = y(:);

            regressor_ids = tree.nodes(leaf_id).ancestral_decision_features;
            X = TF_training(regressor_ids,tree.nodes(leaf_id).sample_ids);
            X = X';
            thisX = repmat(X,size(y,2),1);

            fittedmdl = fitlm(thisX,thisy);
            summaryTable = fittedmdl.Coefficients;

            SST = SST + fittedmdl.SST;
            SSR = SSR + fittedmdl.SSR;

            for j = 1:length(regressor_ids)

                coeff_store(i,j) = table2array(summaryTable(j+1,'Estimate'));
                r2_store(i,j) = fittedmdl.Rsquared.Ordinary;
                pval_store(i,j) = table2array(summaryTable(j+1,'pValue'));

            end
        end
    end
    
    % Plot leaves
    for i = 1:length(leaf_ids)
        leaf_id = leaf_ids(i);
        
        if plotSlope
            M = NaN(length(tree.nodes(leaf_id).ancestral_decision_features), 1);

            % fill M with regression coeff
            regressor_ids = tree.nodes(leaf_id).ancestral_decision_features;

            for j = 1:length(regressor_ids)
                M(j,1) = coeff_store(i,j);
            end
            
        
            this_col_ranges = col_ranges(leaf_id);
            if ismember(leaf_id, cell2mat(tree.lchild_ids.values))
                figure_id = max(node_depths)*figure_ncol+this_col_ranges(1);
            else
                figure_id = max(node_depths)*figure_ncol+this_col_ranges(end);
            end

            subplot(figure_nrow, figure_ncol, figure_id);
            imagesc(M,'AlphaData',~isnan(M));
            colormap(gca, flipud(redgreencmap(64)));
%         if plotSlope
            caxis([-1 1]); % fix colorbar limit: slope more than 1 or -1 are same color
%         else
%             caxis([-8 8]); % mean more than 8 or -8 are same color
%         end
        
            set(gca,'xtick',[]);
            set(gca,'xticklabel',[]);
            set(gca,'ytick',[]);
            set(gca,'yticklabel',[]);

            % Save subplot x position
            pos = get(gca,'Position');
            subplot_lx_pos(leaf_id) = pos(1); % from subplot left edge to figure left boundary
            subplot_rx_pos(leaf_id) = pos(1)+pos(3); % from right edge to left boundary
            subplot_by_pos(leaf_id) = pos(2); % from bottom edge to bottom boundary
            subplot_ty_pos(leaf_id) = pos(2)+pos(4); % from top edge to bottom boundary
        
            % Add regression stats as label on leaves
            for j = 1:length(regressor_ids)
                regressor_id = regressor_ids(j);

                % if coefficient steeper than 1, make background golden
                if abs(coeff_store(i,j)) >= 1
                    bgc = 'y';
                else
                    bgc = 'w';
                end

                text(0.8,... % x coordinate 
                    j,...% y coordinate
                    sprintf('vars:%s\ncoef:%.2f\nR^2  :%.2f\npval:%.2f',TF_name{regressor_id},coeff_store(i,j),r2_store(i,j),pval_store(i,j)),... % text
                    'FontSize',10,'BackgroundColor',bgc);
            end
        
            
        else % Plot violin plot at leaves      
            
            all_gene_level = gene_training(:);
            
            nbins = floor(length(all_gene_level)/20)+5; % set num of bins for hist
            [cnts,ctrs] = hist(all_gene_level, nbins); % bin counts and bin centers
            % cnts = cnts./max(cnts);
            cnts = smooth(cnts);
            cnts = cnts(:)';
            
            % subplot id of this leaf
            this_col_ranges = col_ranges(leaf_id);
            if ismember(leaf_id, cell2mat(tree.lchild_ids.values))
                figure_id = max(node_depths)*figure_ncol+this_col_ranges(1);
            else
                figure_id = max(node_depths)*figure_ncol+this_col_ranges(end);
            end

            subplot(figure_nrow, figure_ncol, figure_id);
            patch([ctrs,ctrs(end:-1:1)],...
                [cnts,-cnts(end:-1:1)],...
                [ctrs,ctrs(end:-1:1)], 'FaceAlpha',0);

            set(gca,'ytick',[]);
            set(gca,'yticklabel',[]);

            hold on
            
            % Add violin plot of this leaf
            thisGene_level = gene_training(:, tree.nodes(leaf_id).sample_ids);
            
            [cnts,ctrs] = hist(thisGene_level(:), nbins); % bin counts and bin centers
            % cnts = cnts./max(cnts);
            cnts = smooth(cnts);
            cnts = cnts(:)';
            
            patch([ctrs,ctrs(end:-1:1)],...
                [cnts,-cnts(end:-1:1)],...
                [ctrs,ctrs(end:-1:1)]);
            
            
            % Add vertical line as mean
            x_coor = mean(thisGene_level,'all');
            xline(x_coor,'--', strcat('mean=',num2str(x_coor,3)),'LineWidth',2,'LabelOrientation','horizontal');
            
            
            % Add leaf ID as xlabel
            xlabel(sprintf('Node %d',leaf_id));
            
            % Save subplot x position
            pos = get(gca,'Position');
            subplot_lx_pos(leaf_id) = pos(1); % from subplot left edge to figure left boundary
            subplot_rx_pos(leaf_id) = pos(1)+pos(3); % from right edge to left boundary
            subplot_by_pos(leaf_id) = pos(2); % from bottom edge to bottom boundary
            subplot_ty_pos(leaf_id) = pos(2)+pos(4); % from top edge to bottom boundary
            
        end
    end
    
    % Add gene names as common xlabel 
    % Also add overall R2
    h=axes(fig,'visible','off');
    h.XLabel.Visible='on';
    
    if plotSlope
        totalR2 = SSR/SST;
    else
        totalR2 = NaN; % total R2 does not apply to MCMC on gaussian learning each leaf
    end
    xlabel(h,sprintf('\nModule Genes: %s\nMedian Factor: %s\nTotal R^2: %.2f',strjoin(gene_name,', '), num2str(1./median_factor,3), totalR2));
    
    h.Title.Visible = 'on';
    title(h,sprintf('Genes Expression Level Distribution at Each Node\n'));
    
    % Draw arrow from parent to children nodes
    for i = 1:length(decision_ids)
        decision_id = decision_ids(i);
        
        lchild_id = tree.lchild_ids(decision_id);
        rchild_id = tree.rchild_ids(decision_id);
        
        % Arrow from decision node bottom to lchild top
        annotation('arrow',[subplot_lx_pos(decision_id)+0.03, subplot_lx_pos(lchild_id)+0.03],...
            [subplot_by_pos(decision_id)-0.01, subplot_ty_pos(lchild_id)+0.01]);
        
        % Arrow from decision node bottom to rchild top
        annotation('arrow',[subplot_rx_pos(decision_id)-0.03, subplot_rx_pos(rchild_id)-0.03],...
            [subplot_by_pos(decision_id)-0.01, subplot_ty_pos(rchild_id)+0.01]);
    end
end


