function main(file)
% Input file name as string: 'data.mat'

% Load data processed locally
loadedFile = load(file);

TF_name = loadedFile.TF_name;
gene_name = loadedFile.gene_name;
TF_level_log = loadedFile.TF_level_log;
gene_level_log = loadedFile.gene_level_log;
weights = loadedFile.weights;


% Genes are normalized to have same median as first gene in TU
median_factor = median(gene_level_log,1);
median_factor = median_factor(1)./median_factor;
gene_level_log = gene_level_log.*repmat(median_factor,size(gene_level_log,1),1);


% split data into training and validation set
all_idx = 1:size(TF_level_log,1);
training_idx = datasample(all_idx, floor(0.8*length(all_idx)), 'Replace',false);
validation_idx = all_idx(~ismember(all_idx, training_idx));

TF_training = TF_level_log(training_idx,:);
TF_validation = TF_level_log(validation_idx,:);

gene_training = gene_level_log(training_idx,:);
gene_validation = gene_level_log(validation_idx,:);

% get hyperparameters and tree
alpha_var = 0.95;
beta_vars = [80];

BLR = false;
if BLR == false
    % single mean variance hyperparam
    mu = mean(gene_level_log(:));
    a = 1;
    nu=1;
    lambda=1;
else
    % bayesian linear regression hyperparam
    mu = 0; % i.e. mu_b
    a = 1; % i.e. V_b
    nu = 1; % i.e. a
    lambda = 1; %i.e. b
end

n = 20;
iter = 5000;

for i = 1:length(beta_vars)
    beta_var = beta_vars(i);

    disp('======================================');
    disp(strcat('alpha = ',num2str(alpha_var)));
    disp(strcat('beta = ',num2str(beta_var)));
    
    % Create tree object
    lchildren = containers.Map([1],[NaN]);
    rchildren = containers.Map([1],[NaN]);

    decision_features = containers.Map([1],[NaN]);
    decision_vals = containers.Map([1],[NaN]);

    T = Tree(lchildren, rchildren, decision_features, decision_vals);

    T.assign_data(TF_training);
    
    [chain, lnprob, acceptance_rate] = MCMC(T, mu, a, nu, lambda, alpha_var, beta_var, n, TF_training, gene_training, iter, BLR, weights);

    % find max a posteriori tree and prediction accuracy
    T_max = chain(lnprob == max(lnprob));
    T_max = T_max(1).copy(); % in case more than 1 highest prob tree
    
    % fill node_mean and node_var on T_max
    T_max.update_leaf_value(gene_training);
    % fill linear_reg_coef on T_max
    T_max.update_leaf_value_linear(TF_training, gene_training)
    
    % CV error score
    if BLR
        predicted_y = T.max.predict_linear(TF_validation);
    else
        predicted_y = T_max.predict(TF_validation);
    end
    Error = sum((predicted_y-mean(gene_validation,2)).^2);
    
    % save chain, lnprob, acceptance_rate
    save(strcat('beta_',num2str(beta_var),'_alpha_',num2str(alpha_var*100),'%','_',date,'_',strtok(file,'.'),'.mat'),'chain','lnprob',...
        'acceptance_rate','T_max','Error','training_idx','validation_idx','TF_level_log','gene_level_log',...
        'TF_name','gene_name','weights','median_factor');
    
    clear lchildren;
    clear rchildren;
    clear decision_features;
    clear decision_vals;
    clear T;
    clear Tree;
    clear chain;
    clear T_max;
    
end
end



