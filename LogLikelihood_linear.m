function L = LogLikelihood_linear(tree, mu_b, V_b, a, b, training_set, Y)
% Calculate log-likelihood of data given a tree.

% Input
% tree: an object of Tree class
% mu_b, V_b, a, b: hyperparameters of linear regression coefficient N(w|mu_b, sigma^2*V_b)IG(sigma^2|a,b)
%   Input mu_b = scalar for scalar*one_vector of length(w); V_b = scalar for scalar*I of size(w)

% training_set: TF expression matrix: nsamples x ngenes
% Y: target gene expression matrix: nsamples x ngenes

% Output
% L: loglikelihood


% if Y is a row vector, convert to col vector
if size(Y,1) == 1 
    Y = Y(:);
end


% Find leaves
allID = cell2mat(tree.lchild_ids.keys);
leafIDs = allID(isnan(cell2mat(tree.lchild_ids.values)));

L = 0;

for i = 1:length(leafIDs)
    
    leafID = leafIDs(i);
    
    dataID = tree.nodes(leafID).sample_ids;
    % target value in this leaf
    y = Y(dataID,:);
    
    if isempty(y)
        continue
    end
    
    features = sort(unique(tree.get_ancestral_decision_features(leafID)));
    if isempty(features) % at root node
        features = 1:size(training_set,2);
    end
    X = training_set(dataID, features);
    
    
    % stack up all col
    thisy = y(:);
    thisX = [ones(length(thisy),1),repmat(X,size(y,2),1)];
    
    % P(data) = multivariate student t pdf
    if isscalar(V_b)
        V_b_autofill = V_b*eye(size(thisX,2));
    end
    
    if isscalar(mu_b)
        mu_b_autofill = repmat(mu_b,size(thisX,2),1);
    end
    nu = 2*a;
    sigma = b/a*(eye(size(thisX,1))+thisX*V_b_autofill*thisX');
    
    % multivariate student t with dof = 2a, cov = b/a(I+X*V_b*X'), evaluate
    % at y-X*mu_b
    L = L + mvtpdf_log((thisy-thisX*mu_b_autofill),sigma,nu);
    
    % If not to use Matlab function
    % n = size(thisX,1);
    % inp = thisy-thisX*mu_b;
    % L = L + gammaln(nu/2+n/2)-(nu/2+n/2)*log(1+1/nu*inp*(sigma\inp))...
    % -gammaln(nu/2)-n/2*log(nu*pi)-0.5*log(det(sigma))

end
