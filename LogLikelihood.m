function L = LogLikelihood(tree, mu, a, nu, lambda, Y)
% Calculate log-likelihood of data given a tree.

% Input
% tree: an object of Tree class
% mu, a, nu, lambda: hyperparamaters. Refer to definition in MCMC function
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
    
    % convert data into vector for easier calculation of sufficient stats
    y = y(:);
    n = length(y);
    
    s = (n-1) * var(y);
    t = n*a/(n+a)*(mean(y) - mu)^2;
    
    
    L = L + ...
        (-n*0.5)*log(pi) + ...
        nu*0.5 * log(lambda*nu) + ...
        0.5*(log(a)-log(n+a)) + ...
        gammaln((n+nu)*0.5) - ... % change to log base 2
        gammaln(nu*0.5) - ...
        ((n+nu)*0.5)*log(s+t+nu*lambda);
    

end