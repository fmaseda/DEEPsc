function [M, L, Y, C] = lmnn2(X)
% LMNN Learns a metric using large-margin nearest neighbor metric learning
%
%   [M, L, Y, C] = lmnn2(X)
%
% The function uses large-margin nearest neighbor (LMNN) metric learning to
% learn a metric on the data set specified by the NxD matrix X, treating
% all points as impostors. The metric is returned in M.
% 
% This is a modification to lmnn() where there is no pull term at all and
% only target for each point is itself; points are moved apart until
% minimum separation distance is sqrt(D)
% ----------
% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    % Initialize some variables
    [N, D] = size(X);
    self = logical(eye(N));
    M = eye(D);
    C = Inf; prev_C = Inf;
    
    % Set learning parameters
    min_iter = 50;          % minimum number of iterations
    max_iter = 10000;       % maximum number of iterations
    eta = .1;               % learning rate
    tol = 1e-3;             % tolerance for convergence (default=1e-3)
    best_C = Inf;           % best error obtained so far
    best_M = M;             % best metric found so far
    
    % Select target neighbors
    slack = zeros(N, N);    % slack(i,j) = error in distance from X(i,:) to X(j,:)
    G = zeros(D, D);        % gradient for learning
    
    % Perform main learning iterations
    iter = 0;
    while (prev_C - C > tol || iter < min_iter) && iter < max_iter
        
        % Compute pairwise (squared) distances under current metric
        XM = X * M;
        DD = sum(XM.*X, 2)+sum(XM.*X, 2)'-2*(XM*X');
        
        % Compute value of slack variables
        old_slack = slack;
        slack = max(0, sqrt(D)-DD);     % slack(i,j)=0 if X(i,:) is more than sqrt(D) away from X(j,:)
        slack(self) = 0;                % don't count self-distance as error
        
        % Compute value of cost function
        prev_C = C;
        C = sum(slack(:));
        
        % Maintain best solution found so far (subgradient method)
        if C < best_C
            best_C = C;
            best_M = M;
        end
        
        % Perform gradient update            
        [r, c] = find(slack > 0 & old_slack == 0);
        G = G - (X(r,:)-X(c,:))'*(X(r,:) - X(c,:)); % add terms for new violations
        [r, c] = find(slack == 0 & old_slack > 0);
        G = G + (X(r,:)-X(c,:))'*(X(r,:)-X(c,:));   % remove terms for resolved violations

        M = M - (eta ./ N) .* G;
        
        % Project metric back onto the PSD cone
        [V, L] = eig(M);
        V = real(V); L = real(L);
        ind = find(diag(L) > 0);
        if isempty(ind)
            warning('Projection onto PSD cone failed. All eigenvalues were negative.'); break
        end
        M = V(:,ind) * L(ind, ind) * V(:,ind)';
        if any(isinf(M(:)))
            warning('Projection onto PSD cone failed. Metric contains Inf values.'); break
        end
        if any(isnan(M(:)))
            warning('Projection onto PSD cone failed. Metric contains NaN values.'); break
        end
        
        % Update learning rate
        if prev_C > C
            eta = eta * 1.01;
        else
            eta = eta * .5;
        end
        
        % Print out progress
        iter = iter + 1;
        no_slack = sum(slack(:) > 0);
        if rem(iter, 100) == 0
            disp(['Iteration ' num2str(iter) ': error is ' num2str(C ./ N) ...
                  ', number of constraints: ' num2str(no_slack)]);
        end
    end
    
    disp('----------')
    disp(['Final iteration ' num2str(iter) ': error is ' num2str(C ./ N) ...
      ', number of constraints: ' num2str(no_slack)]);

    
    % Return best metric and error
    M = best_M;
    C = best_C;
    
    % Compute mapped data
    [L, S, ~] = svd(M);
    L = bsxfun(@times, sqrt(diag(S)), L);
    Y = X * L;
end