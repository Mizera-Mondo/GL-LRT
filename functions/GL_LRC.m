function [X, L, A] = GL_LRC(Y, R, k, options)
%LR-LRC Solve 1/2||D(Y - X)||_F^2 + alpha*Tr{D(X)'*L*D(X)} +
%beta||L||_F^2

arguments
    Y, R, k double
    options.alpha = 0.1;
    options.beta = 0.1;
    options.graphRefineMethod = 'quadprog';
    options.LowRankEst = true;
    options.debug = false;
end
alpha = options.alpha;
beta = options.beta;

% k is the maximum rank of estimated signal matrix X
X = Y;

% Random vaild initial value of A
[n, T] = size(Y);
B = [zeros(T - 1, 1) eye(T - 1); ...
    zeros(1, T)];
D = @(X) X - R*X*B;
A = rand(n, n);
A = A - diag(diag(A));
A = A + A';
A = n*A./sum(A, 'all');
L = diag(sum(A)) - A;

tol = 1e-3;
iter = 1;
maxIter = 100;
isConverge = false;
isMaxIter = false;

while ~isConverge && ~isMaxIter

    if options.debug
        disp("============================================");
        disp("Iter: " + num2str(iter));
    end

    L_old = L;
    X_old = X;

    % Optimizing L a.k.a. A
    DX = D(X);
    M = genM(DX);
    
    if options.debug
        disp("Starting Graph Refinement...");
        tic
        A = SolveSubA(M, alpha, beta,"method","ADMM");
        toc

    else
        A = SolveSubA(M, alpha, beta,"method","ADMM");
    end



    L = diag(sum(A)) - A;

    % Optimizing X
    
    if options.LowRankEst
        if options.debug
        disp("Starting Low Rank Component Estimation...");
        tic
            X = SolveSubX(Y, R, B, L, alpha, k);
        toc
        else
            X = SolveSubX(Y, R, B, L, alpha, k);
        end
    end
    
    isConverge = norm(L_old - L, 'fro')/norm(L_old, 'fro') < tol ...
            && norm(X_old - X, 'fro')/norm(X_old, 'fro') < tol;
    isMaxIter = iter >= maxIter;
    iter = iter + 1;
end
end


