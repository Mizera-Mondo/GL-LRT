function A = ProjectAOntoVaildSet(A, N)
%ProjectAOntoVaildSet
%   Aii = 0
%   Aij >= 0
%   ||A||_1 = N
K = A;
Q = K;
[n, ~] = size(Q);
G = zeros(size(K));
m = 1;
isAdmmConverge = false;
isMaxIter = false;
tol = 1e-3;
iter = 0;
maxIter = 1000;
while ~isAdmmConverge && ~isMaxIter
    iter = iter + 1;
    K_old = K;
    Q_old = Q;

    K = (m*Q - G + 2*A)/(2+m);
    K = K - diag(diag(K));
    K(K < 0) = 0;

    Q = K + 1/m*G;
    r = N - sum(Q, 'all');
    Q = Q + r/n^2;

    G = G + m*(K - Q);

    if norm(K_old - K, 'fro')/norm(K_old, 'fro') < tol && ...
            norm(Q_old - Q, 'fro')/norm(Q_old, 'fro') < tol
        isAdmmConverge = true;
    end
    if iter >= maxIter
        isMaxIter = true;
    end
end
A = K;
end

