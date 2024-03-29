function A = SolveSubA(M, alpha, beta, options)
%solveSubA solve the sub-problem ||A||_F^2 + ||A*1||_2^2 + alpha/(2*beta)*Tr{AM}
arguments
    M, alpha, beta double
    options.method = 'quadprog'
end

[n, ~] = size(M);

if strcmp(options.method, 'quadprog')
    
    Aeq = [];
    beq = [];
    Aie = [];
    bie = [];

    % Construct target function for vectorized A
    U = kron(ones(1, n), eye(n));
    H = eye(n^2) + U'*U;
    f = alpha/beta*mat2vec(M);

    % Construct constraints for vectorized A

    % Equality part: 1'A1 = n, A = A', Aii = 0
    % 1'A1 = n
    Aeq = ones(1, n^2);
    beq = n;
    % A = A'
    for i = 1:n
        for j = 1:i - 1
            aeq = zeros(n);
            aeq(i, j) = 1;
            aeq(j, i) = -1;
            Aeq = [Aeq; (mat2vec(aeq))'];
            beq = [beq; 0];
        end
    end
    % Aii = 0
    for i = 1:n
        aeq = zeros(n);
        aeq(i, i) = 1;
        Aeq = [Aeq; (mat2vec(aeq))'];
        beq = [beq; 0];
    end

    % Inequality part: Aij >= 0, i ~= j
    Aie = -1*eye(n^2);
    bie = zeros(n^2, 1);

    A = quadprog(H, f, Aie, bie, Aeq, beq, [], [], [], optimoptions('quadprog', 'Display','off'));
    A = vec2mat(A, n);

elseif strcmp(options.method, 'CVX')
    on = ones(n, 1);
    ze = zeros(n, 1);
    cvx_begin quiet
        variable A(n, n) symmetric nonnegative
        minimize square_pos(norm(A, 'fro')) + square_pos(norm(A*ones(n, 1), 'fro')) + alpha/(2*beta)*trace(A*M)
        subject to
            diag(A) == ze;
            on'*A*on == n;
    cvx_end
elseif strcmp(options.method, 'ADMM')
    isADMMConverge = false;
    isMaxIter = false;
    % Initialize vars
    A = zeros(n, n);
    K = A;
    Phi = zeros(n, n);
    rho = 1;
    tol = 1e-4;
    iter = 0;
    maxIter = 1000;
    while ~isADMMConverge && ~isMaxIter
        iter = iter + 1;
        A_old = A;
        K_old = K;
        % Update of A
        A = (rho*K - Phi - alpha/beta*M)/((2+rho)*eye(n) + 2*ones(n, n));
        % Update of K
        K = ProjectAOntoVaildSet(A + Phi./rho, n);
        % Update of Phi, rho
        Phi = Phi + rho*(A - K);
        rho = 1.05*rho;
        % Convergence check
        if norm(A_old - A, 'fro')/(norm(A_old, 'fro')+1e-10) < tol && ...
                norm(K_old - K, 'fro')/(norm(K_old, 'fro')+1e-10) < tol
            isADMMConverge = true;
        end
        if iter >= maxIter
            isMaxIter = true;
        end
    end
else
    error('%s is not a vaild solver!', options.method);
end

end