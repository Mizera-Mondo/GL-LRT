function A = SolveSubA_PARA(M, alpha, beta, options)
%solveSubA solve the sub-problem ||A||_F^2 + ||A*1||_2^2 + alpha/(2*beta)*Tr{AM}
arguments
    M, alpha, beta double
    options.method = 'quadprog'
end

if strcmp(options.method, 'quadprog')

    [n, ~] = size(M);
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
    parfor i = 1:n
        for j = 1:i - 1
            aeq = zeros(n);
            aeq(i, j) = 1;
            aeq(j, i) = -1;
            Aeq = [Aeq; (mat2vec(aeq))'];
            beq = [beq; 0];
        end
    end
    % Aii = 0
    parfor i = 1:n
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
    [n, ~] = size(M);
    on = ones(n, 1);
    ze = zeros(n, 1);
    cvx_begin quiet
        variable A(n, n) symmetric nonnegative
        minimize square_pos(norm(A, 'fro')) + square_pos(norm(A*ones(1, n), 'fro')) + alpha/(2*beta)*trace(A*M)
        subject to
            diag(A) == ze;
            on'*A*on == n;
    cvx_end

else
    error('%s is not a vaild solver!', options.method);
end

end