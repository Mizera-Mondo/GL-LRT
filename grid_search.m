% Add current paths
currentPath = pwd;
addpath(genpath(currentPath));

beta_list = [0.1, 0.5, 1, 2, 5, 10, 15, 20, 27, 34, 42, 50];
nodeNum = 25;
signalLength = 100;
eigenPara = 0.5;
times = 1:12;
res_list = zeros(1, 18);
parfor i = times
            usedEigNum = ceil(nodeNum*eigenPara);
            noiseCov = 0.1;
            rPertubation = 0;
            % Signal Generation
            [A, L] = genRandomGraph(nodeNum);
            
            [Y, R] = genRandomSignal(L, usedEigNum, signalLength, noiseCov, rPertubation);
            
            B = zeros(signalLength);
            B(1:end - 1, 2:end) = eye(signalLength - 1);
            D = @(X) X - R*X*B;
            alpha = 0.1;
             for beta = beta_list
                % For Debug: The Target Function
                targ1 = @(X, L) alpha*trace((D(X))'*L*D(X));
                targ2 = @(X, L) beta*norm(L, "fro")^2;
                targ3 = @(X) norm(D(X - Y), 'fro')^2;

                Ar = rand(nodeNum);
                Ar = Ar - diag(diag(Ar));
                Ar = Ar + Ar';
                Ar = Ar./sum(Ar, "all")*nodeNum;
                Lr = diag(sum(Ar)) - Ar;
                
                % Estimation
                tic
                [X, Lest] = GL_LRC(Y, R, usedEigNum, alpha = alpha, beta = beta, debug = false, LowRankEst = true);
                t = toc
                DX = D(X);
                M = genM(DX);
                
                
                % Results
                
                Aest = diag(diag(Lest)) - Lest;

                bestfM = 0;
                bestthreA = 0;
                for threA = 0.01:0.01:0.3
                    [~, ~, ~, fM] = classifierPerformance(A > threA, Aest > threA);
                    if fM >= bestfM
                        bestfM  = fM;
                        bestthreA = threA;
                    end
                     [a, r, p, fM] = classifierPerformance(A > bestthreA, Aest > bestthreA);
                end

                errLap = norm(Lest - L, 'fro')/norm(L, 'fro');
                nmse = 100*errLap;
                res = [beta, nodeNum, signalLength, usedEigNum, t, nmse, i, a, r, p, fM, targ1(X, L), targ2(X, L), targ1(X, Lest), targ2(X, Lest), targ3(X), targ1(X, Lr), targ2(X, Lr)];
                res_list = [res_list; res];
                
            end
end
save(sprintf("data\\nodeNum_%d_signalLength_%d_eigen_%0.1f.mat",nodeNum, signalLength, eigenPara),"res_list");