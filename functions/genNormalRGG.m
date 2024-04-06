function [A, L] = genNormalRGG(nodeNum, options)
%GENRANDOMRGG generate random geometric graph normalized to its node count.
    arguments
        nodeNum double
        options.dimen = 2;
        options.thre = 0.85;
    end
    
    sigma = 0.5;
    reCalcWeight = @(M) exp(-1/2*sigma^-2*M.^2);
    [~, D] = random_geometric_graph(nodeNum, options.dimen, options.thre);
    A = reCalcWeight(D);
    % A is natrually symmetric
    A(A < options.thre) = 0;
    A = A./sum(A, "all")*nodeNum;
    L = diag(sum(A)) - A;
end

