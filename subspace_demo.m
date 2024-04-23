nodeNum = 30;
signalLength = 100;
eigenPara = 0.3;
noiseCov = 0.01;
rPertubation = 0;
% Signal Generation
[A, L] = genNormalRGG(nodeNum);
usedEigNum = ceil(nodeNum*eigenPara);
[Y, R] = genRandomSignal(L, usedEigNum, signalLength, noiseCov, rPertubation);
Y = Y./sqrt(signalLength);
% L2 = L + 0.05*randn(nodeNum);
% L2 = (L2 + L2')./2;
L2 = L;
L2(:, 1) = zeros(nodeNum, 1);
L2(1, :) = zeros(1, nodeNum);
[Y2, R2] = genRandomSignal(L2, usedEigNum, signalLength, noiseCov, rPertubation);
Y2 = Y2./sqrt(signalLength);

B = zeros(signalLength);
B(1:end - 1, 2:end) = eye(signalLength - 1);
D = @(X) X - R*X*B;
V = D(Y);
V2 = D(Y2);
[ve, val] = eig(L);
subBase = ve(:, 1:usedEigNum);
error1 = zeros(1, 100);

error2 = error1;

 for i = 1:100
error1(i) = subspace_error(subBase, V(:, i));
error2(i) = subspace_error(subBase, V2(:, i));
 end

errorExp = zeros(1, signalLength*2);
wd = 20;
alpha = 0.01;
Va = [V V2];
for i = 1:wd - 1
    errorExp(i) = subspace_wmean_error(subBase, Va(:, 1:i), alpha, wd);
end
 for i = wd:signalLength*2
     errorExp(i) = subspace_wmean_error(subBase, Va(:, i - wd + 1:i), alpha, wd);

 end


close all
figure;
subplot(2, 2, 1)
imagesc(L);
title("Original Laplacian");
subplot(2, 2, 2)
imagesc(L2);
title("Perturbed Laplacian");

subplot(2, 2, [3 4])
semilogy([error1 error2], "LineWidth",1);
hold on
semilogy(errorExp, "LineWidth",1);
grid on
legend("Normal", "Weighted"); xlabel("Sample No."); ylabel("Subspace Error");