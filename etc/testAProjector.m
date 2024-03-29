numN = 10;
A = randn(numN);
P = ProjectAOntoVaildSet(A, 10);
dis = norm(A - P, 'fro');
mind = norm(A, 'fro');
for i = 1:1000
    K = P + rand(numN); 
    K = K./sum(K, 'all')*10;
    d = norm(K - A, 'fro');
    if d < mind 
        mind = d; 
    end
end
disp("Projected distance: " + num2str(dis) + ", minimal random distance: " + num2str(mind));