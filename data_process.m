res = zeros(33, 18);
for i = 2:1651
    k = i - 1;
    idx = mod(k, 33);
    idx = idx + (idx == 0)*33;
    res(idx, :) = res(idx, :) + res_list(i, :);
end
res = res./50;