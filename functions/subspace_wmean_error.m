function subError = subspace_wmean_error(subBase, Vw, alpha, len)
%SUBSPACE_WMEAN_ERROR 此处显示有关此函数的摘要
%   此处显示详细说明
[p, l] = size(Vw);
if l<len
    Vw = [zeros(p, len - l) Vw];
end
rectifier = 1/(sum(exp((-1*len+1:0).*alpha)));
subError = 0;
for i = 0:len - 1
    subError = subError + exp((-1*i).*alpha)*subspace_error(subBase, Vw(:, len - i));
end
subError = subError*rectifier;

end

