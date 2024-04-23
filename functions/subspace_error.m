function err = subspace_error(B, v)
%SUBSPACE_ERROR
Pb = B*inv(B'*B)*B';
err = norm(v - Pb*v, "fro")^2;
end

