% ===== PN with the isotropic Gaussian correlation structure =========== %
function [PN] = iso(M, p)
[x] = xp(M, p);
lambda = 0.9;
N = M^p;
PN = zeros(N, N);

for i = 1:N
    for j = 1:N
        d = norm(x(i,:) - x(j,:));
        PN(i, j) = exp(-lambda*(d^2));
    end
end

