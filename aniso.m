% ======= PN with the anisotropic Gaussian correlation structure ======== %
function [PN] = aniso(M, p)
[x] = xp(M, p);
lambda = 0.9;
N = M^p;
PN = zeros(N, N);

for i = 1:N
    for j = 1:N
        dx = x(i,:) - x(j,:);
        d = sqrt((dx*(diag([1 5])))*(transpose(dx)));
        PN(i, j) = exp(-lambda*(d^2));
    end
end


