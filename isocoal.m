% ========= PN for coal-ash data ============================= %
function [PN] = isocoal(x)
N = length(x(:,1));
lambda = 0.027;
PN = zeros(N, N);

for i = 1:N
    for j = 1:N
        d = norm(x(i,:) - x(j,:));
        PN(i, j) = exp(-lambda*(d^2));
    end
end

