% ==== construction of locations x ============================= %
function [x] = xp(M, p)
%% eg: [x] = xp(25, 2)
N = M^p;
x = zeros(N, p);
pi_x = linspace(0, 1, M);
XN = zeros(N, 1);
h = 0;

if p == 1
    x = pi_x'; 
elseif p == 2
    [x1, x2] = meshgrid(pi_x, pi_x);
    X1 = XN;
    X2 = XN;
    for i = 1:M
        for j = 1:M
            h = h + 1;
            X1(h) = x1(i,j);
            X2(h) = x2(i,j);
            x(h, 1) = X1(h);
            x(h, 2) = X2(h);
        end   
    end
elseif p == 3
    [x1, x2, x3] = meshgrid(pi_x, pi_x, pi_x);
    X1 = XN;
    X2 = XN;
    X3 = XN;
    for i = 1:M
        for j = 1:M
            for k = 1:M
                h = h + 1;
                X1(h) = x1(i,j,k);
                X2(h) = x2(i,j,k);
                X3(h) = x3(i,j,k);
                x(h, 1) = X1(h);
                x(h, 2) = X2(h);
                x(h, 3) = X3(h);
            end
        end   
    end
end


