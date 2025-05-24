%% =============== loss function =============================%%
function [lossvalue] = plossnew(N, n, freqP, a, b, c, sigma_epsilon, CN, IN, F)
E = zeros(N, n); %Depends on n
nc = 1;
for i = 1:N
    if freqP(i) > 0
        E(i, nc:(nc + freqP(i) - 1)) = 1;
        nc = nc + freqP(i);
    end
end

In = eye(n); %Depend on n
Cn = ((E'*CN)*E) + ((sigma_epsilon^2).*In);
r1 = rcond(Cn);
if r1 > 1e-16
    V0 = (E/Cn)*E';
    r2 = rcond((F'*V0)*F);
    if r2 > 1e-16
        L0 = ((((IN - CN*V0)*((F/((F'*V0)*F))*F')) + CN)*(E/Cn));
        A0 = IN - (L0*E');
        A = A0*A0';
        L1 = eigs(A, 1);
        L2 = trace(A);
        L3 = trace(n.*((A0*CN)*A0'));
        L4 = trace(n.*(L0*L0'));
    else 
        L1 = 200;
        L2 = 200;
        L3 = 200;
        L4 = 200;
    end 
else 
    L1 = 200;
    L2 = 200;
    L3 = 200;
    L4 = 200;
end
lossvalue = ((1 - a - b - c)*L1) + (a*L2) + (b*L3) + (c*L4);

