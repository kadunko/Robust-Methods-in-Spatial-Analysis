% ============================= twomainsequential ============================= %
function [valmin, minloss] = bfsearch(M, n0, a, b, c, delta, p) 
%% Conditions: sigma_epsilon = delta*sqrt(c/b) [0 < delta < 1]
%% Run the following code, substituting suitable parameter values, 
%% to get the results in Table 1, Figure 1, and Figure 2 of the article.
%% M=5, n0 = 7, and PN with the isotropic Gaussian correlation structure for Figure 1
%% M=6, n0 = 8, and PN with the anisotropic Gaussian correlation structure for Figure 2
%% [valmin, minloss] = bfsearch(6, 8, 0.1, 0.2, 0.5, 0.4, 2)

format short g
close ALL
rng(200); % For reproducibility

sigma0 = 0.3;
sigma_epsilon = delta*sqrt(c/b);

[x] = xp(M, p);

N = M^p;
Index = 1:N;
Ksample = nchoosek(Index,n0);
K = length(Ksample(:,1));
Floss = zeros(K,1);

%PN = iso(M, p); % isotropic Gaussian correlations
PN = aniso(M, p); % anisotropic Gaussian correlations
CN = (sigma0^2).*PN; %For spatial design
IN = eye(N); %Identity matrix

F = zeros(N, p+1);
F(:,1) = 1;

for i = 2:(p+1)
    F(:,i) = x(:,(i-1));
end

for k = 1:K
    Istar = Ksample(k,:);
    freqP = zeros(N,1);
    for i = 1:n0
        freqP(Istar(i)) = 1;
    end  
    Floss(k) = plossnew(N, n0, freqP, a, b, c, sigma_epsilon, CN, IN, F);
end

minloss = min(Floss);
mem = ismember(Floss, minloss);
valmin = sum(mem);
markerSizes = 80;
[~, I] = min(Floss);
Istar = Ksample(I,:);
x1 = x(Istar,1);
x2 = x(Istar,2);
scatter(x1, x2, markerSizes, 'r', 'filled')
yticks([0 0.25 0.5 0.75 1])
xticks([0 0.25 0.5 0.75 1])

