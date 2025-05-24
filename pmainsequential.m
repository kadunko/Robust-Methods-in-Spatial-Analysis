% ========== pmainsequential function ============================= %
function [minloss] = pmainsequential(M, n0, sigma0, delta, a, b, c, p) 
%% Conditions: sigma_epsilon = delta*sqrt(c/b) [0 < delta < 1]
%% Run the following code, substituting suitable parameter values, 
%% to get the results in Table 1, Figure 1, and Figure 2 of the article.
%% M=5, n0 = 7, and PN with the isotropic Gaussian correlation structure for Figure 1
%% M=6, n0 = 8, and PN with the anisotropic Gaussian correlation structure for Figure 2
%% [minloss] = pmainsequential(6, 8, 0.3, 0.4, 0.1, 0.2, 0.5, 2)

format short g
close ALL
rng(202); % For reproducibility

sigma_epsilon = delta*sqrt(c/b);

[x] = xp(M, p);

N = M^p;
Index = 1:N;
K = 1500; %The number of runs
n1 = 4; %The initial sample size

%PN = iso(M, p); % isotropic Gaussian correlations
PN = aniso(M, p); % anisotropic Gaussian correlations
ICN = (sigma0^2).*PN; %For spatial design
IN = eye(N); %Identity matrix

F = zeros(N, p+1);
F(:,1) = 1;

for i = 2:(p+1)
    F(:,i) = x(:,(i-1));
end

Gloss = zeros(1, K);
Gminloss = 1000;
GIndex = zeros(n0, 1);

for k = 1:K
    Istar = randsample(N,n1);
    freqP = zeros(N,1);
    for i = 1:n1
        freqP(Istar(i)) = 1;
    end 
    
    IDesign = freqP; % Initial Design
    UIndex = Istar; %Used Index
    IIndex = UIndex; %Initial Index
    rn = n0 - n1;
    CN = ICN;
    IFloss = zeros(1, rn + 1);
    IFloss(1,1) = plossnew(N, n1, IDesign, a, b, c, sigma_epsilon, CN, IN, F);
    [freqP, UIndex, n, Floss] = sequential(Index, IIndex, n1, rn, ICN, IDesign, IFloss, N, a, b, c, sigma_epsilon, IN, F);

    freqP = freqP - IDesign;
    IDesign = freqP;
    IIndex = setdiff(UIndex, IIndex);
    n1 = n - n1;
    rn = n0 - n1;
    CN = ICN;
    IFloss = zeros(1, rn + 1);
    IFloss(1,1) = plossnew(N, n1, IDesign, a, b, c, sigma_epsilon, CN, IN, F);

    [freqP, UIndex, n, Floss] = sequential(Index, IIndex, n1, rn, ICN, IDesign, IFloss, N, a, b, c, sigma_epsilon, IN, F);
    Gloss(k) = Floss(rn+1);
    if Gminloss > Gloss(k)
        Gminloss = Gloss(k);
        GIndex = UIndex;
    end
end

x1 = x(GIndex,1);
x2 = x(GIndex,2);
markerSizes = 80;
scatter(x1, x2, markerSizes, 'r', 'filled')
yticks([0 0.25 0.5 0.75 1])
xticks([0 0.25 0.5 0.75 1])
minloss = min(Gloss);
