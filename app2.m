% ============================= Application 2 ============================= %
function [mest, mestSE, RMSE] = app2(n0, a, b, c, delta, p) 
%% Conditions: sigma_epsilon = delta*sqrt(c/b) [0 < delta < 1]
%% Run the following code by substituting the appropriate n0 value 
%% to get the results in Table 2 and Figure 5 of the article. 
%% [mest, mestSE, RMSE] = app2(50, 0.3, 0.2, 0.3, 0.79, 2)

%% Run the following code, substituting suitable parameter values, 
%% to get the results in Table 3 and Figure 6 of the article. 
%% [mest, mestSE, RMSE] = app2(30, 0, 0.92, 0.001, 0.79, 2)

format short g
close ALL
rng(702); % For reproducibility
sigma0 = sqrt(0.77);
sigma_epsilon = delta*sqrt(c/b);
[x, y] = coal(); %to get coal-ash data
x1max = max(x(:,1)) + 1;
x2max = max(x(:,2)) + 1;

N = length(y); %The total number of obervations
Index = 1:N;
K = 150; %The number of runs
n1 = round(n0/2,0); %The initial sample size

PN = isocoal(x); % isotropic Gaussian correlations
ICN = (sigma0^2).*PN; %For spatial design
IN = eye(N); %Identity matrix

F = zeros(N, p+1);
F(:,1) = 1;

for i = 2:(p+1)
    F(:,i) = x(:,(i-1));
end

Iloss = 1:K;
Gloss = zeros(1, K);
Gminloss = 100000;
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
    [freqP, UIndex, n, ~] = sequential(Index, IIndex, n1, rn, ICN, IDesign, IFloss, N, a, b, c, sigma_epsilon, IN, F);

    freqP = freqP - IDesign;
    IDesign = freqP;
    IIndex = setdiff(UIndex, IIndex);
    n1 = n - n1;
    rn = n0 - n1;
    CN = ICN;
    IFloss = zeros(1, rn + 1);
    IFloss(1,1) = plossnew(N, n1, IDesign, a, b, c, sigma_epsilon, CN, IN, F);

    [~, UIndex, ~, Floss] = sequential(Index, IIndex, n1, rn, ICN, IDesign, IFloss, N, a, b, c, sigma_epsilon, IN, F);
    Gloss(k) = Floss(rn+1);
    if Gminloss > Gloss(k)
        Gminloss = Gloss(k);
        GIndex = UIndex;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rx1 = x(:,1);
rx2 = x(:,2);
x1 = x(GIndex,1);
x2 = x(GIndex,2);
markerSizes = 20;

%Plot for n0 optimal locations
figure(1)
scatter(rx1, rx2, markerSizes)
hold on;
scatter(x1, x2, markerSizes, 'r', 'filled');
yticks(0:2:(x2max));
xticks(0:2:(x1max));
xlim([0 x1max]);
ylim([0 x2max]);
minloss = min(Gloss);
[~, I] = min(Gloss);

lossmax = round(((max(Gloss) + 1)/1))*1;
lossmin = floor(((min(Gloss))/1))*1;

textCell = arrayfun(@(x1,y1) sprintf('(%d, %3.3f)',x1,y1),I,minloss,'un',0);

%Plot for 150 losses
figure(2)
plot(Iloss, Gloss);
hold on;
scatter(I, minloss, markerSizes, 'r', 'filled')
text(I+.2, minloss-.1,textCell,'FontSize',8) 
yticks(lossmin:1:lossmax);
xticks(0:15:K);
xlim([0 K+1]);
ylim([lossmin lossmax]);

%M-estimates using information in n0 selected optimal locations
[mest, stats] = robustfit(x(GIndex,:), y(GIndex)); 
mest = round(mest, 4); %M-estimates

mestSE = stats.se;
mestSE = round(mestSE, 4); %Standard errors of M-estimates

%[mestF, ~] = robustfit(x, y);
%mestF = round(mestF, 4); %M-estimates using information in the full data
%Edis = round(sqrt(sum((mest - mestF).^2)),4); %Euclidean distance

TF = zeros(N, p+1); %To save test features
TL = zeros(N, 1); %To save test labels
TL(:,1) = y';
TF(:,1) = 1;
TF(:,2) = x(:,1);
TF(:,3) = x(:,1);

yest = TF*mest; %Forecasted y

RMSE = sqrt(mean((TL - yest).^2));


