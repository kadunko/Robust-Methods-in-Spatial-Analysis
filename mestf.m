%% ====== M-estimates using the full data  ======= %%
function [mest, mestSE, RMSE] = mestf(p)
%Run the following code to obtain results for the scenario n0=208 
%in Table 2 of the paper.

%[mest, mestSE, RMSE] = mestf(2)
[x, y] = coal();
[mest, stats] = robustfit(x, y);
mest = round(mest, 4);

mestSE = stats.se;
mestSE = round(mestSE, 4);
N = length(y);
TF = zeros(N, p+1); %To save test features
TL = zeros(N, 1); %To save test labels
TL(:,1) = y';
TF(:,1) = 1;
TF(:,2) = x(:,1);
TF(:,3) = x(:,1);

yest = TF*mest;

RMSE = sqrt(mean((TL - yest).^2));