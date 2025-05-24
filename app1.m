% ============================= Application 1 ============================= %
function [minloss] = app1(r0, n0, a, b, c, delta, p) 
%% Conditions: sigma_epsilon = delta*sqrt(c/b) [0 < delta < 1]
%% Run the following code to get the results in Figure 3 and Figure 4 of the article. 
%% [minloss] = app1(90, 11, 0.4, 0.3, 0.2, 0.6, 2)

format short g
close ALL
rng(202); % For reproducibility

rx = rand(1,r0);
ry = rand(1,r0);
M = round(sqrt(r0),0) + 2;

sigma0 = 1;
sigma_epsilon = delta*sqrt(c/b);

[x] = xp(M, p);

N = M^p;
Index = 1:N;
K = 100; %The number of runs
n1 = round(n0/2,0);

PN = iso(M, p); % isotropic Gaussian correlations
%PN = aniso(M, p); % anisotropic Gaussian correlations
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = x(GIndex,1);
x2 = x(GIndex,2);
markerSizes = 20;

sx = zeros(n0, p); %To collect optimal real locations
dx = zeros(n0, p); %To save grid locations

for i = 1:n0
    dx(i, 1) = x1(i);
    dx(i, 2) = x2(i);
end

gx = zeros(r0, p); %To save permanent locations

for i = 1:r0
    gx(i, 1) = rx(i);
    gx(i, 2) = ry(i);
end

for i = 1:n0
    eds = zeros(r0, 1);
    for j = 1:r0
        eds(j) = norm(dx(i,:) - gx(j,:));        
    end
    [~, J] = min(eds);
    sx(i,:) =  gx(J,:);
    gx(J,:) = [];
    r0 = r0 - 1;
end
%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FIGURE 3(a) in the paper
figure(1)
scatter(rx, ry, markerSizes);
yticks([0 0.25 0.5 0.75 1]);
xticks([0 0.25 0.5 0.75 1]);
xlim([-0.05 1.1]);
ylim([-0.05 1.1]);
ylabel('x_2');
xlabel('x_1');

%FIGURE 3(b)
figure(2)
scatter(x(:,1), x(:,2), markerSizes, 'r');
yticks([0 0.25 0.5 0.75 1]);
xticks([0 0.25 0.5 0.75 1]);
xlim([-0.05 1.1]);
ylim([-0.05 1.1]);
ylabel('x_2');
xlabel('x_1');

%FIGURE 4(a)
figure(3)
scatter(rx, ry, markerSizes);
yticks([0 0.25 0.5 0.75 1]);
xticks([0 0.25 0.5 0.75 1]);
xlim([-0.05 1.3]);
ylim([-0.05 1.1]);
ylabel('x_2');
xlabel('x_1');
hold on;
scatter(sx(:,1), sx(:,2), markerSizes, 'b', 'filled');
hold on;
scatter(x1, x2, markerSizes, 'r', 'filled')
legend('PL', 'POL', 'GOL')

lossmax = round(((max(Gloss) + 0.5)/1))*1;
lossmin = floor(((min(Gloss))/1))*1;
minloss = min(Gloss);
[~, I] = min(Gloss);
textCell = arrayfun(@(x1,y1) sprintf('(%d, %3.3f)',x1,y1),I,minloss,'un',0);

%FIGURE 4(b)
figure(4)
plot(Iloss, Gloss);
hold on;
scatter(I, minloss, markerSizes, 'r', 'filled')
text(I+.2, minloss-.1,textCell,'FontSize',8) 
yticks(lossmin:1:lossmax);
xticks(0:40:K);
xlim([0 K+1]);
ylim([lossmin lossmax]);
