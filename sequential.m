% ============= Sequential Design ========================%
function [freqP, UIndex, n, Floss] = sequential(Index, IIndex, n1, rn, ICN, IDesign, IFloss, N, a, b, c, sigma_epsilon, IN, F)
UIndex = IIndex;
freqP = IDesign;
Floss = IFloss;
n = n1;

for j = 1:rn
    AIndex = setdiff(Index, UIndex); % Active Index
    n = n + 1;
    CN = ICN;
    na = length(AIndex);
    Aloss = zeros(1, na); %To strore losses

    for i = 1:na
        AfreqP = freqP; % To choose the next design
        AfreqP(AIndex(i)) = 1;
        Aloss(i) = plossnew(N, n, AfreqP, a, b, c, sigma_epsilon, CN, IN, F);
    end

    [~, I] = min(Aloss);

    freqP(AIndex(I)) = 1;
    UIndex = [UIndex' AIndex(I)]';
    Floss(1,j+1) = plossnew(N, n, freqP, a, b, c, sigma_epsilon, CN, IN, F);
end

