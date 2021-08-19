function newpopulation = selection(population, fit)
%计算个体选择概率
sumfit = sum(fit);
fitvalue = fit./sumfit;    
%计算累积概率
fitcum = cumsum(fitvalue); 
fiti = 1;
newi = 1;
[m, n] = size(population);   %获取种群规模
newpopulation = zeros(m, n);
randmatrix = sort(rand(m, 1));
while newi<=m
    if randmatrix(newi)<fitcum(fiti)
        newpopulation(newi, :) = population(fiti, :);
        newi = newi+1;
    else
        fiti = fiti+1;
    end
end
end