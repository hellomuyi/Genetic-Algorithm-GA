function newpopulation = selection(population, fit)
%�������ѡ�����
sumfit = sum(fit);
fitvalue = fit./sumfit;    
%�����ۻ�����
fitcum = cumsum(fitvalue); 
fiti = 1;
newi = 1;
[m, n] = size(population);   %��ȡ��Ⱥ��ģ
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