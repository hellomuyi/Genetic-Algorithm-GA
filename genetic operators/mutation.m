%变异操作函数
%population为原种群，pm为变异概率，lu为约束界限（一列向量）
function newpopulation = mutation(population, pm, lu)
[psize, d] = size(population);
newpopulation = population;
for m = 1:psize
    for n = 1:d
        if rand<pm
            newpopulation(m, n) = lu(1)+rand*(lu(2)-lu(1));
        end
    end
end
end
            