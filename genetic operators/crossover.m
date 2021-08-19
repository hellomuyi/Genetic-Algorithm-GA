function newpopulation = crossover(population, pc)
newpopulation = population;
[psize, d] = size(population);
randmatrix = rand(psize, 1);
index = find(randmatrix<pc);   %获取待交叉的个体
for k = 1:2:(length(index)-1)
    location = randi(d-1);%随机获取交叉点的位置，[1, d-1]之间的整数
    newpopulation(index(k), 1:location) = population(index(k+1), 1:location);
    newpopulation(index(k+1), 1:location) = population(index(k), 1:location);
end
end