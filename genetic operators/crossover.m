function newpopulation = crossover(population, pc)
newpopulation = population;
[psize, d] = size(population);
randmatrix = rand(psize, 1);
index = find(randmatrix<pc);   %��ȡ������ĸ���
for k = 1:2:(length(index)-1)
    location = randi(d-1);%�����ȡ������λ�ã�[1, d-1]֮�������
    newpopulation(index(k), 1:location) = population(index(k+1), 1:location);
    newpopulation(index(k+1), 1:location) = population(index(k), 1:location);
end
end