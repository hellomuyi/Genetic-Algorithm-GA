%�����������
%populationΪԭ��Ⱥ��pmΪ������ʣ�luΪԼ�����ޣ�һ��������
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
            