%{
ʹ��ʾ����[x, y] = myGA(@func, 30, [-100; 100]);
*********************************************
                  ����˵��
*********************************************
    f    �����������
    D    ����������ά��
    lu   ����ά��������Ԫ�طֱ�Ϊ�������Ա������½���Ͻ� 
    x    ���������Ž⣨��������
    fbest�����Ž��Ӧ�ĺ���ֵ
%}
function [x, fbest] = myGA(f, D, lu)
%% ���ò���
PN = 50;            %��Ⱥ��ģ
G = 600;            %����������
Pc = 0.6;           %�������
Pm = 0.01;          %�������
%D = 30;            %ά�ȣ���������Ŀ 
L = lu(1, 1);       %�½�
U = lu(2, 1);       %�Ͻ�

%% ��ʼ��
rng('shuffle');
population = L + rand(PN, D)*(U-L); %�����ó�ʼ��Ⱥ
traceBest = zeros(G, 1);            %��¼����������Ŀ�꺯��ֵ
traceMean = zeros(G, 1);            %��¼������ƽ��Ŀ�꺯��ֵ
%% ��������
for gen = 1:G   
    fit = -f(population);  %������ȺĿ�꺯��ֵ��ȡ�෴����
    indexBest = find(max(fit)==fit);    
    pBest = population(indexBest(1), :); %��¼�Ŵ�����ǰ�����Ÿ���

    maxfit = max(fit);
    minfit = min(fit);
    fit = (fit-minfit+1e-7)/(maxfit-minfit+1e-6); %��һ���õ���Ӧ�ȣ�����Ǹ��� 
    
    %ѡ�����
    population = selection(population, fit);
    
    %�������
    population = mutation(population, Pm, lu);
    %�Ľ��������
    %population = new_mutation(population, Pm, lu, gen, G);
    
    %�������
    population = crossover(population, Pc);
    %�Ľ��������
    %population = new_crossover(f, population, Pc);
    
    %��Ӣ����ѡ����� 
    fit = -f(population);
    indexWorst = find(min(fit)==fit);   %��ȡ�������λ��
    population(indexWorst(1), :) = pBest;%�Ŵ�����ǰ�����Ÿ����滻�Ŵ��������������
    
    fit = -f(population);
    traceBest(gen, 1) = -max(fit);         %��¼��ǰ��������Ŀ�꺯��ֵ
    traceMean(gen, 1) = sum(-fit)./PN;     %��¼��ǰ����ƽ��Ŀ�꺯��ֵ
end
fbest = traceBest(G, 1);  %��������ֵ
fit = -f(population);
index = find(max(fit)==fit);
x = population(index(1), :);

%% ********************** ��ͼ ***************************
h = plotyy(1:G, traceMean(:, 1), 1:G, traceBest(:, 1)); 
    grid
    xlabel('��������');
    d1 = get(h(1), 'ylabel');
    set(d1, 'string', 'ƽ��Ŀ�꺯��ֵ');
    d2 = get(h(2), 'ylabel');
    set(d2, 'string', '����Ŀ�꺯��ֵ');
    legend('ƽ��Ŀ�꺯��ֵ', '����Ŀ�꺯��ֵ');
    title('ƽ��Ŀ�꺯��ֵ������Ŀ�꺯��ֵ����������ı仯');
end
%%
%{
****************************************************
            * ����ΪGA�������Ӻ��� *
****************************************************
%}
%% ****************** ѡ����� ***********************
%populationΪԭ��Ⱥ��fitΪ��Ⱥ��һ���õ�����Ӧ��
function newpopulation = selection(population, fit)
%�������ѡ�����
sumfit = sum(fit);
fitvalue = fit/sumfit; 
%�����ۻ�����
fitcum = cumsum(fitvalue);
fiti = 1;
newi = 1;
[m, n] = size(population);   %��ȡ��Ⱥ��ģ
newpopulation = zeros(m, n);
randmatrix = sort(rand(m, 1));
while newi <= m
    if randmatrix(newi)<fitcum(fiti)
        newpopulation(newi, :) = population(fiti, :);
        newi = newi+1;
    else
        fiti = fiti+1;
    end
end
end

%% ****************** ������� *********************
%populationΪԭ��Ⱥ��pcΪ�������
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

%% ****************** ������� **********************
%populationΪԭ��Ⱥ��pmΪ������ʣ�luΪԼ�����ޣ���������
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

%% **************** �Ľ�������� ********************  %�Լ�����ʵ�֣���ǰ����ǿǿ���
%fΪ�����������populationΪԭ��Ⱥ��pcΪ�������
function newpopulation = new_crossover(f, population, pc)
newpopulation = population;
[psize, d] = size(population);
randmatrix = rand(psize, 1);
index = find(randmatrix<pc);   %��ȡ������ĸ���
%����Ϊǿǿ��ϵĸĽ�����
fit = f(population);
fitselected = fit(index);    %����������Ŀ�꺯��ֵ
[~, ptr] = sort(fitselected); %Ŀ�꺯��ֵ����
index = index(ptr);
for k = 1:2:(length(index)-1)
    location = randi(d-1);%�����ȡ������λ�ã�[1, d-1]֮�������
    newpopulation(index(k), 1:location) = population(index(k+1), 1:location);
    newpopulation(index(k+1), 1:location) = population(index(k), 1:location);
end
end

%% ****************** �Ľ�������� ********************** %�Լ�����ʵ�֣���ǰ�������ԷǾ��ȱ���
%populationΪԭ��Ⱥ��pmΪ������ʣ�luΪԼ�����ޣ���������
function newpopulation = new_mutation(population, pm, lu, gen, G)
[psize, d] = size(population);
newpopulation = population;
for m = 1:psize
    for n = 1:d
        if rand<pm
            newpopulation(m, n) = lu(1) + rem(population(m, n)-lu(1)+...
                (-1+2*rand(1))*nonuniform(gen, lu, G)+lu(2)-lu(1), lu(2)-lu(1));
        end
    end
end
end

%% ********* �Ǿ��ȱ��죬���ػ���x�����ƫ��λ�� ***********
%����gen��ǰ������luΪ��Χ����������GΪ����������
function gen = nonuniform(gen, lu, G)
precision = 0.01;   %��⾫��         
k = (lu(2)-lu(1)-precision)/(1-G);    %б��
gen = k*(gen-G) + precision;
end