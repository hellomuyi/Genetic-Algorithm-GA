%{
使用示例：[x, y] = myGA(@func, 30, [-100; 100]);
*********************************************
                  参数说明
*********************************************
    f    ：待求函数句柄
    D    ：待求函数的维度
    lu   ：二维列向量，元素分别为待求函数自变量的下界和上界 
    x    ：所求最优解（行向量）
    fbest：最优解对应的函数值
%}
function [x, fbest] = myGA(f, D, lu)
%% 设置参数
PN = 50;            %种群规模
G = 600;            %最大进化代数
Pc = 0.6;           %交叉概率
Pm = 0.01;          %变异概率
%D = 30;            %维度，即基因数目 
L = lu(1, 1);       %下界
U = lu(2, 1);       %上界

%% 初始化
rng('shuffle');
population = L + rand(PN, D)*(U-L); %随机获得初始种群
traceBest = zeros(G, 1);            %记录历代的最优目标函数值
traceMean = zeros(G, 1);            %记录历代的平均目标函数值
%% 进化过程
for gen = 1:G   
    fit = -f(population);  %计算种群目标函数值（取相反数）
    indexBest = find(max(fit)==fit);    
    pBest = population(indexBest(1), :); %记录遗传操作前的最优个体

    maxfit = max(fit);
    minfit = min(fit);
    fit = (fit-minfit+1e-7)/(maxfit-minfit+1e-6); %归一化得到适应度，满足非负性 
    
    %选择操作
    population = selection(population, fit);
    
    %变异操作
    population = mutation(population, Pm, lu);
    %改进变异操作
    %population = new_mutation(population, Pm, lu, gen, G);
    
    %交叉操作
    population = crossover(population, Pc);
    %改进交叉操作
    %population = new_crossover(f, population, Pc);
    
    %精英个体选择策略 
    fit = -f(population);
    indexWorst = find(min(fit)==fit);   %获取最差个体的位置
    population(indexWorst(1), :) = pBest;%遗传操作前的最优个体替换遗传操作后的最差个体
    
    fit = -f(population);
    traceBest(gen, 1) = -max(fit);         %记录当前代的最优目标函数值
    traceMean(gen, 1) = sum(-fit)./PN;     %记录当前代的平均目标函数值
end
fbest = traceBest(G, 1);  %所得最优值
fit = -f(population);
index = find(max(fit)==fit);
x = population(index(1), :);

%% ********************** 绘图 ***************************
h = plotyy(1:G, traceMean(:, 1), 1:G, traceBest(:, 1)); 
    grid
    xlabel('进化代数');
    d1 = get(h(1), 'ylabel');
    set(d1, 'string', '平均目标函数值');
    d2 = get(h(2), 'ylabel');
    set(d2, 'string', '最优目标函数值');
    legend('平均目标函数值', '最优目标函数值');
    title('平均目标函数值和最优目标函数值随进化代数的变化');
end
%%
%{
****************************************************
            * 以下为GA函数的子函数 *
****************************************************
%}
%% ****************** 选择操作 ***********************
%population为原种群，fit为种群归一化得到的适应度
function newpopulation = selection(population, fit)
%计算个体选择概率
sumfit = sum(fit);
fitvalue = fit/sumfit; 
%计算累积概率
fitcum = cumsum(fitvalue);
fiti = 1;
newi = 1;
[m, n] = size(population);   %获取种群规模
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

%% ****************** 交叉操作 *********************
%population为原种群，pc为交叉概率
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

%% ****************** 变异操作 **********************
%population为原种群，pm为变异概率，lu为约束界限（列向量）
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

%% **************** 改进交叉操作 ********************  %自己创新实现，当前采用强强组合
%f为待求函数句柄，population为原种群，pc为交叉概率
function newpopulation = new_crossover(f, population, pc)
newpopulation = population;
[psize, d] = size(population);
randmatrix = rand(psize, 1);
index = find(randmatrix<pc);   %获取待交叉的个体
%以下为强强组合的改进部分
fit = f(population);
fitselected = fit(index);    %待交叉个体的目标函数值
[~, ptr] = sort(fitselected); %目标函数值排序
index = index(ptr);
for k = 1:2:(length(index)-1)
    location = randi(d-1);%随机获取交叉点的位置，[1, d-1]之间的整数
    newpopulation(index(k), 1:location) = population(index(k+1), 1:location);
    newpopulation(index(k+1), 1:location) = population(index(k), 1:location);
end
end

%% ****************** 改进变异操作 ********************** %自己创新实现，当前采用线性非均匀变异
%population为原种群，pm为变异概率，lu为约束界限（列向量）
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

%% ********* 非均匀变异，返回基因x变异的偏移位置 ***********
%参数gen当前代数，lu为范围的列向量，G为最大进化代数
function gen = nonuniform(gen, lu, G)
precision = 0.01;   %求解精度         
k = (lu(2)-lu(1)-precision)/(1-G);    %斜率
gen = k*(gen-G) + precision;
end