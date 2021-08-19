%{
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
global G;
global Pc;
global Pm;
PN = 50;           %种群规模
G = 1000;          %最大进化代数
Pc = 0.6;          %交叉概率
Pm = 0.01;         %变异概率
%D = 30;           %维度，即基因数目
L = lu(1, 1);      %下界
U = lu(2, 1);      %上界
%% 初始化
rng('default')
rng('shuffle');
population = L + rand(PN, D)*(U-L); %随机获得初始种群
traceBest = zeros(G, 1);            %记录历代的最优目标函数值
traceMean = zeros(G, 1);            %记录历代的平均目标函数值
%% 进化过程
for generation = 1:G
    fit = -f(population);  %计算种群目标函数值（取相反数）
    indexBest = find(max(fit)==fit);
    pBest = population(indexBest(1), :); %记录遗传操作前的最优个体
    
    maxfit = max(fit);
    minfit = min(fit);
    fit = (fit-minfit+1e-7)/(maxfit-minfit+1e-6); %归一化得到适应度，满足非负性
    beforePopulation = population;
    %选择操作
    population = selection(population, fit);
    
    %变异操作
    population = mutation(population, lu);
    %改进变异操作
    %population = new_mutation(population, lu, generation);
    
    %交叉操作
    population =crossover(population);
    %改进交叉操作
    %population = new_crossover(f, population, generation);
    
    %模拟退火
    %population = mySA(population, beforePopulation, generation, f);
    
    %精英个体选择策略
    fit = -f(population);
    indexWorst = find(min(fit)==fit);   %获取最差个体的位置
    population(indexWorst(1), :) = pBest;%遗传操作前的最优个体替换遗传操作后的最差个体
    
    fit = -f(population);
    traceBest(generation, 1) = -max(fit);         %记录当前代的最优目标函数值
    traceMean(generation, 1) = sum(-fit)./PN;     %记录当前代的平均目标函数值
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
while newi<=m
    if randmatrix(newi)<fitcum(fiti)
        newpopulation(newi, :) = population(fiti, :);
        newi = newi+1;
    else
        fiti = fiti+1;
    end
end
end
%% ****************** 交叉操作 *********************
%population为原种群
function newpopulation = crossover(population)
global Pc;
newpopulation = population;
[psize, d] = size(population);
randmatrix = rand(psize, 1);
index = find(randmatrix<Pc);   %获取待交叉的个体
for k = 1:2:(length(index)-1)
    location = randi(d-1);%随机获取交叉点的位置，[1, d-1]之间的整数
    newpopulation(index(k), 1:location) = population(index(k+1), 1:location);
    newpopulation(index(k+1), 1:location) = population(index(k), 1:location);
end
end
%% ****************** 变异操作 **********************
%population为原种群，lu为约束界限（列向量）
function newpopulation = mutation(population, lu)
global Pm;
[psize, d] = size(population);
newpopulation = population;
for m = 1:psize
    for n = 1:d
        if rand<Pm
            newpopulation(m, n) = lu(1)+rand*(lu(2)-lu(1));
        end
    end
end
end
%% **************** 改进交叉操作 ********************
%f为待求函数句柄，population为原种群，generation为当前代数
function newpopulation = new_crossover(f, population, generation)
global Pc;
global G;
if generation < G/4
    PcCur = Pc + (1-Pc)/2;
elseif generation < G*0.75
    PcCur = Pc;
else
    PcCur = Pc - (1-Pc)/2;
end
newpopulation = population;
[psize, d] = size(population);
randmatrix = rand(psize, 1);
index = find(randmatrix<PcCur);   %获取待交叉的个体
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
%% ****************** 改进变异操作 **********************
%population为原种群，lu为约束界限（列向量），generation为当前代数
function newpopulation = new_mutation(population, lu, generation)
global Pm;
global G;
if generation < G/4
    PmCur = Pm - Pm/2;
elseif generation < G*0.75
    PmCur = Pm;
else
    PmCur = Pm + Pm/2;
end
[psize, d] = size(population);
newpopulation = population;
for m = 1:psize
    for n = 1:d
        if rand<PmCur
            newpopulation(m, n) = lu(1) + rem(population(m, n)-lu(1)+...
                (-1+2*rand(1))*nonuniform(generation, lu)+lu(2)-lu(1), lu(2)-lu(1));
        end
    end
end
end
%% ********* 非均匀变异，返回基因x变异的偏移位置 ***********
%参数generation当前代数，lu为范围的列向量
function generation = nonuniform(generation, lu)
global G;
precision = 0.01;   %求解精度
k = (lu(2)-lu(1)-precision)/(1-G);    %斜率
generation = k*(generation-G) + precision;
end
%% ******************** 模拟退火操作 **********************
%newPopulation为当前代遗传操作后的新种群，population为当前代遗传操作前的种群
%generation为当前代数，f为待求函数句柄
function newpopulation = mySA(newPopulation, population, generation, f)
global G;
t = (G-generation+1)/G*150;
newpopulation = newPopulation;
fitNew = f(newPopulation);  %计算函数值
fit = f(population);
[p, ~] = find(fitNew>fit);       %须用准则判断的下标

c = exp((fit(p)-fitNew(p))./t);
randmatrix = rand(length(p), 1);
index = find(randmatrix<c);     %进一步满足准则的下标（对p）
newpopulation(p(index), :) = population(p(index), :);
end