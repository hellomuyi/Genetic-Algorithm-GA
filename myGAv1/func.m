%测试函数Sphere
function fval = func(x)%一行为一个染色体，可能有多行
fval = sum(x .^2,  2);
end