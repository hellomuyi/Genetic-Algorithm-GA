function fval = func(x)%һ��Ϊһ��Ⱦɫ�壬�����ж���
%fval = sum(x .^2, 2);
 fval = sum(x .^ 2 - 10 .* cos(2 .* pi .* x) + 10, 2);
end
%�˺������½��ȡ[-100;100]��ά��20or30

%{
function y=fun_ackley(x)

a = 20;
b = 0.2;
n = size(x,2);

y = a + exp(1) - a*exp(-b*(sum(x.^2, 2)./n).^0.5) - exp(sum(cos(2*pi.*x),2)./n);
y = y';
%}