%How to formulate a Tx1 ( temperature vs x1(composition) ) graph for LLE
clear all;
close all;

x = linspace(0,1,100);

for i = 1:length(x)
    y(i) = -1/ (x(i)*(1-x(i)));
    y2(i) = -4;
end
A = 2;
T_root1 = fsolve(@(T) func(T, A), 100000);
T_root2 = fsolve(@(T) func(T, A), 273.15);

disp(T_root1);
disp(T_root2);

%plot(x,y,x,y2);


function [f] = func(T, A)
    f = ((-1500) / T + 23.9 - (3)*log(T)) - A;
end



