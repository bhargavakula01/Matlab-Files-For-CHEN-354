clear all;
close all;
options = optimoptions('fsolve', 'Display','none');

k = linspace(0,1,1000);
guess = [0.5, 1];
for i = 1:length(k)
    root = fsolve(@(w) sub1(w, k(i)), guess, options);
    guess = root;
end

disp(root);

function f = sub1(w, k)
    p1sat = 15;
    p2sat = 8;

    gamma1_alpha = 1;
    gamma2_beta = 1;

    x1_beta = 0;
    x2_beta = 1- x1_beta;

    x1_alpha = 1;
    x2_alpha = 1 - x1_alpha;

    y1 = w(1);
    y2 = 1 - y1;
    P = w(2);

    f = w*0;
    f(1) = k*((y1*P) - (x1_alpha*gamma1_alpha*p1sat));
    f(2) = k*((y2*P) - (x2_beta*gamma2_beta*p2sat));
end
