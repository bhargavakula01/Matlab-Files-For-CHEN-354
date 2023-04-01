clear all;
close all;
options = optimoptions('fsolve', 'Display','none');


k = linspace(0,1,1000);
guess = [0.5; 0.5];
for i = 1:length(k)
    root = fsolve(@(w) sub6(w,k(i)), guess, options);
    guess = root;
end

disp(root);


function f = sub1(x)
    HCl = (5-4*x) / (6-x);
    v_HCl = -4;
    H2O = (2*x) / (6 - x);
    v_H2O = 2;
    Oxygen = (1 - x) / (6 - x);
    v_oxygen = -1;
    Cl = (2*x) / (6 - x);
    v_Cl = 2;

    f = HCl^(v_HCl) * H2O^(v_H2O) * Oxygen^(v_oxygen) * Cl^(v_Cl) - 2*7.1804;
end

function f = sub2(x, k)
    f = x * 0;

    ep = x(1);
    x1 = x(2);

    numerator = 0.24 * 0.24 * 0.12;

    f(1) = (numerator / ((1 - ep) / (1 + x1 + 2*ep))) - 0.3;
    f(2) = k*(((2*ep) / (1 + x1 + 2*ep)) - 0.24);
end


function f = sub3(w, k)
    y1 = w(1);
    y2 = 1 - y1;
    P = w(2);

    p1sat = 0.0923;
    p2sat = 1.431;

    x1_alpha = 0.0162;
    x2_alpha = 1 - x1_alpha;
    x1_beta = 0.587;
    x2_beta = 1 - x1_beta;

    gamma1_alpha = 100;
    gamma2_alpha = 1;

    f = w*0;
    f(1) = k*((y1*P) - (x1_alpha*gamma1_alpha*p1sat));
    f(2) = k*((y2*P) - (x2_alpha*1*p2sat));
end


function f = sub4(w, k)
    x1_alpha = w(1);
    x2_alpha = 1 - x1_alpha;
    y1 = w(2);
    y2 = 1 - y1;

    p1sat = 0.0923;
    p2sat = 1.431;
    P = 1.5;
    gamma1_alpha = 100;
    gamma2_alpha = 1;
    
    f = w * 0;
    f(1) = k*((y1*P) - (x1_alpha*gamma1_alpha*p1sat));
    f(2) = k*(x1_alpha*gamma1_alpha*p1sat + x2_alpha*gamma2_alpha*p2sat - P);
end


function f = sub5(w, k)
    x1_beta = 0.587;
    x2_beta = 1 - x1_beta;

    p1sat = 0.0923;
    p2sat = 1.431;
    P = 1.5573;

    y1 = 0.0960;
    y2 = 1 - y1;

    gamma1_beta = w(1);
    gamma2_beta = w(2);

    f = w*0;

    f(1) = k*(x1_beta*gamma1_beta*p1sat - y1*P);
    f(2) = k*(x2_beta*gamma2_beta*p2sat - y2*P);


end

function f = sub6(w, k)
    x1_beta = w(1);
    x2_beta = 1 - x1_beta;
    y1 = w(2);
    y2 = 1 - y1;

    p1sat = 0.0923;
    p2sat = 1.431;
    P = 0.5;
    gamma1_beta = 2.7593;
    gamma2_beta = 2.3821;
    
    f = w * 0;
    f(1) = k*((y1*P) - (x1_beta*gamma1_beta*p1sat));
    f(2) = k*(x1_beta*gamma1_beta*p1sat + x2_beta*gamma2_beta*p2sat - P);
end
















