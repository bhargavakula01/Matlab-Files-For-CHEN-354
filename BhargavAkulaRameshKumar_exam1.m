clear all;
close all;
main();


function y = main()
    disp('1: fugacity of pure liquid');
    disp('2: fugacity of pure species');
    disp('3: Rackett Equation');
    disp('4: Fugacity Coefficient of Pure Species');
    disp('5: Fugacity Coefficient of Binary System');
    disp('6: B value');
    disp('7: B_hat value');
    disp('8: T_cij Calculator');
    disp('9: Z_cij calculator');
    disp('10: w_ij calculator');
    disp('11: V_cij Calculator');
    disp('12: P_cij Calculator');
    disp('13: Linear Intrapolation');
    disp('14: Antoines Equation Calculator');
    user_input = input('Enter number for function: ');
    if(user_input == 1)
        P = input('Enter value for pressure: ');
        Pc = input('Enter value for critical pressure: ');
        Psat = input('Enter value for saturation pressure: ');
        T = input("Enter value for temperature: ");
        Tc = input('Enter value for critical temperature: ');
        w = input('Enter value for accentric factor: ');
        Z_c = input('Enter value for Zc: ');
        V_c = input('Enter value of Vc: ');
        R = input('Enter value for R: ');
        fi =fugacity_pure_liquid(P,Pc, Psat, T, Tc, R, w, Z_c, V_c);
        disp('Pure Species (Lquid) value: ' + string(fi));
    elseif(user_input == 2)
        Pc = input('Enter value for critical pressure: ');
        Tc = input('Enter value for critical temperature: ');
        w = input('Enter value for accentric factor: ');
        P = input('Enter value for pressure: ');
        T = input("Enter value for temperature: ");
        fi = fugacity_pure_species(Pc, Tc, w, P, T);
        disp('Fugacity Pure Species: ' + string(fi));
    elseif(user_input == 3)
        Z_c = input('Enter value for Zc: ');
        V_c = input('Enter value of Vc: ');
        Tc = input('Enter value for critical temperature: ');
        T = input("Enter value for temperature: ");
        V_sat = rackett_equation(Z_c, V_c, Tc, T);
        disp('Molar Volume at Saturation Pt: ' + string(V_sat));
    elseif(user_input == 4)
        Pc = input('Enter value for critical pressure: ');
        Tc = input('Enter value for critical temperature: ');
        w = input('Enter value for accentric factor: ');
        P = input('Enter value for pressure: ');
        T = input("Enter value for temperature: ");
        phi_i = fug_coefficient_pure_species(Pc, Tc, w, P, T);
        disp('Fugacity Coefficient: ' + string(phi_i));
    elseif(user_input == 5)
        P = input("Enter value for Pressure: ");
        T = input("Enter value for Absolute Temperature: ");
        R = input('Enter value for R: ');
        T1 = input("Enter critical temperature for comp 1: ");
        T2 = input("Enter critical temperature for comp 2: ");
        w1 = input("Enter accentric factor for comp 1: ");
        w2 = input("Enter accentric factor for comp 2: ");
        Z1 = input("Enter compressibility factor for comp 1: ");
        Z2 = input("Enter compressibility factor for comp 2: ");
        V1 = input("Enter Vc for comp 1: ");
        V2 = input("Enter Vc for comp 2: ");
        k = input('Enter value for constant k: ');
        y1 = input('Enter value for y1: ');
        y2 = input('Enter value for y2: ');
        matrix_ans = fug_coefficient_mixture(P,T,R,T1,T2,w1,w2,Z1,Z2,V1,V2,k,y1,y2);
        disp('Phi_1: ' + string(matrix_ans(1)));
        disp('Phi_2: ' + string(matrix_ans(2)));
        disp('f_1: ' + string(matrix_ans(3)));
        disp('f_2: ' + string(matrix_ans(4)));
    elseif(user_input == 6)
        Pcij = input("Enter value for Pcij: ");
        Tcij = input('Enter value for Tcij: ');
        R = input('Enter value for R: ');
        T = input("Enter value for Absolute Temperature: ");
        w = input('Enter value for accentric factor: ');
        B_ij = b_ij(Pcij, Tcij, R, T, w);
        disp('Value for B_ij: ' + string(B_ij));
    elseif(user_input == 7)
        T = input("Enter value for Absolute Temperature: ");
        Tc = input('Enter value for critical temperature: ');
        w = input('Enter value for accentric factor: ');
        b_hat = B_Hat_ij_calc(T, Tc, w);
        disp('b_hat value: ' + string(b_hat));
    elseif(user_input == 8)
        Tci = input('Enter value for Tci: ');
        Tcj = input('Enter value for Tcj: ');
        k = input('Enter value for k: ');
        T_c = Tcij_calc(Tci, Tcj,k);
        disp('Tcij is: ' + string(T_c));
    elseif(user_input == 9)
        Zci = input('Enter value for Zci: ');
        Zcj = input('Enter value for Zcj: ');
        Z_c = Zcij_calc(Zci, Zcj);
        disp('Z_cij: ' + stirng(Z_c));
    elseif(user_input == 10)
        wi = input('Enter value for wi: ');
        wj = input('Enter value for wj: ');
        w = wij_calc(wi, wj);
        disp('wij: ' + string(w));
    elseif(user_input == 11)
        Vci = input('Enter value for Vci: ');
        Vcj = input('Enter value for Vcj: ');
        V_c = Vcij_calc(Vci, Vcj);
        disp('V_cij: ' + string(V_c));
    elseif(user_input == 12)
        Zcij = ('Enter value for Zcij: ');
        Tcij = ('Enter value for Tcij: ');
        Vcij = ('Enter value for Vcij: ');
        R = input('Enter value for R: ');
        P_c = Pcij_calc(Zcij, Tcij, Vcij, R);
        disp('P_cij value: ' + string(P_c));
    elseif(user_input == 13)
        x1 = input('Enter x1: ');
        y1 = input('Enter y1: ');
        x2 = input('Enter x2: ');
        y2 = input('Enter y2: ');
        x = input('Enter x value: ');
        y = linear_intrapolation(x1,x2, y1,y2, x);
        disp('The value for y based on x is: ' + string(y));
    elseif(user_input == 14)
        a = input('Enter value for A: ');
        b = input('Enter value for B: ');
        c = input('Enter value for C: ');
        t = input('Enter value for temperature (Degrees C): ');
        psat = Antoines_Equ(a, b, c, t);
        disp('The value for Psat using Antoines is: ' + string(psat));
    end
end

function f_i = fugacity_pure_liquid(P,Pc, Psat, T, Tc, R, w, Z_c, V_c)
    phi_sat = fug_coefficient_pure_species(Pc, Tc, w, Psat, T);
    V_sat = rackett_equation(Z_c, V_c, Tc, T);
    power_value = (V_sat*(P - Psat)) / (R*T);
    f_i = (phi_sat)*(Psat)*exp(power_value);
end

function f_i = fugacity_pure_species(Pc, Tc, w, P, T)
    phi = fug_coefficient_pure_species(Pc,Tc,w,P,T);
    f_i = phi * P;
end

function V_sat = rackett_equation(Z_c, V_c, Tc, T)
    reduced_temperature = T/Tc;
    power_val = (1-reduced_temperature)^(2/7);
    V_sat = V_c * (Z_c)^(power_val);
end

function phi_i = fug_coefficient_pure_species(Pc, Tc, w, P, T)
    reduced_pressure = P / Pc;
    reduced_temperature = T/Tc;
    B_zero = 0.083 - (0.422 / (reduced_temperature)^1.6);
    B_one = 0.139 - (0.172 / (reduced_temperature)^4.2);
    phi_i = exp((reduced_pressure/reduced_temperature) *(B_zero + (w)*(B_one)));
end

function phi_f = fug_coefficient_mixture(P,T,R,T1,T2,w1,w2,Z1,Z2,V1,V2,k,y1,y2)
    for i = 1:3
        if(i == 1)
            w_list(i) = wij_calc(w1,w1);
            T_list(i) = Tcij_calc(T1, T1, k);
            Z_list(i) = Zcij_calc(Z1, Z1);
            V_list(i) = Vcij_calc(V1, V1);
            P_list(i) = Pcij_calc(Z_list(i), T_list(i), V_list(i), R);
        elseif(i == 2)
            w_list(i) = wij_calc(w1,w2);
            T_list(i) = Tcij_calc(T1, T2, k);
            Z_list(i) = Zcij_calc(Z1, Z2);
            V_list(i) = Vcij_calc(V1, V2);
            P_list(i) = Pcij_calc(Z_list(i), T_list(i), V_list(i), R);
        elseif(i == 3)
            w_list(i) = wij_calc(w2,w2);
            T_list(i) = Tcij_calc(T2, T2, k);
            Z_list(i) = Zcij_calc(Z2, Z2);
            V_list(i) = Vcij_calc(V2, V2);
            P_list(i) = Pcij_calc(Z_list(i), T_list(i), V_list(i), R);
        end
    end

    for j = 1:3
        B(j) = b_ij(P_list(j), T_list(j),R, T, w_list(j));
    end
    
    sigma_ij = 2*(B(2)) - B(1) - B(3);
    phi_1 = exp((P/(R*T))*((B(1))+(y2*y2)*sigma_ij));
    phi_2 = exp((P/(R*T))*((B(3))+(y1*y1)*sigma_ij));

    f1 = P * phi_1 * y1;
    f2 = P * phi_2 * y2;

    phi_f = [phi_1, phi_2, f1,f2];

end

function B_ij = b_ij(Pcij, Tcij, R, T, w)
    B_Hat_ij = B_Hat_ij_calc(T, Tcij, w);
    B_ij = (B_Hat_ij * R * Tcij) / Pcij;
end

function b_hat = B_Hat_ij_calc(T, Tc, w)
    reduced_temperature = T/Tc;
    B_zero = 0.083 - (0.422 / (reduced_temperature)^1.6);
    B_one = 0.139 - (0.172 / (reduced_temperature)^4.2);
    b_hat = B_zero + (w)*(B_one);
end

function T_c = Tcij_calc(Tci, Tcj,k)
    p_val = 1/2;
    T_c = ((Tci *Tcj)^(p_val)) *(1-k);
end

function Z_c = Zcij_calc(Zci, Zcj)
    Z_c = (Zci + Zcj) / 2;
end

function w = wij_calc(wi, wj)
    w = (wi + wj) / 2;
end

function V_c = Vcij_calc(Vci, Vcj)
    power_frac = 1/3;
    Vci_cal = (Vci)^(power_frac);
    Vcj_cal = (Vcj)^(power_frac);
    V_c = ((Vci_cal + Vcj_cal) / 2)^3;
end

function P_c = Pcij_calc(Zcij, Tcij, Vcij, R)
    P_c = (Zcij*R*Tcij) / (Vcij);
end

function y = linear_intrapolation(x1,x2, y1,y2, x)
    y = y1 + (((x - x1) * (y2-y1)) / (x2 - x1));
end

function psat = Antoines_Equ(a,b,c, t)
    psat = exp(a - (b / (t + c)));
end