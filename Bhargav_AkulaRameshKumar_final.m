clear all;
close all;
options = optimoptions('fsolve', 'Display','none');

disp('1 - Finding K using -G divided by RT method');
disp('2 - Finding K using K0*K1*K2 method');
disp('3 - Determining temperature using K value');
disp('4 - Determining the Reaction Coordinate value(General Gas Phase)');
disp('5 - Determining the Reaction Coordinate value(General liquid phase)');
disp('6 - Determining the Reaction Coordinate value(Ideal liquid phase)');
disp('7 - Determining temperature for SLE case 2');
disp('8 - Determining Eutectic Temperature and composition (SLE)');
user_choice = input('What function do you want to perform?: ');

%Percent conversion is 
%(how much was used up / initial amount of component) * 100!!!

if(user_choice == 1)
    delt_A = 0;
    delt_B = 0;
    delt_C = 0;
    delt_D = 0;
    delt_Gibbs = 0;
    delt_Enthalpy = 0;
    number_components = input('How many components do you have?: ');
    T = input('Enter value for temperature (Kelvin): ');
    T_ref = input('Enter value for reference temperature (Kelvin): ');
    R = input('Enter value for the gas constant: ');
    
    for i = 1:number_components
        comp_type = input('Select Reactant (0) or Product (1): ');
        stoich_coeff = input('Type abs value for stoichiometric coefficient: ');
        A = input('value for A in table C.1 for component: ');
        B = input('value for B in table C.1 for component: ');
        C = input('value for C in table C.1 for component: ');
        D = input('value for D in table C.1 for component: ');
        G_ref = input('Enter value for Gibbs at 298 K (Table C.4): ');
        H_ref = input('Enter value for Enthalpy at 298 K (Table C.4): ');
    
        if(comp_type == 0)
            delt_A = delt_A - stoich_coeff*A;
            delt_B = delt_B - stoich_coeff*B;
            delt_C = delt_C - stoich_coeff*C;
            delt_D = delt_D - stoich_coeff*D;
            delt_Gibbs = delt_Gibbs - stoich_coeff*G_ref;
            delt_Enthalpy = delt_Enthalpy - stoich_coeff*H_ref;
        elseif(comp_type == 1)
            delt_A = delt_A + stoich_coeff*A;
            delt_B = delt_B + stoich_coeff*B;
            delt_C = delt_C + stoich_coeff*C;
            delt_D = delt_D + stoich_coeff*D;   
            delt_Gibbs = delt_Gibbs + stoich_coeff*G_ref;
            delt_Enthalpy = delt_Enthalpy + stoich_coeff*H_ref;
        end
    end
    [DeltGRT, K] = determine_K_method_1(delt_A, delt_B, delt_C, delt_D,delt_Gibbs, delt_Enthalpy, T, T_ref, R);
    disp('Delta G Divided By RT: ' + string(DeltGRT));
    disp('Equilibrium constant: ' + string(K));
elseif(user_choice == 2)
    delt_A = 0;
    delt_B = 0;
    delt_C = 0;
    delt_D = 0;
    delt_Gibbs = 0;
    delt_Enthalpy = 0;
    number_components = input('How many components do you have?: ');
    T = input('Enter value for temperature (Kelvin): ');
    T_ref = input('Enter value for reference temperature (Kelvin): ');
    R = input('Enter value for the gas constant: ');
    
    for i = 1:number_components
        comp_type = input('Select Reactant (0) or Product (1): ');
        stoich_coeff = input('Type abs value for stoichiometric coefficient: ');
        A = input('value for A in table C.1 for component: ');
        B = input('value for B in table C.1 for component: ');
        C = input('value for C in table C.1 for component: ');
        D = input('value for D in table C.1 for component: ');
        G_ref = input('Enter value for Gibbs at 298 K (Table C.4): ');
        H_ref = input('Enter value for Enthalpy at 298 K (Table C.4): ');
    
        if(comp_type == 0)
            delt_A = delt_A - stoich_coeff*A;
            delt_B = delt_B - stoich_coeff*B;
            delt_C = delt_C - stoich_coeff*C;
            delt_D = delt_D - stoich_coeff*D;
            delt_Gibbs = delt_Gibbs - stoich_coeff*G_ref;
            delt_Enthalpy = delt_Enthalpy - stoich_coeff*H_ref;
        elseif(comp_type == 1)
            delt_A = delt_A + stoich_coeff*A;
            delt_B = delt_B + stoich_coeff*B;
            delt_C = delt_C + stoich_coeff*C;
            delt_D = delt_D + stoich_coeff*D;   
            delt_Gibbs = delt_Gibbs + stoich_coeff*G_ref;
            delt_Enthalpy = delt_Enthalpy + stoich_coeff*H_ref;
        end
    end
    K = three_K_method(delt_A, delt_B, delt_C, delt_D,delt_Gibbs, delt_Enthalpy, T, T_ref, R);
    disp('Equilibrium constant: ' + string(K));
elseif(user_choice == 3)
    delt_A = 0;
    delt_B = 0;
    delt_C = 0;
    delt_D = 0;
    delt_Gibbs = 0;
    delt_Enthalpy = 0;
    number_components = input('How many components do you have?: ');
    T_ref = input('Enter value for reference temperature (Kelvin): ');
    R = input('Enter value for the gas constant: ');
    K = input('What is the Equilibrium constant for this reaction?: ');
    
    for i = 1:number_components
        comp_type = input('Select Reactant (0) or Product (1): ');
        stoich_coeff = input('Type abs value for stoichiometric coefficient: ');
        A = input('value for A in table C.1 for component: ');
        B = input('value for B in table C.1 for component: ');
        C = input('value for C in table C.1 for component: ');
        D = input('value for D in table C.1 for component: ');
        G_ref = input('Enter value for Gibbs at 298 K (Table C.4): ');
        H_ref = input('Enter value for Enthalpy at 298 K (Table C.4): ');
    
        if(comp_type == 0)
            delt_A = delt_A - stoich_coeff*A;
            delt_B = delt_B - stoich_coeff*B;
            delt_C = delt_C - stoich_coeff*C;
            delt_D = delt_D - stoich_coeff*D;
            delt_Gibbs = delt_Gibbs - stoich_coeff*G_ref;
            delt_Enthalpy = delt_Enthalpy - stoich_coeff*H_ref;
        elseif(comp_type == 1)
            delt_A = delt_A + stoich_coeff*A;
            delt_B = delt_B + stoich_coeff*B;
            delt_C = delt_C + stoich_coeff*C;
            delt_D = delt_D + stoich_coeff*D;   
            delt_Gibbs = delt_Gibbs + stoich_coeff*G_ref;
            delt_Enthalpy = delt_Enthalpy + stoich_coeff*H_ref;
        end
    end
    initial_guess = 1000;
    k = linspace(0,1,1000);
    for i = 1:length(k)
        temp_root = fzero(@(T) determine_Temp(delt_A, delt_B, delt_C, delt_D,delt_Gibbs, delt_Enthalpy, T, T_ref, R, K, k(i)), initial_guess);
        initial_guess = temp_root;
    end
    disp('Temperature to get K value is: ' + string(temp_root));
elseif(user_choice == 4)
    K = input('What is the Equilibrium constant for this reaction?: ');
    P = input('What is the pressure of this reaction: ');
    P_ref = 1; %bar
    v = input('What is the overall stoichiometric coefficient of the rxn?: ');
    total_moles = input('Total num moles for rxn: ');

    number_components = input('How many components do you have?: ');
    component_matrix = zeros(number_components, 3);
    for i = 1:number_components
        comp_type = input('Select Reactant (0) or Product (1): ');
        stoich_coeff = input('Type abs value for stoichiometric coefficient: ');
        component_matrix(i, 2) = input('How many moles did you start with for this component?: ');
        component_matrix(i, 3) = input('What is the phi value for this component');%phi_hat = 1 for ideal gas rxns!!
        if(comp_type == 0)
            component_matrix(i, 1) = -stoich_coeff;
        else
            component_matrix(i, 1) = stoich_coeff;
        end
    end
    epsilon = fzero(@(x) determine_epsilon_gas(K, P, P_ref, v, total_moles, component_matrix, x, number_components), 0.5);
    disp('Epsilon value for the rxn: ' + string(epsilon));
elseif(user_choice == 5)
    K = input('What is the Equilibrium constant for this reaction?: ');
    P = input('What is the pressure of this reaction: ');
    P_ref = 1; %bar
    T = input('Enter value for temperature (Kelvin): ');
    v = input('What is the overall stoichiometric coefficient of the rxn?: ');
    total_moles = input('Total num moles for rxn: ');

    number_components = input('How many components do you have?: ');
    component_matrix = zeros(number_components, 4);
    for i = 1:number_components
        comp_type = input('Select Reactant (0) or Product (1): ');
        stoich_coeff = input('Type abs value for stoichiometric coefficient: ');
        component_matrix(i, 2) = input('How many moles did you start with for this component?: ');
        component_matrix(i, 3) = input('What is the gamma value for this component');
        component_matrix(i, 4) = input('What is the sat volume?: ');
        if(comp_type == 0)
            component_matrix(i, 1) = -stoich_coeff;
        else
            component_matrix(i, 1) = stoich_coeff;
        end
    end
    epsilon = fzero(@(x) determine_epsilon_liquid(K, P, P_ref,T, v, total_moles, component_matrix, x, number_components), 0.5);
    disp('Epsilon value for the rxn: ' + string(epsilon));
elseif(user_choice == 6)
    K = input('What is the Equilibrium constant for this reaction?: ');
    v = input('What is the overall stoichiometric coefficient of the rxn?: ');
    total_moles = input('Total num moles for rxn: ');

    number_components = input('How many components do you have?: ');
    component_matrix = zeros(number_components, 4);
    for i = 1:number_components
        comp_type = input('Select Reactant (0) or Product (1): ');
        stoich_coeff = input('Type abs value for stoichiometric coefficient: ');
        component_matrix(i, 2) = input('How many moles did you start with for this component?: ');
        component_matrix(i, 3) = input('What is the gamma value for this component'); %gamma = 1 for ideal soln
        component_matrix(i, 4) = input('What is the sat volume?: ');
        if(comp_type == 0)
            component_matrix(i, 1) = -stoich_coeff;
        else
            component_matrix(i, 1) = stoich_coeff;
        end
    end
    epsilon = fzero(@(x) determine_epsilon_liquid_ideal(K, v, total_moles, component_matrix, x, number_components), 0.5);
    disp('Epsilon value for the rxn: ' + string(epsilon));
elseif(user_choice == 7)
    x1 = input('What is the value of x1 = phi_1: ');
    Tm = input('What is the melting point temperature (K): ');
    delH = input('What is the heat of fusion value?: ');
    R = input('What is the rate constant that you will be using: ');
    
    T_root = fzero(@(T) determineT_SLE(x1, Tm, delH, R, T), 1000);
    disp('Temperature (K): ' + string(T_root));
elseif(user_choice == 8)
    Tm1 = input('What is the melting point temperature of comp 1 (K): ');
    Tm2 = input('What is the melting point temperature of comp 2 (K): ');
    delH1 = input('What is the heat of fusion value (comp 1)?: ');
    delH2 = input('What is the heat of fusion value (comp 2)?: ');
    R = input('What is the rate constant that you will be using: ');

    k = linspace(0,1,1000);
    T_guess = (Tm1 + Tm2) / 2;
    guess = [0.5 ; 1000];
    for i = 1:length(k)
        root = fsolve(@(w) SLE_eutectic_calculations(w,k(i),Tm1, Tm2, delH1, delH2, R), guess, options);
        guess = root;
    end
    disp('x1 value: ' + string(root(1)));
    disp('Temperature: ' + string(root(2)));
end



function [DeltGRT, K] = determine_K_method_1(delt_A, delt_B, delt_C, delt_D,delt_Gibbs, delt_Enthalpy, T, T_ref, R)
    IDCPS_1 = delt_A*log(T/T_ref);
    IDCPS_2 = ( delt_B + (delt_C + ( delt_D / ( (T^2)*(T_ref^2) ) ))*((T+T_ref)/2))*(T-T_ref); 
    IDCPS = IDCPS_1 + IDCPS_2;
    IDCPH = (delt_A*(T-T_ref)+(1/2)*delt_B*((T^2)-(T_ref^2))+(1/3)*delt_C*((T^3)-(T_ref^3))-delt_D*((T^(-1))-(T_ref^(-1))));
    
    comp1_ans = (delt_Gibbs - delt_Enthalpy) / (R*T_ref);
    comp2_ans = (delt_Enthalpy) / (R*T);
    DeltGRT = comp1_ans + comp2_ans + (1/T)*IDCPH - IDCPS;
    K = exp(-DeltGRT);
end

function K = three_K_method(delt_A, delt_B, delt_C, delt_D,delt_Gibbs, delt_Enthalpy, T, T_ref, R)
    IDCPS_1 = delt_A*log(T/T_ref);
    IDCPS_2 = ( delt_B + (delt_C + ( delt_D / ( (T^2)*(T_ref^2) ) ))*((T+T_ref)/2))*(T-T_ref); 
    
    IDCPS = IDCPS_1 + IDCPS_2;
    IDCPH = (delt_A*(T-T_ref)+(1/2)*delt_B*((T^2)-(T_ref^2))+(1/3)*delt_C*((T^3)-(T_ref^3))-delt_D*((T^(-1))-(T_ref^(-1))));

    K0 = exp((-delt_Gibbs) / (R*T_ref));
    K1 = exp(((delt_Enthalpy) / (R*T_ref)) * (1 - (T_ref/T)));
    K2 = exp(((-1) / T)*IDCPH + IDCPS);

    K = K0*K1*K2;
end

function f = determine_Temp(delt_A, delt_B, delt_C, delt_D,delt_Gibbs, delt_Enthalpy, T, T_ref, R, K, k)
    IDCPS_1 = delt_A*log(T/T_ref);
    IDCPS_2 = ( delt_B + (delt_C + ( delt_D / ( (T^2)*(T_ref^2) ) ))*((T+T_ref)/2))*(T-T_ref); 
    IDCPS = IDCPS_1 + IDCPS_2;
    IDCPH = (delt_A*(T-T_ref)+(1/2)*delt_B*((T^2)-(T_ref^2))+(1/3)*delt_C*((T^3)-(T_ref^3))-delt_D*((T^(-1))-(T_ref^(-1))));
    
    comp1_ans = (delt_Gibbs - delt_Enthalpy) / (R*T_ref);
    comp2_ans = (delt_Enthalpy) / (R*T);
    DeltGRT = comp1_ans + comp2_ans + (1/T)*IDCPH - IDCPS;

    f = k*(-DeltGRT - log(K));
end

function f = determine_epsilon_gas(K, P, P_ref, v, total_moles, component_matrix, x, N)
    yi = 1;
    for i = 1:N
        vi = component_matrix(i, 1);
        phi = component_matrix(i,3);
        yi = yi * ((((component_matrix(i, 2) + vi*x) / (total_moles + v*x)) * phi)^(vi));
    end
    right_side = ((P / P_ref)^(-v)) * K;
    f = yi - right_side;
end

function f = determine_epsilon_liquid(K, P, P_ref,T, v, total_moles, component_matrix, x, N)
    R = 8.314;
    xi = 1;
    sum_volume = 0;
    for i = 1:N
        vi = component_matrix(i, 1);
        gamma = component_matrix(i,3);
        xi = xi * ((((component_matrix(i, 2) + vi*x) / (total_moles + v*x)) * gamma)^(vi));

        sum_volume = sum_volume + vi*component_matrix(i,4);
    end
    right_side = exp(((P - P_ref) / (R*T)) * sum_volume) * K;
    f = xi - right_side;
end

function f = determine_epsilon_liquid_ideal(K, v, total_moles, component_matrix, x, N)
    R = 8.314;
    xi = 1;
    for i = 1:N
        vi = component_matrix(i, 1);
        gamma = component_matrix(i,3);
        xi = xi * ((((component_matrix(i, 2) + vi*x) / (total_moles + v*x)) * gamma)^(vi));
    end
    f = xi - K;
end

function f = determineT_SLE(x1, Tm, delH, R, T)
    f = exp(((delH) / (R*Tm)) * ((T - Tm) / (T))) - x1;
end

function f = SLE_eutectic_calculations(w,k,Tm1, Tm2, delH1, delH2, R)
    x1 = w(1);
    x2 = 1 - x1;
    T = w(2);

    f = w*0;

    f(1) = k*(exp(((delH1) / (R*Tm1)) * ((T - Tm1) / (T))) - x1);
    f(2) = k*(exp(((delH2) / (R*Tm2)) * ((T - Tm2) / (T))) - x2);

end
























