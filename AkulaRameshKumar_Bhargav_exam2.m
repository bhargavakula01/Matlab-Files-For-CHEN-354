user = input('Enter value: ');
options = optimoptions('fsolve','Display','off');

%WILL NEED TO CHANGE GAMMA EQUATIONS BASED ON THE MODEL BEING USED!

if(user == 1)
    disp('Running Bubble P calculations');
    %Bubble P
    %Given x1 and T and need to find y1 and P
    p1sat = input('Enter value for comp 1 sat pressure (kPa): ');
    p2sat = input('Enter value for comp 2 sat pressure (kPa): ');
    pavg = (p1sat + p2sat) / 2;
    x1 = input('Enter value for x1: ');
    gamma1 = input('Enter value for gamma1: ');
    gamma2 = input('Enter value for gamma2: ');
    k = linspace(0,1,1000);
    w_initial = [pavg ; 0.5];
    for i = 1:length(k)
        root = fsolve(@(w) BubbleP(w,k(i), p1sat, p2sat, x1, gamma1, gamma2), w_initial, options);
        w_initial = root;
    end
    disp('Bubble Pressure (kPa): ' + string(w_initial(1)));
    disp('y1: ' + string(w_initial(2)));
elseif(user == 2)
    disp('Running Dew P calculations');
    %Dew P
    %Given x1 and T and need to find y1 and P
    p1sat = input('Enter value for comp 1 sat pressure (kPa): ');
    p2sat = input('Enter value for comp 2 sat pressure (kPa): ');
    pavg = (p1sat + p2sat) / 2;
    y1 = input('Enter value for y1: ');
    k = linspace(0,1,1000);
    w_initial = [pavg ; 0.5];
    for i = 1:length(k)
        root = fsolve(@(w) DewP(w,k(i), p1sat, p2sat, y1), w_initial, options);
        w_initial = root;
    end
    disp('Dew Pressure (kPa): ' + string(w_initial(1)));
    disp('x1: ' + string(w_initial(2)));
elseif(user == 3)
    disp('Running Bubble T calculations');
    % Using T1sat and T2sat to determine initial T value
    T1sat = input('Enter value for comp 1 sat temperature (degrees C): ');
    T2sat = input('Enter value for comp 2 sat temperature (degrees C): ');
    T_initial = ((T1sat + T2sat) /2) + 273.15;
    P = input('Enter value for overall pressure (kPa): ');
    x1 = input('Enter value for x1: ');
    A1 = input('Enter value from Antoines for comp 1 for A: ');
    B1 = input('Enter value from Antoines for comp 1 for B: ');
    C1 = input('Enter value from Antoines for comp 1 for C: ');
    A2 = input('Enter value from Antoines for comp 2 for A: ');
    B2 = input('Enter value from Antoines for comp 2 for B: ');
    C2 = input('Enter value from Antoines for comp 2 for C: ');
    %solving for the bubble point temperature
    T_bubble = fzero(@(T) BubbleT(T,A1,B1, C1, A2, B2, C2, x1, P), T_initial);
    %displaying answer
    disp('Bubble Temperature: ' + string(T_bubble) + " K");
    p1sat = exp(A1 - (B1 / ((T_bubble - 273.15) + C1)));
    gamma1 = exp((2.771-0.00523*T_bubble)*(x2)^2);
    disp('y1: ' + string(BubbleT_calcy1(p1sat, gamma1, P, x1)));
elseif(user == 4)
    disp('Running Dew T calculations');
    % Using T1sat and T2sat to determine initial T value
    T1sat = input('Enter value for comp 1 sat temperature (degrees C): ');
    T2sat = input('Enter value for comp 2 sat temperature (degrees C): ');
    T_initial = ((T1sat + T2sat) /2) + 273.15;
    P = input('Enter value for overall pressure (kPa): ');
    y1 = input('Enter value for y1: ');
    A1 = input('Enter value from Antoines for comp 1 for A: ');
    B1 = input('Enter value from Antoines for comp 1 for B: ');
    C1 = input('Enter value from Antoines for comp 1 for C: ');
    A2 = input('Enter value from Antoines for comp 2 for A: ');
    B2 = input('Enter value from Antoines for comp 2 for B: ');
    C2 = input('Enter value from Antoines for comp 2 for C: ');
    
    w_initial = [0.5; T_initial]; %[x1, T]
    k = linspace(0,1,1000); 
    for i = 1:length(k)
        %Solving function below and updating guess
        DewT_root = fsolve(@(w) DewT(w,k(i),A1, B1, C1, A2, B2, C2, y1, P), w_initial, options);
        w_initial = DewT_root;
    end
    %displaying the values
    disp('Dew Point Temperature: ' + string(DewT_root(2)) + " K");
    disp('x1: ' + string(DewT_root(1)));
elseif(user == 5)
    %Azeotropic pressure and composition calculations
    p1sat = input('Enter value for comp 1 sat pressure (kPa): ');
    p2sat = input('Enter value for comp 2 sat pressure (kPa): ');
    T = input('Enter value for Temperature (Kelvin): ');
    pavg = (p1sat + p2sat) / 2;
    w_initial = [0.5; pavg];
    k = linspace(0,1,1000);
    for i = 1:length(k)
        root = fsolve(@(w) AzeotropeP(w,k(i),p1sat, p2sat, T), w_initial, options);
        w_initial = root;
    end
    disp('Azeotrope Pressure: ' + string(w_initial(2)) + " kPa");
    disp('x1: ' + string(w_initial(1)));
elseif(user == 6)
    %Azeotrpoic Temperature and composition calculations
    A1 = input('Enter value from Antoines for comp 1 for A: ');
    B1 = input('Enter value from Antoines for comp 1 for B: ');
    C1 = input('Enter value from Antoines for comp 1 for C: ');
    A2 = input('Enter value from Antoines for comp 2 for A: ');
    B2 = input('Enter value from Antoines for comp 2 for B: ');
    C2 = input('Enter value from Antoines for comp 2 for C: ');
    T1sat = input('Enter value for comp 1 sat temperature (degrees C): ');
    T2sat = input('Enter value for comp 2 sat temperature (degrees C): ');
    T_initial = ((T1sat + T2sat) /2) + 273.15;
    P = input('Enter value for overall pressure (kPa): ');
    w_initial = [0.5;T_initial];
    k = linspace(0,1,1000);
    for i = 1:length(k)
        root = fsolve(@(w) AzeotropeT(w,k(i),A1, B1, C1, A2, B2, C2), w_initial, options);
        w_initial = root;
    end
    disp('Azeotrope Temperature: ' + string(w_initial(2)) + " K");
    disp('x1: ' + string(w_initial(1)));
end



function f = BubbleP(w,k, p1sat, p2sat,x1, gamma1, gamma2)
    f = w*0;
    y1 = w(2);
    y2 = 1 - y1;
    P = w(1);
    x2 = 1 - x1;
    f(1) = k*(x1*gamma1*p1sat + x2*gamma2*p2sat - P);
    f(2) = k*(y1*P - x1*gamma1*p1sat);
end

function f = DewP(w,k, p1sat, p2sat,y1)
    A = 0.95;
    f = w*0;
    x1 = w(2);
    x2 = 1 - x1;
    y2 = 1 - y1;
    P = w(1);
    gamma1 = exp(A*(x2)^2);
    gamma2 = exp(A*(x1)^2);
    f(1) = k*(x1*gamma1*p1sat - y1*P);
    f(2) = k*(x2*gamma2*p2sat - y2*P);
end

function f = BubbleT(T,A1, B1, C1, A2, B2, C2, x1, P)
    %initializing constants
    A = (2.771-0.00523*T);
    x2 = 1 - x1;
    p1sat = exp(A1 - (B1 / ((T - 273.15) + C1)));
    p2sat = exp(A2 - (B2 / ((T - 273.15) + C2)));

    gamma1 = exp(A*(x2)^2);
    gamma2 = exp(A*(x1)^2);
    % using the modified raoult's law to determine bubble point temperature
    f = x1*gamma1*p1sat + x2*gamma2*p2sat - P;
end

function y1 = BubbleT_calcy1(p1sat, gamma1, P, x1)
    y1 = (x1*p1sat*gamma1) / P;
end

function f = DewT(w,k, A1, B1, C1, A2, B2, C2, y1, P)
    %intializing constants
    A = (2.771-0.00523*T);
    y2 = 1 - y1;
    %x1, x2, and T from vector "w"
    x1 = w(1);
    x2 = 1 - w(1);
    T = w(2);
    %determininng p1sat and p2sat
    p1sat = exp(A1 - (B1 / ((T - 273.15) + C1)));
    p2sat = exp(A2 - (B2 / ((T - 273.15) + C2)));

    gamma1 = exp(A*(x2)^2);
    gamma2 = exp(A*(x1)^2);
    %set of non linear functions that we are using to solve for x1 and T
    f = w * 0;
    f(1) = k*(x1*gamma1*p1sat) - y1*P;
    f(2) = k*(x2*gamma2*p2sat) - y2*P;

end

function f = AzeotropeP(w,k,p1sat, p2sat, T)
    x1 = w(1);
    x2 = 1 - x1;
    P = w(2);
    A = 0.95;
    gamma1 = exp(A*(x2)^2);
    gamma2 = exp(A*(x1)^2);
    
    f = w*0;
    f(1) = k*(gamma1*p1sat) - P;
    f(2) = k*(gamma2*p2sat) - P;
end

function f = AzeotropeT(w,k,A1, B1, C1, A2, B2, C2)
    x1 = w(1);
    x2 = 1 - x1;
    T = w(2);
    A = (2.771-0.00523*T);
    gamma1 = exp(A*(x2)^2);
    gamma2 = exp(A*(x1)^2);
    p1sat = exp(A1 - (B1 / ((T - 273.15) + C1)));
    p2sat = exp(A2 - (B2 / ((T - 273.15) + C2)));
    
    f = w*0;
    f(1) = k*(gamma1*p1sat) - P;
    f(2) = k*(gamma2*p2sat) - P;
end



