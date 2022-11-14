clc
clear
close all
options = optimset('Display','off');
R = 8.314;
syms f(x) ;
syms g(x) ;
T_m_Bi = 544;
T_m_Cd = 594;
delta_H0_Bi = 10900;
delta_H0_Cd = 6400;
delta_S0_Bi = 20 ;
delta_S0_Cd = 10.77;
f(x) =   1.2 - (16.45e-3 * x  ) + (21.1e5 .* x^-2 );
g(x) =  7.5 - 12.3e-3 .* x;
G1 = g(x) / x;
F1(x) = f(x) / x;
F = int(f,T_m_Bi,x,'IgnoreAnalyticConstraints',true);
G = int(g,T_m_Cd,x,'IgnoreAnalyticConstraints',true);
F1_prime = int(F1,x,T_m_Bi,x,'IgnoreAnalyticConstraints',true);
G1_prime = int(G1,T_m_Cd,x,'IgnoreAnalyticConstraints',true);
delta_G_Bi(x) = delta_H0_Bi + F - x * (delta_S0_Bi + F1_prime);
delta_G_Cd(x) = delta_H0_Cd + G - x * (delta_S0_Cd + G1_prime);
X_Bi(x) = exp(-delta_G_Bi/(R * x));
X_Cd(x) = exp(-delta_G_Cd/(R * x));
X1 = matlabFunction(X_Bi);
X2 = matlabFunction(X_Cd);
Y = @(x) X1(x) + X2(x) - 1;
S = fsolve(Y,406,options);
T = S : 1 : 600;
P = 0 : 0.1 : 1;
x1 = [0 1];
y = [S S];
x = input('Enter the amount of X_Cd : ');
t = input('Enter the temperature : ');
if t > S
    if x < 1-X1(t)
        disp('LIQUID + SOLID Bi')
    elseif x > 1-X1(t) && x < X2(t)
        disp('LIQUID')
    else
        disp('LIQUID + SOLID Cd')
    end
else
    disp('SOLID Bi +SOLID Cd')
end
plot(1-X1(T),T,'g--')
hold on
plot(X2(T),T,'--')
line(x1,y);
axis([0,1,380,600]);
title('Cd-Bi phase diagram','color','red')
text(0.05,460,'liquid + solid Bi')
text(0.75,460,'liquid + solid Cd')
text(0.5,500,'liquid')
text(0.38,395,'solid Bi + solid Cd')
xlabel('X_C_d');
ylabel('T, K');
fprintf('The eutectic point is %7.2f \n',S)











