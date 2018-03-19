% Kinetic submodel tray bioreactor: Degradation of a mixture of dyes
% Sara Jimenez, Maura Guérios and David Mitchell
% Research groups SIRYTCOR and LTEF

close all; clear all; clc;

%% Enzyme production 

% dLdt = alpha.*r_CO2(t)

% Experimental data

    % Laccase activity from the day 3 to 13 in U/kg_ds (0.5 g)

    L = 10^3*[9.3646; 8.3516; 7.1141; 10.1183; 13.9929; 15.3791; ...
         20.3626; 15.5069; 29.5100; 29.4214; 32.0675];

    % Time in hours

    t = 24*[3; 4; 5; 6; 7; 8; 9; 10; 11; 12; 13];

    % Prelocation

    dL = zeros(length(L)-1,1);

    r_CO2 = zeros(length(t)-1,1);

    tl = zeros(length(t)-1,1);
    
    for i=1:10

        dL(i) = (L(i+1)-L(i))/(t(i+1)-t(i));

        tl(i) = (t(i+1)+t(i))/2;

        r_CO2(i) = -0.0001.*(tl(i))+0.0549;

    end
    
% Fitting to a one degree polynomial

    lineL = fitlm(r_CO2, dL,'Intercept',false);
    alpha=lineL.Coefficients.Estimate;
    alpha=1615;

figure(1)
plot(lineL), title('Lacasse production'),xlabel('(molCO_2)/(kgds h)'),ylabel('[(U/kgds)/h]'); hold on;



%% Solution of the ODE for the enzyme production

ti=[0 384];

%Li=10^3*9.3646;

Ki=[0, 1416.18];

[t,K]=ode45(@Enzyme,ti,Ki, [], alpha);

L=K(:,1);
D=K(:,2);

t=t/24;

% Laccase activity from the day 2 to 13 in U/kg_ds (0.5 g)
texp = [2; 4; 6; 8; 14; 16];
Lexp = 10^3*[0.2559 3.6438 8.0834 9.3111 26.0593 29.4135];

% Dye concentration
tdye = xlsread('DyeDegradation.xlsx','Raw data','A3:A19'); 
TA = xlsread('DyeDegradation.xlsx','Raw data','N3:N19'); 
AR = xlsread('DyeDegradation.xlsx','Raw data','O3:O19');
Mix = xlsread('DyeDegradation.xlsx','Raw data','P3:P19');

figure(2)
plot(t, L, 'k', texp, Lexp, 'bx'), title('Lacasse production'), xlabel('t [days]'), ylabel('[(U/kgds)]');
figure(3)
plot(t, D, 'k', tdye, Mix, 'k*', tdye, TA, 'kx', tdye, AR, 'k+'), title('Dye degradation'), xlabel('t [days]'), ylabel('Dye concentration [(mg_d_y_e/kg_d_s)]');



%% 

function dK = Enzyme(t,K, alpha)

L=K(1);
D=K(2);

dL=alpha.*(-0.0001.*t+0.0549);

dD=-0.0017.*L.*D./(2499+D);

dK=[dL;dD];

end


