% Model of a tray bioreactor operating under SSF conditions
% Sara Jimenez, Maura Guerios and David Mitchell
% Research groups SIRYTCOR and LTEF

% Main program 

close all; clear all; clc;

tic

%% Version 2.0: The second version will only consider the mass balances in the bed, treated as a isothermic system, the production of enzyme and the dye degradation

% Abbreviations 
    % a:    air
    % s:    solids
    % b:    bed
    % g:    gas
    % ds:   dry solids
    % da:   dry air
    
% Geometry: rectangular

% Parameters given by the user

    Input=xlsread('Input.xlsx','Sheet1','C2:C14');

    % Pressure $Pa$
    par.P=Input(1);
    % Temperature $ºC$
    par.T=Input(2);
    % Weight of dry solids to ferment  $kg$
    par.Wds=Input(3);
    % Width of the bioreactor $m$
    par.D=Input(4);
    % Length of the bioreactor $m$
    par.L=Input(5);
    % Initial moisture content in wet basis $\frac{kg_{H_2O}}{kg_{ws}}$
    par.phii_wb=Input(6);
    % Porosity $\frac{m^3_a}{m^3_b} = \frac{m^2_a}{m^2_b}$
    par.epsilon=Input(7);
    % Density of dry solids $\frac{kg_{ds}}{m^3_s}$
    par.rho_s=Input(8); 
    % Dry-weight loss factor
    par.f_dwl=Input(9);
    % Airflow $\frac{L}{s}$
    par.G=Input(10);
    % Relative humidity of the headspace $%$
    par.RH=Input(11); 
    % Oxygen concentration in the headspace $%v/v$
    par.O_air=Input(12); 
    % Carbon dioxide concentration in the headspace $%v/v$
    par.C_air=Input(13); 

% Parameters associated with the geometry (Rectangular)
    
    % Volume of the bed $m^3$
    par.Vb=par.Wds*(1/par.rho_s)*(1/(1-par.epsilon));  
    % Height of the bed $m$
    par.Height=par.Vb/(par.D*par.L);
    % Cross-sectional area of the flow channel $m^2$
    Acs=(pi*par.R^2)-(par.R^2/2*(theta-sin(theta)));
    % Perimeter of the flow channel $m$
    Pfc=S_h+c;
    
    % Number of layers
    par.nl=input('In how many layers do you want to divide the bed?');
    
    % Height of each layer in the bed $m$
    par.dz=par.Height/par.nl;
    
    % Pre-location 
    h=zeros(par.nl,1); auxr=zeros(par.nl,1); par.At=zeros(par.nl,1); par.Vt=zeros(par.nl,1);
    
    h(1)=par.Height;
    
    for i=1:par.nl
        
        % Area of transference of each layer of the bed $m^2$ 
        auxr(i)=sqrt(2*par.R*h(i)-h(i)^2);
        par.At(i)=2*auxr(i)*par.L;
        % Volumen of each layer of the bed $m^3$ 
        par.Vt(i)=par.At(i)*par.dz;
        h(i+1)=h(i)-par.dz;
      
    end
    
% Saturation vapor pressure of water (Antoine equation) $Pa$
    par.Psat=Antoine(par.T);   
    
% Parameters associated with the air 
    
    % Viscosity of air $\frac{Pa}{s}$
    par.miu_a=1.83e-5;
    % Density of dry air $\frac{kg_{da}}{m^3_a}$
    par.rho_a=1.14;
    % Water activity in the bed if it were in equilibrium with the headspace 
    par.aw_eq=par.RH/100; 
    % Humidity of the headspace $\frac{kg_{H_2O}}{kg_{da}}$ 
    par.H_h=definition(par.aw_eq,par);
    % Oxygen concentration in the headspace $\frac{molO_2}{m^3_{a}}$
    par.O_h=((par.O_air*par.rho_a)/(100*0.0288));
    % Carbon dioxide concentration in the headspace $\frac{molCO_2}{m^3_{a}}$
    par.C_h=((par.C_air*par.rho_a)/(100*0.0288));
    
% Diffusion coefficients

    % Oxygen in air $\frac{m^2_{a}}{h}$
    par.D_O=0.0634;
    % Carbon dioxide in air $\frac{m^2_{a}}{h}$
    par.D_C=0.05;
    % Water vapor in air $\frac{m^2_{a}}{h}$
    par.D_W=0.1015;
    
% Thickness of the stagnant layer

    delta=zeros(1,par.L*100); delta(1)=0;

        for j=1:par.L*100 % each centimeter

            % Reynolds number for flat plates
            Re=((j/100)*par.rho_a*(par.G/(1000*Acs)))/par.miu_a;

            if Re <= 500000
                delta(j+1)=5.0*(j/100)*Re^(-0.5); % Laminar flow
            else
                delta(j+1)=0.37*(j/100)*Re^(-0.2); % Turbulent flow
            end

        end

    j=(0:0.01:par.L);

%     figure(1)
%     plot(j*100,delta*100,'k'),title('Thickness of the stagnant layer'),xlabel('x [cm]'),ylabel('[cm]');

    delta=mean(delta);
     
% Heat capacities 

    % Dry air $\frac{J}{kg_{da}ºC}$
    Cp_da=(28.94+(0.4147e-2)*par.T+(0.3191e-5)*(par.T^2)-(1.965e-9)*(par.T^3))*(1/0.0288);
    % Nitrogen $\frac{J}{kg_{N2}ºC}$
    Cp_N=(29+(0.2199e-2)*par.T+(0.5723e-5)*(par.T^2)-(2.871e-9)*(par.T^3))*(1000/28.02);
    % Oxygen $\frac{J}{kg_{O2}ºC}$
    Cp_O=(29.1+(1.158e-2)*par.T-(0.6076e-5)*(par.T^2)+(1.311e-9)*(par.T^3))*(1000/32);
    % Carbon dioxide $\frac{J}{kg_{CO2}ºC}$
    Cp_C=(36.11+(4.233e-2)*par.T-(2.887e-5)*(par.T^2)+(7.464e-9)*(par.T^3))*(1000/44.01);
    % Water vapor $\frac{J}{kg_{W}ºC}$
    Cp_wv=(33.46+(0.6880e-2)*par.T+(0.7604e-5)*(par.T^2)-(3.593e-9)*(par.T^3))*(1000/18.016);
    % Liquid water $\frac{J}{kg_{W}ºC}$
    Cp_lw=(18.2964+(47.212e-2)*(par.T+273.15)-(133.88e-5)*((par.T+273.15)^2)+(1314.2e-9)*((par.T+273.15)^3))*(1000/18.016);
    % Gas phase in the headspace $\frac{J}{kg_{da}ºC}$
    Cp_gh=Cp_da+par.H_h*Cp_wv;
    
% Overall mass transfer coefficients 

    % Oxygen due to diffusion $\frac{m_{a}}{h}$
    par.Ky_O=par.D_O/delta; 
    % Carbon dioxide due to diffusion $\frac{m_{a}}{h}$
    par.Ky_C=par.D_C/delta; 
    % Water due to diffusion $\frac{m_{a}}{h}$
    par.Ky_dW=par.D_W/delta;
    % Water due to evaporation $\frac{kg_{H_2O}}{h m^2_{b} driving force}$
        % Gas mass flux $\frac{kg}{m^2 s}$
        J=par.G*par.rho_a/(1000*Acs);
        % Characteristic dimension of the system
        Dc=4*Acs/Pfc;
        % Heat transfer coefficient $\frac{W}{m^2 ºC}$
        h=8.8*J^(0.8)/Dc^(0.2);
    par.Ky_eW=((h/Cp_gh)*(3600))*(0.62413/((par.P/par.Psat)-1));
     
% Initial conditions 
    
    trun=input('How many hours do you want to run the simulation?');

    % Pre-location 
        Mi=zeros(1,par.nl); auxt=zeros(1,trun); 
        auxrho_s=zeros(trun,1); 
        auxO=zeros(trun,par.nl); 
        auxC=zeros(trun,par.nl); 
        auxW=zeros(trun,par.nl); 
        aw=zeros(trun,par.nl); 
        H=zeros(trun,par.nl); 
        phi=zeros(trun,par.nl);
        auxLac=zeros(trun,1);
        auxDye=zeros(trun,1);
    
    for j=1:par.nl  
        
    % Initial moisture content in dry basis $\frac{kg_{H_2O}}{kg_{ds}}$
        phi(1,j)=(par.phii_wb/(100-par.phii_wb))*100;
        par.phi=phi(1,1);
    % Initial water activity of the solids (aws)
        aw(1,j)=isotherm(phi(1),par);   
    % Initial humidity in the interparticle space air $\frac{kg_{H_2O}}{kg_{da}}$
        H(1,j)=definition(aw(1),par);   
    % Initial total amount of water $\frac{kg_{H_2O}}{m^3_{b}}$
        auxW(1,j)=(1-par.epsilon)*par.rho_s*phi(1)+par.epsilon*par.rho_a*H(1);
    % Initial oxygen concentration $\frac{molO_2}{m^3_{a}}$ - how affect the humidity just see the maximum error?
        auxO(1,j)=par.O_h;
    % Initial carbon dioxide concentration $\frac{molCO_2}{m^3_{a}}$
        auxC(1,j)=par.C_h;        
    % Initial density of dry-solids $\frac{kg_{ds}}{m^3_{b}}$
        auxrho_s(1,j)=par.rho_s;
    % Initial Laccase $\frac{{U}{kg_{ds}}}$
        auxLac(1,j)=0;
    % Initial dye concentration $\frac{{mg}{kg_{ds}}}$
        auxDye(1,j)=1416.18;
        
    end

        
    for i=1:trun % Hours
        
        Mi(1:par.nl)=auxO(i,:); 
        Mi((par.nl+1):(2*par.nl))=auxC(i,:); 
        Mi((2*par.nl+1):(3*par.nl))=auxW(i,:); 
        Mi((3*par.nl+1):(4*par.nl))=auxrho_s(i,:);
        Mi((4*par.nl+1):(5*par.nl))=auxLac(i,:);
        Mi((5*par.nl+1):(6*par.nl))=auxDye(i,:);

        auxt(1)=0; 

        ti=[auxt(i):1/60:i]; % the ode is gonna run for 1 hour with steps of 1 min

        [t,M]=ode15s(@TrayFunction,ti,Mi,[],par);
        
        auxt(i+1)=auxt(i)+1;

        auxO(i+1,:)=M(end,(1:par.nl)); auxC(i+1,:)=M(end,((par.nl+1):(2*par.nl))); 
        auxW(i+1,:)=M(end,((2*par.nl+1):(3*par.nl))); auxrho_s(i+1,:)=M(end,((3*par.nl+1):(4*par.nl)));
        auxLac(i+1,:)=M(end,((4*par.nl+1):(5*par.nl))); auxDye(i+1,:)=M(end,((5*par.nl+1):(6*par.nl)));
        
        O((i*60+1-60):(i*60+1),:)=M(:,1:par.nl); 
        C((i*60+1-60):(i*60+1),:)=M(:,(par.nl+1):(2*par.nl)); 
        W((i*60+1-60):(i*60+1),:)=M(:,(2*par.nl+1):(3*par.nl)); 
        rho_s((i*60+1-60):(i*60+1),:)=M(:,(3*par.nl+1):(4*par.nl));
        Lac((i*60+1-60):(i*60+1),:)=M(:,(4*par.nl+1):(5*par.nl));
        Dye((i*60+1-60):(i*60+1),:)=M(:,(5*par.nl+1):(6*par.nl));
        
    for k=1:par.nl
        
        eqn2=@(phi_v)(1-par.epsilon)*rho_s(end,k).*phi_v+par.epsilon*par.rho_a*((par.Psat/(0.0288*par.P))*(1-exp(-0.0004373*(par.T+57.926).*phi_v^(1.2626))))/(55.5556-((55.5556*par.Psat)/par.P)*(1-exp(-0.0004373*(par.T+57.926).*phi_v^(1.2626))))-W(end,k);

        phi(i+1,k)=fsolve(eqn2, par.phi); phi(i+1,k)=double(phi(i+1,k));

        aw(i+1,k)=isotherm(phi(i+1,k),par);

        H(i+1,k)=definition(aw(i+1,k),par);
      
    end
    
    end
      
    t=(0:1:i); t1=(0:1:i*60);
    t=t/(24); t1=t1/(24*60);
    
    for i=1:par.nl
               
%         figure(2)
%         subplot(2,2,1),plot(t1,rho_s(:,par.nl),'k'),title('Density of dry solids'),xlabel('t [days]'),ylabel('[kg_d_s/m^3_b]');
%         subplot(2,2,2),plot(t1,O(:,i),'k'),title('Oxygen'),xlabel('t [days]'),ylabel('[molO_2/m^3_a]'); hold on;
%         subplot(2,2,3),plot(t1,C(:,i),'k'),title('Carbon dioxide'),xlabel('t [days]'),ylabel('[molCO_2/m^3_a]'); hold on;
%         subplot(2,2,4),plot(t1,W(:,i),'k'),title('Water'),xlabel('t [days]'),ylabel('[kg_H_2_O/m^3_b]'); hold on;
%             
%         figure(3)
%         subplot(1,3,1),plot(t,phi(:,i),'k'),title('Moisture content'),xlabel('t [days]'),ylabel('[kg_H_2_O/m^3_s]'); hold on;
%         subplot(1,3,2),plot(t,aw(:,i),'k'),title('Water activity'),xlabel('t [days]'); hold on;
%         subplot(1,3,3),plot(t,H(:,i),'k'),title('Humidity'),xlabel('t [days]'),ylabel('[kg_H_2_O/m^3_a]'); hold on;

        figure(2)
        plot(t1,rho_s(:,par.nl),'k'),title('Density of dry solids'),xlabel('t [days]'),ylabel('[kg_d_s/m^3_b]');
        figure(3)
        plot(t1,Lac(:,par.nl),'k'),title('Laccase production'),xlabel('t [days]'),ylabel('[U/kg_d_s]');
        figure(4)
        plot(t1,Dye(:,par.nl),'k'),title('Dye degradation'),xlabel('t [days]'),ylabel('[mg_D_y_e/kg_d_s]');
        figure(5)
        plot(t1,O(:,i),'k'),title('Oxygen'),xlabel('t [days]'),ylabel('[molO_2/m^3_a]'); hold on;
        figure(6)
        plot(t1,C(:,i),'k'),title('Carbon dioxide'),xlabel('t [days]'),ylabel('[molCO_2/m^3_a]'); hold on;
        figure(7)
        plot(t1,W(:,i),'k'),title('Water'),xlabel('t [days]'),ylabel('[kg_H_2_O/m^3_b]'); hold on;  
        figure(8)
        plot(t,phi(:,i),'k'),title('Moisture content'),xlabel('t [days]'),ylabel('[kg_H_2_O/m^3_s]'); hold on;
        figure(9)
        plot(t,aw(:,i),'k'),title('Water activity'),xlabel('t [days]'); hold on;
        figure(10)
        plot(t,H(:,i),'k'),title('Humidity'),xlabel('t [days]'),ylabel('[kg_H_2_O/m^3_a]'); hold on;


    end
    
    toc
%% Auxiliary functions
    
    function aw = isotherm(phi,par)
    
    % Modified Henderson equilibrium moisture equation for corn cob - White et al. 1985
    
    aw=(1-exp(-0.0004373*(par.T+57.926)*phi^(1.2626)));
    
    end 
    
    function Psat = Antoine(T)
    
    % Saturation vapor pressure of pure water (Antoine equation) $Pa$
    
    Psat=133.322*exp(18.3036-(3816.44/(T+227.02)));
    
    end
    
    function H = definition(aw,par)
    
    % From the assumption of pseudohomogeneous phase aws = awg and awg = Pw/Psat. 
    % Pw can be expressed as (nw/nt)*P, in such expression the only unknown variable 
    % is the humidity content in the interparticle space air (H)
    
    % Rearranging from the definition of awg - Initial humidity in the interparticle space air $\frac{kg_{H_2O}}{kg_{da}}$
    
    H=((1/0.0288)*(par.Psat/par.P)*aw)/(55.5556-55.5556*aw*(par.Psat/par.P));
    
    end
    
%% Mass balances and kinetic sub-model as a function of time

    function dM = TrayFunction(t, M, par)
        
    % Components
    
        O=M(1:par.nl);                      % Oxygen concentration $\frac{molO_2}{m^3_{a}}$
        C=M((par.nl+1):(2*par.nl));         % Carbon dioxide concentration $\frac{molCO_2}{m^3_{a}}$
        W=M((2*par.nl+1):(3*par.nl));       % Total amount of water in a control volume $\frac{kg_{H_2O}}{m^3_{b}}$
        rho_s=M((3*par.nl+1):(4*par.nl));   % Density of dry solids $\frac{{kg_{ds}}{m^3_{s}}}$
        Lac=M((4*par.nl+1):(5*par.nl));     % Laccase $\frac{U/kg_{ds}}$
        Dye=M((5*par.nl+1):(6*par.nl));     % Dye $\frac{mg/kg_{ds}}$
        
    % Pre-location 
        
        dO=zeros(par.nl,1); 
        dC=zeros(par.nl,1); 
        dW=zeros(par.nl,1); 
        drho_s=zeros(par.nl,1);
        dLac=zeros(par.nl,1);
        dDye=zeros(par.nl,1);
        phi_v=zeros(par.nl,1); 
        aw_v=zeros(par.nl,1); 
        H_v=zeros(par.nl,1); 
        WV=zeros(par.nl,1);
        
    % Respirometric

        if t<12.0000
            
            % Oxygen uptake rate $\frac{molO_2}{kg_{solids}h}$
            r_O=0;
            % Carbon dioxide production rate $\frac{molCO_2}{kg_{solids}h}$
            r_C=0; 

        elseif (12.0000<=t)&&(t<=72.0000) % the ode runs each hour

            % Oxygen uptake rate $\frac{molO_2}{kg_{solids}h}$
            r_O=(0.0534/60)*t-(0.6408/60);
            % Carbon dioxide production rate $\frac{molCO_2}{kg_{solids}h}$
            r_C=(0.0477/60)*t-(0.5724/60);

        else

            % Oxygen uptake rate $\frac{molO_2}{kg_{solids}h}$
            r_O=-0.0002*t+0.0678;
            % Carbon dioxide production rate $\frac{molCO_2}{kg_{solids}h}$
            r_C=-0.0001*t+0.0549;  

        end
    

    % Auxiliary equations - parameters
    
    for i=1:par.nl
        
        eqn2=@(phi_v)(1-par.epsilon)*rho_s(i).*phi_v+par.epsilon*par.rho_a*((par.Psat/(0.0288*par.P))*(1-exp(-0.0004373*(par.T+57.926).*phi_v^(1.2626))))/(55.5556-((55.5556*par.Psat)/par.P)*(1-exp(-0.0004373*(par.T+57.926).*phi_v^(1.2626))))-W(i);

        % Moisture content $\frac{kg_{H2O}}{kg_{ds}}$
        phi_v(i)=fsolve(eqn2, par.phi); phi_v(i)=double(phi_v(i));

        % Water activity
        aw_v(i)=isotherm(phi_v(i),par);

        % Humidity $\frac{kg_{H2O}}{kg_{da}}$
        H_v(i)=definition(aw_v(i),par);
            
        % Water vapour $\frac{kg_{H2O}}{kg_b}$
        WV(i)=par.epsilon*par.rho_a*H_v(i);
        
        % Density of solids 
        drho_s(i)=(-par.f_dwl.*r_C.*par.Wds)./(par.Vb*(1-par.epsilon));
        
        % Laccase production 
        dLac(i)=(2277.*r_C.*par.Wds)./(rho_s(i).*par.Vb*(1-par.epsilon));
        
        % Dye degradation 
        dDye(i)=(-0.001.*Lac(i).*Dye(i)./(2000+Dye(i))).*par.Wds./(rho_s(i).*par.Vb*(1-par.epsilon));
        
    end

    % Top layer balances
    
        dO(1)=(1/(par.epsilon.*par.Vt(1)))*(par.Ky_O*par.epsilon.*par.At(1).*(par.O_h-O(1))-(1-par.epsilon)*par.rho_s.*par.Vt(1)*r_O-par.epsilon.*par.At(2)*par.D_O.*(O(1)-O(2))/par.dz); 
        dC(1)=(1/(par.epsilon.*par.Vt(1)))*(par.epsilon.*par.At(2)*par.D_C*(C(2)-C(1))/par.dz+(1-par.epsilon)*par.rho_s.*par.Vt(1)*r_C-par.Ky_C*par.epsilon.*par.At(1)*(C(1)-par.C_h));
        dW(1)=(1./par.Vt(1))*(par.epsilon.*par.At(2)*par.D_W*(WV(2)-WV(1))*(1/par.epsilon)/par.dz-par.Ky_eW.*par.At(1)*(aw_v(1)-par.aw_eq)-par.Ky_dW*par.epsilon*par.rho_a.*par.At(1)*(H_v(1)-par.H_h));

        
    % Middle layers balances
    
    for i=2:(par.nl-1)        
 
        dO(i)=(1/(par.epsilon.*par.Vt(i)))*(par.epsilon.*par.At(i)*par.D_O.*(O(i-1)-O(i))/par.dz-(1-par.epsilon)*par.rho_s.*par.Vt(i)*r_O-par.epsilon.*par.At(i+1)*par.D_O.*(O(i)-O(i+1))/par.dz);
        dC(i)=(1/(par.epsilon.*par.Vt(i)))*(par.epsilon.*par.At(i+1)*par.D_C*(C(i+1)-C(i))/par.dz+(1-par.epsilon)*par.rho_s.*par.Vt(i)*r_C-par.epsilon.*par.At(i)*par.D_C*(C(i)-C(i-1))/par.dz);
        dW(i)=(1./par.Vt(i))*(par.epsilon.*par.At(i+1)*par.D_W*(WV(i+1)-WV(i))*(1/par.epsilon)/par.dz-par.epsilon.*par.At(i)*par.D_W*(WV(i)-WV(i-1))*(1/par.epsilon)/par.dz);
        
    end
        
    % Bottom layer balances

        dO(par.nl)=(1/(par.epsilon.*par.Vt(par.nl)))*(par.epsilon.*par.At(par.nl)*par.D_O*(O(par.nl-1)-O(par.nl))/par.dz-(1-par.epsilon)*par.rho_s.*par.Vt(par.nl)*r_O);
        dC(par.nl)=(1/(par.epsilon.*par.Vt(par.nl)))*((1-par.epsilon)*par.rho_s.*par.Vt(par.nl)*r_C-par.epsilon.*par.At(par.nl)*par.D_C*(C(par.nl)-C(par.nl-1))/par.dz);              
        dW(par.nl)=(1./par.Vt(par.nl))*(-par.epsilon.*par.At(par.nl)*par.D_W*(WV(par.nl)-WV(par.nl-1))*(1/par.epsilon)/par.dz);
        
      
        dM=[dO; dC; dW; drho_s; dLac; dDye];

    end 