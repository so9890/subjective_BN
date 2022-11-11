%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%COMET Model Supplement: Abatement Cost Function Coefficient Estimation %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Author: Lint Barrage, Brown University (lint.barrage@gmail.com)
%%Supplement to: "Optimal Dynamic Carbon Taxes in a Climate-Economy Model
%%with Distortionary Fiscal Policy" in The Review of Economic Studies (2019)

%SUMMARY:
%This m-file produces the COMET model's abatement cost function coefficients
%Inputs: 2010 DICE Model (Nordhaus, 2010) vectors of forecast business-as-usual
%(i)   emissions intensity "DICE_2010_Sigma"
%(ii)  global GDP "DICE_2010_YBAU"
%(iii) industrial carbon emissions "DICE_Emiss_BAU"
%M-file inputs: objfun_fitlogit_Abatement_Cost_Curve.m
%Output: COMET_Abatement_Coefficients.xlsx - vector of abatement cost function coefficients

clear
cd 'C:\Users\lintb\Desktop\COMET\Supplementary\COMET Calibration\COMET Abatement Costs'

%%% Set Parameters %%%
%%%%%%%%%%%%%%%%%%%%%%
T = 24;                 %Direct optimization period time horizon, 1 t = 10 years
periods = 10;           %Simulation periods after T
y = (1:1:T);            %Time in calendar years
for i = 2:1:T+periods;
    y(i) = 2015+(i*10);
end
y(1) = 2015;
%Discount rate:
r = 0.05;

%%% Load DICE Inputs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
Pc0 = 1.26;                             %DICE base year price of backstop technology, thousands of dollars per mtC
Pc = zeros((periods+T),1);
YearComp = 2250;                        %DICE year when the backstop technology becomes cost competitive (no additional cost for clean energy)
CostlyYears = YearComp - 2010;
CostlyPeriods = CostlyYears/10;
vau1 = 2;                               
vau2 = 0.05;
for i = 1:1:(CostlyPeriods);
    Pc(i) = Pc0*(vau1-1+exp((-vau2)*(i)))*(1/vau1);
end
phi2 = 2.8;  %DICE Abatement cost function exponent
  load('DICE_2010_Sigma','DICE_2010_Sigma')
  Sigma = DICE_2010_Sigma(1:T+1);
  load('DICE_2010_YBAU','DICE_2010_YBAU')
  DICE_2010_YBAU = DICE_2010_YBAU(1:T+1);
  load('DICE_Emiss_BAU','DICE_Emiss_BAU')
  %Compute DICE abatement costs function coefficients:
   phi1_t = zeros(T,1);
   phi1t_hat = zeros(T,1);
  for i = 1:1:T;
      phi1_t(i) = Pc(i)*(Sigma(i+1)/phi2);  %DICE abatement cost coefficient (for fraction abated)
      phi1t_hat(i) = phi1_t(i)*(DICE_2010_YBAU(i+1)*10*1000)*(DICE_Emiss_BAU(i)*10)^(-phi2);    %Implied abatement cost coefficient per TON of clean energy
  end

%%% Vector of Clean Energy / Abatement LEVELS %%%
EVec = (0:10:250);

%%% Compute DICE-Implied Total Costs %%%
TC_DICE = zeros(T,(length(EVec)));
for i = 1:1:T;
    for j = 1:1:length(EVec);
        if EVec(j)<=(DICE_Emiss_BAU(i)*10)
             TC_DICE(i,j) = phi1t_hat(i)*EVec(j)^phi2;             
        else
            TC_DICE(i,j) = (phi1t_hat(i)*(DICE_Emiss_BAU(i)*10)^phi2)+Pc(i)*1000*(EVec(j)-(DICE_Emiss_BAU(i)*10));
        end
    end
end

%%% Choose Parameters to Minimize Sum of Squared Errors %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Initial guess:
x0 = zeros(7,1);
x0(1) = 1;  
x0(2) = 1;    
x0(3) = 50;
x0(4) = 20;
x0(5) = 1;
x0(6) = 1;
x0(7) = 5;
x0(8) = -1;

%Test:
f = objfun_fitlogit_Abatement_Cost_Curve(x0,TC_DICE,EVec,T,Pc,r)

%Bounds:
lb = (-Inf)*ones(8,1);
ub = (Inf)*ones(8,1);

%Minimize:
options = optimset('Display','iter','MaxIter',10000,'MaxFunEvals',320000,'TolFun',1e-9);
[x,fval] = fminunc(@(x)objfun_fitlogit_Abatement_Cost_Curve(x,TC_DICE,EVec,T,Pc,r),x0,options)

%Display:
gamma = x(1)
b3 = x(2)
a1 = x(3)
a2 = x(4)
a3 = x(5)
a4 = x(6)
b1 = x(7)
b2 = x(8)

%Save:
M = {'COMET Parameter', 'Value', 'Alt. Notation';'gamma', gamma, '"alpha-bar"'; 'a1', a1, '"alpha1"'; 'a4', a4, '"alpha2"'; 'a2', a2, '"b01"'; 'a3', a3, '"b02"'; 'b1', b1, '"b11"'; 'b2', b2, '"b22"'; 'b3', b3, '"b2"'}
filename = 'COMET_Abatement_Coefficients.xlsx'
xlswrite('COMET_Abatement_Coefficients.xlsx',M)


