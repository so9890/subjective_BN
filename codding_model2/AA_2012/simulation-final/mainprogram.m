%******* This is the main program AcemogluaAghion etc 2012
%1. Fix the parameters in the parameter section
%2. Choose the options (DTC or not, delay and 2 instruments or only the carbon tax in the option section)
%3. Choose an initial guess for the allocation of researchers and the input tax.
% The program then computes the following vectors:
% - xx: the share of scientists in the clean sector
% - Acc and Add: the productivity of the clean and dirty sector respectively
% - Qq: the minimum clean research subsidy necessary to implement the optimal policy, when we assume that when multiple equilibria are possible for the allocation of scientists, the interior one is chosen.
% - tauu: the carbon tax
% - Cc: the consumption
% - Ss: environmental quality
% - Tt: the increase in temperature since preindustrial times
% - Ratio: the ration of clean to dirty inputs in production
% - Util: is the total utility at time 0
% - Utilinit: the utility with the initial guess is provided for the case of optimization using the carbon tax only as a check: with a carbon tax only the economy evolves very abruptly and the program sometimes fail to improve on the initial guess.
% Moreover the program plots the share of scientists, the path of Acc, Add,and Cc, the input tax and the increase in temperature.

%clear all; % uncomment to clear memory

%********************************************************************
%********************************************************************
%1. initial parameters
global dt rho sigma psi alpha gamma eta_d eta_c qsi epsilon delta numsim phi S_bar lambda t_disaster max_t
dt=5;  % number of years in a period
rho = 0.001*dt; % discount rate
epsilon =10;  % elasticity of substitution, must be greater than 1 for the program to work properly.
sigma = 2; % Coefficient of relative risk aversion 
alpha = 1/3;% 1/3; % share of machines in production
psi= alpha^2; % cost of machines
gamma = 1; % size of innovation
eta_d = 0.02*dt; % probability of success in the dirty sector
eta_c = 0.02*dt; % probability of success in the the clean sector
numsim=80; % number of periods
% initial values for the clean and dirty inputs, source table 11.1 in the Annual Energy Review 2008
Yc0=307.77; % production of non fossil fuel energy in the world primary supply of energy from 2002 to 2006 source in Quadrillion of Btu 
Yd0=1893.25; % production of fossil fuel energy in the world primary supply of energy from 2002 to 2006 source in Quadrillion of Btu
% initial emissions, source table 11.19 Annual Energy Review 2008 
emission0=17.48251782; % world emissions in ppm from 2002 to 2006 (converted from original data in GtCO2)
t_disaster=6;% disaster temperature
S_bar = 280*(2^(t_disaster/3)-1); % corresponding value for S_bar
S0 = S_bar - 99; % 99pm increase in CO2 concentration since preindustrial times
max_t=3; % temperature up to which the damage function is matched with Nordhaus's one
lambda = fmincon(@(l)damage_calibS(l),0.35,[],[],[],[],0.00001,0.999999,[],optimset('Tolfun',1e-11));
% the formula for Nordhaus's model is used in the function damage_calibS and can be found at http://nordhaus.econ.yale.edu/RICEmodels.htm
delta =.5*emission0/S0; % half of CO2 emissions are absorbed at current atmospheric levels


%********************************************************************
%********************************************************************
%2. options
mode = 1; % 1: for the model with directed technical change and both instruments, 
          % 0.5: for a model with DTC but only the carbon tax,
          % 0: for a model without dtc
% the program for the case with a carbon tax only is correct only if 
% epsilon>(2-alpha)/(1-alpha), moreover if this is the case,
% it assumes that when several allocation of scientists are possible
% equilibria, the equilibrium chosen is the interior one.
delay = 0; % delay before applying optimal policy (possible only for mode = 1)
display_iter=0;    %=1 if you want to see iterations, ~1 if otherwise
display_diag=1;    %=1 if you want to see diagnostics, ~1 if otherwise


%********************************************************************
%********************************************************************
%3. initial guesses
if mode == 1
    sce10rho01=0.4;
    taue10rho01=2;
    sc0=sce10rho01; % initial guess for the share of scientists in clean research
    t0=taue10rho01; % initial guess for the input tax
    x0=[sc0; t0]; % vector stacking the guess for the share of scientists in clean research and the guess for the input tax
    ub=[ones(numsim,1); 100000*ones(numsim,1)]; % upper bound for the optimization
    lb=[zeros(numsim,1); zeros(numsim,1)]; % lower boud for the optimization
else
    t0=taue10rho01taxonly; % initial guess for the input tax
    x0=t0;
    ub=100000*ones(numsim,1);
    lb=zeros(numsim,1);
end

%**** No intervention in the code is necessary beyond this point

%********************************************************************
%********************************************************************
% 4. Warnings
if epsilon <= 1
    display ('*********the program is not stable for that range******************')
end
if  mode ==.5 && espilon <= (2-alpha)/(1-alpha)
    display('*********epsilon is too low for the code to make sense ******************')
end
if  mode <1 && delay>0
    display('*********delay is a possible option only with DTC and two instruments******************')
end


%********************************************************************
%********************************************************************
% 5. preliminary computations
phi=(1-alpha)*(1-epsilon);
qsi = emission0/Yd0;
Ac0=(alpha/psi)^(-alpha/(1-alpha))*Yc0*(1+(Yc0/Yd0)^(1/epsilon-1))^(alpha/phi+1); % corresponding initial value for Ac0 (assuming that the optimal subsidy to the use of all machines is implemented)
Ad0 =(alpha/psi)^(-alpha/(1-alpha))*Yd0*(1+(Yd0/Yc0)^(1/epsilon-1))^(alpha/phi+1); % and for Ad0
numsimleft=numsim-delay;
StockAc=[];
StockAd=[];
StockC=[];
StockQ=[];
Stocktau=[];
Stocksc=[];
StockS=[];
utcomp=0;
% evolution of the economy during the delay period
if delay>0
    numsim = delay;
    comp=[zeros(delay,1);zeros(delay,1)];
    % comp stacks the allocation of scientists to the clean sector and the initial carbon tax
    RespC=mysimenvtaxnew2(comp, Ac0, Ad0, S0); % compute the relevant parameters for the economy when the optimal subsidies to all machines is used, and comp is a vector stacking the share of scientists allocated to the clean sector and the input tax
    StockAc=RespC.Ac;
    StockAd=RespC.Ad;
    StockC=RespC.C;
    StockQ=zeros(delay,1);
    Stocktau=zeros(delay,1);
    Stocksc=zeros(delay,1);
    StockS=RespC.S;
    utcomp=-sum(RespC.Util(1:end));
    S0 = min(max(0.00000000000000001,-qsi*(alpha/psi)^(alpha/(1-alpha))*(RespC.Ac(delay)^(((1-alpha)*(1-epsilon))+alpha)*RespC.Ad(delay))/((RespC.Ac(delay)^((1-alpha)*(1-epsilon))+RespC.Ad(delay)^((1-alpha)*(1-epsilon)))^(alpha/((1-alpha)*(1-epsilon)))*(RespC.Ac(delay)^((1-alpha)*(1-epsilon))+RespC.Ad(delay)^((1-alpha)*(1-epsilon)))) + (1+delta)*RespC.S(delay)),S_bar); %notice that S should not be very close to zero 
    Ac0=StockAc(delay);
    Ad0=StockAd(delay);
end
numsim=numsimleft;

%********************************************************************
%********************************************************************
%6. optimization
if display_iter==1
    if display_diag==1
    options = optimset('Display','iter', 'Diagnostics','on','FunValCheck','on','TolFun',1e-11,'TolX',1e-9);
    else
    options = optimset('Display','iter', 'Diagnostics','off', 'FunValCheck','on','TolFun',1e-11,'TolX',1e-9);
    end
else
    if display_diag==1
    options = optimset('Display','off', 'Diagnostics','on','FunValCheck','on','TolFun',1e-11,'TolX',1e-9);
    else 
    options = optimset('Display','off', 'Diagnostics','off', 'FunValCheck','on','TolFun',1e-11,'TolX',1e-9);
    end
end    
if mode == 1
    [x,fval,exitflag] = fmincon(@(x)mysimopttaxnew2(x, Ac0, Ad0, S0),x0,[],[],[],[],lb,ub,[],options); % optimization with DTC and 2 instruments
    % mysimopttaxnew2 computes the utility given x - a vector stacking the path for the share of scientists in clean research and the input tax - and the initial values of clean and dirty productivities and quality of the environment
else
    if mode == 0
        [x,fval,exitflag] = fmincon(@(x)mysimopttaxnew2noDTC(x, Ac0, Ad0, S0),x0,[],[],[],[],lb,ub,[],options); % optimization without DTC
        % similarly mysimopttaxnew2noDTC computes the utility given a vector of input tax and assuming that there is no DTC
        x=[eta_d/(eta_c+eta_d)*ones(numsim,1);x]; % x now combines the (constant) allocation of scientists with the tax rate
    else
        [x,fval,exitflag] = fmincon(@(x)mysimopttaxnew2onlytau(x, Ac0, Ad0, S0),x0,[],[],[],[],lb,ub,[],options); % optimization with DTC and carbon tax only
        % mysimopttaxnew2onlytau computes the utility given a vector of input tax x and deriving the corresponding allocation for scientists
        % if several allocation of scientists are equilibria, the interior one is chosen
    end
end
if exitflag<=0
    display('*********problem, optimization not converged******************')
end


%********************************************************************
%********************************************************************
%7. take results
Util=utcomp-1/(1+rho)^delay*fval;
if mode == 0.5
    Utilinit=-mysimopttaxnew2([sce3rho15;M2],Ac0,Ad0,S0);
    Resp = mysimenvtaxnew2onlytau(x, Ac0, Ad0, S0);
    % mysimenvtaxnew2onlytau computes all the relevant parameters of the economy given a path of input taxes: if there are several solution for the allocation of scientists, the interior one is picked
    xx=[Stocksc;Resp.S_c];
    Acc=[StockAc;Resp.Ac];
    Add=[StockAd;Resp.Ad];
    tauu=[Stocktau;Resp.tau];
    Cc=[StockC;Resp.C];
    Ss=[StockS;Resp.S];
else
    Resp = mysimenvtaxnew2(x, Ac0, Ad0, S0);
    % mysimenvtaxnew2 computes all the relevant parameters of the economy given the vector x which stacks the share of scientists in clean research and the input tax
    xx=[Stocksc;x(1:numsim)];
    Acc=[StockAc;Resp.Ac];
    Add=[StockAd;Resp.Ad];
    Qq=[StockQ;Resp.Q]; % Qq is the clean research subsidy, it makes sense only if mode=1
    tauu=[Stocktau;Resp.tau];
    Cc=[StockC;Resp.C];
    Ss=[StockS;Resp.S];
end
Tt=zeros(numsim+delay,1);
Ratio=zeros(numsim+delay,1);
for i=1:(numsim+delay)
    Tt(i)=3/log(2)*log((280*2^(t_disaster/3)-Ss(i))/280);
    Ratio(i)=1/(1+(Add(i)/Acc(i))^(epsilon*(1-alpha))*(1+tauu(i))^(-epsilon));
end


%********************************************************************
%********************************************************************
%8. plot results

numsimc=numsim+delay;
subplot(2,2,1)
plot(linspace(1,numsimc,numsimc),[xx(1:numsimc),ones(numsimc,1)-xx(1:numsimc)]);
xlabel('period')
ylabel('Fraction of Scientists')
legend('S_c','S_d')
subplot(2,2,2)
plot(linspace(1,numsimc,numsimc),[Acc(1:numsimc),Add(1:numsimc),Cc(1:numsimc)]);%,
xlabel('period')
ylabel('economy')
legend('Ac','Ad','C')%
subplot(2,2,3)
plot(linspace(1,numsimc,numsimc),tauu(1:numsimc));
xlabel('period')
ylabel('input tax')
legend('tau')
subplot(2,2,4)
plot(linspace(1,numsimc,numsimc),Tt(1:numsimc));
xlabel('period')
ylabel('temp')
legend('delta')

