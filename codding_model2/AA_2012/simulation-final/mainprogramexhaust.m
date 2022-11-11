%******* This is the main program in the exhaustible case
%1. Fix the parameters in the parameter section
%2. Choose an initial guess for the allocation of researchers, the input tax and the resource tax.
% The program then computes the following vectors:
% - xx: the share of scientists in the clean sector
% - Acc and Add: the productivity of the clean and dirty sector respectively
% - cleansubss: the minimum clean research subsidy necessary to implement the optimal policy, when we assume that when multiple equilibria are possible for the allocation of scientists, the interior one is chosen.
% - tauu: the carbon tax
% - Cc: the consumption
% - Ss: environmental quality
% - thetaa: the resource tax
% - Qq: the stock of resource 
% - Tt: the increase in temperature since preindustrial times
% - Ratio: the ration of clean to dirty inputs in production
% - Util: is the total utility at time 0
% Moreover the program plots the share of scientists, the path of Acc, Add,and Cc, the input tax, the increase in temperature, the resource tax and the evolution of the stock of the resource
% *********CAREFUL*********: Qq is now the stock of resource, the clean research subsidy is now given by cleansubss



%clear all; % uncomment to clear memory


%********************************************************************
%********************************************************************
%1. parameters
global dt rho sigma psi alpha gamma eta_d eta_c qsi epsilon delta numsim phi S_bar lambda t_disaster max_t kappa alpha1 alpha2 k1 k2 phi1 pf
dt=5;  % number of years in a period
rho = 0.001*dt; % discount rate
epsilon =10;  % elasticity of substitution, must be greater than 1 for the program to work properly.
sigma = 2; % Coefficient of relative risk aversion 
alpha = 1/3;% 1/3; % share of machines in production
gamma = 1; % size of innovation
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
delta =.5*emission0/S0; % half of CO2 emissions are absorbed at current atmospheric levels
% Parameters linked to the resource: see file calibrationexhaustible.tex
R00=375.8307453; % resource production in 2002 - 2006 (in billions of barrels)
Q00=4178.829103; % stock of resource in 2002 (in billion of barrels)
Q0=Q00-R00; % stock of resource in 2006 (in billion of barrels)
proil0 = 1000000*(-394.726*Q00 +0.0421224*Q00^2 + 942115); % price of resource in 2002 in 2000$
p0c0= 1.86605E+14; % GDP in 2002 - 2006 in 2000$
spendRgdp= (proil0*R00)/p0c0; % ratio spending on resource over GDP
alpha2= (1-alpha)*(1+(Yd0/Yc0)^(1/epsilon-1))*spendRgdp;
alpha1=alpha-alpha2;
eta_d = (1-alpha)/(1-alpha1)*0.02*dt; % probability of success in the dirty sector
eta_c = 0.02*dt; % probability of success in the the clean sector

%********************************************************************
%********************************************************************
%2. options
display_iter=0;    %=1 if you want to see iterations, ~1 if otherwise
display_diag=1;    %=1 if you want to see diagnostics, ~1 if otherwise

%********************************************************************
%********************************************************************
%3. initial guesses
sc0 =exsce10rho01; % initial guess for the share of scientists in clean research
t0=extaue10rho01; % initial guess for the input tax
theta0=exthetae10rho01; % initial guess for the resource tax
x0=[sc0; t0; theta0]; % vector stacking the guesses for the share of scientists in clean research, the input tax and the resource tax.
ub=[ones(numsim,1); 100*ones(numsim,1);1000*ones(numsim,1)]; % upper bound
lb=[zeros(numsim,1); zeros(numsim,1);zeros(numsim,1)]; % lower bound

%**** No intervention in the code is necessary beyond this point

%********************************************************************
%********************************************************************
% 4. Warning
if epsilon <= 1
    display ('*********the program is not stable for epsilon <= 1******************')
end

%********************************************************************
%********************************************************************
% 5. preliminary computations
psi= alpha^2; % cost of machines (the expressions for Ad0 and Ac0 assume this relationship)
phi=(1-alpha)*(1-epsilon);
phi1=(1-alpha1)*(1-epsilon);
kappa = (1-alpha)/(1-alpha1)*(psi^alpha2*alpha1^alpha1*alpha2^alpha2/alpha^alpha)^(1-epsilon);
k1 = ((psi^alpha2)*(alpha1^alpha1)*(alpha2^alpha2))^(1-epsilon);
k2 = (psi^(-alpha1/(1-alpha)))*(alpha^(alpha/(1-alpha)))*(alpha1^(alpha1/(1-alpha)))*(alpha2^(alpha2/(1-alpha)));
pf=alpha2^(-1)*proil0*R00/Yd0*(1+(Yd0/Yc0)^(1/epsilon-1))^(1/(1-epsilon));
Ad0 = (alpha1^alpha1*alpha2^alpha2/(alpha^(2*alpha1)))^(-1/(1-alpha1))*Yd0^((1-alpha)/(1-alpha1))*(1+(Yd0/Yc0)^(1/epsilon-1))^((alpha+phi)/phi1)*cost(Q00)^(alpha2/(1-alpha1));
Ac0=(alpha1^alpha1*alpha2^alpha2/alpha^(alpha1-alpha2)*cost(Q00)^(-alpha2)*Ad0^(1-alpha1)*(Yc0/Yd0)^(1/epsilon))^(1/(1-alpha));
qsi = emission0/Yd0;
StockAc=[];
StockAd=[];
StockC=[];
Stocktau=[];
Stocksc=[];
StockS=[];
StockQ=[];
Stockprofitq=[];
Stocktheta=[];
utcomp=0;


%********************************************************************
%********************************************************************
%4. optimization
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
[x,fval,exitflag] = fmincon(@(x)mysimopttaxnewexhaust2(x, Ac0, Ad0, S0, Q0),x0,[],[],[],[],lb,ub,[],options);
 % mysimopttaxnew2 computes the utility given x - a vector stacking the path for the share of scientists in clean research, the input tax and the resource tax- and the initial values of clean and dirty productivities and quality of the environment
if exitflag<=0
    display('*********problem, optimization not converged******************')
end


%********************************************************************
%********************************************************************
%5. take results
Util=utcomp-1/(1+rho)^delay*fval;
Resp = mysimenvtaxnewexhaust2(x, Ac0, Ad0, S0, Q0);
% mysimenvtaxnewexhaust2 computes all the relevant parameters of the economy given the vector x which stacks the share of scientists in clean research, the input tax and the resource tax
xx=[Stocksc;x(1:numsim)];
Acc=[StockAc;Resp.Ac];
Add=[StockAd;Resp.Ad];
Qq=[StockQ;Resp.Q];
cleansubss=[Stockprofitq;Resp.cleansubs];
tauu=[Stocktau;Resp.tau];
thetaa=[Stocktheta;Resp.theta];
Cc=[StockC;Resp.C];
Ss=[StockS;Resp.S];
Tt=zeros(numsim+delay,1);
Ratio=Tt;
for i=1:(numsim+delay)
    Tt(i)=3/log(2)*log((280*2^(t_disaster/3)-Ss(i))/280);
    Ratio(i)=1/(1+(Add(i)/Acc(i))^(epsilon*(1-alpha))*(1+tauu(i))^(-epsilon));
end

%********************************************************************
%********************************************************************
%6. plot results
subplot(3,2,1)
plot(linspace(1,numsim,numsim),[xx(1:numsim),ones(numsim,1)-xx(1:numsim)]);
xlabel('period')
ylabel('Fraction of Scientists')
legend('S_c','S_d')
subplot(3,2,2)
plot(linspace(1,numsim,numsim),[Acc(1:numsim),Add(1:numsim),Cc(1:numsim)]);%,
xlabel('period')
ylabel('economy')
legend('Ac','Ad','C')
subplot(3,2,3)
plot(linspace(1,numsim,numsim),tauu(1:numsim));
xlabel('period')
ylabel('input tax')
legend('tau')
subplot(3,2,4)
plot(linspace(1,numsim,numsim),Tt(1:numsim));
xlabel('period')
ylabel('temp')
legend('delta')
subplot(3,2,5)
plot(linspace(1,numsim,numsim),thetaa(1:numsim));%);
xlabel('period')
ylabel('resource taxe')
legend('theta')
subplot(3,2,6)
plot(linspace(1,numsim,numsim),Qq(1:numsim));%,
xlabel('period')
ylabel('stock')
legend('Q')