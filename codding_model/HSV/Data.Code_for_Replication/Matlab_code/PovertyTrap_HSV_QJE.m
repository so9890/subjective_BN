% %----------------------------------------------------
% POVERTY TRAP MODEL with transition
% - ASSUMES \gamma=beta
% This file plots social welfare as a function of \tau
%-----------------------------------------------------
clear;


%------------------
%Basic parameters
%------------------

eta     = 1.0;          % parameter of the exponential for cost distribution (normalization?)

%Externally set values
UStau   = 0.181;        % US progressivity rate 2000-2006: 0.151 FROM PSID, 0.20 from CBO with adjustment; 
sigma   = 2.0;          % curvature of disutility of work (from framework) 
rstar   = 0.04;         % Interest rate
g       = 0.189;        % fraction of output devoted to (valued) expenditures;

%From Framework estimates (average 2000-2006)
veps    = 0.164;        % cross-sectional variance of insurable wage shocks 
theta   = 3.124  %3.402 % elasticity of substitution in production & 1/theta^2 is cons. ineq. at age 0
vphi    = 0.036;        % cross-sectional variance of preference heterogeneity    


%Cross-sectional residual uninsurable inequality
valpha  = 0.098;        % valpha = valpha0+delta/(1-delta)*vomega;
valpha0  = 0;           % Initial residual uninsurable inequality = zero       
delta   = 0.971;        % survival rate

psi = 0.65;                % New parameter defining elasticity of skill investment

%------------------
%Derivative parameters
%------------------

taum1 = UStau;          % Define the tau inherited from the past  
chi   = g/(1-g);        % weight on utility from G (used to be 0.25)
                        % chi = 0.057;   value for chi that rationalizes tau = 0.151

vomega = (valpha-valpha0)*(1-delta)/delta;

rho   = rstar +(1-UStau)*(2-UStau)*vomega/2;    % Rate of time preference
beta  = 1/(1+rho);                              % Discount factor


%%%%%%%%%%%
%% Parameters specific to the poverty trap model
gamma = beta;                                      % Weight the planner uses to discount future generations

% Set the median low-skilled wage equal to 50% of the median wage of 
% skilled workers in the economy with tau=UStau
% - equivalent to the median earnings of high school dropouts being 50% of
% the median earnings of workers with at least high school
w0 = 0.54 * exp( 1/(theta*(1+psi)) + ...
     (1/(theta-1))*log(theta/(theta-1))+ ...
     (psi/((theta-1)*(1+psi)))*log((1- UStau)/theta)-(1/((theta-1)*(1+psi)))*log(eta) );

% Probability of unskilled parents having unskilled children
PrUU = 0.355;
% Probability of skilled parents having unskilled children
PrSU = 0.135;
 
% The implied steady-state share of unskilled workers
xi_US = PrSU / (1+PrSU-PrUU);
	
% Define some useful functions: 
% Base-price for skills
pi0 = @(tau) (1/(theta-1))*log(theta/(theta-1))+(psi/((theta-1)*(1+psi)))*log((1-tau)/theta)-(1/((theta-1)*(1+psi)))*log(eta);

% log of lambda(tau,xi0)
loglambda = @(tau,xi0) log(1-chi/(1+chi)) + tau/(1+sigma)*log(1-tau) ...
              + tau*(1-tau)/(sigma+tau)*(2+sigma+(1+sigma)/(sigma+tau))*veps/2 ...
              + tau*(1-tau)*vphi/2 ...
              + tau*(1-tau)*valpha0/2+log( ( 1-delta*exp( -tau*(1-tau)*vomega/2 ) )/(1-delta) )...
              + log(xi0*w0 + (1-xi0)*theta/(theta-1) *exp(pi0(tau))) ...
              - log(xi0*(w0)^(1-tau) + (1-xi0)*theta/(theta-1+tau) *exp((1-tau) * pi0(tau)));
   
% skilled/unskilled children treshold for stochastic variable alpha-phi + eta/theta * kappa
zeta = @(tau,LL,c_bar) 1/(1-tau) * (log(c_bar) -(LL + (1-tau)/(1+sigma)*log(1-tau) ...
              + (1-tau)*(1-2*tau-sigma*tau)/(sigma+tau) *veps/2));


% find the consumption threshold and the share X (which always are unskilled)
% consistent with Pr(unskilled|unskilled) and Pr(skilled|unskilled) given UStau 
    
stopcrit = 1;
%initial guess for c_bar
c_bar=w0;
c_hi = 10*w0;
c_lo = 0.001*w0;
while stopcrit > 1e-7
   c_bar = (c_lo + c_hi)/2;

   % some useful variables:
   LL = loglambda(UStau,xi_US);
   u = theta * (zeta(UStau,LL,c_bar) - pi0(UStau) + (valpha+vphi)/2);
   v = theta * sqrt(valpha+vphi);
        
   % compute P(c<c_bar | skilled) using the forumla for the cdf of an
   % exponentially modified normal distribution:
   PrCskilled = normcdf(u,0,v) - exp(-u + v^2/2 + log(normcdf(u,v^2,v)));
   
   % compute P(c<c_bar | UNskilled) using the cdf of a normal
   PrC_UNskilled = normcdf(zeta(UStau,LL,c_bar) - log(w0), ...
           -(valpha+vphi)/2, sqrt(valpha+vphi));

   stopcrit = abs( (1-PrC_UNskilled)/(1-PrCskilled) - (1-PrUU)/(1-PrSU));
   % update guess and bisect
   if ((1-PrC_UNskilled)/(1-PrCskilled) > (1-PrUU)/(1-PrSU))
       c_lo = c_bar;
   else
      c_hi = c_bar;
   end
   %[c_lo c_bar c_hi PrC_UNskilled/PrCskilled]
end   
% calculate the implied X
X= 1 - (1-PrSU)/(1-PrCskilled)


%Decomposition of SWF, given X and c_bar
N=2000;  % number of different tau(i) to evaluate
T=2000;  % number of periods to simulate
Mscript  =zeros(N,1);
sigmahat =zeros(N,1);
tau      =zeros(N,1);
xi       =zeros(N,T);
W_SKILL1 =zeros(1,N);
W_SKILL2 =zeros(1,N);
W_EDU    =zeros(1,N);
W_SKILL1_utopic = zeros(1,N);
W_SKILL2_utopic = zeros(1,N);
W_EDU_utopic    = zeros(1,N);
for i=1:N
    % calculate implied long-run share of unskilled, given different tau
    tau(i)= -0.4 + 1/N *i;
	
    sigmahat(i)=(sigma+tau(i))/(1-tau(i));

    %M term in allocations
    Mscript(i)=( (1-tau(i))*(1-tau(i)*(1+sigmahat(i)))/sigmahat(i) )*(veps/2);
    
  
    % Simulate a long series of periods given tau(i) and
    % starting from an initial guess xi_US
    xi(i,1) = xi_US;  
    xi_utopic(i,1) = xi_US;  
    % education costs of all old individuals in the initial period
    W_EDU(i) = -(1-xi(i,1)) *(1-gamma)/gamma * beta*delta/(1-beta*delta) *psi/(1+psi) *(1-tau(i))/theta;
    W_EDU_utopic(i) = W_EDU(i);
    for t=1:T 
        % calculate lambda for this period
        %LL = loglambda(tau(i),xi(i,t));
        
        % some useful variables
        u = theta * (zeta(tau(i),loglambda(tau(i),xi(i,t)),c_bar) ...
                     - pi0(tau(i)) + (valpha+vphi)/2);
        v = theta * sqrt(valpha+vphi);
        
        % compute P(c<c_bar | skilled) using the forumla for the cdf of an
        % exponentially modified normal distribution:
        PrCskilled = normcdf( u,0,v) - exp(-u + v^2/2 + log(normcdf(u,v^2,v))) ;
        
        % compute P(c<c_bar | UNskilled) using the cdf of a normal
	    PrC_UNskilled =normcdf( zeta(tau(i),loglambda(tau(i),xi(i,t)),c_bar) -log(w0),...
                                 -(valpha+vphi)/2, sqrt(valpha+vphi));

		% calculate the next-period share of unskiled, xi_{t+1}
		%xi(i) = PrCskilled / (1 - (PrC_UNskilled - PrCskilled));
		xi(i,t+1) = delta*xi(i,t) ...
              + (1-delta)*xi(i,t)    *(X +(1-X)*PrC_UNskilled) ...
              + (1-delta)*(1-xi(i,t))*(X +(1-X)*PrCskilled);
        
        % Calculate SWF piece associated to productivity gain from skill investment, 
        % when xi is ENDOGENOUS (add up the discounted gain)
        Temp  = (1+chi)* log( xi(i,t) *w0 + (1-xi(i,t))*theta/(theta-1) *exp(pi0(tau(i))));
        W_SKILL1(i) = W_SKILL1(i) + (1-gamma)*gamma^(t-1)*Temp;
    
        % Calculate SWF piece associated to consumption inequality when xi is ENDOGENOUS:
        Temp  =  - log( xi(i,t) *w0^(1-tau(i)) ...
                          + (1-xi(i,t))* exp((1-tau(i))*pi0(tau(i))) *theta/(theta-1+tau(i)) ) ...
            + (1-xi(i,t))*(1-tau(i))*( pi0(tau(i)) + 1/theta) + xi(i,t)*(1-tau(i))*log(w0);
        W_SKILL2(i) = W_SKILL2(i) + (1-gamma)*gamma^(t-1)*Temp;

        % Calculate SWF piece associated to skill investment costs when xi is ENDOGENOUS:
        Temp  = -(1-xi(i,t)) *(gamma-beta*delta)/(gamma *(1-beta*delta)) *psi/(1+psi) *(1-tau(i))/theta;
        W_EDU(i) = W_EDU(i) + (1-gamma)*gamma^(t-1)*Temp;
        
        
        %%% AND NOW THE UTOPIC CASE WITH c_bar=0
        % calculate next-period share of unskilled in UTOPIC case: c_bar=0
        % and xi_utopic is ENDOGENOUS
        xi_utopic(i,t+1) = delta*xi_utopic(i,t) + (1-delta)* X;
        % Calculate SWF piece associated to productivity gain from skill investment, 
        Temp  = (1+chi)* log( xi_utopic(i,t) *w0 + (1-xi_utopic(i,t))*theta/(theta-1) *exp(pi0(tau(i))));
        W_SKILL1_utopic(i) = W_SKILL1_utopic(i) + (1-gamma)*gamma^(t-1)*Temp;    
        % Calculate SWF piece associated to consumption inequality:
        Temp  =  - log( xi_utopic(i,t) *w0^(1-tau(i)) ...
                          + (1-xi_utopic(i,t))* exp((1-tau(i))*pi0(tau(i))) *theta/(theta-1+tau(i)) ) ...
            + (1-xi_utopic(i,t))*(1-tau(i))*( pi0(tau(i)) + 1/theta) + xi_utopic(i,t)*(1-tau(i))*log(w0);
        W_SKILL2_utopic(i) = W_SKILL2_utopic(i) + (1-gamma)*gamma^(t-1)*Temp;
        % Calculate SWF piece associated to skill investment costs:
        Temp  = -(1-xi_utopic(i,t)) *(gamma-beta*delta)/(gamma *(1-beta*delta)) *psi/(1+psi) *(1-tau(i))/theta;
        W_EDU_utopic(i) = W_EDU_utopic(i) + (1-gamma)*gamma^(t-1)*Temp;
    end
    
    % calculate the social welfare function
    % SWF of the representative agent with chi=0 (and g=0)
    W_RA_NOG(i)    = ( 1/(1+sigma) )*( log(1-tau(i)) - (1-tau(i)) );
    
    % SWF of the representative agent with valued G
    % Note that log(1+g)+chi log g = chi*log(chi) - (1+chi)*log(1+chi)
    W_RA(i)   = W_RA_NOG(i) + chi*log(chi) - (1+chi)*log(1+chi) + chi*( log(1-tau(i))/(1+sigma) );

    %SWF piece associated to preference heterogeneity (second line for right-wing planner)
    W_PREF(i)      = -(1-tau(i))^2*vphi/2;
    
    %SWF piece associated with within skill-group consumption inequality
    W_OMEGA(i)      = log( ( 1-delta*exp( -tau(i)*(1-tau(i))*vomega/2 ) )/(1-delta) ) ...
                         - (1-tau(i))*beta*delta/(gamma-beta*delta)*vomega/2;
    
    %SWF piece associated to initial variance of uninsuranle shocks
    W_ALPHA0(i) = -(1-tau(i))^2*valpha0/2;

    %SWF piece associated to insurable shocks (second line for chi=0 (and g=0))
    W_EPS(i)  = (1+chi)*(1-tau(i))/(sigma+tau(i))*veps - (1+chi)*sigma*( (1-tau(i))/(sigma+tau(i)) )^2*veps/2 ;
    W_EPS_NOG(i)  = (1-tau(i))/(sigma+tau(i))*veps - sigma*( (1-tau(i))/(sigma+tau(i)) )^2*veps/2 ;
    
    %SWC piece associated to insurable shock for progressive consumption tax
    W_CONSTAX_EPS(i)  = ((1+chi)/sigma)*veps - (1+chi)*sigma*(1/sigma^2)*(veps/2) ;

    % SWF piece associated to productivity gain from skill investment when xi=EXOGENOUS:
    W_SKILL1_EX(i) =(1+chi)* log( xi_US *w0 + (1-xi_US)*theta/(theta-1) ...
                   * exp(pi0(tau(i))));

    % SWF piece associated to consumptoin inequality and skill investment costs
    % when xi=xi_US is EXOGENOUS:
    W_SKILL2_EX(i)  =  - log( xi_US *w0^(1-tau(i)) + (1-xi_US)* exp((1-tau(i))*pi0(tau(i))) *theta/(theta-1+tau(i)) ) ...
         + (1-xi_US)*(1-tau(i))*( pi0(tau(i)) + 1/((1+psi)*theta)) + xi_US*(1-tau(i))*log(w0);
end

%-----------------------------------------
%Aggregation of components of SWF
%and calculation of optimal taus
%----------------------------------------

%Piece together terms involving skill
W_SKILL        = W_SKILL1    + W_SKILL2 + W_EDU;
W_SKILL_EX     = W_SKILL1_EX + W_SKILL2_EX;
W_SKILL_utopic = W_SKILL1_utopic + W_SKILL2_utopic + W_EDU_utopic;
 
%Piece together initial dispersion and lifetime shocks
W_ALPHA = W_ALPHA0+W_OMEGA;

%Piece together all components (second line is approximation of lifetime cumulation of Vomega)
W_TOT        = W_RA + W_SKILL        + W_PREF + W_ALPHA + W_EPS;
W_TOT_EX     = W_RA + W_SKILL_EX     + W_PREF + W_ALPHA + W_EPS;
W_TOT_utopic = W_RA + W_SKILL_utopic + W_PREF + W_ALPHA + W_EPS;



%SWF for progressive consumption tax
W_CONSTAX_TOT    = W_RA + W_SKILL    + W_PREF +W_ALPHA +W_CONSTAX_EPS;
W_CONSTAX_TOT_EX = W_RA + W_SKILL_EX + W_PREF +W_ALPHA +W_CONSTAX_EPS;

%SWF for case when chi=0 (and g=0)
%W_TOT_NOG       = W_RA_NOG + W_SKILL1_NOG + W_SKILL2 + W_PREF + W_ALPHA + W_EPS_NOG;


%-----------------------------
%Compute implied tax system
%-----------------------------

% Representative Agent
disp('RA taustar');
[c,i]=max((W_RA)');
optau_RA = tau(i)

% Piece associated to skills investment
disp('Skill Investment taustar');
[c,i]          = max((W_SKILL)');
optau_SKILL    = tau(i)
[c,i]          = max((W_SKILL_EX)');
optau_SKILL_EX = tau(i)

tau_iUS = sum(tau<=UStau);
% Optimal tau when xi is ENDOGENOUS
disp('Exact taustar ENDOGENOUS xi');
[c,i]=max(W_TOT');
optau_EXACT = tau(i)
disp('Welfare gain from tax reform');
exp(max(W_TOT')-W_TOT(tau_iUS))-1

% Optimal tau when xi is EXOGENOUS
disp('Exact taustar EXOGENOUS xi');
[c,i]=max(W_TOT_EX');
optau_EXACT_EX = tau(i)
disp('Welfare gain from tax reform');
exp(max(W_TOT_EX')-W_TOT_EX(tau_iUS))-1

% Optimal tau when xi is based on c_bar=0 (UTOPIC case)
disp('Exact taustar UTOPIC xi');
[c,i]=max(W_TOT_utopic');
optau_EXACT_utopic = tau(i)
disp('Welfare gain from tax reform in utopic xi case (c_bar=0)');
exp(max(W_TOT_utopic)-W_TOT_utopic(tau_iUS))-1
disp('Welfare gain from obtaining utopic xi with benchmark tau_US');
exp(W_TOT_utopic(tau_iUS)-W_TOT(tau_iUS))-1



disp('Exact taustar RA+SKILL');
[c,i]=max((W_RA+W_SKILL)');
optau2_EXACT = tau(i)

disp('Exact taustar RA+SKILL+PREF');
[c,i]=max((W_RA+W_SKILL+W_PREF)');
optau3_EXACT = tau(i)

disp('Exact taustar RA+SKILL+PREF+ALPHA');
[c,i]=max((W_RA+W_SKILL+W_PREF+W_ALPHA)');
optau4_EXACT = tau(i)

% Consumption tax
disp('Consumption tax taustar');
[c,i]=max(W_CONSTAX_TOT');
optau_CONS = tau(i)

exp(max(W_CONSTAX_TOT')-W_TOT(tau_iUS))-1



%-----------------------------
%Plots of SWF
%-----------------------------

%Max of welfare function, needed for normalization
welf_max          = max(W_TOT);
welf_max_CONSTAX  = max(W_CONSTAX_TOT);
%welf_max_NOG      = max(W_TOT_NOG);
%welf_max_XG2      = max(W_TOT_XG2(1:1500));
%welf_max_XG3      = max(W_TOT_XG3(1:1500));


save poverty_trap_transition
keyboard

%%Plot RA component of Social Welfare Function
%figure(2);
%plot(tau,W_RA_NOG,tau,W_RA,'Linewidth',2.5);
%grid on;
%xlim([-0.25 0.6]);
%legend('RA G non-valued','RA G valued','Location','SouthEast','fontsize',18);
%xlabel('Progressivity rate (\tau)', 'fontsize', 18);
%ylabel('Social Welfare', 'fontsize', 18);
%title('Representative Agent Economies', 'fontsize', 18);

%%Plot Pref. het. and uninsurable risk component of Social Welfare Function
%figure(4);
%plot(tau,100*(exp(W_PREF)-1),tau,100*(exp(W_ALPHA)-1),tau,100*(exp(W_EPS)-1),'Linewidth',2.5);
%legend('Pref Het','Unins. Shocks','Insurable Shocks','Location','SouthEast', 'fontsize', 18);
%grid on;xlim([-0.25 0.6]);
%xlabel('Progressivity rate (\tau)', 'fontsize', 18);ylabel('welfare contrib. (% of cons.)', 'fontsize', 18);
%title('Pref heterogeneity, Efficiency Shocks', 'fontsize', 18);

%Plot components of Social Welfare Function adding sequentially
figure(5);
plot(tau,100*(exp(W_RA-welf_max)-1),tau,100*(exp(W_RA+W_SKILL-welf_max)-1),tau,100*(exp(W_RA+W_SKILL+W_PREF-welf_max)-1),...
     tau,100*(exp(W_RA+W_SKILL+W_PREF+W_ALPHA-welf_max)-1),tau,100*(exp(W_TOT-welf_max)-1),'Linewidth',2.5);
legend('Rep. Agent','+Skill Het.','+Pref Het','+Unins. Shocks','+Ins. Shocks = Baseline','Location','Best')
grid on;xlim([-0.3 0.6]); ylim([-30,2]);
xlabel('Progressivity rate (\tau)', 'fontsize', 14);
ylabel('welf change rel. to baseline optimum (% of cons.)', 'fontsize', 14);
title('Social Welfare Function', 'fontsize', 16);

%Plot final Social Welfare Function, including when xi is exogenous
figure(6);
plot(tau,100*(exp(W_TOT-welf_max)-1),tau,100*(exp(W_TOT_EX-welf_max)-1),'Linewidth',2.5);
legend('Baseline (endog. \xi)','Exogenous \xi','Location','Best');
grid on;xlim([-0.25 0.6]); ylim([-30,2]);
xlabel('Progressivity rate (\tau)', 'fontsize', 14);
ylabel('welf change rel. to optimum (% of cons.)', 'fontsize', 14);
title('Social Welfare Function', 'fontsize', 16);

%Progressive consumption tax
figure(11);
plot(tau,100*(exp(W_TOT-welf_max)-1),...
     tau,100*(exp(W_CONSTAX_TOT-welf_max)-1),'Linewidth',2.5);
legend('Baseline Model','Consumption Tax','Location','Best');
grid on;xlim([0 0.4]); ylim([-4,1]);
xlabel('Progressivity rate (\tau)', 'fontsize', 14);
ylabel('welf change rel. to optimum (% of cons.)', 'fontsize', 14);
title('Social Welfare Function', 'fontsize', 16);

% Fixed points for share of unskilled
figure(12)
hold off
plot(tau,xi(:,T),'Linewidth',2.5);
grid on;
xlim([0.05 0.6]);ylim([0.13 0.2]);
hold on
[c,j2]=min((tau-UStau).^2);
plot(tau(j2),xi(j2,T),'*r');
[c,i]=max(W_TOT');
plot(tau(i),xi(i,T),'*r');
xlabel('Progressivity rate (\tau)', 'fontsize', 14);
ylabel('steady state \xi', 'fontsize', 14);
title('Steady states for share of unskilled (\xi)', 'fontsize', 16);

% Dynamics of xi
figure(2)
hold off
[c,i]=max(W_TOT');
plot(xi(i,1:200),'Linewidth',2.5);
grid on;
%xlim([0.05 0.5]);ylim([0.125 0.23]);
%hold on
xlabel('Years', 'fontsize', 14);
ylabel('Share of unskilled, \xi_t', 'fontsize', 14);
title('Dynamics of share of unskilled (\xi_t)', 'fontsize', 16);

figure(1)
N=1000000;
bins=-2:0.02:2;
[c,i]=min((tau-UStau).^2);
lnC = loglambda(tau(i),xi(i)) + (1-tau(i))/(1+sigma) * log(1-tau(i)) + ...
    (1-tau(i))*(1-2*tau(i)-sigma*tau(i))/(sigma+tau(i))*veps/2;
% distribution of log consumption for the unskilled
TempU=  lnC + (1-tau(i))* ...
    (log(w0) -(valpha+vphi)/2 + sqrt(valpha+vphi)*randn(N,1));
TempUlnYmed = log(w0) -(valpha+vphi)/2 - (sigma+1)/(sigma+UStau) *veps/2 ...
           + log(1-UStau)/(1+sigma) ...
           - (1-UStau)*(1-2*UStau-sigma*UStau)/(sigma+UStau)^2 *veps/2;
distU=hist(TempU,bins);
% distribution of log uninsurable income for the skilled
TempS=  lnC + (1-tau(i))* (pi0(tau(i)) ...
    -(valpha+vphi)/2 + sqrt(valpha+vphi)*randn(N,1) ...
    + exprnd(1/theta,N,1));
distS=hist(TempS,bins);

plot(bins,distU,bins,distS,'Linewidth',2.5);
legend('Unskilled','Skilled','Location','Best');
hold on
plot(log(c_bar)*[1 1],[0 max(max(distU),max(distS))+500],':k')
xlim([-1.98 1.98]); 
ylim([-10 (max(max(distU),max(distS))+100)]);
xlabel('log consumption', 'fontsize', 14);
ylabel('density', 'fontsize', 14);
title('Skill-specific consumption distribution for \tau=0.181', 'fontsize', 14);

figure(3)
subplot(2,1,1)
plot(bins,distU,'Linewidth',2.5);
hold on
plot(bins,distS,'--r','Linewidth',2.5);
legend('Unskilled','Skilled','Location','Best');
hold on
plot(log(c_bar)*[1 1],[0 max(max(distU),max(distS))+500],':k')
xlim([-1.98 1.98]); 
ylim([-10 (max(max(distU),max(distS))+1000)]);
xlabel('log consumption', 'fontsize', 14);
ylabel('density', 'fontsize', 14);
title('Skill-specific consumption distribution for \tau=0.181', 'fontsize', 14);

% Dynamics of xi
%subplot(2,1,2)
%hold off
%[c,i]=max(W_TOT');
%plot(xi(i,1:200),'Linewidth',2.5);
%grid on;
%xlabel('Years', 'fontsize', 12);
%ylabel('Share of unskilled, \xi_t', 'fontsize', 12);
%title('Dynamics for the share of unskilled (\xi_t)', 'fontsize', 14);



% generate empirical consumption distribution
ConsDist=[TempU(1:ceil(xi_US*N)); TempS(ceil(xi_US*N)+1:N)];
[c_bar/median(exp(ConsDist)) c_bar/mean(exp(ConsDist))]

% print stuff
%figure(1)
%print -depsc cons_distn
%figure(5) 
%print -depsc PT_SWF_dcmp
%figure(6) 
%print -depsc PT_SWF_xi
%figure(11) 
%print -depsc PT_SWF_cons
%figure(12)
%print -depsc PT_FP_zoom
%figure(2)
%print -depsc PT_xi_t
%figure(3)
%print -depsc PT_2panelB

