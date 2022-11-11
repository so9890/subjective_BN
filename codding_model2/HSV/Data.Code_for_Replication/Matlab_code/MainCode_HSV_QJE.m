% % ----------------3.12------------------------------------
% THIS IS THE MAIN FILE FOR PRODUCING THE FIGURES AND NUMERICAL RESULTS IN
% THE PAPER
% THE ONLY VERSION OF THE MODEL NOT COMPUTED HERE IS THE POVERTY TRAP
% VERSION WITH INVESTMENT CONSTRAINTS
%-----------------------------------------------------
 clear;

%pause on;

%-------------------------------------------------------
%Model Parameterization
%-------------------------------------------------------

%Indicator for Cases
% THE POLITICAL ECONOMY SPECIFICATION, SECTION 6F IN THE PAPER 
Polecon         = 'N';
% THE ALTERNATIVE CALIBRATION IN THE SECOND LINE OF TABLE 1 
Alternative_Cal = 'N';
% THE SPECIFICATION WITH FIXED INVESTMENT, SECTION 6C IN THE PAPER
Irreversible    = 'N';

%------------------
%Basic parameters
%------------------

eta     = 1.0;          % parameter of the exponential for cost distribution (normalization?)

%Externally set values
UStau   = 0.181;        % US progressivity rate 2000-2006: 0.181 FROM PSID, 0.20 from CBO with adjustment; 
sigma   = 2.0;          % curvature of disutility of work (from framework) 
rstar   = 0.04;         % Interest rate
g       = 0.189;        % fraction of output devoted to (valued) expenditures;

%From Framework estimates (average 2000-2006)
veps    = 0.164;        % cross-sectional variance of insurable wage shocks 
theta   = 3.124;        % elasticity of substitution in production & 1/theta^2 is cons. ineq. at age 0
vphi    = 0.036;        % cross-sectional variance of preference heterogeneity    


%Cross-sectional residual uninsurable inequality
valpha  = 0.098;        % valpha = valpha0+delta/(1-delta)*vomega;
valpha0  = 0;           % Initial residual uninsurable inequality = zero       
delta   = 0.971;        % survival rate

psi = 0.65;                % New parameter defining elasticity of skill investment

nul = 0.5;                 % Low value for inequality aversion parameter
nuh = 2;                   % High value for inequality aversion parameter

%-----------------------------
%Original draft values for basic parameters
%-----------------------------
%eta     = 1.0;          % parameter of the exponential for cost distribution (normalization?)
%UStau   = 0.151;        % US progressivity rate 2000-2006: 0.151 FROM PSID, 0.20 from CBO with adjustment; 
%sigma   = 2.165;        % curvature of disutility of work (from framework) 
%rstar   = 0.04;         % Interest rate
%g       = 0.189;        % fraction of output devoted to (valued) expenditures;
%veps    = 0.166;        % cross-sectional variance of insurable wage shocks 
%theta   = 3.144;        % elasticity of substitution in production & 1/theta^2 is cons. ineq. at age 0
%vphi    = 0.035;        % cross-sectional variance of preference heterogeneity    
%valpha  = 0.097;        % valpha = valpha0+delta/(1-delta)*vomega;
%valpha0  = 0;           % Initial residual uninsurable inequality = zero       
%delta   = 0.971;        % survival rate
%psi = 1;                % New parameter defining elasticity of skill investment

%------------------
%Derivative parameters
%------------------

taum1 = UStau;          % Define the tau inherited from the past  


chi   = g/(1-g);        % weight on utility from G (used to be 0.25)
                        % chi = 0.057;   value for chi that rationalizes tau = 0.151

vomega = (valpha-valpha0)*(1-delta)/delta;

rho   = rstar +(1-UStau)*(2-UStau)*vomega/2;    % Rate of time preference
beta  = 1/(1+rho);                              % Discount factor
gamma = beta;                                   % Weight the planner uses to discount future generations


%valpha=0;
%valpha0=0;
%vomega=0;
%vphi=0;
%veps=0;

% Following using to compute G in the US given g and UStau

%UStau=0.084;


sigmahat_UStau=(sigma+UStau)/(1-UStau);

logYUStau = (theta/(theta-1))*log(theta/(theta-1))+((psi/((theta-1)*(1+psi)))+(1/(1+sigma)))*log(1-UStau)-(psi/((theta-1)*(1+psi)))*log(theta)-(1/((theta-1)*(1+psi)))*log(eta)+...
    ((1/sigmahat_UStau^2)*(UStau+sigmahat_UStau+sigmahat_UStau*UStau)*(veps/2));

Y_UStau = exp(logYUStau);
           
agg_G   = g*Y_UStau;        % Needed to compute average output in the baseline optimum



%Alternative calibration with theta = 2 based on Pareto tail from PSID
if Alternative_Cal == 'Y'

    theta   = 2.000;        % elasticity of substitution in production & 1/theta^2 is cons. ineq. at age 0
    vphi    = 0.021;        % cross-sectional variance of preference heterogeneity    
    veps    = 0.139;        % cross-sectional variance of insurable wage shocks 
    
    vunins0 = (1/theta)^2;  % cross-sectional variance of uninsurable initial conditions 
    valpha   = 0.0;
    valpha0  = 0; 
    vomega  = 0.000;        % variance of uninsurable life-cycle shocks 

    rho   = rstar +(1-UStau)*(2-UStau)*vomega/2;
    beta = 1/(1+rho);
    gamma = beta;           % Define the weight the planner uses to discount future generations

    psi = 1;                % New parameter defining elasticity of skill investment
    
    logYUStau = (theta/(theta-1))*log(theta/(theta-1))+((psi/((theta-1)*(1+psi)))+(1/(1+sigma)))*log(1-UStau)-(psi/((theta-1)*(1+psi)))*log(theta)-(1/((theta-1)*(1+psi)))*log(eta)+...
    ((1/sigmahat_UStau^2)*(UStau+sigmahat_UStau+sigmahat_UStau*UStau)*(veps/2));

    Y_UStau = exp(logYUStau);

    agg_G   = g*Y_UStau;     % need to recompute average output in the baseline optimum  
 
end
    


%Decomposition of SWF
for i=1:1999;
    
    tau(i)=-1+0.001*i; % tax progressivity

    sigmahat(i)=(sigma+tau(i))/(1-tau(i));
    
    %Aggregate Output
 
    logY = (theta/(theta-1))*log(theta/(theta-1))+((psi/((theta-1)*(1+psi)))+(1/(1+sigma)))*log(1-tau(i))-(psi/((theta-1)*(1+psi)))*log(theta)-(1/((theta-1)*(1+psi)))*log(eta)+...
    ((1/sigmahat(i)^2)*(tau(i)+sigmahat(i)+sigmahat(i)*tau(i))*(veps/2));

    Y(i) = exp(logY);
    
    %M term in allocations
    Mscript(i)=( (1-tau(i))*(1-tau(i)*(1+sigmahat(i)))/sigmahat(i) )*(veps/2);
    
    %Base-price for skills
    pi0(i) = (1/(theta-1))*log(theta/(theta-1))+(psi/((theta-1)*(1+psi)))*log((1-tau(i))/theta)-(1/((theta-1)*(1+psi)))*log(eta);

    %Integral of y_{i}^{1-tau} across agents (w/o lambda term)
    Yd_1(i)=(1-tau(i))^((1-tau(i))/(1+sigma))*exp(-(1/sigmahat(i))*Mscript(i))*exp(((1-tau(i))*(1+sigmahat(i))/sigmahat(i))*( ((1-tau(i))*(1+sigmahat(i))/sigmahat(i)) - 1 )*(veps/2)); 
    Yd_2(i)= exp((1-tau(i))*pi0(i))*exp(-tau(i)*(1-tau(i))*vphi/2);
    Yd_3(i)=((1-delta)*exp(-tau(i)*(1-tau(i))*(valpha0/2)) / (1-delta*exp((-tau(i)*(1-tau(i))/2)*vomega)))*(theta/(-1+tau(i)+theta));
    
    Yd(i)=Yd_1(i)*Yd_2(i)*Yd_3(i);
    
    %Lambda from aggregate resource which is used in the case where G is
    %non-valued and aggregate G is fixed (case 3)
    lambda3(i) = (Y(i)-agg_G)/Yd(i);
    
    %SWF with chi=0 and an exogenous amount G wasted (case 3)
    W_TOT_XG3(i) =     log(lambda3(i))+((1-tau(i))/(1+sigma))*log(1-tau(i))+Mscript(i)-(1-tau(i))/(1+sigma)+...                   % this line has terms in log(c) (except unins. wage) and disutility from hours  
      (1-tau(i))*pi0(i)+(1-tau(i))/theta-(1-tau(i))*(valpha0/2)-(1-tau(i))*(delta/(1-delta))*(vomega/2)-(1-tau(i))*(vphi/2)-...   % this line has E[(1-tau)*(log(p)+alpha-phi)]
      (psi/(1+psi))*(1-tau(i))/theta+(((1-beta)*delta)/(1-beta*delta))*(psi/(1+psi))*(1/theta)*(1-taum1); ;                      % this line has education cost  

    %SWF with chi=0 and an exogenous amount G wasted for  right-wing planner (case 3)
    
    % SWF of the representative agent with chi=0 (and g=0)
    W_RA_NOG(i)    = ( 1/(1+sigma) )*( log(1-tau(i)) - (1-tau(i)) );
    
    % SWF of the representative agent with valued G
    % Note that log(1+g)+chi log g = chi*log(chi) - (1+chi)*log(1+chi)
    W_RA(i)   = W_RA_NOG(i) + chi*log(chi) - (1+chi)*log(1+chi) + chi*( log(1-tau(i))/(1+sigma) );
    
    %SWF piece associated to productivity gain from skill investment
    %Second line is for the case chi=0 (and g=0)
    W_EDU1(i)      = (1+chi)*( (psi/(1+psi))*(1/(theta-1))*log(1-tau(i)) + (1/((theta-1)*(1+psi)))*log( (1/(eta*theta^psi))*(theta/(theta-1))^(theta*(1+psi)) ) );
    W_EDU1_NOG(i) = ( (psi/(1+psi))*(1/(theta-1))*log(1-tau(i)) + (1/((theta-1)*(1+psi)))*log( (1/(eta*theta^psi))*(theta/(theta-1))^(theta*(1+psi)) ) );
    
    %SWF piece associated to between skill-group consumption inequality
    %(second line for right-wing planner)
    W_EDU2(i)      =  (1-tau(i))/theta  + log(1-((1-tau(i))/theta));

    W_EDU2_NUL(i)      =  log(1-((1-tau(i))/theta)) - (nul/(nul-1))*log(((1-nul)/nul)*(1/theta)*(1-tau(i))+1);
    W_EDU2_NUH(i)      =  log(1-((1-tau(i))/theta)) - (nuh/(nuh-1))*log(((1-nuh)/nuh)*(1/theta)*(1-tau(i))+1);
    W_EDU2_NUINF(i)    =  0;

    %SWF piece associated to average education cost (second line for right-wing planner)
    W_EDU3(i)      =     -(psi/(1+psi))*(1/theta)*(1-tau(i)) + (((1-beta)*delta)/(1-beta*delta))*(psi/(1+psi))*(1/theta)*(1-taum1);
        
    %SWF piece associated to preference heterogeneity (second line for right-wing planner)
    W_PREF(i)      = -(1-tau(i))^2*vphi/2;

    W_PREF_NUL(i)      = -(1-tau(i))^2*(1/nul)*vphi/2;
    W_PREF_NUH(i)      = -(1-tau(i))^2*(1/nuh)*vphi/2;
    W_PREF_NUINF(i)    = 0;
    
    %SWF piece associated with within skill-group consumption inequality
    W_OMEGA(i)      = log( ( 1-delta*exp( -tau(i)*(1-tau(i))*vomega/2 ) )/(1-delta) ) - (1-tau(i))*delta/(1-delta)*vomega/2;
%   W_OMEGA_APPX(i) = log( ( 1-delta*exp( -tau(i)*(1-tau(i))*vomega/2 ) )/(1-delta) ) - tau(i)*(1-tau(i))*delta/(1-delta)*vomega/2;
    W_ALPHA_APPROX(i) = -(1-tau(i))^2*valpha/2;
    
 
    W_OMEGA_NUCOMMON(i) = -(1-tau(i))^2*((beta*delta)/(1-beta*delta))*(vomega/2); 

    W_OMEGA_NUL(i) = -(1-tau(i))^2*(1/nul)*((1-beta)/(1-beta*delta))*(delta/(1-delta))*(vomega/2); 
    W_OMEGA_NUH(i) = -(1-tau(i))^2*(1/nuh)*((1-beta)/(1-beta*delta))*(delta/(1-delta))*(vomega/2); 
    W_OMEGA_NUINF(i) = 0; 
    
    %SWF piece associated to initial variance of uninsuranle shocks
    W_ALPHA0(i) = -(1-tau(i))^2*valpha0/2;

    %SWF piece associated to insurable shocks (second line for chi=0 (and g=0))
    W_EPS(i)  = (1+chi)*(1-tau(i))/(sigma+tau(i))*veps - (1+chi)*sigma*( (1-tau(i))/(sigma+tau(i)) )^2*veps/2 ;
    W_EPS_NOG(i)  = (1-tau(i))/(sigma+tau(i))*veps - sigma*( (1-tau(i))/(sigma+tau(i)) )^2*veps/2 ;
    
    %SWC piece associated to insurable shock for progressive consumption tax
    W_CONSTAX_EPS(i)  = ((1+chi)/sigma)*veps - (1+chi)*sigma*(1/sigma^2)*(veps/2) ;

    %Additional terms for SWF defined as ex-ante expected welfare of agent
    %with median and average values of (kappa, alpha0, phi) before you draw your epsilon
    
    W_ADJ(i) =           (1-tau(i))*(vphi/2) + (1-tau(i))*(delta/(1-delta))*(vomega/2) - (1-tau(i))/theta + (psi/(1+psi))*((1-tau(i))/theta - (1-taum1)/theta) - (beta*delta/(1-beta*delta))*(1-tau(i))*(vomega/2);  % Here we undo the expectations of the idiosyncratic terms in welfare
    
    W_TOT_MEAN_ADD(i) = -(1-tau(i))*(vphi/2) - (1-tau(i))*(delta/(1-delta))*(vomega/2) + (1-tau(i))/theta - (psi/(1+psi))*((1-tau(i))/theta - (1-taum1)/theta);
    W_TOT_MED_ADD(i) =  -(1-tau(i))*(vphi/2) - (1-tau(i))*(delta/(1-delta))*(vomega/2) + (1-tau(i))*log(2)/theta - (psi/(1+psi))*((1-tau(i))*log(2)/theta - (1-taum1)*log(2)/theta);
    
       
end;

%-----------------------------------------
%Aggregation of components of SWF
%and calculation of optimal taus
%----------------------------------------

%Piece together edu terms (second line is right-wing planner)
W_EDU   = W_EDU1 + W_EDU2 + W_EDU3;

W_EDU_NUL = W_EDU1 + W_EDU2_NUL + W_EDU3;
W_EDU_NUH = W_EDU1 + W_EDU2_NUH + W_EDU3;
W_EDU_NUINF = W_EDU1 + W_EDU2_NUINF + W_EDU3;

%Piece together edu terms in case chi=0 (and g=0) (second line is right-wing planner)
W_EDU_NOG   = W_EDU1_NOG + W_EDU2 + W_EDU3;

%Piece together initial dispersion and lifetime shocks
W_ALPHA = W_ALPHA0+W_OMEGA;

%Piece together all components (second line is approximation of lifetime cumulation of Vomega)
W_TOT = W_RA + W_EDU + W_PREF + W_ALPHA + W_EPS;
W_TOT_APPROX = W_RA+W_EDU+W_PREF+W_ALPHA_APPROX+W_EPS;

W_TOT_NUL = W_RA + W_EDU_NUL + W_PREF_NUL + W_OMEGA_NUCOMMON + W_OMEGA_NUL + W_EPS;
W_TOT_NUH = W_RA + W_EDU_NUH + W_PREF_NUH + W_OMEGA_NUCOMMON + W_OMEGA_NUH + W_EPS;
W_TOT_NUINF = W_RA + W_EDU_NUINF + W_PREF_NUINF + W_OMEGA_NUCOMMON + W_OMEGA_NUINF + W_EPS;

%SWF in the case where skills are exogenous
W_TOT_NOSKILL2 = W_TOT - W_EDU + W_EDU2; 

%SWF for progressive consumption tax
W_CONSTAX_TOT = W_RA+W_EDU+W_PREF+W_ALPHA+W_CONSTAX_EPS;

%Ex-ante SWF for agent with median and mean (kappa, phi,alpha)
W_TOT_MED = W_TOT + W_ADJ + W_TOT_MED_ADD;
W_TOT_MEAN = W_TOT + W_ADJ + W_TOT_MEAN_ADD;

%Various SWF for case when chi=0 (and g=0)
W_TOT_NOG       = W_RA_NOG + W_EDU_NOG + W_PREF + W_ALPHA + W_EPS_NOG;
%SWF for case where chi=0 and an exogenous fraction of output g is wasted 
W_TOT_XG2     = W_TOT + log(1-g) - chi*W_EPS/(1+chi) - chi*W_EDU1/(1+chi) + (W_RA_NOG - W_RA);


%-----------------------------
%Compute implied tax system
%-----------------------------

% Representative Agent
disp('RA taustar');
[c,i]=max((W_RA)');
optau_RA = tau(i)

% Piece associated to skills investment
disp('Skill Investment taustar');
[c,i]=max((W_EDU)');
optau_EDU = tau(i)

% Optimal tau
disp('Exact taustar');
[c,i]=max(W_TOT');
optau_EXACT = tau(i)
disp('Welfare gain from tax reform');
100*(exp(max(W_TOT')-W_TOT(1181))-1)

disp('Exact taustar RA+EDU');
[c,i]=max((W_RA+W_EDU)');
optau2_EXACT = tau(i)

disp('Exact taustar RA+EDU+PREF');
[c,i]=max((W_RA+W_EDU+W_PREF)');
optau3_EXACT = tau(i)

disp('Exact taustar RA+EDU+PREF+ALPHA');
[c,i]=max((W_RA+W_EDU+W_PREF+W_ALPHA)');
optau4_EXACT = tau(i)

disp('Approximated SWF taustar');
[c,i]=max(W_TOT_APPROX');
optau_APPROX = tau(i)

disp('NUL taustar');
[c,i]=max(W_TOT_NUL');
optau_NUL = tau(i)

disp('NUH taustar');
[c,i]=max(W_TOT_NUH');
optau_NUH = tau(i)

disp('NUINF taustar');
[c,i]=max(W_TOT_NUINF');
optau_NUINF = tau(i)


% Consumption tax
disp('Consumption tax taustar');
[c,i]=max(W_CONSTAX_TOT');
optau_CONS = tau(i)
disp('Welfare gain from tax reform');
100*(exp(max(W_CONSTAX_TOT')-W_TOT(1181))-1)


% Ex-ante welfare of median agent 
% This is not the true median voter model. Here we take the median of
% alpha, phi and kappa separately, rather than the median of the
% combination
disp('Ex-ante welfare of median agent taustar');
[c,i]=max(W_TOT_MED');
optau_MED = tau(i)

% Ex-ante welfare of mean agent 
disp('Ex-ante welfare of mean agent taustar');
[c,i]=max(W_TOT_MEAN');
optau_MEAN = tau(i)

% Exogenous skills taustar 
disp('Exogenous skills taustar');
[c,i]=max(W_TOT_NOSKILL2');
optau_NOSKILL2 = tau(i)

% chi=0 taustar 
disp('Chi=0 taustar');
[c,i]=max(W_TOT_NOG');
optau_NOG = tau(i)
disp('Welfare gain from tax reform');
100*(exp(max(W_TOT_NOG')-W_TOT_NOG(1181))-1)

% chi=0 and fraction of output g wasted
disp('Chi=0 and fraction of output g wasted -- taustar');
[c,i]=max(W_TOT_XG2(1:1500)');
optau_XG2 = tau(i)
disp('Welfare gain from tax reform');
100*(exp(max(W_TOT_XG2(1:1500)')-W_TOT_XG2(1181))-1)

% chi=0 and amount G wasted
disp('Chi=0 and amount G wasted -- taustar');
[c,i]=max(W_TOT_XG3(1:1500)');
optau_XG3 = tau(i)
disp('Welfare gain from tax reform');
100*(exp(max(W_TOT_XG3(1:1500)')-W_TOT_XG3(1181))-1)

%------------------------------------------
%HOW TAU CHANGES WITH THETA AND PSI
%------------------------------------------

for j=1:200;
    
    psil = 0.1;
    psih = 5;
    psi0 = 0;
    
    thetavec(j)=1.001+0.05*(j-1);
     
    % General case with sigma finite and chi>0 (but other sources of
    % heterogeneity besides kappa shut down) See
    % hsv_tex_revision_clean_v2_temp.tex
    
    psiv(1)=0;
    psiv(2)=psi;
    psiv(3)=5;
    
    for jj=1:3;
        
    psit=psiv(jj);    
    
    aq(j)=(1/((1+psit)*thetavec(j)))-(1/(1+sigma));  
    bq(j)=((1+chi+thetavec(j))/(1+sigma))+(psit/(1+psit))*((chi+thetavec(j))/(thetavec(j)-1));
    cq(j)=-((1/(1+sigma))+(psit/((1+psit)*(thetavec(j)-1))))*(1+chi)*thetavec(j);

    optau_EDUVEC1(j,jj) = 1-(-bq(j)+(bq(j)^2-4*aq(j)*cq(j))^0.5)/(2*aq(j));

    bq_nochi(j)=((1+thetavec(j))/(1+sigma))+(psit/(1+psit))*((thetavec(j))/(thetavec(j)-1));
    cq_nochi(j)=-((1/(1+sigma))+(psit/((1+psit)*(thetavec(j)-1))))*thetavec(j);
      
    optau_EDUVEC1_nochi(j,jj) = 1-(-bq_nochi(j)+(bq_nochi(j)^2-4*aq(j)*cq_nochi(j))^0.5)/(2*aq(j));
    
    end;
    
    % Case for paper plot with sigma=infinity and chi = 0 (to shut down other channels) 
    optau_EDUVEC2(j) = 1 - (-thetavec(j)+(thetavec(j)^2+4*(thetavec(j)-1)/psi)^0.5)/(2*(thetavec(j)-1)/(psi*thetavec(j)));

    optau_EDUVEC3(j) = 1 - (-thetavec(j)+(thetavec(j)^2+4*(thetavec(j)-1)/psil)^0.5)/(2*(thetavec(j)-1)/(psil*thetavec(j)));
    optau_EDUVEC4(j) = 1 - (-thetavec(j)+(thetavec(j)^2+4*(thetavec(j)-1)/psih)^0.5)/(2*(thetavec(j)-1)/(psih*thetavec(j)));    
    optau_EDUVEC5(j) = 1 - (-thetavec(j)+(thetavec(j)^2+4*(thetavec(j)-1)/psi0)^0.5)/(2*(thetavec(j)-1)/(psi0*thetavec(j))); 
    % Now compute the optimal tau with exogenous skills

    optau_EDUVEC_XS(j) = 1 - (1/(thetavec(j)-1-sigma))*(thetavec(j)/2)*(1+thetavec(j)-(5-2*thetavec(j)+thetavec(j)^2+4*sigma)^0.5);
end;


figure(20);
plot(thetavec,optau_EDUVEC4,thetavec,optau_EDUVEC3,thetavec,optau_EDUVEC2,'Linewidth',2.5);
grid on; legend('\psi=0.1','\psi=0.65 (baseline)','\psi=5');xlim([1 10]);
xlabel('\theta', 'fontsize', 14);ylabel('\tau^{\ast}', 'fontsize', 14);
    


% THIS IS FIGURE II IN THE PAPER
figure(21);
plot(thetavec,optau_EDUVEC1_nochi(:,1),thetavec,optau_EDUVEC1_nochi(:,2),thetavec,optau_EDUVEC1_nochi(:,3),'Linewidth',2.5);
grid on; legend('\psi=0.0 (exogenous skills)','\psi=0.65 (baseline)','\psi=2.0');xlim([1 10]);
xlabel('\theta', 'fontsize', 14);ylabel('\tau^{\ast}', 'fontsize', 14);

%-----------------------------
%Plots of SWF
%-----------------------------

%Max of welfare function, needed for normalization
welf_max          = max(W_TOT);
welf_max_CONSTAX  = max(W_CONSTAX_TOT);
welf_max_NOG      = max(W_TOT_NOG);
welf_max_XG2      = max(W_TOT_XG2(1:1500));
welf_max_XG3      = max(W_TOT_XG3(1:1500));



%Plot RA component of Social Welfare Function
figure(1);
plot(tau,W_RA_NOG,tau,W_RA,'Linewidth',2.5);
grid on;
xlim([-0.9 0.9]);
legend('RA G non-valued','RA G valued');
xlabel('Progressivity rate (\tau)', 'fontsize', 18);
ylabel('Social Welfare', 'fontsize', 18);
title('Representative Agent Economies', 'fontsize', 18);

%Plot skill component of Social Welfare Function
figure(2);
plot(tau,100*(exp(W_EDU1+W_EDU3)-1),tau,100*(exp(W_EDU2)-1),tau,100*(exp(W_EDU)-1),'Linewidth',2.5);
legend('(A) Prod Gain - Edu Cost','(B) Btw-Skill Cons Ineq','(A)+(B) Net Effect');
grid on;xlim([-0.9 0.9]);
xlabel('Progressivity rate (\tau)', 'fontsize', 14); ylabel('welfare change rel. to optimum (% of cons.)', 'fontsize', 14);
title('Skill Accumulation', 'fontsize', 18);

%Plot skill component of Social Welfare Function even more decomposed
figure(3);
plot(tau,100*(exp(W_EDU1)-1),tau,100*(exp(W_EDU2)-1),tau,100*(exp(W_EDU3)-1),tau,100*(exp(W_EDU)-1),'Linewidth',2.5);
legend('(A) Prod Gain ','(B) Btw-Skill Cons Ineq','(C) Edu Cost' ,'(A)+(B)+(C) Net Effect');
grid on;xlim([-0.9 0.9]);
xlabel('Progressivity rate (\tau)', 'fontsize', 14); ylabel('welfare change rel. to optimum (% of cons.)', 'fontsize', 14);
title('Skill Accumulation', 'fontsize', 18);


%Plot Pref. het. and uninsurable risk component of Social Welfare Function
figure(4);
plot(tau,100*(exp(W_PREF)-1),tau,100*(exp(W_ALPHA)-1),tau,100*(exp(W_EPS)-1),'Linewidth',2.5);
legend('Pref Het','Unins. Shocks','Insurable Shocks');
grid on;xlim([-0.9 0.9]);
xlabel('Progressivity rate (\tau)', 'fontsize', 18);ylabel('welfare contrib. (% of cons.)', 'fontsize', 18);
title('Pref heterogeneity, Efficiency Shocks', 'fontsize', 18);

%Plot components of Social Welfare Function adding sequentially
figure(5);
plot(tau,100*(exp(W_RA-welf_max)-1),tau,100*(exp(W_RA+W_EDU-welf_max)-1),tau,100*(exp(W_RA+W_EDU+W_PREF-welf_max)-1),...
     tau,100*(exp(W_RA+W_EDU+W_PREF+W_ALPHA-welf_max)-1),tau,100*(exp(W_TOT-welf_max)-1),'Linewidth',2.5);
legend('Rep. Agent','+Skill Het.','+Pref Het','+Unins. Shocks','+Ins. Shocks = Baseline');
grid on;xlim([-.5 .5]); ylim([-20,10]);
xlabel('Progressivity rate (\tau)', 'fontsize', 18);ylabel('welf change rel. to baseline optimum (% of cons.)', 'fontsize', 18);
title('Social Welfare Function', 'fontsize', 18);

%Plot final Social Welfare Function
figure(6);
plot(tau,100*(exp(W_TOT-welf_max)-1),'Linewidth',2.5);
legend('Baseline Model');
grid on;xlim([-.5 .5]); ylim([-20,10]);
xlabel('Progressivity rate (\tau)', 'fontsize', 18);ylabel('welf change rel. to optimum (% of cons.)', 'fontsize', 18);
title('Social Welfare Function', 'fontsize', 18);

%Plot final Social Welfare Function corresponding to the three assumptions on G
figure(7);
plot(tau,100*(exp(W_TOT-welf_max)-1),tau,100*(exp(W_TOT_XG2-welf_max)-1),tau,100*(exp(W_TOT_XG3-welf_max)-1),'Linewidth',2.5);
legend('Baseline Model','Exog G/Y','Exog G');
grid on;xlim([-.5 .5]); ylim([-20,120]);
xlabel('Progressivity rate (\tau)', 'fontsize', 18);ylabel('welf change rel. to optimum (% of cons.)', 'fontsize', 18);
title('Social Welfare Function', 'fontsize', 18);

%Progressive consumption tax
figure(8);
plot(tau,100*(exp(W_TOT-welf_max)-1),tau,100*(exp(W_CONSTAX_TOT-welf_max)-1),'Linewidth',2.5);
legend('Baseline Model','Consumption Tax');
grid on;xlim([-.5 .5]); ylim([-10,10]);
xlabel('Progressivity rate (\tau)', 'fontsize', 18);ylabel('welf change rel. to optimum (% of cons.)', 'fontsize', 18);
title('Social Welfare Function', 'fontsize', 18);

%-----------------------------
%Political Economy
%-----------------------------
if Polecon == 'Y'
   
   S = 100000;
   alphavar = -valpha/2 + sqrt(valpha)*randn(S,1);
   phivar   = vphi/2 + sqrt(vphi)*randn(S,1);
   kappavar = exprnd(1/eta,S,1);
   zvar = alphavar - phivar + (eta/theta)*kappavar - (psi/(1+psi))*(eta/theta)*kappavar;   
%  zvar = alphavar - phivar;         
   zmed = median(zvar);
         
%Decomposition of SWF
for i=1:1999;
    
    tau(i)=-1+0.001*i;

    W_ADJ(i) =  (1-tau(i))*(vphi/2) + (1-tau(i))*(delta/(1-delta))*(vomega/2) - (1-tau(i))/theta + (psi/(1+psi))*((1-tau(i))/theta - (1-taum1)/theta) - (beta*delta/(1-beta*delta))*(1-tau(i))*(vomega/2);  % Here we undo the expectations of the idiosyncratic terms in welfare      
    W_TOT_VOTER(i)  = W_TOT(i) + W_ADJ(i) + (1-tau(i))*zmed ;
 
end

welf_max_VOTER   = max(W_TOT_VOTER);

% Political economy tau
disp('Political economy taustar');
[c,i]=max(W_TOT_VOTER');
optau_VOTER = tau(i)



%Plot final Social Welfare Function
figure(9);
plot(tau,100*(exp(W_TOT-welf_max)-1),tau,100*(exp(W_TOT_VOTER-welf_max_VOTER)-1),'Linewidth',2.5);
legend('Baseline Model','Median Voter');
grid on;xlim([-.5 .5]); ylim([-20,10]);
xlabel('Progressivity rate (\tau)', 'fontsize', 18);ylabel('welf change rel. to optimum (% of cons.)', 'fontsize', 18);
title('Social Welfare Function', 'fontsize', 18);
end

%Sufficient condition for progressivity
chithresh = ( psi/((1+psi)*(theta-1)) + 1/(sigma+1) )^(-1)*(1/(1+psi))*( (1/(theta-1))-(1/theta)+vphi+valpha );
%-----------------------------
%Plots of tax_systems
%-----------------------------

if Polecon == 'Y'
    tauvec=[UStau optau_EXACT optau_NOG optau_VOTER];
else
    tauvec=[UStau optau_EXACT optau_NOG optau_NUL optau_NUINF];
end

for j=1:size(tauvec,2);
    
    %Select case \chi = 0
    if j == 3
        chi = 0;
    else
        chi = g/(1-g);
    end
    
    tau_temp = tauvec(j);
    
    sigmahat_tau_temp = (sigma+tau_temp)/(1-tau_temp);
    
    %Lambda
    loglambda = log(1-chi/(1+chi)) + tau_temp/(1+sigma)*log(1-tau_temp) +...
                tau_temp*(1-tau_temp)/(sigma+tau_temp)*(2+sigma+(1+sigma)/(sigma+tau_temp))*veps/2 +...
                tau_temp*(1-tau_temp)*vphi/2 +...
                tau_temp*(1-tau_temp)*valpha0/2+log( ( 1-delta*exp( -tau_temp*(1-tau_temp)*vomega/2 ) )/(1-delta) )+...
                (psi/(1+psi))*tau_temp/(theta-1)*log((1-tau_temp)/theta)-(1/(1+psi))*tau_temp/(theta-1)*log(eta)+...
                (theta-1+tau_temp)/(theta-1)*log(theta/(theta-1))+log((theta-1+tau_temp)/theta);
                
    lambda=exp(loglambda);
    
    lambdav(j)=lambda;

    logybar = (theta/(theta-1))*log(theta/(theta-1))+((psi/((theta-1)*(1+psi)))+(1/(1+sigma)))*log(1-tau_temp)-(psi/((theta-1)*(1+psi)))*log(theta)-(1/((theta-1)*(1+psi)))*log(eta)+...
    ((1/sigmahat_tau_temp^2)*(tau_temp+sigmahat_tau_temp+sigmahat_tau_temp*tau_temp)*(veps/2));

    %Computes ybar (average Y) for each different value of tau
    %ybar = exp(logybar);
    
    %Compute ybar for US tau
    ybar = Y_UStau;
    
    for i=1:50000;
        yv(i,j)=0.1*ybar*0.001*i;
        mtax(i,j) = 1 - lambda*(1-tau_temp)*yv(i)^(-tau_temp);
        atax(i,j) = 1 - lambda*yv(i)^(-tau_temp);

	    yv(i,j) = yv(i,j)/ybar;
    end;
end;


figure(10);
plot(yv(:,1),mtax(:,1),yv(:,1),atax(:,1),'Linewidth',2.5);
legend('US Marginal (\tau^{US}=0.181)','US Average (\tau^{US}=0.181)');
grid on;xlim([0.1 max(yv(:,1))]);
xlabel('Income (1 = Average Income)', 'fontsize', 18);ylabel('Marginal and average tax rates', 'fontsize', 18);

% THIS IS FIGURE D2 IN THE APPENDIX
%Plot marginal tax rate--baseline
figure(11);
subplot(2,2,1);
plot(yv(:,1),mtax(:,1),'-b',yv(:,2),mtax(:,2),'--r',yv(:,3),mtax(:,3),':k','Linewidth',2.5);
legend('Actual US: \tau^{US} = 0.181', 'Utilitarian: \tau^{\ast} = 0.084','\chi = 0: \tau^{\ast} = 0.200' ,'Location','SouthEast', 'fontsize', 18,'Interpreter','latex');
grid on;xlim([0.0 max(yv(:,1))]);ylim([-0.2 0.501]);
xlabel('Income (1 = average income at \tau^{US})', 'fontsize', 14);ylabel('Marginal tax rate', 'fontsize', 14);title('Utilitarian Planner', 'fontsize', 14);
%Plot marginal tax rate-- Inequality Aversion
subplot(2,2,2);
plot(yv(:,2),mtax(:,2),'--r',yv(:,4),mtax(:,4),':k',yv(:,5),mtax(:,5),'-.g','Linewidth',2.5);
legend('\nu=1 (utilit.): \tau^{\ast} = 0.084', '\nu=2: \tau^{\ast} = 0.190','\nu=0: \tau^{\ast} = -0.159' ,'Location','SouthEast', 'fontsize', 18,'Interpreter','latex');
grid on;xlim([0.0 max(yv(:,1))]);ylim([-0.2 0.501]);
xlabel('Income (1 = average income at \tau^{US})', 'fontsize', 14);ylabel('Marginal tax rate', 'fontsize', 14);title('Inequality Averse Planner', 'fontsize', 14);
%Plot average tax rate--baseline
subplot(2,2,3);
plot(yv(:,1),atax(:,1),'-b',yv(:,2),atax(:,2),'--r',yv(:,3),atax(:,3),':k','Linewidth',2.5);
legend('Actual US: \tau^{US} = 0.181', 'Utilitarian: \tau^{\ast} = 0.084', '\chi = 0: \tau^{\ast} = 0.20' ,'Location','SouthEast', 'fontsize', 18,'Interpreter','latex');
grid on;xlim([0.0 max(yv(:,1))]);ylim([-0.401 0.4501]);
xlabel('Income (1 = average income at \tau^{US})', 'fontsize', 14);ylabel('Average tax rate', 'fontsize', 14);title('Utilitarian Planner', 'fontsize', 14);
%Plot average tax rate-- Inequality Aversion
subplot(2,2,4);
plot(yv(:,2),atax(:,2),'--r',yv(:,4),atax(:,4),':k',yv(:,5),atax(:,5),'-.g','Linewidth',2.5);
legend('\nu=1 (utilit.): \tau^{\ast} = 0.084', '\nu=2: \tau^{\ast} = 0.190', '\nu=0: \tau^{\ast} = -0.159' ,'Location','SouthEast', 'fontsize', 18,'Interpreter','latex');
grid on;xlim([0.0 max(yv(:,1))]);ylim([-0.401 0.4501]);
xlabel('Income (1 = average income at \tau^{US})', 'fontsize', 14);ylabel('Average tax rate', 'fontsize', 14);title('Inequality Averse Planner', 'fontsize', 14);

%return

if Polecon == 'Y'
    
%Compare median voter and utilitarian
figure(12);
plot(yv(:,2),mtax(:,2),yv(:,4),mtax(:,4),'Linewidth',2.5);
legend( 'Utilitarian \tau^{\ast}=0.084', 'Median voter \tau^{\ast}=0.108');
grid on;xlim([0.1 max(yv(:,1))]);
xlabel('Income (1 = Average Income)', 'fontsize', 18);ylabel('Marginal tax rate', 'fontsize', 18);

end


if Irreversible == 'Y'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%THIS PART OF CODE COMPUTES THE MODEL WITH IRREVERSIBLE INVESTMENT
%SEE THE FILE TRANSITION.TEX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ngpt = 599; % number of tau values to consider 
ngpg = 20;  % number of gamma values to consider
simlength = 10000;       % length of simulation for computing welfare 

zeta = (-1/((1+psi)*(theta-1)))*(log(eta)+psi*log(theta))+(1/(theta-1))*log(theta/(theta-1));

% Construct grid for different values for gamma 
% Note that we need (1-beta*delta/gamma) to be larger than one, which
% implies gamma > beta*delta = 0.9302

gama_min = 0.94;
gama_max = 0.999;


gamav(1)=gama_min;

for j=2:ngpg;
    
    gamav(j)=gamav(j-1)+(gama_max-gama_min)/(ngpg-1);
    
end;


Yt = zeros(simlength,ngpt);
Ydt = zeros(simlength,ngpt);
lambdat = zeros(simlength,ngpt);
DPV = zeros(ngpt,ngpg);

taut=zeros(ngpt,1);
Mscript=zeros(ngpt,1);
X=zeros(ngpt,1);
Z=zeros(ngpt,1);
omterm=zeros(ngpt,1);

% First compute sequence for lambda for each possible tau

for i=1:ngpt;

  
    taut(i)=0.0+0.001*i;
    
    Mscript(i) = ((1-taut(i))/(sigma+taut(i)))*(1-2*taut(i)-sigma*taut(i))*veps/2;
    
    X(i) = exp(Mscript(i))*(1-taut(i))^((1-taut(i))/(1+sigma))*exp((1-taut(i))*zeta)*exp(-taut(i)*(1-taut(i))*vphi/2)*(theta/(theta-(1-taut(i))));
    
    Z(i) = exp( ((1-taut(i))/(sigma+taut(i))^2)*(sigma+2*taut(i)+sigma*taut(i))*veps/2 )*exp(zeta)*(theta/(theta-1))*(1-taut(i))^(1/(1+sigma));
    
    omterm(i) = exp(-(1-taut(i))*taut(i)*vomega/2);
         
 
 %  Piece of code to construct welfare in baseline model with reversible
 %  investment model (only difference relative to above is that part of the
 %  omega term now varies with gamma
 
    W_RA_NOG(i) = ( 1/(1+sigma) )*( log(1-taut(i)) - (1-taut(i)) );
    W_RA(i)     = W_RA_NOG(i) + chi*log(chi) - (1+chi)*log(1+chi) + chi*( log(1-taut(i))/(1+sigma) );       
 
    
    W_EDU1(i)     = (1+chi)*( (psi/(1+psi))*(1/(theta-1))*log(1-taut(i)) + (1/((theta-1)*(1+psi)))*log( (1/(eta*theta^psi))*(theta/(theta-1))^(theta*(1+psi)) ) );
    W_EDU2(i)   =  (1-taut(i))/theta  + log(1-((1-taut(i))/theta));
    W_EDU3(i)   =  -(psi/(1+psi))*(1/theta)*(1-taut(i)) + (((1-beta)*delta)/(1-beta*delta))*(psi/(1+psi))*(1/theta)*(1-taum1);

    W_EDU(i)    = W_EDU1(i) + W_EDU2(i) + W_EDU3(i);
    
    W_PREF(i)   = -(1-taut(i))^2*vphi/2;
   
    W_ALPHA0(i) = -(1-taut(i))^2*valpha0/2;
   
    W_EPS(i)    = (1+chi)*(1-taut(i))/(sigma+taut(i))*veps - (1+chi)*sigma*( (1-taut(i))/(sigma+taut(i)) )^2*veps/2 ;

    W_OMEGA1(i) = log( ( 1-delta*exp( -taut(i)*(1-taut(i))*vomega/2 ) )/(1-delta) );
    
    W_OMEGA2(i) = -(1-taut(i))*vomega/2;
        
    W_TOT_XOM(i) = W_RA(i) + W_EDU(i) + W_PREF(i) + W_ALPHA0(i) + W_EPS(i) + W_OMEGA1(i);

  
    
    for t=1:simlength;                   
        
                
        Yt(t,i) = (delta^t*(1-taum1)^(psi/((1+psi)*(theta-1)))+(1-delta^t)*(1-taut(i))^(psi/((1+psi)*(theta-1))))*Z(i);
      
        Ydt(t,i) = ( (1-taum1)^(psi*(1-taut(i))/((1+psi)*(theta-1)))*delta^t*omterm(i)^t/(1-delta*omterm(i)) + (1-taut(i))^(psi*(1-taut(i))/((1+psi)*(theta-1)))*(1-delta^t*omterm(i)^t)/(1-delta*omterm(i)) )*X(i)*(1-delta);
        
        lambdat(t,i) = (1/(1+chi))*Yt(t,i)/Ydt(t,i);

    end;    
    
    % Now for each tau, construct social welfare for alternative values for
    % gamma
  
    
    for j=1:ngpg;
        
        gamma = gamav(j);

        W_TOT_FIN(i,j) = W_TOT_XOM(i) + (beta*delta/(gamma-beta*delta))*W_OMEGA2(i);
        
        DPV(i,j) = 0;
 
        U = zeros(simlength);

        for t=1:simlength;
            
            U(t) = log(lambdat(t,i)) + ((1-taut(i))/(1+sigma))*log(1-taut(i)) - (1-taut(i))*vphi/2 + (1-taut(i))*zeta + (1-taut(i))/theta + Mscript(i) + ...
            (1-taut(i))*(-vomega/2)*( (beta*delta/gamma) /(1-(beta*delta/gamma)) )   + ...
            (1- (beta*delta/gamma)^t )*(1-taut(i))*psi*log(1-taut(i))/((1+psi)*(theta-1)) + ((beta*delta/gamma)^t)*(1-taut(i))*psi*log(1-taum1)/((1+psi)*(theta-1))  + ...
            chi*log(chi/(1+chi)) + chi*log(Yt(t,i)) - (1-taut(i))/(1+sigma) - ...
            ((1-(beta*delta/gamma))/(1-beta*delta))*(1-taut(i))*psi/((1+psi)*theta) ;
        
            DPV(i,j) = DPV(i,j) + gamma^(t-1)*U(t);
        
        end;
        
        DPV(i,j) = DPV(i,j)*(1-gamma)/(1-gamma^simlength);
        
    end;
    
    i
    
end;

for j=1:ngpg;

    [c,i]=max(DPV(:,j));
    optau_trans(j) = taut(i);

    [c,i]=max(W_TOT_FIN(:,j));
    optau_base(j) = taut(i);

    
    j
end;

% THIS IS FIGURE D3 IN THE APPENDIX
%Plot final Social Welfare Function

figure(13);
hold off;
plot(gamav,optau_trans,'Linewidth',2.5,'LineStyle','-');
hold on;
plot(gamav,optau_base,'Linewidth',2.5, 'LineStyle','--');
legend('Fixed Investment','Flexible Investment','Interpreter','latex');
hold on;
x=[beta,beta];
y=[0,0.3];
text(0.96, 0.01, '\beta', 'Color', 'b','fontsize', 18);
plot(x,y);
grid on; xlim([0.94 1]); ylim([0,0.3]);
xlabel('Generational Weight \gamma', 'fontsize', 18);
ylabel('\tau^{\ast}', 'fontsize', 18);
%title('Progressivity versus Generational Weight', 'fontsize', 16);

qsp = spline(gamav,optau_trans);
ts=ppval(qsp,beta)



end;
