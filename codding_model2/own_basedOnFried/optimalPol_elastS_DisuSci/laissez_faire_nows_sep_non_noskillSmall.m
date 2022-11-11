function f=laissez_faire_nows_sep_non_noskillSmall(x, params, list, pol, laggs, indic)
% Model
% equilibrium for one period!
% takes policy as given
list.choice=list.choice_small;

% starts to solve the model in 2020-2024;
% i.e. Aj_laggs  refer to 2015-2019 period
% A_lag is except for the initial condition is 
%- read in policy and parameters
read_in_params;
read_in_pol;
if indic.util==1 || thetaa~=1
    error('wrong utility specification')
end

% choice variables
%- transform variables directly instead of in code

 F      = exp(x(list.choice=='F'));
 pf     = exp(x(list.choice=='pf'));
 gammalh = x(list.choice=='gammalh')^2;

 if indic.xgrowth==0
     Af     = exp(x(list.choice=='Af'));
     Ag     = exp(x(list.choice=='Ag'));
     minn   =indic.minn;
     sff     = upbarH*minn/(1+exp(x(list.choice=='sff')));%exp(x(list.choice=='S')); % total labour supply
     sg      = upbarH*minn/(1+exp(x(list.choice=='sg')));%exp(x(list.choice=='S')); % total labour supply
     wsf     = exp(x(list.choice=='wsf'));
     wsg     = exp(x(list.choice=='wsg'));
      gammasg = x(list.choice=='gammasg')^2;
     gammasf = x(list.choice=='gammasf')^2;
 end
 
%% - read in auxiliary equations
%- initial condition
Af_lag=laggs(list.init=='Af0');
Ag_lag=laggs(list.init=='Ag0');

if indic.xgrowth==1
    Ag=(1+vg)*Ag_lag;
    Af=(1+vf)*Af_lag;
end
A_lag   = (rhof*Af_lag+rhog*Ag_lag)/(rhof+rhog);
h       = ((1-taul)/chii)^(1/(1+sigmaa)); 
w       = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; % labour demand fossil

lambdaa =(w.*h+tauf.*pf.*F)./(w.*h).^(1-taul); 
SGov    = (w.*h-lambdaa.*(w.*h).^(1-taul))...
         +tauf.*pf.*F; 
C       = lambdaa.*(w.*h).^(1-taul)+SGov;

Lf      = F./(((1-tauf).*alphaf.*pf).^(alphaf./(1-alphaf)).*Af); 
pg      = (w/((1-alphag)*alphag^(alphag/(1-alphag))*Ag))^(1-alphag);
G       = F.*(pf./pg).^eppse; % green demand
Lg      = G./(Ag*(pg.*alphag).^(alphag./(1-alphag)));

muu   = C.^(-thetaa); % same equation in case thetaa == 1
E     = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));


%% model equations
q=0;

% market clearing labour
q=q+1;
f(q) = Lg +Lf - h;

q=q+1;
f(q) = 1 - (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); %definition

%- definitions prices
%15
q=q+1;
f(q)= gammalh.*(upbarH-h);

if indic.xgrowth==0
    %6- demand green scientists
    q=q+1;
    f(q)= wsf - (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*(1-tauf).*F.*(1-alphaf).*Af_lag)./(rhof^etaa.*Af); 

    %7
    q=q+1;
    f(q)= wsg - (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G.*(1-alphag).*Ag_lag)./(rhog^etaa.*(1-taus)*Ag);

    %8
    q=q+1;
    f(q) = Af-Af_lag.*(1+gammaa.*(sff./rhof).^etaa*(A_lag./Af_lag).^phii);
    %9
    q=q+1;
    f(q) = Ag-Ag_lag.*(1+gammaa.*(sg./rhog).^etaa*(A_lag./Ag_lag).^phii);

    % optimality scientists
    q=q+1;
    f(q)= (chiis)*sff^sigmaas-(wsf-gammasf); % scientist hours supply
    q=q+1;
    f(q)= (chiis)*sg^sigmaas-((wsg-gammasg));
    % Kuhn tucker scientists
    q=q+1;
    f(q)= gammasf.*(upbarH-sff);
    q=q+1;
    f(q)= gammasg.*(upbarH-sg);
end
% fprintf('number equations: %d; number variables %d', q, length(list.choice));
end