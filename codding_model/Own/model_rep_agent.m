%% Representative agent model

% competitive equilibrium and objective function gov. 
% in symbolic variables

% with two skill types if zetaa!=1

% firm specific variables already substituted by optimality conditions
% machines, price of machines, and machine-specific productivity drop
%% Variables and parameters
% today and tomorrow

syms c cp ...           % consumption 
    muu muup ...          % lagrange muultiplier budget
    lambdaa lambdaap ...% shifter government revenues
    tauul tauulp ...    % progressivity parameter
    wl wlp ...          % wage rate low skill
    wh whp ...          % wage rate high skill
    H Hp ...            % disutility weighted hours supplied
    hl hlp ...          % low skill hours supplied
    hh hhp ...          % high skill hours supplied
    ...
    ... % production
    pc pcp ...          % price cleaner good
    pd pdp ...          % price dirty good
    yc ycp ...          % clean sector output
    yd ydp ...          % dirty sector output
    p pp ...            % aggregate price level
    Y Yp ...            % aggregate output
    Lc Lcp ...          % labour input clean sector
    Ld Ldp ...          % labour input dirty sector
    pcL pcLp ...        % price of clean labour input good
    pdL pdLp ...        % price of dirty labour input good
    lhc lhcp ...        % high skill input clean sector
    llc llcp ...        % low skill input clean sector
    lhd lhdp ...        % high skill input dirty sector
    lld lldp ...        % low skill input dirty sector
    xd xdp ...          % machines dirty sector
    xc xcp ...          % machines clean sector
    ...
    ... % productivity
    Ac Acp ...          % productivity clean sector
    Ad Adp ...          % productivity dirty sector
    vdd vddp ...          % productivity growth dirty sector - gov choice
    vc vcp ...          % productivity growth clean sector - gov choice
    real
%    ...
%     ... % inflation
%     pic picp ...        % inflation in the clean sector = pcp/pc
%     pid pidp ...        % inflation in the dirty sector = pdp/pd
    

% parameters
syms sigmaa...      % 1/sigmaa = Frisch elasticity of labour
     zetaa ...      % disutility increase from high skill labour
     eppsilon ...   % elasticity of substitution clean and dirty product
     alphaa ...     % 'capital' share production
     psii ...       % cost machine production 
     thetac ...     % high skill weight clean sector
     thetad ...     % high skill weight dirty sector
     Uppsilon ...   % aggregate growth
     betaa ...      % time preference household
     gammaa ...     % coefficient of relative risk aversion
     etaa ...       % disutility from labour
     E ...          % vector of net emission targets
     deltaa ...     % regeneration rate nature
     kappaa ...     % emission share of dirty output
     G ...          % exogenous gov. revenues
     real 
 
symms.params = [sigmaa, zetaa, eppsilon, alphaa, psii, thetac, thetad, Uppsilon, betaa, gammaa, etaa, G];     
symms.targets = [deltaa, kappaa];

%% Model f(yp, y, xp, x)=0 

%-- auxiliary stuff
 
% aggregate price level
p  = (pc^(1-eppsilon)+pd^(1-eppsilon))^(1/(1-eppsilon));
pp = (pcp^(1-eppsilon)+pdp^(1-eppsilon))^(1/(1-eppsilon));

% (disposable) income (tax system)

I    = hh*wh+hl*wl;
DI   = lambdaa*(I)^(1-tauul);
DIhh = lambdaa*(1-tauul)*(I)^(-tauul)*wh;
DIhl = lambdaa*(1-tauul)*(I)^(-tauul)*wl;

% Utility function 
if indic.util == 0
    % log utility, KPR
    Muc = c^(-gammaa);
    Muhh = -H^sigmaa*zetaa;
    Muhl = -H^sigmaa;
    
elseif indic.util==1
    % more general BGP
    Muc  = c^(-gammaa)*(H)^(-etaa*(1-gammaa));
    Muhh = -etaa*zetaa*H^(-1)*c*Muc;
    Muhl = -etaa*H^(-1)*c*Muc;

elseif indic.util==2
    % MaCurdy (as discussed in Boppart, Krusell)
    % sigmaa=1/theta
    % gammaa=sigmaa
    Muc  = c^(-gammaa);
    Muhh = -zetaa*H^(sigmaa);
    Muhl = -H^(sigmaa);
end


%-- model equations
q= 0; % to count model equations

%-- household
% budget
q=q+1;
f(q)= p*c-DI;

% consumption foc
q=q+1;
f(q) = Muc-p*muu;

% labour supply focs = more generally to be able to use different utility
% functions
% q=q+1;
% f(q)= H-(1-tauul)^(1/(1+sigmaa));
%q=q+1;
%f(q)= wh/wl-zetaa;

% high skill supply
q=q+1;
f(q) = -Muhh-muu*DIhh;

% low skill supply
q=q+1;
f(q) = -Muhl-muu*DIhl;

% definition total disutility hours
q=q+1;
f(q)= zetaa*hh+hl- H;

%-- final good production
% production 
q=q+1;
f(q)= Y-(yc^((eppsilon-1)/eppsilon)+yd^((eppsilon-1)/eppsilon))^(eppsilon/(eppsilon-1));

% profit max
q=q+1;
f(q)= yd-(pc/pd)^eppsilon*yc;

%-- clean sector
% production
q=q+1;
f(q)= yc-(alphaa/psii)^(alphaa/(1-alphaa))*pc^(alphaa/(1-alphaa))*Ac*Lc;

% labour demand
q=q+1;
f(q)= pcL-(1-alphaa)*(alphaa/psii)^(alphaa/(1-alphaa))*pc^(1/(1-alphaa))*Ac;

% machine demand: already substituted!

%-- dirty sector
%production
q=q+1;
f(q)= yd-(alphaa/psii)^(alphaa/(1-alphaa))*pd^(alphaa/(1-alphaa))*Ad*Ld;

% labour demand
q=q+1;
f(q)= pdL-(1-alphaa)*(alphaa/psii)^(alphaa/(1-alphaa))*pd^(1/(1-alphaa))*Ad;

% machine demand: already substituted!

%-- labour sector
% clean labour input good
q=q+1;
f(q)= Lc-lhc^thetac*llc^(1-thetac);

% dirty labour input good
q=q+1;
f(q)= Ld-lhd^thetad*lld^(1-thetad);

% skill demand
q=q+1;
f(q)= lhc-(pcL/wh*(thetac))^(1/(1-thetac))*llc;

q=q+1;
f(q)= llc-(pcL/wl*(1-thetac))^(1/(thetac))*lhc;

q=q+1;
f(q)= lhd-(pdL/wh*(thetad))^(1/(1-thetad))*lld;

q=q+1;
f(q)=lld-(pdL/wl*(1-thetad))^(1/(thetad))*lhd;

%-- productivity
% clean
q=q+1;
f(q)= Acp-(1+vc)*Ac;

%dirty
q=q+1;
f(q)= Adp-(1+vdd)*Ad;

% aggregate growth: write as inequality condition in gov problem! 
%c=c+1;
%f(c)=vc+vdd- Uppsilon;

%-- market clearing
% numeraire (aggregate goods market clears by walras' law)
q=q+1;
if indic.fullDisposal==0
    f(q) = p-1;
else
    f(q) = Y-c-psii*(xd+xc);
end
% high skill market
q=q+1;
f(q)= lhc+lhd-hh;

% low skill market
q=q+1;
f(q)= llc+lld-hl;

% additional variables to be tracked
%machines dirty
q=q+1;
f(q)= xd-(alphaa/psii*pd)^(1/(1-alphaa))*Ad*Ld;

% machines clean
q=q+1;
f(q)= xc-(alphaa/psii*pc)^(1/(1-alphaa))*Ac*Lc;

% government budget
q=q+1;
f(q)= G-(I-DI);

fprintf('number model equations: %d', q);

%% summarising variables

% exogenous variables
x  = [Ac Ad];
xp = [Acp Adp]; 

% endogenous variables
y  =[c wl wh hl hh pc pd yc yd Y Lc Ld pcL pdL lhc llc lhd lld H xc xd muu lambdaa];
yp =[cp wlp whp hlp hhp pcp pdp ycp ydp Yp Lcp Ldp pcLp pdLp lhcp llcp lhdp lldp Hp xcp xdp muup lambdaap];

% policy variables: in laissez faire as if parameters
pol  = [tauul vc vdd];
polp = [tauulp vcp vddp];

% marginal utilities
symms.marginals= [Muc, Muhh, Muhl];

% save lists of variables
list.x=string(x);
list.y=string(y);
list.xp=string(xp);
list.yp=string(yp);

list.pol     = string(pol);
list.params  = string(symms.params);
list.targets = string(symms.targets);
list.marginals = string(["Muc" "Muhh" "Muhl"]);

%% Government problem
% competitive equilibrium solution (symbolic)
solution_eqbm_syms;

%% lagrange muultiplier govs
syms muu_target real % exogenous emission target

% vector of symbolic variables for which to solve problem

if indic.withtarget==1
    symms.optim = sort([muu_target, tauul]);
elseif indic.withtarget==0
    symms.optim = sort(tauul);
end
list.optim  = string(symms.optim);

%% social welfare function and constraints

U=log(c)-(hl+zetaa*hh)^(1+sigmaa)/(1+sigmaa);

W        = U;                                     % value function 
target   = yd-(deltaa+E)/kappaa;
%growthAc = Acp - (1+vc)*Ac;
%growthAd = Adp - (1+vdd)*Ad;
%budget  = (H*wl-lambdaa*(H*wl)^(1-tauul))-G;


Obj_ram = W-muu_target*indic.withtarget*target; %...
           %- muu_growthAc*growthAc- muu_growthAd * growthAd;%-muu_budget*budget; %... % no beta as static; solution to static problem same as to infinite sum

    
Obj_sp = 0; %W + muuLr*lambdaa*(L-lr)+muuLp*(1-lambdaa)*(L-lp)...
%     +muun*(An*hn-lambdaa*cnr-(1-lambdaa)*cnp)...
%     +muus*(As*hs-lambdaa*csr-(1-lambdaa)*csp)...
%     +muul*(-lab_m);
%end




















