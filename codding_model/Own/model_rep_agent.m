%% Representative agent model

% laissez faire economy
% model written with symbolic variables

% with two skill types if zetaa!=1

% firm specific variables already substituted by optimality conditions
% machines, price of machines, and machine-specific productivity drop
%% Variables and parameters
% today and tomorrow

syms c cp ...           % consumption 
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
    ...
    ... % productivity
    Ac Acp ...          % productivity clean sector
    Ad Adp ...          % productivity dirty sector
    vd vdp ...          % productivity growth dirty sector - gov choice
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
     real 
 
symsparams = [sigmaa, zetaa, eppsilon, alphaa, psii, thetac, thetad, Uppsilon, betaa];     

%% Model f(yp, y, xp, x)=0 

% auxiliary stuff
 
p=(pc^(1-eppsilon)+pd^(1-eppsilon))^(1/(1-eppsilon));
pp=(pcp^(1-eppsilon)+pdp^(1-eppsilon))^(1/(1-eppsilon));

q= 0; % to count model equations

%-- household
% consumption foc
q=q+1;
f(q)= c- lambdaa*(H*wl)^(1-tauul);

% labour supply focs
q=q+1;
f(q)= H-(1-tauul)^(1/(1+sigmaa));

q=q+1;
f(q)= wh/wl-zetaa;

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
f(q)= Adp-(1+vd)*Ad;

% aggregate growth: write as inequality condition in gov problem! 
%c=c+1;
%f(c)=vc+vd- Uppsilon;

%-- market clearing
% numeraire (aggregate goods market clears by walras' law)
q=q+1;
f(q)= p-1;

% high skill market
q=q+1;
f(q)= lhc+lhd-hh;

% low skill market
q=q+1;
f(q)= llc+lld-hl;


fprintf('number model equations: %d', q);

%% summarising variables

% exogenous variables
x  = [Ac Ad];
xp = [Acp Adp]; 

% endogenous variables
y  =[c wl wh hl hh pc pd yc yd Y Lc Ld pcL pdL lhc llc lhd lld H];
yp =[cp wlp whp hlp hhp pcp pdp ycp ydp Yp Lcp Ldp pcLp pdLp lhcp llcp lhdp lldp Hp];

% policy variables: in laissez faire as if parameters
pol  = [lambdaa tauul vc vd];
polp = [lambdaap tauulp vcp vdp];

% save lists of variables
list.x=string(x);
list.y=string(y);
list.xp=string(xp);
list.yp=string(yp);

list.pol=string(pol);
list.params=string(symsparams);
