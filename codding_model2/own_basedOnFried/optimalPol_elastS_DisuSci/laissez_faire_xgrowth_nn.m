function f=laissez_faire_xgrowth_nn(x, params, list, pol, laggs, indic)
% Model
% equilibrium for one period!
% takes policy as given

% starts to solve the model in 2020-2024;
% i.e. Aj_laggs  refer to 2015-2019 period
% A_lag is except for the initial condition is 
%- read in policy and parameters
read_in_params;
read_in_pol;

%- initial condition
Af_lag=laggs(list.init=='Af0');
Ag_lag=laggs(list.init=='Ag0');

% choice variables
%- transform variables directly instead of in code
if indic.noskill==0
     hhf    = exp(x(list.choice=='hhf'));
     hhg    = exp(x(list.choice=='hhg'));
     hlf    = exp(x(list.choice=='hlf'));
     hlg    = exp(x(list.choice=='hlg'));
     hl     = upbarH/(1+exp(x(list.choice=='hl')));
     hh     = upbarH/(1+exp(x(list.choice=='hh')));
     gammall = x(list.choice=='gammall')^2;
     wh     = exp(x(list.choice=='wh'));
     wl     = exp(x(list.choice=='wl'));
else
    h      = upbarH/(1+exp(x(list.choice=='h')));
    w      = exp(x(list.choice=='w'));
    Lg     = exp(x(list.choice=='Lg'));
    Lf     = exp(x(list.choice=='Lf'));
end

if indic.ineq==0
    if indic.BN==0
      C      = exp(x(list.choice=='C'));
    else
      C      =B/(1+exp(x(list.choice=='C')));
    end
elseif indic.ineq==1
    if indic.BN==0
     Ch      = exp(x(list.choice=='Ch'));
     Cl      = exp(x(list.choice=='Cl'));
    else
      Ch     =Bh./(1+exp(x(list.choice=='Ch')));
      Cl     =Bl./(1+exp(x(list.choice=='Cl')));
    end
end
 F      = exp(x(list.choice=='F'));
 G      = exp(x(list.choice=='G'));
 gammalh = x(list.choice=='gammalh')^2;
 pg     = exp(x(list.choice=='pg'));
 pee     = exp(x(list.choice=='pee'));
 pf     = exp(x(list.choice=='pf'));
 lambdaa  = exp(x(list.choice=='lambdaa'));

%% - read in auxiliary equations
% A_lag   = (rhof*Af_lag+rhon*An_lag+rhog*Ag_lag)/(rhof+rhon+rhog);
%- technology
Ag=(1+vg)*Ag_lag;
Af=(1+vf)*Af_lag;

if indic.noskill==0
    Lg      = hhg.^thetag.*hlg.^(1-thetag);
    Lf      = hhf.^thetaf.*hlf.^(1-thetaf);

    SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
                +(1-zh)*(wl.*hl-lambdaa.*(wl*hl).^(1-taul))...
                +tauf.*pf.*F; % with scientists belonging to household sector, scientists sallary dont cancel from profits
                % subsidies and profits and wages scientists cancel
else
    SGov    = (w.*h-lambdaa.*(w.*h).^(1-taul))...
            +tauf.*pf.*F;  
end

if indic.ineq==0
    if indic.BN==0
        muu   = C.^(-thetaa); % same equation in case thetaa == 1
    else
        muu   = -(C-B).^(zetaa-1);
    end
elseif indic.ineq==1
    if indic.BN==0
        muuh   = Ch.^(-thetaa); % same equation in case thetaa == 1
        muul   = Cl.^(-thetaa); % same equation in case thetaa == 1

    else
        muul   = -(Cl-Bl).^(zetaa-1);
        muuh    = -(Ch-Bh).^(zetaa-1);
    end

end
E       = (F.^((eppse-1)/eppse)+G.^((eppse-1)/eppse)).^(eppse/(eppse-1));

Y       = E; % (deltay^(1/eppsy).*E.^((eppsy-1)/eppsy)+(1-deltay)^(1/eppsy).*N.^((eppsy-1)/eppsy)).^(eppsy/(eppsy-1)); % production function Y 

wlg     = pg.^(1/(1-alphag)).*(1-alphag).*alphag.^(alphag/(1-alphag)).*Ag;
wlf     = (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; 

%% model equations
q=0;

%1- household optimality (muu auxiliary variable determined above)
if indic.noskill==0
    
    if indic.ineq==0
        q=q+1;
        f(q)= chii*hh.^(sigmaa+taul)- ((muu.*lambdaa.*(1-taul).*(wh).^(1-taul))-gammalh./zh.*hh.^taul); %=> determines hh
        %2
        q=q+1;
        f(q)= chii*hl.^(sigmaa+taul) - ((muu.*lambdaa.*(1-taul).*(wl).^(1-taul))-gammall./(1-zh).*hl.^taul); %=> determines hl
        %3- budget
        q=q+1;
        f(q) = zh.*lambdaa.*(wh.*hh).^(1-taul)+(1-zh).*lambdaa.*(wl.*hl).^(1-taul)+SGov-C; %=> determines C

    else   
        q=q+1;
        f(q)= chii*hh.^(sigmaa+taul)- ((muuh.*lambdaa.*(1-taul).*(wh).^(1-taul))-gammalh.*hh.^taul); %=> determines hh
        %2
        q=q+1;
        f(q)= chii*hl.^(sigmaa+taul) - ((muul.*lambdaa.*(1-taul).*(wl).^(1-taul))-gammall.*hl.^taul); %=> determines hl
    
        %3- budget
        q=q+1;
        f(q) = lambdaa.*(wl.*hl).^(1-taul)+SGov-Cl; %=> determines C
        q=q+1;
        f(q) = lambdaa.*(wh.*hh).^(1-taul)+SGov-Ch;
    end
else
    q=q+1;
    f(q)= chii*h.^(sigmaa+taul)- ((muu.*lambdaa.*(1-taul).*(w).^(1-taul))-gammalh.*h.^taul); %=> determines hh
   %3- budget
    q=q+1;
    f(q) = lambdaa.*(w.*h).^(1-taul)+SGov-C; %=> determines C
end

%4- output fossil
q=q+1;
f(q) = ((1-tauf).*alphaf.*pf).^(alphaf./(1-alphaf)).*Af.*Lf -F; 

%6- output green
q=q+1;
f(q)=  G-Ag.*Lg.*(pg.*alphag).^(alphag./(1-alphag));

%13- optimality labour input producers
if indic.noskill==0

    %14
    q=q+1;
    f(q)= thetag*Lg.*wlg-wh.*hhg;
    %15
    q=q+1;
    f(q)=(1-thetag)*Lg.*wlg-wl.*hlg;
    %18- demand skill
    q=q+1;
    f(q) = wh - thetaf*Lf./hhf.*wlf; % from optimality labour input producers fossil, and demand labour fossil
    %19
    q=q+1;
    f(q) = wl-(1-thetaf)*Lf./hlf.*wlf;

else
    q=q+1; % labour demand fossil
    f(q) =  w  - (1-alphaf)*alphaf^(alphaf/(1-alphaf)).*((1-tauf).*pf).^(1/(1-alphaf)).*Af; % labour demand fossil
    q=q+1;
    f(q) = Lg - pg.*(1-alphag).*G./w;
end

% prices and wages
%17- optimality energy producers
q=q+1;
f(q) = pf.*F.^(1/eppse)- (G).^(1/eppse).*pg; 


%- definitions prices
%20
q=q+1;
f(q) = pee - (pf.^(1-eppse)+pg.^(1-eppse)).^(1/(1-eppse)); %definition
%21
q=q+1;
f(q) =  1-pee;

%22- market clearing (consumption good=> numeraire)
if indic.noskill==0

    q=q+1;
    f(q) = zh.*hh-(hhf+hhg); % high skill market clearing
    %23
    q=q+1;
    f(q) = (1-zh).*hl-(hlf+hlg); % low skill market clearing

    %13- Kuhn Tucker Labour supply and scientists
    %25
    q=q+1;
    f(q)= gammalh.*(upbarH-hh);
    %26
    q=q+1;
    f(q)= gammall.*(upbarH-hl);
else
    q=q+1;
    f(q) = Lf+Lg-h;
    q=q+1;
    f(q)= gammalh.*(h-upbarH);
end

% balanced budget
q=q+1;
f(q)= SGov;
% fprintf('number equations: %d; number variables %d', q, length(list.choice));
end