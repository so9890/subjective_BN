%% solution to competitive eqbm 

% can be read in as numeric or symbolic version 
% depending on how variables are written.

gammad = (thetad/(zetaa*(1-thetad)))^thetad;
gammac = (thetac/(zetaa*(1-thetac)))^thetac;
chii   = (1-thetad)*(1-thetac)/(thetac-thetad);

chitilde =  (thetac^thetac*thetad^(-thetad))^((1-alphaa)*(1-eppsilon))...
            *(1-thetac)^(-thetac-(1-thetac)*(alphaa+eppsilon*(1-alphaa)))...
            *(1-thetad)^(thetad+(1-thetad)*(alphaa+eppsilon*(1-alphaa)));
        
helpper  = (Ac/Ad)^((1-alphaa)*(1-eppsilon))...
            *zetaa^(-(thetac-thetad)*(1-alphaa)*(1-eppsilon))*chitilde;
        
kappatilde = ((1-thetac)*(1-thetad)*(helpper+1))...
            /((1-thetad)+(1-thetac)*(helpper)); 
zd      = thetad^thetad*(1-thetad)^(1-thetad);
zc      = thetac^thetac*(1-thetac)^(1-thetac);

%% solution equilibrium

% skills
H  = (1-tauul)^(1/(1+sigmaa)); 
hl = kappatilde*H; 
hh = 1/zetaa*(H-hl);

% labour input good
Lc= gammac*chii*(1-kappatilde/(1-thetad))*H;
Ld= gammad*chii*(kappatilde/(1-thetac)-1)*H;

% sector good prices
pd = ((Ad/Ac*(zd/zc)* zetaa*(thetac-thetad))^((1-alphaa)*(1-eppsilon))+1)^(-1/(1-eppsilon));
pc = pd*(Ad/Ac*zd/zc*zetaa^(thetac-thetad))^(1-alphaa);

% labour input prices
pcL = (1-alphaa)*(alphaa/psii)^(alphaa/(1-alphaa))*pc^(1/(1-alphaa))*Ac;
pdL = (1-alphaa)*(alphaa/psii)^(alphaa/(1-alphaa))*pd^(1/(1-alphaa))*Ad;

% skill specific wage
wl = pdL*zetaa^(-thetad)*thetad^thetad*(1-thetad)^(1-thetad);
wh = zetaa*wl;

% consumption
c = lambdaa*(H*wl)^(1-tauul);

% skill inputs by sector
llc = Lc/gammac;
lld = Ld/gammad;

lhc = gammac^(1/thetac)*llc;
lhd = gammad^(1/thetad)*lld;

% sector output
yc = (alphaa/psii*pc)^(alphaa/(1-alphaa))*Ac*Lc;

yd= (pc/pd)^eppsilon*yc;

% technology in next period
Acp = (1+vc)*Ac;
Adp = (1+vd)*Ad; 

% machines
xd = (alphaa/psii*pd)^(1/(1-alphaa))*Ad*Ld;
xc = (alphaa/psii*pc)^(1/(1-alphaa))*Ac*Lc;

% gov budget
G = (H*wl-lambdaa*(H*wl)^(1-tauul));

% goods market clearing
Y = c+psii*(xd+xc)+G;

% shadow value income
mu= 1/c;