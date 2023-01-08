function [RESall]=add_vars(RES, list, params, indic, varlist, symms, MOM)

read_in_params;

%- new container for additional variables
RESall=RES;

% function to calculate additional variables for graphs
for k = keys(RES)
    kk =string(k);
    %- extract variable matrix
    varrs= RES(kk);
    
    %- read in variables
    C=varrs(varlist=='C',:)';
    hh = varrs(varlist=='hh',:)';
    hl = varrs(varlist=='hl',:)';
    sg = varrs(varlist=='sg',:)';
    sff = varrs(varlist=='sff',:)';
    sn = varrs(varlist=='sn',:)';
    S = varrs(varlist=='S',:)';
    Ag = varrs(varlist=='Ag',:)';
    Af = varrs(varlist=='Af',:)';
    An = varrs(varlist=='An',:)';
    A = varrs(varlist=='A',:)';
    Lg = varrs(varlist=='Lg',:)';
    Lf = varrs(varlist=='Lf',:)';
    G = varrs(varlist=='G',:)';
    F = varrs(varlist=='F',:)';
    wh = varrs(varlist=='wh',:)';
    wl = varrs(varlist=='wl',:)';
    ws = varrs(varlist=='ws',:)';

    E = varrs(varlist=='E',:)';
    Y = varrs(varlist=='Y',:)';
    tauf = varrs(varlist=='tauf',:)';
    taul = varrs(varlist=='taul',:)';
    lambdaa = varrs(varlist=='lambdaa',:)';

    pf = varrs(varlist=='pf',:)';
    pee = varrs(varlist=='pee',:)';
    pg = varrs(varlist=='pg',:)';
    pn = varrs(varlist=='pn',:)';

% welfare measures

if thetaa~=1
    Utilcon = (C.^(1-thetaa))./(1-thetaa);
elseif thetaa==1
    Utilcon = log(C);
end

Utillab = chii.*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);

if indic.sep==0
      Utilsci = chiis*S.^(1+sigmaas)./(1+sigmaas);
 else
      Utilsci = chiis*sff.^(1+sigmaas)./(1+sigmaas)+chiis*sg.^(1+sigmaas)./(1+sigmaas)+chiis*sn.^(1+sigmaas)./(1+sigmaas);
 end
%  SWFtt = Utilcon-Utillab-Utilsci;
% continuation value
% T= length(Y);
gammay = Y(end)/Y(end-1)-1;
PVconsump= 1/(1-betaa*(1+gammay)^(1-thetaa))*Utilcon(end);
PVwork =indic.PVwork* 1/(1-betaa)*(Utillab(end)+Utilsci(end));
PV= ones(size(Y)).*betaa^length(Y)*(PVconsump-PVwork); % continuation value in period 0


% ratios
AgAf=Ag./Af;
if max(sff~=zeros(size(sff)))
    sgsff= sg./sff;
    snS= sn./S;
    sffS =sff./S;
    sgS =sg./S; 
    sffsg=1./sgsff;
else
    sgsff=zeros(size(AgAf));
    snS=zeros(size(AgAf)); 
    sgS=snS;
    sffS=snS;
    sffsg=sgsff; 
end


GFF = G./F;
EY= E./Y;
LgLf = Lg./Lf;
pgpftf =pg./(pf+tauf);
pepn = pee./pn;
CY = C./Y;
hhhl = hh./hl;
whwl = wh./wl;
Hagg = zh.*hh+(1-zh).*hl;
%- growth rates
gAg = [(Ag(2:end)-Ag(1:end-1))./Ag(1:end-1);0]*100;
gAagg = [(A(2:end)-A(1:end-1))./A(1:end-1);0]*100;

gAn = [(An(2:end)-An(1:end-1))./An(1:end-1); 0]*100;
gAf = [(Af(2:end)-Af(1:end-1))./Af(1:end-1);0]*100;

% analytical measure of taul in integrated policy regime
analyTaul = tauf.*pf.*F./Y;

%- average tax rates 
% dTaulHH dTaulHl: marginal tax rate
% dTaulAv: average income weighted marginal tax rate

dTaulHh = (1-(1-taul).*lambdaa.*(wh.*hh).^(-taul)).*100;
dTaulHl = (1-(1-taul).*lambdaa.*(wl.*hl).^(-taul)).*100;
dTaulS = (1-(1-taul).*lambdaa.*(ws.*S).^(-taul)).*100;

dTaulAv = 1./((1-zh)*wl.*hl+zh*wh.*hh).*(zh*wh.*hh.*dTaulHh+(1-zh).*wl.*hl.*dTaulHl);
dTaulAvS = 1./((1-zh)*wl.*hl+zh*wh.*hh+ws.*S).*(zh*wh.*hh.*dTaulHh+(1-zh).*wl.*hl.*dTaulHl+ws.*S.*dTaulS);
%- tauf in tons per C02

tauf_CO2=tauf./omegaa;
tauf_perton2019 = tauf_CO2*(MOM.GDP1519MILLION*1e6)./(1e9); % denominator to go from gigaton to ton in 2019 prices
Tauf=tauf_perton2019*1.12; % to have it in 2022 prices
%       
%- update variables and varlist to include additional variables
jj= eval(symms.plotsvarsAdd); 

RESall(kk)=[varrs; jj'];

end % keys
end % function
