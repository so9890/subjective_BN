function [RESall]=add_vars(RES, list, params, indic, varlist, symms)

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
    E = varrs(varlist=='E',:)';
    Y = varrs(varlist=='Y',:)';
    tauf = varrs(varlist=='tauf',:)';
    pf = varrs(varlist=='pf',:)';
    
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
if indic.xgrowth==0
    sgsff= sg./sff;
else
    sgsff=zeros(size(AgAf));
end
GFF = G./F;
EY= E./Y;
LgLf = Lg./Lf;

CY = C./Y;
hhhl = hh./hl;
whwl = wh./wl;

%- growth rates
gAg = [(Ag(2:end)-Ag(1:end-1))./Ag(1:end-1);0]*100;
gAagg = [(A(2:end)-A(1:end-1))./A(1:end-1);0]*100;

gAn = [(An(2:end)-An(1:end-1))./An(1:end-1); 0]*100;
gAf = [(Af(2:end)-Af(1:end-1))./Af(1:end-1);0]*100;

% analytical measure of taul in integrated policy regime
analyTaul = tauf.*pf.*F./Y;

%- update variables and varlist to include additional variables
jj= eval(symms.plotsvarsAdd); 

RESall(kk)=[varrs; jj'];

end % keys
end % function
