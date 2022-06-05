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
    if indic.ineq==0
        C=varrs(varlist=='C',:)';
    else
        Ch=varrs(varlist=='Ch',:)';
        Cl=varrs(varlist=='Cl',:)';
    end
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
    
% welfare measures
if indic.ineq==0
    if indic.BN==0
        if thetaa~=1
            Utilcon = (C.^(1-thetaa))./(1-thetaa);
        elseif thetaa==1
            Utilcon = log(C);
        end
    else
        Utilcon=-(C-B).^(zetaa)./(zetaa); 
    end
else
    if indic.BN==0
        if thetaa~=1
            Utilcon = zh.*(Ch.^(1-thetaa))./(1-thetaa)+(1-zh).*(Cl.^(1-thetaa))./(1-thetaa);
        elseif thetaa==1
            Utilcon = zh.*log(Ch)+(1-zh).*log(Cl);
        end
    else
        Utilcon=zh.*(-(Ch-Bh).^(zetaa)./(zetaa))+(1-zh).*(-(Cl-Bl).^(zetaa)./(zetaa)); 
    end

end

 Utillab = -chii.*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);

if indic.sep==0
      Utilsci =- chiis*S.^(1+sigmaas)./(1+sigmaas);
 else
      Utilsci = -(chiis*sff.^(1+sigmaas)./(1+sigmaas)+chiis*sg.^(1+sigmaas)./(1+sigmaas)+chiis*sn.^(1+sigmaas)./(1+sigmaas));
end

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

if indic.ineq==0
    CY = C./Y;
else
    CY = (zh*Ch+(1-zh).*Cl)./Y;
end
hhhl = hh./hl;
whwl = wh./wl;

%- growth rates
gAg = [(Ag(2:end)-Ag(1:end-1))./Ag(1:end-1);0]*100;
gAagg = [(A(2:end)-A(1:end-1))./A(1:end-1);0]*100;

gAn = [(An(2:end)-An(1:end-1))./An(1:end-1); 0]*100;
gAf = [(Af(2:end)-Af(1:end-1))./Af(1:end-1);0]*100;

%- update variables and varlist to include additional variables
jj= eval(symms.plotsvarsAdd); 

RESall(kk)=[varrs; jj'];

end % keys
end % function
