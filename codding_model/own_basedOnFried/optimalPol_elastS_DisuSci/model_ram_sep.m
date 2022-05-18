% get objective function of ramsey planner
function [OB_RAM, list, symms, Ftarget]= model_ram_sep( list, params, T, init201519, indic, Ems, symms)

% prepare variables
syms    mu_IMP  mu_NProd ...
        mu_OPThhn mu_OPThln mu_OPThlg mu_OPThhg  mu_OPThh...
        mu_LOMAg mu_LOMAn mu_LOMAf mu_demSff mu_demSn mu_demSg... %KT_hh KT_hl KT_S...
        hhf hhg hlf hlg C F G Af Ag An ...
        HL HH sff sg sn mu_target real
    
% symms.optsym=[ mu_IMP mu_MarketS  mu_NProd mu_OPThhn mu_OPThh ... %KT_hh KT_hl KT_S ...
%         mu_OPThhg mu_OPThln mu_OPThlg mu_LOMAg mu_LOMAn mu_LOMAf...
%         hhf hhg hlf hlg C F G Af Ag An HL HH S]; % dropped: mu_wageG
%       

symms.optsym=[ mu_IMP   mu_NProd mu_OPThh ... %KT_hh KT_hl KT_S ...
          mu_OPThhn mu_OPThln mu_OPThlg mu_OPThhg mu_LOMAg mu_LOMAf mu_LOMAn mu_demSff mu_demSn mu_demSg...
        hhf hhg hlf hlg C F G Af Ag An HL HH sff sg sn]; % dropped: mu_wageG
      
if indic.target== 1
    symms.optsym=[symms.optsym mu_target];
end
list.optsym=string(symms.optsym);
    
% create symbolic variables
vecs=sym('a',[T,length([list.optsym])]);
    for s = [list.optsym]  % loop over list entries
        vecs(:,[list.optsym]==s)=sym(sprintf('%s%d',s),  [T,1]); 
    end
% stack into one column vector
x=vecs(:); 
list.optALL=string(x);
symms.optALL=x;

%% read in stuff
%- variables to be passed as x
read_in_SYMMODEL_sep;
%% social welfare

%- create discount vector
 disc=repmat(betaa, 1,T);
 expp=0:T-1;
 vec_discount= disc.^expp;

% %- vector of utilities
% if thetaa~=1
%  Utilcon = (C.^(1-thetaa))./(1-thetaa);
% elseif thetaa==1
%  Utilcon = log(C);
% end
%  Utillab = chii.*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);
%  Utilsci = chiis*S.^(1+sigmaas)./(1+sigmaas);
% 
%  SWF = Utilcon-Utillab-Utilsci;

 %% constraints
%IMP     = C-zh.*lambdaa.*(wh.*hh).^(1-taul)-(1-zh).*lambdaa.*(wl.*hl).^(1-taul)-SGov;
% replace by FOCS skill supply
 IMP     = C-zh.*lambdaa.*(wh.*hh).^(1-taul)-(1-zh).*lambdaa.*(wl.*hl).^(1-taul)-SGov;
%  MarketS = (sff+sg+sn)-S;  % LOM neutral technology 
 LOMAg   = Ag-Ag_lag.*(1+gammaa.*(sg./rhog).^etaa.*(A_lag./Ag_lag).^phii);
 LOMAn   = An-An_lag.*(1+gammaa.*(sn./rhon).^etaa.*(A_lag./An_lag).^phii);
 LOMAf   = Af-Af_lag.*(1+gammaa.*(sff./rhof).^etaa.*(A_lag./Af_lag).^phii);
 NProd   = N-(An.*Ln).*(pn.*alphan).^(alphan./(1-alphan)); % from production function neutral good
 OPThhn  = thetan*Ln.*wln-wh.*hhn; % optimality labour good producers neutral high skills
 OPThhg  = thetag*Lg.*wlg-wh.*hhg; % optimality labour good producers green high
 OPThln  = (1-thetan)*Ln.*wln-wl.*hln; % optimality labour good producers neutral low
 OPThlg  = (1-thetag)*Lg.*wlg-wl.*hlg; % optimality labour good producers green low
 OPThh   = chii*hh.^(sigmaa+taul)-(muu.*lambdaa.*(1-taul).*(wh).^(1-taul));
demSff   = wsf.*(rhof^etaa.*Af) - (gammaa*etaa*(A_lag./Af_lag).^phii.*sff.^(etaa-1).*pf.*(1-tauf).*F.*(1-alphaf).*Af_lag); 
demSg    = wsg.*(rhog^etaa.*Ag) - (gammaa*etaa*(A_lag./Ag_lag).^phii.*sg.^(etaa-1).*pg.*G.*(1-alphag).*Ag_lag);
demSn    = wsn.*(rhon^etaa.*An) - (gammaa*etaa*(A_lag./An_lag).^phii.*sn.^(etaa-1).*pn.*N.*(1-alphan).*An_lag);

Ftarget = (Ems(1:T)'+deltaa)/omegaa;
Target  = Ftarget-F; % emission target only enters constraints if indic.target==1 see below
% KTHH    = upbarH-hh;
% KTHL    = upbarH-hl;
% KTS     = upbarH-S;
 %% objective function 
 OB_RAM=vec_discount*(SWF)...
         + vec_discount*(-mu_IMP.*IMP- mu_demSff.*demSff- mu_demSn.*demSn- mu_demSg.*demSg...
         - mu_LOMAf.*LOMAf- mu_LOMAn.*LOMAn- mu_LOMAg.*LOMAg...
         - mu_NProd.*NProd-mu_OPThhn.*OPThhn...
         - mu_OPThhg.*OPThhg -mu_OPThln.*OPThln...
         - mu_OPThlg.*OPThlg- mu_OPThh.*OPThh);% - KT_hl.*KTHL- KT_hh.*KTHH- KT_S.*KTS);
% OB_RAM=vec_discount*(SWF)...
%         + vec_discount*(-mu_IMP.*IMP- mu_MarketS.*MarketS...
%         - mu_LOMAf.*LOMAf- mu_LOMAg.*LOMAg- mu_LOMAn.*LOMAn...
%         - mu_NProd.*NProd...
%         - mu_OPTLn.*OPTLn- mu_OPTLg.*OPTLg...
%         - mu_OPThh.*OPThh);% - KT_hl.*KTHL- KT_hh.*KTHH- KT_S.*KTS);

if indic.target==1
    OB_RAM=OB_RAM-sum(mu_target.*Target);
end
end