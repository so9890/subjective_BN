function [x0LF, SL, SP, SR, Sall, ...
    Sinit201014, init201014 , Sinit201519, init201519]=fsolution(symms, trProd, trLab, resSci, paramss, list, poll, MOM) 

read_in_pars_calib;
%soltions
%- as structure
cell_par=arrayfun(@char, symms.calib, 'uniform', 0);
SL=cell2struct(num2cell(trLab), cell_par, 2);
cell_par=arrayfun(@char, symms.prod, 'uniform', 0);
SP=cell2struct(num2cell(trProd), cell_par, 2);
% for science: calib3 contains phis from inital run and the rest from final
% run with phis fixed
cell_par=arrayfun(@char, symms.calib3, 'uniform', 0);
SR=cell2struct(num2cell(resSci), cell_par, 2);

% - division into parameters and variables
hhn = SL.hhn;
hhg = SL.hhg;
hhf = SL.hhf;
gammalh = SL.gammalh;
gammall = SL.gammall;
hl    = SL.hl;
hh    = SL.hh;
wh = SL.wh;
wl = SL.wl;
pn=SP.pn;
pg=SP.pg;

% parameters
deltay = SP.deltay;
omegaa = SP.omegaa;
thetan = SL.thetan;
thetaf = SL.thetaf;
thetag = SL.thetag;
zh     = SL.zh;
chii   = SL.chii;
lambdaa = SL.lambdaa;
phii    = SR.phii;
chiis  = SR.chiis;
etaa  = SR.etaa;
% rhog  = SR.rhog;
% rhon  = SR.rhon;
% rhof  = SR.rhof;
gammaa  = SR.gammaa;


% Remaining variables
hln = hhn*(1-thetan)/(thetan)*MOM.whwl; % hln
hlf = hhf*(1-thetaf)/(thetaf)*MOM.whwl; % hlf
hlg = hhg*(1-thetag)/(thetag)*MOM.whwl; % hlg 

Lg = hhg.^thetag.*hlg.^(1-thetag);
Ln = hhn.^thetan.*hln.^(1-thetan);
Lf = hhf.^thetaf.*hlf.^(1-thetaf);

[ C, Lnwln, Lgwlg, Lfwlf, pf, F, pn, pg, pee, E, Y, N, G, xn, xg, xf, ...
            AfLf, AgLg, AnLn, omegaa, deltay]=resProd_GE(list, trProd, MOM, ...
            paramss, poll, 'calib');

muu = C^(-thetaa);

Af   = AfLf/Lf; % production 
Ag   = AgLg/Lg; % production 
An   = AnLn/Ln; % production 
% to save as initial values
Af0=Af; An0=An; Ag0=Ag;
init201519 = eval(symms.init);
clearvars Af0 Ag0 An0

wln   = Lnwln/Ln; % price labour input neutral sector
wlg   = Lgwlg/Lg;
wlf   = Lfwlf/Lf; 

sff = SR.sff;
sg = SR.sg;
sn = SR.sn;
S = sff+sg+sn; % determined by chii in science problem

% wse = SR.wse;
wsf = SR.wsf;
wsg = SR.wsg;
wsn = SR.wsn;
% chiisT=ws/(C^thetaa*S^sigmaas);
Af0 = SR.Af_lag; % 2010-2014
Ag0 = SR.Ag_lag;
An0 = SR.An_lag;

A  = (rhof*Af+rhon*An+rhog*Ag)/(rhof+rhon+rhog);

A0  = (rhof*Af0+rhon*An0+rhog*Ag0)/(rhof+rhon+rhog);

Agtest= Ag0*(1+gammaa*(sg/rhog)^etaa*(A0/Ag0)^phii); 
 
 if abs(Ag-Agtest)>1e-10
     error('growth rate off')
 else
     fprintf('growth rate works!!')
 end
 
SGov    = zh*(wh.*hh-lambdaa.*(wh.*hh).^(1-taul))...
            +(1-zh)*(wl.*hl-lambdaa.*(wl*hl).^(1-taul))...
            +tauf.*pf.*F;
Emnet     = omegaa*F-deltaa; % net emissions

if thetaa~=1
 Utilcon = (C.^(1-thetaa))./(1-thetaa);
elseif thetaa==1
 Utilcon = log(C);
end
 Utillab = chii.*(zh.*hh.^(1+sigmaa)+(1-zh).*hl.^(1+sigmaa))./(1+sigmaa);
 SWF = Utilcon-Utillab;

% save
cell_par=arrayfun(@char, symms.allvars, 'uniform', 0);
Sall=cell2struct(num2cell(eval(symms.allvars)), cell_par, 2);

init201014 = eval(symms.init);
cell_par=arrayfun(@char, symms.init, 'uniform', 0);
Sinit201014=cell2struct(num2cell(init201014), cell_par, 2);
Sinit201519=cell2struct(num2cell(init201519), cell_par, 2);



%- save for laissez faire
gammasf =0;
gammasg =0;
gammasn =0;
x0LF= eval(symms.choice);
% test
zh=paramss(list.paramsdir=='zh');
tauf = poll(list.poldir=='tauf');
Cincome= zh*hh*wh+(1-zh)*wl*hl+tauf*pf*F;

if abs(Cincome-C)>1e-10
    error('goods market does not clear!')
else
    fprintf('goods market clears!')
end

if abs(hh*zh-hhn-hhf-hhg)>1e-10
    error('high skill market does not clear!')
else
    fprintf('high skill market clears!')
end
if abs((1-zh)*hl-(hlg+hlf+hln))>1e-10
    error('low skill market does not clear!')
else
    fprintf('low skill market clears!')
end
end