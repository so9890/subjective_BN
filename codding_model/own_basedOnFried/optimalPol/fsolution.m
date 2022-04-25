function [x0LF, SL, SP, SR, Sall, Sinit, init]=fsolution(symms, trProd, trLab, trR, paramss, list, poll, MOM) 

read_in_pars_calib;
%soltions
%- as structure
cell_par=arrayfun(@char, symms.calib, 'uniform', 0);
SL=cell2struct(num2cell(trLab), cell_par, 2);
cell_par=arrayfun(@char, symms.prod, 'uniform', 0);
SP=cell2struct(num2cell(trProd), cell_par, 2);
cell_par=arrayfun(@char, symms.calib3, 'uniform', 0);
SR=cell2struct(num2cell(trR), cell_par, 2);

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
deltay=SP.deltay;
omegaa=SP.omegaa;
thetan = SL.thetan;
thetaf = SL.thetaf;
thetag = SL.thetag;
el = SL.el;
eh = SL.eh;
chii = SL.chii;
lambdaa = SL.lambdaa;

% Remaining variables
hln = hhn*(1-thetan)/(thetan)*MOM.whwl; % hln
hlf = hhf*(1-thetaf)/(thetaf)*MOM.whwl; % hlf
hlg = hhg*(1-thetag)/(thetag)*MOM.whwl; % hlg 

Lg = hhg.^thetag.*hlg.^(1-thetag);
Ln = hhn.^thetan.*hln.^(1-thetan);
Lf = hhf.^thetaf.*hlf.^(1-thetaf);

[ C, Lnwln, Lgwlg, Lfwlf, pf, F, pn, pg, pee, E, Y, N, G, xn, xg, xf, ...
            AfLf, AgLg, AnLn, omegaa, deltay]=resProd(list, trProd, MOM, paramss, poll);

muu = C^(-thetaa);

Af   = AfLf/Lf; % production 
Ag   = AgLg/Lg; % production 
An   = AnLn/Ln; % production 

wln   = Lnwln/Ln; % price labour input neutral sector
wlg   = Lgwlg/Lg;
wlf   = Lfwlf/Lf; 

sff = SR.sff;
sg = SR.sg;
sn = SR.sn;
Af0 = SR.Af_lag;
Ag0 = SR.Ag_lag;
An0 = SR.An_lag;
ws = SR.ws;

 A0  = (rhof*Af0+rhon*An0+rhog*Ag0)/(rhof+rhon+rhog);
 Agtest= Ag0*(1+gammaa*(sg/rhog)^etaa*(A0/Ag0)^phii); 
 
 if abs(Ag-Agtest)>1e-10
     error('growth rate off')
 else
     fprintf('growth rate works')
 end

% save
cell_par=arrayfun(@char, symms.allvars, 'uniform', 0);
Sall=cell2struct(num2cell(eval(symms.allvars)), cell_par, 2);

init = eval(symms.init);
cell_par=arrayfun(@char, symms.init, 'uniform', 0);
Sinit=cell2struct(num2cell(init), cell_par, 2);


%- save for laissez faire
x0LF= eval(symms.choice);
% test
zh=paramss(list.paramsdir=='zh');
zl=paramss(list.paramsdir=='zl');
tauf = poll(list.poldir=='tauf');
Cincome= zh*eh*hh*wh+zl*wl*el*hl+tauf*pf*F;

if abs(Cincome-C)>1e-10
    error('goods market does not clear!')
else
    fprintf('goods market clears!')
end

if abs(hh*zh*eh-hhn-hhf-hhg)>1e-10
    error('high skill market does not clear!')
else
    fprintf('high skill market clears!')
end
if abs(zl*el*hl-(hlg+hlf+hln))>1e-10
    error('low skill market does not clear!')
else
    fprintf('low skill market clears!')
end
end