function [Sparams, Spol, params, pol]=parsSol_GOOD(symms,trProd, trLab, resSci, paramss, list, poll, phis) 
% parameters
read_in_pars_calib;

%soltions
cell_par=arrayfun(@char, symms.calib, 'uniform', 0);
SL=cell2struct(num2cell(trLab), cell_par, 2);
cell_par=arrayfun(@char, symms.prod, 'uniform', 0);
SP=cell2struct(num2cell(trProd), cell_par, 2);
cell_par=arrayfun(@char, symms.calib3, 'uniform', 0);
SciL=cell2struct(num2cell(resSci), cell_par, 2);

% parameters
deltay=SP.deltay;
omegaa=SP.omegaa;
thetan = SL.thetan;
thetaf = SL.thetaf;
thetag = SL.thetag;
zh = SL.zh;
chii = SL.chii;
lambdaa = SL.lambdaa;
chiis=SciL.chiis;
gammaa=SciL.gammaa;
% rhog=SciL.rhog;
% rhon=SciL.rhon;
% rhof=SciL.rhof;

params = eval(subs(symms.params));
pol    = eval(symms.pol);

cell_par=arrayfun(@char, symms.params, 'uniform', 0);
Sparams=cell2struct(num2cell(params), cell_par, 2);

cell_par=arrayfun(@char, symms.pol, 'uniform', 0);
Spol=cell2struct(num2cell(pol), cell_par, 2);

end