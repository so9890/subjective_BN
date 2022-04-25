function [Sparams, Spol, Starg, params, pol, targets]=parsSol(symms,trProd, trLab, paramss, list, poll, tarr) 
% parameters
read_in_pars_calib;
deltaa= tarr(list.tardir=='deltaa');

%soltions
cell_par=arrayfun(@char, symms.calib, 'uniform', 0);
SL=cell2struct(num2cell(trLab), cell_par, 2);
cell_par=arrayfun(@char, symms.prod, 'uniform', 0);
SP=cell2struct(num2cell(trProd), cell_par, 2);

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

params= eval(symms.params);
pol    = eval(symms.pol);
targets = eval(symms.targets);

cell_par=arrayfun(@char, symms.params, 'uniform', 0);
Sparams=cell2struct(num2cell(params), cell_par, 2);

cell_par=arrayfun(@char, symms.pol, 'uniform', 0);
Spol=cell2struct(num2cell(pol), cell_par, 2);

cell_par=arrayfun(@char, symms.targets, 'uniform', 0);
Starg=cell2struct(num2cell(targets), cell_par, 2);

end