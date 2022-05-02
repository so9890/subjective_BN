function [Sparams, Spol, params, pol]=parsSol(symms,trProd, trLab, paramss, list, poll) 
% parameters
read_in_pars_calib;

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
zh = SL.zh;
chii = SL.chii;
lambdaa = SL.lambdaa;


params = eval(symms.params);
pol    = eval(symms.pol);

cell_par=arrayfun(@char, symms.params, 'uniform', 0);
Sparams=cell2struct(num2cell(params), cell_par, 2);

cell_par=arrayfun(@char, symms.pol, 'uniform', 0);
Spol=cell2struct(num2cell(pol), cell_par, 2);

end