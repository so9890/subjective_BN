%% Model file of the simple, static model
% no innovation 
% no monopolistic competition
% no labour sector choice
clear, clc

cd('/home/sonja/Documents/projects/Overconsumption/codding_model')
% variables
syms cnh cnl csh csl Ch Cl hh hl yn ys lls lln lhs lhn ...
     psL pnL ps pn wh wl Ls Ln real

%% vector of variables for which to solve model
 symvariables=[cnh, cnl, csh, csl, Ch, Cl, hh, hl, yn,ys, lls, lln, lhs, lhn, ...
     psL, pnL, ps, pn, wh, wl, Ls, Ln];
% parameters
syms eppsilons eppsilonn omegas alphaa psii B lambdaa As An real

%-- model equations
c=0;
% HHs

%high skill
c=c+1;
f(c)=csh-wh^2*(omegas/ps)^(1+omegas)*((1-omegas)/pn)^(1-omegas)/(1+wh^2*(omegas/ps)^(2*omegas)*((1-omegas)/pn)^(2*(1-omegas)))*B;
c=c+1;
f(c)=cnh-wh^2*(omegas/ps)^(omegas)*((1-omegas)/pn)^(2-omegas)/(1+wh^2*(omegas/ps)^(2*omegas)*((1-omegas)/pn)^(2*(1-omegas)))*B;
c=c+1;
f(c)=hh- wh*(omegas/ps)^(omegas)*((1-omegas)/pn)^(1-omegas)/(1+wh^2*(omegas/ps)^(2*omegas)*((1-omegas)/pn)^(2*(1-omegas)))*B;
c=c+1;
f(c)=Ch-csh^omegas*cnh^(1-omegas);

% low skill
c=c+1;
f(c)=csl-wl^2*(omegas/ps)^(1+omegas)*((1-omegas)/pn)^(1-omegas)/(1+wl^2*(omegas/ps)^(2*omegas)*((1-omegas)/pn)^(2*(1-omegas)))*B;
c=c+1;
f(c)=cnl-wl^2*(omegas/ps)^(omegas)*((1-omegas)/pn)^(2-omegas)/(1+wl^2*(omegas/ps)^(2*omegas)*((1-omegas)/pn)^(2*(1-omegas)))*B;
c=c+1;
f(c)=hl- wl*(omegas/ps)^(omegas)*((1-omegas)/pn)^(1-omegas)/(1+wl^2*(omegas/ps)^(2*omegas)*((1-omegas)/pn)^(2*(1-omegas)))*B;
c=c+1;
f(c)=Cl-csl^omegas*cnl^(1-omegas);

% Production 

% production labour input good
c=c+1;
f(c)=Ls-lhs^eppsilons*lls^(1-eppsilons);
c=c+1;
f(c)=Ln-lhn^eppsilonn*lln^(1-eppsilonn);

% optimality labour producing firm
c=c+1;
f(c)=lhs-(psL*eppsilons/wh)^(1/(1-eppsilons))*lls;
c=c+1;
f(c)=lhn-(pnL*eppsilonn/wh)^(1/(1-eppsilonn))*lln;

c=c+1;
f(c)=lls-(psL*(1-eppsilons)/wl)^(1/(eppsilons))*lhs;
c=c+1;
f(c)=lln-(pnL*(1-eppsilonn)/wl)^(1/(eppsilonn))*lhn;

% production final good
c=c+1;
f(c)=ys-(alphaa*ps/psii)^(alphaa/(1-alphaa))*As*Ls;
c=c+1;
f(c)=yn-(alphaa*pn/psii)^(alphaa/(1-alphaa))*An*Ln;

% labour input demand => price for labour good
c=c+1;
f(c)=psL-(1-alphaa)*(alphaa/psii)^(alphaa/(1-alphaa))*ps^(1/(1-alphaa))*As;
c=c+1;
f(c)=pnL- (1-alphaa)*(alphaa/psii)^(alphaa/(1-alphaa))*pn^(1/(1-alphaa))*An;

% machine demand
% => already replaced!
%xis=(alphaa*ps/pis)^(1/(1-alphaa))*Ls*Ais

% optimality machine producers=> already replaced: pis=pin=psii

% Market clearing

% sust final good
c=c+1;
f(c)=ys-( lambdaa*csh+(1-lambdaa)*csl);

% unsust final good => walras' law
c=c+1;
f(c)=pn-1;

% high skill labour market
c=c+1;
f(c)=hh-(lhs+lhn);
% low skill labour market
c=c+1;
f(c)=hl-(lls+lln);


%--end model block
fprintf('model equations %d', length(f))
fprintf('variables %d', length(symvariables))

%-- Solve symbolically
S=solve(f, symvariables);