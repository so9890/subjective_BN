function [c, ceq] = constraints(y, T, params, init, list, Ems, indic)
% function to read in constraints on government problem

% pars
read_in_params;

% transform x: all are exponentially transformed
 x=exp(y);

% except for hours
x((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T) = upbarH./(1+exp(y((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)));
x((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T) = upbarH./(1+exp(y((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)));
if indic.target == 0
    x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T) = upbS./(1+exp(y((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T)));
else
    x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T) = (y((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T)).^2;
end
x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = (y((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)).^2;

if indic.target==1
    Ftarget=(Ems'+deltaa)./omegaa; 
    x((find(list.opt=='F')-1)*T+1:find(list.opt=='F')*T)   = Ftarget./(1+exp(y((find(list.opt=='F')-1)*T+1:find(list.opt=='F')*T)));
end

%- auxiliary variables
if indic.taus==1 % with taus
[hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, hl, hh, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, ws, taus, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, wsgtil, S]= OPT_aux_vars(x, list, params, T, init, indic);
else
    [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, hl, hh, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, ws, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, S]= OPT_aux_vars_notaus(x, list, params, T, init, indic);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%    Inequality Constraints    %%%
 % only for direct periods
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % time period specific constraints! 
c = []; %  periods and 2 additional ones: 
c(1:T)=S-upbS;
if indic.target==1
    c(T+1:T*2) = F-Ftarget;
end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%    Equality Constraints    %%%
 % include missing equations here %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ceq = [];

 ceq(1:T)       = Ag-Ag_lag.*(1+gammaa.*(sg./rhog).^etaa.*(A_lag./Ag_lag).^phii);
 ceq(T+1:T*2)   = An-An_lag.*(1+gammaa.*(sn./rhon).^etaa.*(A_lag./An_lag).^phii);
 ceq(T*2+1:T*3) = N-(An.*Ln).*(pn.*alphan).^(alphan./(1-alphan)); % from production function neutral good
 % optimality skills (for fossil used to determine wage rates)
 ceq(T*3+1:T*4) = thetan*Ln.*wln-wh.*hhn; % optimality labour good producers neutral high skills
 ceq(T*4+1:T*5) = thetag*Lg.*wlg-wh.*hhg; % optimality labour good producers green high
 ceq(T*5+1:T*6) = (1-thetan)*Ln.*wln-wl.*hln; % optimality labour good producers neutral low
 ceq(T*6+1:T*7) = (1-thetag)*Lg.*wlg-wl.*hlg; % optimality labour good producers green low
 ceq(T*7+1:T*8) = Af-Af_lag.*(1+gammaa.*(sff./rhof).^etaa.*(A_lag./Af_lag).^phii);
 ceq(T*8+1:T*9) = C-zh.*lambdaa.*(wh.*hh).^(1-taul)-(1-zh).*lambdaa.*(wl.*hl).^(1-taul)-SGov;
% add foc for one skill type (ratio respected in taul derivation) 
 ceq(T*9+1:T*10) = chii*hh.^(sigmaa+taul)-(muu.*lambdaa.*(1-taul).*(wh).^(1-taul));
 ceq(T*10+1:T*11) = sg+sff+sn-S;
ceq = ceq';
end