function [c, ceq] = constraints(y, T, params, init, list, Ems, indic)
% function to read in constraints on government problem

% pars
read_in_params;

% transform x: all are exponentially transformed
 x=exp(y);

% except for hours
if indic.noskill==0 %version with skill
    x((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T) = upbarH./(1+exp(y((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)));
    x((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T) = upbarH./(1+exp(y((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)));
else
    x((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T) = upbarH./(1+exp(y((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T)));
end

if indic.target == 0
    x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T) = upbarH./(1+exp(y((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T)));
else
    x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T) = (y((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T)).^2;
end

if indic.taus==1
    x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = (y((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)).^2;
end

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
    if indic.noskill==0
    [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, hl, hh, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, ws, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, S, gammac]= OPT_aux_vars_notaus(x, list, params, T, init, indic);
    else
%         fprintf('using noskill aux vars')
        [xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, h, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee,  ws,  tauf, taul, lambdaa,...
            w, SWF, S]= OPT_aux_vars_notaus_skillHom(x, list, params, T, init, indic);
    end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%    Inequality Constraints    %%%
 % only for direct periods
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % time period specific constraints! 
c = []; %  periods and 2 additional ones: 

c(1:T)=S-upbarH;
% gammac is the bgp growth rate should be at least as high as average
% growth in direct periods
% c(T+1:2*T) = Af(T)./Af_lag(T)-1-gammac; % the growth rate in T+1 should at least be as big as from T-1 to T
% 
% if indic.target==1
%     c(2*T+1:T*3) = F-Ftarget;
% end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%    Equality Constraints    %%%
 % include missing equations here %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ceq = [];

 if indic.noskill==0
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
         if indic.notaul==1
            ceq(T*11+1:T*12) = chii*hl.^(sigmaa+taul)-(muu.*lambdaa.*(1-taul).*(wl).^(1-taul));
         end

 elseif indic.noskill==1
     
        ceq(1:T)       = Ag-Ag_lag.*(1+gammaa.*(sg./rhog).^etaa.*(A_lag./Ag_lag).^phii);
        ceq(T+1:T*2)   = An-An_lag.*(1+gammaa.*(sn./rhon).^etaa.*(A_lag./An_lag).^phii);
        ceq(T*2+1:T*3) =  sg+sff+sn-S;%chii*h.^(sigmaa+taul)-(muu.*lambdaa.*(1-taul).*(w).^(1-taul));
        ceq(T*3+1:T*4) = Ln -pn.*(1-alphan).*N./w; % labour demand by sector 
        ceq(T*4+1:T*5) = Lg - pg.*(1-alphag).*G./w;
        ceq(T*5+1:T*6) = Lf- h./(1+Ln./Lf+Lg./Lf); % labour market clearing 
        ceq(T*6+1:T*7) = Af-Af_lag.*(1+gammaa.*(sff./rhof).^etaa.*(A_lag./Af_lag).^phii);
        ceq(T*7+1:T*8) = C-lambdaa.*(w.*h).^(1-taul)-SGov;
        if indic.notaul==1
            ceq(T*8+1:T*9)= chii*h.^(sigmaa+taul)-(muu.*lambdaa.*(1-taul).*(w).^(1-taul));
        end
 end
 %
ceq = ceq';
end