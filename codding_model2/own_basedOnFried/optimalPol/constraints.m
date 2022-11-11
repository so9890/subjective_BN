function [c, ceq] = constraints(y, T, params, init, list, Ems, indic)
% function to read in constraints on government problem

% pars
read_in_params;

% transform x: all are exponentially transformed
 x=exp(y);

% except for taus
x(T*(find(list.opt=='taus')-1)+1:T*(find(list.opt=='taus')))=y(T*(find(list.opt=='taus')-1)+1:T*(find(list.opt=='taus'))) ;

% except for hours
x((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T) = upbarH./(1+exp(y((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)));
x((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T) = upbarH./(1+exp(y((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)));
if indic.target==1
    Ftarget=(Ems'+deltaa)./omegaa; 
    x((find(list.opt=='F')-1)*T+1:find(list.opt=='F')*T)   = Ftarget./(1+exp(y((find(list.opt=='F')-1)*T+1:find(list.opt=='F')*T)));
end

%- auxiliary variables
[hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, hl, hh, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, wsf, wsn, wsg, taus, tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, wsgtil]= OPT_aux_vars(x, list, params, T, init, indic);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%    Inequality Constraints    %%%
 % only for direct periods
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % time period specific constraints! 
c = zeros(T,1); %  periods and 2 additional ones: 
             %  IMP; skill time endowment

             % the constraints valid for T periods are:

             % 1) budget constraint HH (IMP)
                         

% 4) emission targets coded as linear constraint!
%%% 1. Implementability constraint %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- create discount vector
      
% rather without savings technology: one binding constraint per period
% c(1:T)= C-zh.*wh.*hh-(1-zh).*wl.*hl-tauf.*pf.*F;
% c(1:T)= C-zh.*wh.*hh-(1-zh).*wl.*hl-tauf.*pf.*F;

%%% 2./3. Time endowment labour: Dropped due to 
%%% transformation of variables! %%%
% => can drop kuhn tuckker constraint! in Household problem, the government
% always chooses an interior solution! 
% c(T*1+1:T*2)=hh-upbarH;
% c(T*2+1:T*3)=hl-upbarH;

%%% 4. Emission constraint %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if indic.target==1
% c(T+1:T*2) = F-Ftarget;
% end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%    Equality Constraints    %%%
 % include missing equations here %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ceq = [];

 ceq(1:T)       = S-(sff+sg+sn);  % LOM neutral technology 
 ceq(T+1:T*2)   = wsf-wsgtil./(1-taus); % wage scientists green/ taus
 ceq(T*2+1:T*3) = N-(An.*Ln).*(pn.*alphan).^(alphan./(1-alphan)); % from production function neutral good
 % optimality skills (for fossil used to determine wage rates)
 ceq(T*3+1:T*4) = thetan*Ln.*wln-wh.*hhn; % optimality labour good producers neutral high skills
 ceq(T*4+1:T*5) = thetag*Lg.*wlg-wh.*hhg; % optimality labour good producers green high
 ceq(T*5+1:T*6) = (1-thetan)*Ln.*wln-wl.*hln; % optimality labour good producers neutral low
 ceq(T*6+1:T*7) = (1-thetag)*Lg.*wlg-wl.*hlg; % optimality labour good producers green low
 ceq(T*7+1:T*8) = wsf-wsn; % wage scientists neutral
 ceq(T*8+1:T*9) = C-zh.*lambdaa.*(wh.*hh).^(1-taul)-(1-zh).*lambdaa.*(wl.*hl).^(1-taul)-SGov;
% add foc for one skill type (ratio respected in taul derivation) 
 ceq(T*9+1:T*10) = chii*hh.^(sigmaa+taul)-(muu.*lambdaa.*(1-taul).*(wh).^(1-taul));

if indic.target==1 % scientists are coded as a choice variable    
 ceq(T*10+1:T*11)   = sff - ((Af./Af_lag-1).*rhof^etaa/gammaa.*(Af_lag./A_lag).^phii).^(1/etaa);
 ceq(T*11+1:T*12)  = sg  - ((Ag./Ag_lag-1).*rhog^etaa/gammaa.*(Ag_lag./A_lag).^phii).^(1/etaa);
 ceq(T*12+1:T*13)  = sn  - ((An./An_lag-1).*rhon^etaa/gammaa.*(An_lag./A_lag).^phii).^(1/etaa);
end

ceq = ceq';
end