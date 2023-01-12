function [c, ceq] = constraints_flexetaa(y, T, params, init, list, Ems, indic, MOM, percon, taulFixed)
% function to read in constraints on government problem

% pars
read_in_params;

% Ftarg_20s=(MOM.US_Budget20_30+3*deltaa)./omegaa; 

% transform x: all are exponentially transformed
 x=exp(y);

% except for hours
if indic.noskill==0 %version with skill
    x((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T) = upbarH./(1+exp(y((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)));
    x((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T) = upbarH./(1+exp(y((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)));
else
    x((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T) = upbarH./(1+exp(y((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T)));
end

% if indic.sep==1
    x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = upbarH./(1+exp(y((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)));
    x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = upbarH./(1+exp(y((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)));
    x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = upbarH./(1+exp(y((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)));
if indic.sep==0
%      x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = (y((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)).^2;
%      x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = (y((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)).^2;
%      x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = (y((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)).^2;
     x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T) = upbarS./(1+exp(y((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T)));
end


if indic.target==1
    Ftarget=(Ems'+deltaa)./omegaa; 
    x((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)   = Ftarget./(1+exp(y((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)));
end


%- auxiliary variables

if indic.noskill==0
      [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, Ch, Cl, muuh, muul, hl, hh, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, wsf, wsg, wsn, ws,  tauf, taul,taus, lambdaa,...
            wln, wlg, wlf, SWF, S, GovCon, Tls, Tlsall, PV,PVSWF, objF]= OPT_aux_vars_notaus_flex_newTauf(x, list, params, T, init, indic, MOM, taulFixed);
else

     [xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, h, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee,  ws, wsf, wsn, wsg,  tauf, taul, lambdaa,...
            w, SWF, S, GovCon, Tls, PV,PVSWF, objF]= OPT_aux_vars_notaus_skillHom(x, list, params, T, init, indic, MOM);

end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%    Inequality Constraints    %%%
 % only for direct periods
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % time period specific constraints! 
c = []; %  periods and 2 additional ones: 
% c(1)= sum(F(1:percon))-Ftarg_20s; % when percon=0 then this constraint is not active
c(1:T) = -GovCon; % government budget cannot be negative
c(T+1:T*2) = -Tls; % lump-sum transfers cannot be negative

if indic.xgrowth==0
    if indic.sep==0
        %c(T*2+1:T*3)   =S- upbarS;
    else
        c(T*2+1:T*3)   = sg-upbarH;
        c(3*T+1:4*T)   = sff-upbarH;
        c(4*T+1:5*T)   = sn-upbarH;
    end
end

if indic.notaul>=7
    if indic.sep==0 
        c(T*2+1:T*3)   =-tauf; % tauf has to be positive
    elseif indic.sep~=0 && indic.xgrowth== 0 
        c(T*5+1:T*6)   =-tauf; % tauf has to be positive
    end
end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%    Equality Constraints    %%%
 % include missing equations here %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 ceq = [];
 if indic.noskill==0
     if indic.xgrowth==0

         ceq(1:T)       = chii*hh.^(sigmaa+taul)-(muu.*lambdaa.*(1-taul).*(wh).^(1-taul)); % labor supply
         ceq(T*1+1:T*2) = N-(An.*Ln).*(pn.*alphan).^(alphan./(1-alphan)); % from production function neutral good
         % optimality skills (for fossil used to determine wage rates)
         ceq(T*2+1:T*3) = thetan*Ln.*wln-wh.*hhn; % optimality labour good producers neutral high skills
         ceq(T*3+1:T*4) = thetag*Lg.*wlg-wh.*hhg; % optimality labour good producers green high
         ceq(T*4+1:T*5) = (1-thetan)*Ln.*wln-wl.*hln; % optimality labour good producers neutral low
         ceq(T*5+1:T*6) = (1-thetag)*Lg.*wlg-wl.*hlg; % optimality labour good producers green low
           if indic.sep==0
               if indic.Sun~=2
                  ceq(T*6+1:T*7) = ws-wsf; % wage clearing
               elseif indic.Sun==2
                  ceq(T*6+1:T*7) = chiis.*S.^(sigmaas+taul)-muu.*lambdaa.*(1-taul).*ws.^(1-taul) ; % when scientists are taxed than include optimal supply scientists here
               end
            if indic.notaul<7 % otherwise used to determine taus
                ceq(T*7+1:T*8) = ws-wsg;
             else
                 ceq(T*7+1:T*8) =  zeros(size(T*7+1:T*8));
             end
            ceq(T*8+1:T*9) = ws-wsn;
            
            ceq(T*9+1:T*10)= zh*hh-(hhf+hhg+hhn);
         ceq(T*10+1:T*11)= (1-zh)*hl-(hlf+hlg + hln );
          
         % hh budget different for with and without inequality version
         if indic.Sun==2
            ceq(T*11+1:T*12)   = C-zh.*lambdaa.*(wh.*hh).^(1-taul)-(1-zh).*lambdaa.*(wl.*hl).^(1-taul)-lambdaa.*(ws.*S).^(1-taul)-Tlsall;
         else
            ceq(T*11+1:T*12)   = C-zh.*lambdaa.*(wh.*hh).^(1-taul)-(1-zh).*lambdaa.*(wl.*hl).^(1-taul)-ws.*S-Tlsall;
         end
            ceq(T*12+1:T*13) = S-sn-sff-sg;
            
            if indic.notaul==1 || indic.notaul == 2 ||  indic.notaul == 5 || indic.notaul==8 || indic.notaul==9% when no taul is available
                ceq(T*13+1:T*14) = chii*hl.^(sigmaa+taul)-(muu.*lambdaa.*(1-taul).*(wl).^(1-taul));
            end
         
           elseif indic.sep>=1
             ceq(T*6+1:T*7) = (chiis)*sff.^sigmaas-muu.*wsf; % scientist hours supply
             ceq(T*7+1:T*8) = (chiis)*sg.^sigmaas-muu.*wsg;
             ceq(T*8+1:T*9) = (chiis)*sn.^sigmaas-muu.*wsn;
             ceq(T*9+1:T*10)= zh*hh-(hhf+hhg+hhn);
            ceq(T*10+1:T*11)= (1-zh)*hl-(hlf+hlg + hln );
         
         % hh budget different for with and without inequality version
         if indic.Sun==2
            ceq(T*11+1:T*12)   = C-zh.*lambdaa.*(wh.*hh).^(1-taul)-(1-zh).*lambdaa.*(wl.*hl).^(1-taul)-lambdaa.*(ws.*S).^(1-taul)-Tlsall;
         else
            ceq(T*11+1:T*12)   = C-zh.*lambdaa.*(wh.*hh).^(1-taul)-(1-zh).*lambdaa.*(wl.*hl).^(1-taul)-ws.*S-Tlsall;
         end
            if indic.notaul==1 || indic.notaul == 2 ||  indic.notaul == 5 || indic.notaul==8 || indic.notaul==9 % when no taul is available
                ceq(T*12+1:T*13) = chii*hl.^(sigmaa+taul)-(muu.*lambdaa.*(1-taul).*(wl).^(1-taul));
            end
         
           end
           
         
     elseif indic.xgrowth==1
         ceq(1:T)       = chii*hh.^(sigmaa+taul)-(muu.*lambdaa.*(1-taul).*(wh).^(1-taul));
         ceq(T*1+1:T*2) = N-(An.*Ln).*(pn.*alphan).^(alphan./(1-alphan)); % from production function neutral good
         % optimality skills (for fossil used to determine wage rates)
         ceq(T*2+1:T*3) = thetan*Ln.*wln-wh.*hhn; % optimality labour good producers neutral high skills
         ceq(T*3+1:T*4) = thetag*Lg.*wlg-wh.*hhg; % optimality labour good producers green high
         ceq(T*4+1:T*5) = (1-thetan)*Ln.*wln-wl.*hln; % optimality labour good producers neutral low
         ceq(T*5+1:T*6) = (1-thetag)*Lg.*wlg-wl.*hlg; % optimality labour good producers green low
         ceq(T*6+1:T*7) = zh*hh-(hhf+hhg+hhn);
         ceq(T*7+1:T*8) = (1-zh)*hl-(hlf+hlg + hln );
         ceq(T*8+1:T*9) = C-zh.*lambdaa.*(wh.*hh).^(1-taul)-(1-zh).*lambdaa.*(wl.*hl).^(1-taul)-Tls;
        
         if indic.notaul==1 || indic.notaul == 2 ||  indic.notaul == 5 || indic.notaul==8 || indic.notaul==9 % when no taul is available
                ceq(T*9+1:T*10) = chii*hl.^(sigmaa+taul)-(muu.*lambdaa.*(1-taul).*(wl).^(1-taul));
         end
     end
       
 elseif indic.noskill==1
     if indic.xgrowth==0
            ceq(T*0+1:T*1) = Ln -pn.*(1-alphan).*N./w; % labour demand by sector 
            ceq(T*1+1:T*2) = Lg - pg.*(1-alphag).*G./w;
            ceq(T*2+1:T*3) = Lf- h./(1+Ln./Lf+Lg./Lf); % labour market clearing 
            if indic.Sun~=2
                ceq(T*3+1:T*4) = C-lambdaa.*(w.*h).^(1-taul)-ws.*S-Tlsall;
            else
                ceq(T*3+1:T*4) = C-lambdaa.*(w.*h).^(1-taul)-lambdaa.*(ws.*S).^(1-taul)-Tlsall;
            end
              if indic.sep==0
                ceq(T*4+1:T*5) = ws-wsf; % wage clearing
                if indic.notaul>7
                    ceq(T*5+1:T*6) = ws-wsg;
                end
                ceq(T*6+1:T*7) = ws-wsn;
                ceq(T*7+1:T*8) = S-sn-sg-sff;
                if indic.notaul==1 || indic.notaul == 2 ||  indic.notaul == 5||  indic.notaul == 8 || indic.notaul==9% when no taul is available
                ceq(T*8+1:T*9)= chii*h.^(sigmaa+taul)-(muu.*lambdaa.*(1-taul).*(w).^(1-taul));
                end
              else
                 ceq(T*4+1:T*5) = (chiis)*sff.^sigmaas-muu.*wsf; % scientist hours supply
                 ceq(T*5+1:T*6) = (chiis)*sg.^sigmaas-muu.*wsg;
                 ceq(T*6+1:T*7) = (chiis)*sn.^sigmaas-muu.*wsn;
                if indic.notaul==1 || indic.notaul == 2 ||  indic.notaul == 5 ||  indic.notaul == 8 || indic.notaul==9% when no taul is available
                ceq(T*7+1:T*8)= chii*h.^(sigmaa+taul)-(muu.*lambdaa.*(1-taul).*(w).^(1-taul));
                end
               end

     elseif indic.xgrowth==1
            ceq(T*0+1:T*1) = Ln -pn.*(1-alphan).*N./w; % labour demand by sector 
            ceq(T*1+1:T*2) = Lg - pg.*(1-alphag).*G./w;
            ceq(T*2+1:T*3) = Lf- h./(1+Ln./Lf+Lg./Lf); % labour market clearing 
            ceq(T*3+1:T*4) = C-lambdaa.*(w.*h).^(1-taul)-Tls;

            if indic.notaul==1 || indic.notaul == 2 ||  indic.notaul == 5 ||  indic.notaul == 8 || indic.notaul==9 % when no taul is available
                ceq(T*4+1:T*5)= chii*h.^(sigmaa+taul)-(muu.*lambdaa.*(1-taul).*(w).^(1-taul));
            end

     end
 end

 %
ceq = ceq';
end