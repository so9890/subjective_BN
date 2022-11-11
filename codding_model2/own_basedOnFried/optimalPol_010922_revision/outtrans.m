function [opta] =outtrans(x, list, T, Ftarget, indic, params , symms,init201519,MOM, percon)

read_in_params;
out_trans=exp(x);
            if indic.noskill==0
                out_trans((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)=upbarH./(1+exp(x((find(list.opt=='hl')-1)*T+1:find(list.opt=='hl')*T)));
                out_trans((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)=upbarH./(1+exp(x((find(list.opt=='hh')-1)*T+1:find(list.opt=='hh')*T)));
            else
                out_trans((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T)=upbarH./(1+exp(x((find(list.opt=='h')-1)*T+1:find(list.opt=='h')*T)));
            end

            if indic.sep==1
                out_trans((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = upbarH./(1+exp(x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)));
                out_trans((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = upbarH./(1+exp(x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)));
                out_trans((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = upbarH./(1+exp(x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)));
            elseif indic.sep==0
                out_trans((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T) = (x((find(list.opt=='sn')-1)*T+1:find(list.opt=='sn')*T)).^2;
                out_trans((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T) = (x((find(list.opt=='sg')-1)*T+1:find(list.opt=='sg')*T)).^2;
                out_trans((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T) = (x((find(list.opt=='sff')-1)*T+1:find(list.opt=='sff')*T)).^2;
                out_trans((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T) = upbarH./(1+exp(x((find(list.opt=='S')-1)*T+1:find(list.opt=='S')*T)));
            end    
        if indic.target==1
            out_trans((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)=Ftarget./(1+exp(x((find(list.opt=='F')-1)*T+1+percon:find(list.opt=='F')*T)));
        end
        
            [hhf, hhg, hhn, hlg, hlf, hln, xn,xf,xg,Ag, An, Af,...
            Lg, Ln, Lf, Af_lag, An_lag, Ag_lag,sff, sn, sg,  ...
            F, N, G, E, Y, C, Ch, Cl, muuh, muul, hl, hh, A_lag, SGov, Emnet, A,muu,...
            pn, pg, pf, pee, wh, wl, wsf, wsg, wsn, ws,  tauf, taul, lambdaa,...
            wln, wlg, wlf, SWF, S, GovCon, Tls, PV,PVSWF, objF]...
            = OPT_aux_vars_notaus_flex_newTauf(out_trans, list, params, T, init201519, indic, MOM);
            taus = zeros(size(pn));
            gammall = zeros(size(pn));
            gammalh = zeros(size(pn));
            gammasn = zeros(size(pn));

            opta=eval(symms.allvars);
%         cell_par=arrayfun(@char, symms.allvars, 'uniform', 0);
%         SL=cell2struct(num2cell(eval(symms.allvars)), cell_par, 2);
end