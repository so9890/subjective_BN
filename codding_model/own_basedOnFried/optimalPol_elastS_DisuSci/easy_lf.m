% simple model to compare social planner and 
% competitive equilibrium



%- model for numeric solution as function of pg and Lg
function f = easy_lf(x, params, list, pol,  init201519, indic)

    pg = exp(x(1));
    Lg = exp(x(2));

    % auxiliary

    read_in_params;
    read_in_pol;
    % exogenous variables: tauf, taul
    % endogenous variables: pg, pf, w, h, Lf, Lg, C, lambdaa, G, F, Y

    Af=(1+vf)*init201519(list.init=='Af0');
    Ag=(1+vg)*init201519(list.init=='Ag0');

    G = Ag*Lg;
    w = pg*Ag;
    pf = w/((1-tauf)*Af); % wage clearing 
    Lf = pg/pf*Ag/Af*eppsy/(1-eppsy)*Lg;
    F = Af*Lf;
%   Y = (F)^(eppsy)*(G)^(1-eppsy);
    hdem = Lf+Lg;
    if indic.taxsch==0
        lambdaa = (w*hdem+tauf*pf*F)/((w*hdem)^(1-taul));
        hsup =  (lambdaa^(1-thetaa)*(1-taul)*w^((1-taul)*(1-thetaa))/chii)^(1/(sigmaa+taul+thetaa*(1-taul)));
    %    C  = lambdaa*(w*hdem)^(1-taul);
    elseif indic.taxsch==1 % linear tax system with env tax revenues in gov budgte
        lambdaa = 0;
        hsup = (((w+tauf*pf*F/hdem)^(-thetaa)*w*(1-taul))/(chii))^(1/(sigmaa+thetaa));
    elseif indic.taxsch==2 % linear tax but no transfers; gov consumes
        lambdaa=0; 
        Gov = tauf*pf*F;
        hsup = (w^(1-thetaa)*(1-taul)/chii)^(1/(thetaa+sigmaa));
    elseif indic.taxsch==3 % non linear and no transfers
        lambdaa=w*h/(w*h)^(1-taul);
        Gov=tauf*pf*F;
        hsup=(w^(1-thetaa)*(1-taul)/chii)^(1/(thetaa+sigmaa));
    end
        
    % 
    f(1)= 1-(pf/eppsy)^eppsy*(pg/(1-eppsy))^(1-eppsy); % determines pg; could have alternatively 
    f(2)= hdem-hsup; % labour market clearing determines Lg; (goods market clears by walras law)
end