% simple model to compare social planner and 
% competitive equilibrium



%- model for numeric solution as function of pg and Lg
function f = easy_opt(x, params, list, pol,  init201519, indic)

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
    pf = Ag/((1+tauf)*Af);
    Lf = pg/pf*Ag/Af*eppsy/(1-eppsy)*Lg;
    F = Af*Lf;
    Y = (F)^(eppsy)*(G)^(1-eppsy);
    hdem = Lf+Lg;
    lambdaa = (w*hdem+tauf*pf*F)/((w*hdem)^(1-tauf));
    hsup =  (lambdaa^(1-thetaa)*(1-taul)*w^((1-taul)*(1-thetaa))/chii)^(1/(sigmaa+taul+thetaa*(1-taul)));

    % 
    f(1)= 1- (eppsy*pf+(1-eppsy)*pg); % determines pg
    f(2)= hdem-hsup ; % labour market clearing determines Lg; (goods market clears by walras law)

end