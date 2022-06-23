function f = easy_opt(x, params, list,  init201519, indic)
% solves for optimal allocation by choosing s and h, where s=Lf/h 

  read_in_params;
    % exogenous variables: tauf, taul
    % endogenous variables: pg, pf, w, h, Lf, Lg, C, lambdaa, G, F, Y

    Af=(1+vf)*init201519(list.init=='Af0');
    Ag=(1+vg)*init201519(list.init=='Ag0');

    
if indic.notaul==0
    s = exp(x(1));
    h = exp(x(2));
    
    % auxiliary

  
    tauf = 1-((1-eppsy)/eppsy)*s/(1-s); % tauf determines s
    pg = eppsy^eppsy*(1-eppsy)^(1-eppsy)*((1-tauf)*Af/Ag)^eppsy;
    w = pg*Ag;
    pf = w/((1-tauf)*Af); 
    % labour market clearing: 
    Lg = (1-s)*h; 
    Lf = s*h;
    % production
    G = Ag*Lg;
    F = Af*Lf;
    % income tax/ gov budget
    if indic.taxsch==0
        taul= 1-h^(thetaa+sigmaa)*chii*(w+tauf*pf*Af*s)^(thetaa-1);
        lambdaa = (w*h+tauf*pf*F)/((w*h)^(1-taul));
    elseif indic.taxsch==1
        taul = 1-h^(thetaa+sigmaa)*chii*(w+tauf*pf*Af*s)^(thetaa)/w;
        T = w*taul*h+pf*tauf*F;
    end
    % good market clearing and final production 
    C = (F)^(eppsy)*(G)^(1-eppsy); 
   
    % derivatives
    dCdh = (Af*s)^(eppsy)*(Ag*(1-s))^(1-eppsy);
    dCds = dCdh*h*(eppsy*(1-s)-s*(1-eppsy))/(s*(1-s));
    Uc = C^(-thetaa);
    Uh = -chii*h^sigmaa;
    Uf = -weightext*extexpp*(omegaa)^extexpp*F.^(extexpp-1);
    dFdh = Af*s;
    dFds = Af*h;
    
    % Optimality conditions planner 
    f(1) = Uc*dCdh + Uh +indic.extern* Uf*dFdh; % optimality wrt h => taul
    f(2) = Uc*dCds + indic.extern*Uf*dFdh; % optimality government wrt s => tauf

elseif indic.notaul==1
    
    s = exp(x(1));

    tauf = 1-((1-eppsy)/eppsy)*s/(1-s); % tauf determines s
    pg = eppsy^eppsy*(1-eppsy)^(1-eppsy)*((1-tauf)*Af/Ag)^eppsy;
    w = pg*Ag;
    pf = w/((1-tauf)*Af); 
    % labour supply
    if indic.taxsch==0
        h = ((w+tauf*pf*Af*s)^(1-thetaa)/chii)^(1/(sigmaa+thetaa));
    else
        % with linear tax: 
        h = ((w+tauf*pf*Af*s)^(-thetaa)*w/chii)^(1/(sigmaa+thetaa));
    end
    % labour market clearing: 
    Lg = (1-s)*h; 
    Lf = s*h;
    % production
    G = Ag*Lg;
    F = Af*Lf;
    % gov budget: without income tax prog lambdaa s.t. hh consume all
    lambdaa = (w*h+tauf*pf*F)/((w*h));
    % good market clearing and final production 
    C = (F)^(eppsy)*(G)^(1-eppsy); 
   
    % derivatives
    dCdh = (Af*s)^(eppsy)*(Ag*(1-s))^(1-eppsy);
    dCds = dCdh*h*(eppsy*(1-s)-s*(1-eppsy))/(s*(1-s));
    Uc = C^(-thetaa);
    Uh = -chii*h^sigmaa;
    Uf = -weightext*extexpp*(omegaa)^extexpp*F.^(extexpp-1);
    dFdh = Af*s;
    dFds = Af*h;
    
    % Optimality conditions planner 
    %f(1) = Uc*dCdh + Uh +indic.extern* Uf*dFdh; % optimality wrt h => taul
    f(1) = Uc*dCds + indic.extern*Uf*dFdh; % optimality government wrt s => tauf

end
end