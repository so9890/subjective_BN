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
    if indic.taxsch==0 % non-linear and no lump sum redistribution
        taul= 1-h^(thetaa+sigmaa)*chii*(w+tauf*pf*Af*s)^(thetaa-1);
        lambdaa = (w*h+tauf*pf*F)/((w*h)^(1-taul));
    elseif indic.taxsch==1
        taul = 1-h^(thetaa+sigmaa)*chii*(w+tauf*pf*Af*s)^(thetaa)/w;
        T  = w*taul*h+pf*tauf*F;
    elseif indic.taxsch==2 %linear but no transfers
        taul = 1-(h^(thetaa+sigmaa)*chii)/(w)^(1-thetaa);
        Gov = tauf*pf*F;
    end
    % good market clearing and final production 
    Y = (F)^(eppsy)*(G)^(1-eppsy); 
   
    if indic.taxsch<=1
        C=Y;
    else
        C=Y-Gov;
    end
    
    % derivatives
    dFdh = Af*s;
    dFds = Af*h;
    dCdh = (Af*s)^(eppsy)*(Ag*(1-s))^(1-eppsy);
    dCds = dCdh*h*(eppsy/s-(1-eppsy)/(1-s));
    
    if indic.taxsch>1
        dCdh =dCdh-tauf*pf*dFdh;
        dtaufds = -(1-eppsy)/eppsy*(1/(1-s)^2);
        dpfds =-pf*(eppsy-1)/(1-tauf)*dtaufds;
        dGovds = pf*F*dtaufds+tauf*F*dpfds+tauf*pf*dFds;
        dCds = dCds-dGovds;
    end
    
    Uc = C^(-thetaa);
    Uh = -chii*h^sigmaa;
    Uf = -weightext*extexpp*(omegaa)^extexpp*F.^(extexpp-1);
    
    % Optimality conditions planner 
    f(1) = Uc*dCdh + Uh +indic.extern* Uf*dFdh; % optimality wrt h => taul
    f(2) = Uc*dCds + indic.extern*Uf*dFds; % optimality government wrt s => tauf

elseif indic.notaul==1
    s = exp(x(1));

    tauf = 1-((1-eppsy)/eppsy)*s/(1-s); % tauf determines s
    pg = eppsy^eppsy*(1-eppsy)^(1-eppsy)*((1-tauf)*Af/Ag)^eppsy;
    w = pg*Ag;
    pf = w/((1-tauf)*Af); 
    % labour supply
    if indic.taxsch==0
        h = ((w+tauf*pf*Af*s)^(1-thetaa)/chii)^(1/(sigmaa+thetaa));
    elseif indic.taxsch==1
        % with linear tax: 
        h = ((w+tauf*pf*Af*s)^(-thetaa)*w/chii)^(1/(sigmaa+thetaa));
    elseif indic.taxsch == 2
        h = (w^(1-thetaa)/chii)^(1/(sigmaa+thetaa));
        Gov = tauf*pf*F;
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
    Y = (F)^(eppsy)*(G)^(1-eppsy); 
    
    if indic.taxsch<=1
        C=Y;
    else
        C=Y-Gov;
    end
    
    % derivatives
    dFdh = Af*s;
    dFds = Af*h;
    dCdh = (Af*s)^(eppsy)*(Ag*(1-s))^(1-eppsy);
    dCds = dCdh*h*(eppsy*(1-s)-s*(1-eppsy))/(s*(1-s));
    
    if indic.taxsch>1
%         dCdh =dCdh-tauf*pf*dFdh;
        dtaufds = -(1-eppsy)/eppsy*(1/(1-s)^2);
        dpfds =-pf*(eppsy-1)/(1-tauf)*dtaufds;
        dGovds = pf*F*dtaufds+tauf*F*dpfds+tauf*pf*dFds;
        dCds = dCds-dGovds;
    end
    
        
    Uc = C^(-thetaa);
%     Uh = -chii*h^sigmaa;
    Uf = -weightext*extexpp*(omegaa)^extexpp*F.^(extexpp-1);
    
    % Optimality conditions planner 
    %f(1) = Uc*dCdh + Uh +indic.extern* Uf*dFdh; % optimality wrt h => taul
    f(1) = Uc*dCds + indic.extern*Uf*dFds; % optimality government wrt s => tauf

end
end