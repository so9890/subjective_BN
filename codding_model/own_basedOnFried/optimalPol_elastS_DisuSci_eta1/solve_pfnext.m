function pfT1=solve_pfnext(x0LF, laggs, indic, indexx)

% solves competitive equilibrium with constant policy and


if indic.noskill==1
    x0=x0(list.choice~='hl'&list.choice~='hh'& list.choice~='hln'& list.choice~='hlf'&list.choice~='hlg'& list.choice~='gammall'&...
                list.choice~='hhn'& list.choice~='hhg'& list.choice~='hhf'& list.choice~='wl'& list.choice~='wh');
    syms h lambdaa w Lf Lg Ln real
    symms.choice= [symms.choice(list.choice~='hl'&list.choice~='hh'&list.choice~='hhn'&list.choice~='hhg'...
                &list.choice~='hhf'&list.choice~='hln'&list.choice~='hlg'&list.choice~='hlf'&list.choice~='wh'...
                &list.choice~='wl'&list.choice~='gammall'), h, w, Lf, Lg, Ln];
    list.choice=string(symms.choice);
    x0=[x0,Sall.hh, Sall.wh, Sall.Lf,Sall.Lg,Sall.Ln]; % order has to match how variables are added to symms.choice!

        
    indexxLF.lab = boolean(zeros(size(list.choice)));
    indexxLF.exp = boolean(zeros(size(list.choice)));
    indexxLF.sqr = boolean(zeros(size(list.choice)));
    indexxLF.oneab = boolean(zeros(size(list.choice)));

    indexxLF.lab( list.choice=='h' | list.choice=='S')=1;
    indexxLF.exp(list.choice~='h'&list.choice~='S'& list.choice~='gammalh'& list.choice~='gammas' )=1;
    indexxLF.sqr(list.choice=='gammalh'| list.choice=='gammas' )=1;
    indexx('LF_noskill')=indexxLF;
    
    if indic.noskill==0
        guess_trans=trans_guess(indexx('LF'), x0, params, list.params);
    else
        guess_trans=trans_guess(indexx('LF_noskill'), x0, params, list.params);
    end
    
    modFF = @(x)laissez_faire_nows(x, params, list, pol, laggs, indic);

    options = optimoptions('fsolve', 'TolFun', 10e-12, 'MaxFunEvals',8e3, 'MaxIter', 3e5,  'Algorithm', 'levenberg-marquardt');%, );%, );%, 'Display', 'Iter', );
    [sol2, fval, exitf] = fsolve(modFF, sol3, options);

end
