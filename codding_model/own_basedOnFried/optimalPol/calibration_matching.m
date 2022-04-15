%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: calibrate.m
% By: Stephie Fried (adpoted from routine 10.4 in  Numerical Recipies in C 
% by Press, Teukolsky, Vetterling, and Flannery))
% Date: Spring 2017
% Purpose: Chooses parameter values to minimize distance between empirical
% and model values of the moments. Uses Nelder-Meade algorithm  for
% minimization. Same as minDist.m, but wrapped in a function. 
% Inputs: 
% Empirical values of moments: momD
% Initial guess for parameters: p0
% Directly calibrated parameters and model and data values: alphaf, alpham, deltaef,
%    epse, epsf, epsy, L,phi ,rhonf, rhong,S, theta0, deltaPfStarD
% Initial guesses: guess0, guessBGP
% Options for fsolve: options
% Number of periods: Ttarg
% Bootstrap repitation: rep
% Outputs:
% Sum squared distance betwen model and empircal values of moments: dist
% Vector of calibrated parameter values: params
% Calls programs: momModel.m, amotry.m, bounds.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dist, params] = calibrate(momD, p0, alphaf, alpham, deltaef, epse, epsf, epsy, L,phi ,rhonf, rhong,...
    S, theta0, deltaPfStarD, guess0, guessBgp, options, Ttarg, rep)

n0 = (1.02)^5 -1; %5 year grwoth rate

    
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iniatilize Nelder-Meade Routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%COMPUTATIONAL PARAMETERS
tiny = 1e-10;
ftol  = 1e-10;

%DEFINE INITIAL GUESS1
pp=p0; 
ndim = length(pp); %number of parameters
step = .05*pp; %step size is 5%
p = [pp; repmat(pp, ndim, 1)+diag(step)];
%Each row is a vertex of the starting simplex
counter =0;
yHigh = 1e8;
y = zeros(length(p0) +1, 1);
a0=1;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nelder-Meade Routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(p0) +1
    i
   if i >1
        yHigh = max(y(1:i-1)); 
   end 
    
   if i ==5 
       y(i) = y(1);
   else
       
       p2 = bounds(p(i, :));
       [y(i), solve, guessBgp] =  momModel(p2,a0,alphaf, alpham, deltaef, epse, epsf, epsy, L,phi ,rhonf, rhong,...
    S, theta0, deltaPfStarD, n0, guess0, guessBgp, momD, options,p0, Ttarg,  yHigh);
      
      
   end

    if i ==1
        guess0 = guessBgp;
    end
end

done =0;
psum = sum(p, 1);

while(done ==0)
    %FIND THE HIGHEST, SECOND HIGHEST AND LOWEST POINTS ON THE SIMPLEX
    ihi = find(y == max(y), 1, 'first'); %highest 
    y1 = y(y ~= max(y));
    inhi = find(y == max(y1), 1, 'first'); %second highest
    ilo = find(y == min(y), 1, 'first');  %lowest
    
    pTest = bounds(p(ilo,:)); 
    if abs(max(p(ilo,:) - pTest)) > 1e-12
        disp('somethings wrong');
    end 
 
    rtol = 2*(y(ihi) - y(ilo))/(y(ihi) + y(ilo) + tiny); %stopping criteria
    if rtol < ftol
        done=1; 
    end 
    
    %BEGIN A NEW ITERATION
    %reflect the simplex from ihi
    
    fac =-1;
    [ytry, ptry, solve, guessBgp, y, p, psum] = amotry(p,y, psum, ndim, ihi, fac, ...
        a0,alphaf, alpham, deltaef, epse, epsf, epsy, L,phi ,rhonf, rhong,  S, theta0,...
        deltaPfStarD, n0, guess0, guessBgp, momD, options,p0, Ttarg,  y(ihi));
    
    if solve > 1e-12
        solve
        disp('did not solve');
    end
    if ytry == 1e10 disp('error in refletion'); end
    
    %if ytry < y(ihi), replace y(ihi) with ytry, adjust psum accordingly
    
    %if ytry is better than the best point, try an expolation by a factor of
    %two about ihi (which is the new point).  (expansion)
    
    if ytry <= y(ilo)
        %  disp('Expansion!');
        fac =2;
        [ytry,ptry,  solve, guessBgp, y,p, psum] = amotry(p,y, psum, ndim, ihi, fac, ...
            a0,alphaf, alpham, deltaef, epse,  epsf, epsy, L, phi ,rhonf, rhong,  S, theta0,...
        deltaPfStarD, n0, guess0, guessBgp, momD, options,p0, Ttarg,  y(ihi));
    
        if ytry == 1e10 disp('error in expansion'); end
        if solve > 1e-12
            solve
            disp('did not solve');
        end
        
        %if ytry is worse than the second highest point,do a contraction
        %(and it is not because ytry would not solve)
    elseif ytry >= y(inhi)
        % disp('1D Contraction')
        fac = .5;
        ysave = y(ihi);
        [ytry, ptry, solve, guessBgp, y, p, psum] =  amotry(p,y, psum, ndim, ihi, fac, ...
            a0,alphaf, alpham, deltaef, epse, epsf, epsy, L, phi, rhonf, rhong, S, theta0,...
        deltaPfStarD, n0, guess0, guessBgp, momD, options,p0, Ttarg,  y(ihi));
        if ytry == 1e8 disp('error in 1D contraction'); end
        if solve > 1e-12
            solve
            disp('did not solve');
        end
        %If the new point is worse than ysave, contract around the lowest
        
        if ytry >= ysave
            disp('Many D contraction')
            for i =1:ndim+1
                %if (i ~= ilo)
                for j =1:ndim
                    p(i,j) = .5*(p(i,j) + p(ilo, j));
                end
                p(i,:) = bounds(p(i, :));
                yHigh = max(y);
                [y(i), solve, guessBgp] = momModel(p(i, :),a0,alphaf, alpham, deltaef, epse, epsf, epsy, L, phi , rhonf, rhong,...
                    S, theta0,...
                    deltaPfStarD, n0, guess0, guessBgp, momD, options,p0, Ttarg,  yHigh);
                
                if solve > 1e-12
                    solve
                    disp('did not solve');
                end
                if y(i) == 1e10 disp('error in multiD Contraction'); end
            end
            psum = sum(p,1);
        end
    end

    [ytry, y(ilo)];
    [p(ilo,1), p(ilo, 2)];
    [p(ilo,4), p(ilo, 6)];
    [p(ilo,3)];

  [rep, rtol]
    counter = counter+1;  
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate squaed distance between model and data values of moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

pp =bounds(p(ilo, :)); %vector of chosen parameter values
% 
 alphag = pp(1);  eta = pp(2); deltaye = pp(3); 
 epsf = pp(4); deltaff = pp(5); nu1 = pp(6);


 pp = [alphag, eta, deltaye, epsf, deltaff, nu1]; %redundant if ran previous cell

[dist, solve, guessBgp,  momM, gamma]  =  momModel(pp,a0,alphaf, alpham, deltaef, epse,epsf, epsy, L,phi ,rhonf, rhong,...
    S, theta0, deltaPfStarD, n0, guess0, guessBgp, momD, options,p0, Ttarg,  100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store parameter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

 params = [alphag, eta, deltaye, epsf, deltaff, nu1, gamma]; 





     