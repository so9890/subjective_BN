%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program: amotry.m
% By: Stephie Fried (adpoted from routine 10.4 in  Numerical Recipies in C 
% by Press, Teukolsky, Vetterling, and Flannery)
% Date: Spring 2017
% Purpose: "Extrapolates by a factor of fac through the face of the simplex
% across from the high point, tries it, and replaces the high point if the
% new point is better" ((Press et. al) pg 421)
% Inputs: 
% Simplex features: p, y, psum, ndim, ihi, fac
% Directly calibrated parameters: alphaf, alpham, deltaef, epse, epsf, epsy,
% L, phi ,rhonf, rhong, S
% Percent change in the foriegn oil price: deltaPfStarD
% Empirical ratio between foreign oil price and domestic fossil energy
% price on BGP: theta0
% BGP growth rate: n0
% Empirical values of the moments: momD
% Initial guesses: guess0, guessBGP
% fsolve options: options
% Number of time periods: Ttarg
% Maximum value of the simplex from minDist: yHigh
% Aggregate technology on the BGP: a0
% Outputs:
% Trial point: ptry
% Value of trial point: ytry
% New simplex: pNew
% Values of the new simplex: yNew
% New value of psum: psumNew
% Check that fsolve worked: solve
% New initial guess: guessBGP
% Calls program momModel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [ytry,ptry, solve, guessBgp, yNew, pNew, psumNew] = amotry(p,y, psum, ndim, ihi, fac, ...
    a0,alphaf, alpham, deltaef, epse, epsf, epsy, L,phi, rhonf, rhong, S, theta0, deltaPfStarD, ...
    n0, guess0, guessBgp, momD, options,p0, Ttarg,  yHigh)


fac1 = (1-fac)/ndim;
fac2 = fac1 - fac;
ptry = zeros(1, ndim);
%form the new point

for (j = 1:ndim)
    ptry(j) = psum(j)*fac1 - p(ihi, j)*fac2;
end
%Stay in bounds
ptry = bounds(ptry);
%evaluate the function at the trial point

[ytry, solve, guessBgp] =  momModel(ptry,a0, alphaf, alpham,deltaef, epse, epsf, epsy, L, phi, rhonf, rhong,...
      S, theta0, deltaPfStarD, n0, guess0, guessBgp, momD, options,p0, Ttarg,  yHigh);

%if ytry is better than the highest point, replace the hightest point
%with ytry

%in the case of the expansion, replace reflected point with the
%expanded point if the expanded point is better.

if ytry < y(ihi)
    y(ihi) = ytry;
    p(ihi, :) = ptry;
end
psumNew = sum(p,1);
yNew = y;
pNew =p;




