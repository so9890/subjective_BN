function [symms, Obj_ramPA ]=primal_problem_dynamic(y, x, list, symms, E)

% periods for which to solve directly
P=30; 
%% generate vectors for each variable; save in matrix
list.joint=[list.y, list.x, list.ramsey_mu];
vecs=sym('a',[length([list.joint]),P]);
    for s = [list.joint]  % loop over list entries
        vecs([list.joint]==s,:)=sym(sprintf('%s%d',s),  [1,P]); 
    end

%% replace indexed variables in static problem
% note that script "primal_problem" reads in the sum of swf and constraints
% for one period => replace variables with dynamic ones to get sum of welfare functions
% and add up to get dynamic problem
% remains to connect variables=> use growth rates; should be linked automatically through  

% can now take derivatives if static, or do value function iteration if constant tax rate!
% if dynamic policy => look at sum over next 30 years

% vector of objective function by period, already substituted technology,
% important to variables linked across periods: which they are via A

vec_Obj_ram=sym('a',[P,1]);

for i=1:P
vec_Obj_ram(i)=subs(Obj_ramPA,[]);

vec_discount=[betaa betaa^2 betaa^3 betaa^4 betaa^5 betaa^6 betaa^7 betaa^8 betaa^9 betaa^10 betaa^11 betaa^12 betaa^13 betaa^14 betaa^15 betaa^16 betaa^17 betaa^18 betaa^19 betaa^20 betaa^21 betaa^22 betaa^23 betaa^24 betaa^25 betaa^26 betaa^27 betaa^28 betaa^29 betaa^30];
Obj_ramPA_dynamic= vec_discount*vec_SWF'
