clc; clear; close all;
restoredefaultpath;

addpath(genpath('../SOSTOOLS'))
addpath(genpath('../../../mosek'))

x = mpvar('x',2,1); % dummy x 

h = [x(1)^2-1;
     x(2)^2-1];
g = [monomials(x,0);
     x(1)];
p = 4*x(1)^3 + 6*x(2);
prog = sosprogram(x);

gam = dpvar('gam');
prog = sosdecvar(prog,gam);

kappa = 2;

lhs = p - gam;
lams = {};
for i = 1:length(h)
    hi = h(i);
    deg_hi = hi.maxdeg;
    deg_lam = 2*kappa - deg_hi;
    [prog,lam] = sospolyvar(prog,monomials(x,0:deg_lam));
    lhs = lhs - lam*hi;
    lams{end+1} = lam;
end
sigs = {};
for i = 1:length(g)
    gi = g(i);
    deg_gi = gi.maxdeg;
    if deg_gi <= 2*kappa
        deg_sig = floor((2*kappa - deg_gi)/2);
        [prog,sig] = sossosvar(prog,monomials(x,0:deg_sig));
        lhs = lhs - sig*gi;
        sigs{end+1} = sig;
    end
end

prog = soseq(prog,lhs);
prog = sossetobj(prog,-gam);
options.solver = 'mosek';
prog = sossolve(prog,options);
gam_val = sosgetsol(prog,gam);



