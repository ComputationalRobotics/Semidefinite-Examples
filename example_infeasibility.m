clc; clear; close all;
restoredefaultpath;

addpath(genpath('../SOSTOOLS'))
addpath(genpath('../../../mosek'))

x = mpvar('x',2,1); % dummy x 
h = [x'*x - 1];
g = [3*x(2)-x(1)^3-2;
     x(1)-8*x(2)^3];
g = [g; g(1)*g(2)];

prog = sosprogram(x);

kappa = 4;

lhs = 0;
lams = {};
for i = 1:length(h)
    hi = h(i);
    deg_hi = hi.maxdeg;
    deg_lam = 2*kappa - deg_hi;
    [prog,lam] = sospolyvar(prog,monomials(x,0:deg_lam));
    lhs = lhs + lam*hi;
    lams{end+1} = lam;
end
sigs = {};
for i = 1:length(g)
    gi = g(i);
    deg_gi = gi.maxdeg;
    if deg_gi <= 2*kappa
        deg_sig = ceil((2*kappa - deg_gi)/2);
        [prog,sig] = sossosvar(prog,monomials(x,0:deg_sig));
        lhs = lhs + sig*gi;
        sigs{end+1} = sig;
    end
end

prog = soseq(prog,lhs+1);
options.solver = 'mosek';
prog = sossolve(prog,options);

lams_val = {};
for i = 1:length(lams)
    lams_val{end+1} = cleanpoly(sosgetsol(prog,lams{i}),1e-6);
end
sigs_val = {};
for i = 1:length(sigs)
    sigs_val{end+1} = cleanpoly(sosgetsol(prog,sigs{i}),1e-6);
end
