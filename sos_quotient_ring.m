clc; clear; close all;
restoredefaultpath;

addpath(genpath('../SOSTOOLS'))
addpath(genpath('../../../mosek'))

x = mpvar('x',2,1); % dummy x 
p = 10 - x(1)^2 - x(2);
h = x(1)^2 + x(2)^2 - 1;
prog = sosprogram(x);
kappa = 2; % bound of degree
[prog, sig0] = sossosvar(prog,monomials(x,0:kappa));
[prog, lam] = sospolyvar(prog,monomials(x,0:(2*kappa-2)));
prog = soseq(prog,p-sig0-lam*h);
options.solver = 'mosek';
prog = sossolve(prog,options);

sig0_val = cleanpoly(sosgetsol(prog,sig0),1e-6);
lam_val = cleanpoly(sosgetsol(prog,lam),1e-6);

