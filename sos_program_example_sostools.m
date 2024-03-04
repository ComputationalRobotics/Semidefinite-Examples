clc; clear; close all;
restoredefaultpath;

addpath(genpath('../SOSTOOLS'))
addpath(genpath('../../../mosek'))

dim_y = 2;
dim_x = 1;
x = mpvar('x',1); % dummy x 
prog = sosprogram(x);
y_name = {}; % decision variable y
for i = 1:dim_y
    y_name{i} = sprintf('y_%d',i);
end
y = dpvar(y_name);
y = y(:);
prog = sosdecvar(prog,y);
% two polynomials
p1 = x^4 + y(1)*x + (2 + y(2));
p2 = (y(1)-y(2)+1)*x^2 + y(2)*x + 1;
% two SOS polynomials
[prog, sig1] = sossosvar(prog,monomials(x,0:2));
[prog, sig2] = sossosvar(prog,monomials(x,0:1));
% matching coefficients
prog = soseq(prog,p1-sig1);
prog = soseq(prog,p2-sig2);
% objective
prog = sossetobj(prog,-(y(1)+y(2)));
% choose solver
options.solver = 'mosek';
prog = sossolve(prog,options);
% get solution
ystar = double(sosgetsol(prog,y));
Q1star = reshape(prog.solinfo.x(2+1:2+9),3,3);
Q2star = reshape(prog.solinfo.x(2+9+1:end),2,2);



