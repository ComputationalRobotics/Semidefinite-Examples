clc; clear; close all;
restoredefaultpath;

addpath(genpath('../YALMIP'))
addpath(genpath('../../../mosek'))

% define all variables
x = sdpvar(1);
y = sdpvar(2,1);
Q1 = sdpvar(3,3);
Q2 = sdpvar(2,2);
% define polynomials
p1 = x^4 + y(1)*x + (2 + y(2));
p2 = (y(1)-y(2)+1)*x^2 + y(2)*x + 1;
% SOS polynomials
sig1 = monolist(x,2)' * Q1 * monolist(x,2);
sig2 = monolist(x,1)' * Q2 * monolist(x,1);
% define all constraints
F = [Q1>=0, Q2>=0,...
    coefficients(p1-sig1,x)==0,...
    coefficients(p2-sig2,x)==0];
% objective
obj = -(y(1)+y(2));
% solve
options = sdpsettings('solver','mosek');
optimize(F,obj,options);
% extract solution
ystar = value(y);
Q1star = value(Q1);
Q2star = value(Q2);