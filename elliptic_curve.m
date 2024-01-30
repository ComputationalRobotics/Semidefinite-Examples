clc; clear; close all;

x = -2:0.01:2; 
y = -2:0.01:2; 
[X,Y] = meshgrid(x,y);

ineq = (-X - 5 <= 0) & ...
    (-X.^2 + 2*X - Y.^2 + 7 >=0) & ...
    (3 + X - X.^3 - 3*X.^2 - 2*Y.^2 >= 0);

h = pcolor(X,Y,double(ineq)) ;
h.EdgeColor = 'none' ;


