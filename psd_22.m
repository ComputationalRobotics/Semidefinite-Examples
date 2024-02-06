clc; clear; close all;

syms x y z lam
A = [x, y;
     y, z];
p = det(lam*eye(2) - A);
p = collect(p,lam);


xx = -2:0.02:2; 
yy = -2:0.02:2;
zz = -2:0.02:2;
[X,Y,Z] = meshgrid(xx,yy,zz);

ineq = (-X - Z <= 0) & ...
    (- Y.^2 + X.*Z >=0);

Xin = X(ineq);
Yin = Y(ineq);
Zin = Z(ineq);

cloud = [Xin,Yin,Zin];
k = boundary(cloud);
trisurf(k,cloud(:,1),cloud(:,2),cloud(:,3))
axis equal