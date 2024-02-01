clc; clear; close all;

syms x y z lam
A = [1, x, y;
     x, 1, z;
     y, z, 1];
p = det(lam*eye(3) - A);
p = collect(p,lam);


xx = -2:0.02:2; 
yy = -2:0.02:2;
zz = -2:0.02:2;
[X,Y,Z] = meshgrid(xx,yy,zz);

ineq = (-X.^2 - Y.^2 - Z.^2 + 3 >= 0) & ...
    (X.^2 - 2*X.*Y.*Z + Y.^2 + Z.^2 -1 <= 0);

Xin = X(ineq);
Yin = Y(ineq);
Zin = Z(ineq);

cloud = [Xin,Yin,Zin];
k = boundary(cloud);
trisurf(k,cloud(:,1),cloud(:,2),cloud(:,3))
axis equal