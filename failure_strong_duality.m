%% Example 2.3: failure of strong duality
%% Example 2.4: strong duality holds
%% Observe how the solvers behave strangly when strong duality fails, and 
%% how they solve the problem really fast when strong duality holds

clc; clear; close all;
addpath(genpath("../sedumi"))
addpath(genpath("../../../mosek"))


%% Example 1
alpha = 1;
C = [alpha, 0,  0;
     0, 0, 0;
     0, 0, 0];
c = C(:);
b = [0;1];
A1 = [0, 0, 0;
      0, 1, 0;
      0, 0, 0];
A2 = [1, 0, 0;
      0, 0, 1;
      0, 1, 0];
A = [A1(:),A2(:)]';
K.s = 3;
[xopt,yopt,info] = sedumi(A,b,c,K);
prob = convert_sedumi2mosek(A, b, c, K);
[~,res] = mosekopt('minimize info',prob);

%% Example 2
C = [0,0;0,0];
c = C(:);
b = [0;2];
A1 = [1, 0; 0, 0];
A2 = [0,1;1,0];
A = [A1(:),A2(:)]';
K.s = 2;
[xopt,yopt,info] = sedumi(A,b,c,K);
prob = convert_sedumi2mosek(A, b, c, K);
[~,res] = mosekopt('minimize info',prob);

%% Example 3
C = [0,1;1,0];
c = C(:);
b = [-1;0];
A1 = [-1,0;0,0];
A2 = [0,0;0,-1];
A = [A1(:),A2(:)]';
K.s = 2;
[xopt,yopt,info] = sedumi(A,b,c,K);
prob = convert_sedumi2mosek(A, b, c, K);
[~,res] = mosekopt('minimize info',prob);

%% Strong Duality holds: Example 2.4
C = [2,1;1,0];
c = C(:);
b = [1];
A1 = [1,0;0,1];
A = [A1(:)]';
K.s = 2;
[xopt,yopt,info] = sedumi(A,b,c,K);
prob = convert_sedumi2mosek(A, b, c, K);
[~,res] = mosekopt('minimize info',prob);





