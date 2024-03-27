clc; clear; close all;

addpath(genpath("../../../mosek"))
addpath(genpath("./spotless"))
addpath(genpath("./SOSprograms"))

x = msspoly('x',2);
p = -( x(1)-1 )^2 - ( x(1)-x(2) )^2 - ( x(2)-3 )^2;
g = [1 - ( x(1)-1 )^2;
     1 - ( x(1)-x(2) )^2;
     1 - ( x(2)-3 )^2];

problem.vars            = x;
problem.objective       = p;
problem.equality        = []; 
problem.inequality      = g;
kappa                   = 2; % relaxation order
[SDP,info]              = dense_sdp_relax(problem,kappa);

prob       = convert_sedumi2mosek(SDP.sedumi.At,...
                                  SDP.sedumi.b,...
                                  SDP.sedumi.c,...
                                  SDP.sedumi.K);
[~,res]    = mosekopt('minimize info',prob);
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);

Xmom = Xopt{1};

if rank(Xmom,1e-6) == rank(Xmom(1:3,1:3),1e-6)
    fprintf("Flat extension condition satisfied!\n")
else
    fprintf("Flat extension condition not satisfied!\n")
end

[V,D] = sorteig(Xmom);

ids = (diag(D) > 1e-6);
V = V(:,ids);
U = rref(V');
U(abs(U) < 1e-6) = 0;

