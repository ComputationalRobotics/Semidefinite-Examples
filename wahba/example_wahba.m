clc; clear; close all
restoredefaultpath;

mosekpath   = '../../../../mosek';
addpath(genpath(pwd));
addpath(genpath('../spotless'));
addpath(genpath('../SOSprograms'));
addpath(genpath(mosekpath))

%% generate random problem
N       = 10; % number of measurements
outrate = 0.6; % outlier rate
problem.N               = N;
problem.Covariance      = 0*eye(3); % no noise on inliers
problem.v1_distribution = 'uniform';
problem.nrOutliers      = round(N*outrate);

[a, b, R_gt, problem]   = createWahbaProblem(problem);
betasq = 0.1;

%% Part 3: Implement the brute-force search algorithm


%% Part 4: Implement the dense moment-SOS hierarchy
n = 4 + N;
q = msspoly('q',4);
Rq = Rofq(q);
th = msspoly('th',N);
% equality constraints

% objective


%% Part 5: flat extension and extract solutions


%% Part 6: adversarial outliers


%% Part 7: correlative sparsity


%% Part 8: term sparsity
SDP = QUASAR_Problem(a,b,betasq);
At = sparsevec(SDP.blk,SDP.Acell);
c = sparsevec(SDP.blk,SDP.C);
K.s = SDP.blk{1,2};
prob = convert_sedumi2mosek(At, SDP.b, c, K);
[~,res]    = mosekopt('minimize info',prob);
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);


