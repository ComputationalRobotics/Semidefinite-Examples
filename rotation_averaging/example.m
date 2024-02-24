clc; clear; close all;

%% number of nodes
n = 50;
noise_bound = 120;

%% random groundtruth rotations at each node
R_gt = [];
for i = 1:n
    if i == 1
        R_gt = [R_gt; eye(3)];
    else
        R_gt = [R_gt; rand_rotation];
    end
end

%% relative rotations at each edge, plus random noise
edges = [];
R_edges = [];
for i = 1:n
    for j = (i+1):n
        Ri = R_gt(blkIndices(i,3),:);
        Rj = R_gt(blkIndices(j,3),:);
        Rij = Ri'*Rj;
        Rij_noisy = Rij * rand_rotation('RotationBound',noise_bound/180*pi);
        edges = [edges; i, j];
        R_edges = [R_edges;Rij_noisy];
    end
end

%% Form Rtld
Rtld = zeros(3*n,3*n);
for edge_id = 1:size(edges,1)
    edge = edges(edge_id,:);
    i = edge(1); 
    j = edge(2);

    Rtld(blkIndices(i,3),blkIndices(j,3)) = R_edges(blkIndices(edge_id,3),:);
    Rtld(blkIndices(j,3),blkIndices(i,3)) = R_edges(blkIndices(edge_id,3),:)';
end

%% call YALMIP to solve
addpath(genpath('../../YALMIP'))
addpath(genpath('../../../../mosek'))

X = sdpvar(3*n,3*n);
F = [X >= 0];
for i = 1:n
    F = [F,
        X(blkIndices(i,3),blkIndices(i,3)) == eye(3)];
end
obj = trace(-Rtld*X);
optimize(F,obj);
Xval = value(X);
f_sdp = value(obj);

%% extract solutions
R_est = [];
R_errs = [];
for i = 1:n
    if i == 1
        Ri = eye(3);
    else
        Ri = project2SO3(Xval(1:3,blkIndices(i,3)));
    end
    R_errs = [R_errs, getAngularError(Ri, R_gt(blkIndices(i,3),:))];
    R_est = [R_est, Ri];
end
X_est = R_est'*R_est;
f_est = trace(-Rtld*X_est);

gap = abs(f_est - f_sdp) / (1 + abs(f_est) + abs(f_sdp));

figure;
plot(R_errs,'LineWidth',3);
xlabel('Node ID','FontSize',24)
ylabel('Rotation Error (Deg)','FontSize',24);
ax = gca; ax.FontSize = 20; grid on;

fprintf('Relative Suboptimality Gap: %3.2e.\n',gap);








