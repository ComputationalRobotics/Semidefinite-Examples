%% Primal-dual path following
%% with the AHO newton direction

clc; clear; close all;
restoredefaultpath;
% Add SDPT3 to path as we will need helper functions there
addpath(genpath('../../SDPT3'))

%% Example problem
C = [2,1;1,0];
b = [1];
A = [1,0;0,1];
blk = {'s', 2};
n = blk{1,2};

c = svec(blk,{C});
c = c{1};
At = svec(blk,{A});
At = At{1};
A = At';

%% initialization
X = eye(2); x = svec(blk,X); % identity X
Z = eye(2); z = svec(blk,Z); % identity Z
y = lsqr(At,c - z); % solve the equation At*y + z = c
svecI = svec(blk,eye(n));

%% Primal-dual path following
itr     = 0;
max_itr = 100;
sigma   = 0.5;
tau     = 0.9;

% print header
breakyes = 0;
fprintf(" ITR |    r_d   |    r_p   |    r_c   |  max_kkt |    mu    |    pobj   |    dobj   |\n")

while true
    itr = itr + 1;
    X = smat(blk,x);
    Z = smat(blk,z);
    P = skron(blk,Z,eye(n));
    Q = skron(blk,X,eye(n));
    if itr > max_itr
        breakyes = 1;
        msg = "Maximum number of iterations reached.";
    end

    %% Step 1: computer mu
    mu = error("Please compute mu.");

    % primal, dual, and gap residuals
    r_d = c - At*y - z;
    r_p = b - A*x;
    r_c = mu*svecI - P*x;
    r_d_norm = norm(r_d);
    r_p_norm = norm(r_p);
    r_c_norm = norm(r_c);
    max_kkt = max([r_d_norm,r_p_norm,r_c_norm]);
    if max_kkt < 1e-6
        breakyes = 2;
        msg = "Max KKT residual below 1e-6.";
    end

    %% Step 2: form and solve linear system
    % Option 1: you can directly form the linear system in (2.23) and use
    % Matlab's "\" to solve the linear system directly
    % Option 2: you can also follow (2.24)-(2.26)
   
    direction = error("Please compute Newton direction.");

    dx = direction(1:length(x));
    dy = direction(1+length(x):length(y)+length(x));
    dz = direction(length(y)+length(x)+1:end);
    dX = smat(blk,dx);
    dZ = smat(blk,dz);

    %% Step 3: Compute step sizes
    % Complete the line_search function below
    alpha_hat = line_search(X,dX);
    beta_hat = line_search(Z,dZ);

    alpha = min(1,tau*alpha_hat);
    beta = min(1,tau*beta_hat);
    
    %% Step along Newton direction
    x = x + alpha*dx;
    y = y + beta*dy;
    z = z + beta*dz;
    
    %% Print info and exit if necessary
    fprintf("  %02d | %3.2e | %3.2e | %3.2e | %3.2e | %3.2e | %3.2e | %3.2e |\n",...
        itr,r_d_norm,r_p_norm,r_c_norm,max_kkt,mu,c'*x,b'*y);
    if breakyes
        fprintf("%s\n",msg);
        break;
    end
end



function alpha_hat = line_search(X,dX)
error("Please complete the line_search function.")
alpha_hat = 1;
end


