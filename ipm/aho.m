%% Primal-dual path following
%% with the AHO newton direction
%% Newton direction computed directly using Matlab backslash

clc; clear; close all;
addpath(genpath(pwd))

%% test example: 
C = [2,1;1,0];
b = [1];
A = [1,0;0,1];
blk = {'s', 2};

c = svec(blk,{C});
At = svec(blk,{A});
