clear; close all; clc;

%% Tatsis, Orfeas Emmanouil
%% Fernando Cruz Ceravalls
%% Yuechen Chen

%% SESSION_03
%  TUM - Ass. Professorship for Thermo Fluid Dynamics
%  WS022-023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to solve the 2D steady heat equation in a non-Cartesian Grid by
% the Finite Volumes Method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialize variables
InitFVM

%% Set up the mesh
[X, Y] = setUpMesh(dimY, dimX, l, formfunction);

%% Fill matrix A and vector B. Solve the linear system.
T = solveFVM(dimY, dimX, X, Y, boundary, TD, lamda, alpha, Tinf);

%% Make some plots
postprocess;