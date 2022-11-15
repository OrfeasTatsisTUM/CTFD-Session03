function T = solveFVM(dimY, dimX, X, Y, boundary, TD, lamda, alpha, Tinf)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File solveFVM.m
%
% This routine set up the linear system and solve it
%
% input
% M         Spatial Matrix M
% X         Matrix x coordinates 
% Y         Matrix y coordinates
% boundary  String vector. Boundary types.
% TD        Temperature for each boundary (if Dirichlet)
% alpha     convective heat transfer coefficient
% Tinf      Temperature of the surrouding fluid 
%
% output
% T         Temperature field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Index maps the node position to the correct linear equation

index = @(ii, jj) ii + (jj-1) * dimY;


% B is the right-hand side of the linear system
B = zeros(1,dimY*dimX);

%% Set boundary conditions
% in case of Neumann BC we have qdot=0 so there is nothing to be added to B

% North
if strcmp(boundary.north, 'Dirichlet')
    for i = 1:dimX
%         if i ~= 1 && i ~= dimX
            B(index(1,i)) = TD.north;
%         else
%             B(index(1,i)) = B(index(1,i)) + TD.north/2;
%         end
    end
elseif strcmp(boundary.north, 'Robin')
    for i = 1:dimX
        B(index(1,i)) = -alpha/lamda(1,i) * Tinf;
    end
end

% South
if strcmp(boundary.south, 'Dirichlet')
    for i = 1:dimX
%         if i ~= 1 && i ~= dimX
            B(index(dimY,i)) = TD.south;
%         else
%             B(index(dimY,i)) = B(index(dimY,i)) + TD.south/2;
%         end
    end
elseif strcmp(boundary.south, 'Robin')
    for i = 1:dimX
        B(index(dimY,i)) = -alpha/lamda(dimY,i) * Tinf;
    end
end

% East
if strcmp(boundary.east, 'Dirichlet')
    for i = 1:dimY
%         if i~= 1 && i~=dimY
            B(index(i,dimX)) = TD.east;
%         else
%             B(index(i,dimX)) = B(index(i,dimX)) + TD.east/2;
%         end
    end
elseif strcmp(boundary.east, 'Robin')
    for i = 1:dimX
        B(index(i,dimX)) = -alpha/lamda(i,dimX) * Tinf;
    end
end

% West
if strcmp(boundary.west, 'Dirichlet')
    for i = 1:dimY
%         if i~= 1 && i~=dimY
            B(index(i,1)) = TD.west;
%         else
%             B(index(i,1)) = B(index(i,1)) + TD.west/2;
%         end
    end
elseif strcmp(boundary.west, 'Robin')
    for i = 1:dimX
        B(index(i,1)) = -alpha/lamda(i,1) * Tinf;
    end
end

%% Set up the system matrix A
A = zeros(dimY*dimX);

for i = 1:dimY
    for j = 1:dimX
        % Fill the system matrix and the right-hand side for node (i,j)
        [A(index(i,j), :)] =  stamp(i, j, X, Y, lamda, alpha, Tinf, boundary);
    end
end


%% Solve the linear system
T1(:) = A \ B';

% Convert solution vector into matrix
T = zeros(dimY,dimX);
for j = 1:dimX
    for i = 1:dimY
        T(i,j) = T(i,j) + T1(index(i,j));
    end
end
