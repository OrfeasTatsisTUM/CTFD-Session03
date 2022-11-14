function [stencil] = stamp(i, j, X, Y, lamda, alpha, Tinf, boundary)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%stencil calculates the linear equation for node (i, j)

%  Input:
%      i         node number in x direction
%      j         node number in y direction
%      X         x position of the nodes
%      Y         y position of the nodes
%      b         right-hand side value for node (i,j)
%      alpha     alpha
%      Tinf      Tinf for Robin BC
%      boundary  defines the boundary conditions
%      verbose   verbosity level
%
%  Output:
%      stencil   linear equation for node (i,j)
%      b         new right-hand side value for node (i,j)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Init

n = size(X, 1);
m = size(X, 2);
stencil = zeros(1, n*m);
index=@(ii, jj) ii + (jj-1)*n;

% Determine the node positon
if j == 1
    nodePosition = 'West';
elseif j == size(X,2)
    nodePosition = 'East';
elseif i== size(X,1)
    nodePosition = 'South';
elseif i == 1
    nodePosition = 'North';
else
    nodePosition = 'inner Node';
end


% Calculate the equation for the correct node position
switch nodePosition
%% Inner
    case 'inner Node'
        
        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        data_inner
        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
        %$$$$$$$$$$$$$$$$$$$$$ Stencil $$$$$$$$$$$$$$$$$$$
        build_inner
        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
              
        % P
        stencil(index(i, j))     = lamda(i,j)      * D0;

        % East
        stencil(index(i, j+1))   = lamda(i,j+1)    * D3;

        % West
        stencil(index(i, j-1))   = lamda(i,j-1)    * D_3;

        % South
        stencil(index(i+1, j))   = lamda(i+1,j)    * D1;
        
        % North
        stencil(index(i-1, j))   = lamda(i-1,j)    * D_1;
        
        % NW
        stencil(index(i-1, j-1)) = lamda(i-1,j-1)  * D_4;
        
        % NE
        stencil(index(i-1, j+1)) = lamda(i-1,j+1)  * D2;
        
        % SW
        stencil(index(i+1, j-1)) = lamda(i+1, j-1) * D_2;
        
        % SE
        stencil(index(i+1, j+1)) = lamda(i+1, j+1) * D4;
        
%% South
    case 'South'

        if strcmp(boundary.south, 'Dirichlet')
            stencil(index(i, j))     = 1;
        else

            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            data_south
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

            bc_control = strcmp(boundary.south, 'Robin'); % factor that includes T_P in 3.16 (A.14)
            
            %$$$$$$$$$$$$$$$$$$$$$ Stencil $$$$$$$$$$$$$$$$$$$
            build_south
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

            % P
            stencil(index(i, j))     = lamda(i,j)      * D0;

            % East
            stencil(index(i, j+1))   = lamda(i,j+1)    * D3;

            % West
            stencil(index(i, j-1))   = lamda(i,j-1)    * D_3;

            % North
            stencil(index(i-1, j))   = lamda(i-1,j)    * D_1;

            % NW
            stencil(index(i-1, j-1)) = lamda(i-1,j-1)  * D_4;

            % NE
            stencil(index(i-1, j+1)) = lamda(i-1,j+1)  * D2;
        end
      
%% North
    case 'North'

        if strcmp(boundary.north, 'Dirichlet')
            stencil(index(i, j))     = 1;
        else
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            data_north
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            
            bc_control = strcmp(boundary.north, 'Robin'); % factor that includes T_P in 3.16 (A.14)

            %$$$$$$$$$$$$$$$$$$$$$ Stencil $$$$$$$$$$$$$$$$$$$
            build_north
            %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

            % P
            stencil(index(i, j))     = lamda(i,j)      * D0;

            % East
            stencil(index(i, j+1))   = lamda(i,j+1)    * D3;

            % West
            stencil(index(i, j-1))   = lamda(i,j-1)    * D_3;

            % South
            stencil(index(i+1, j))   = lamda(i+1,j)    * D1;

            % SW
            stencil(index(i+1, j-1)) = lamda(i+1, j-1) * D_2;
            
            % SE
            stencil(index(i+1, j+1)) = lamda(i+1, j+1) * D4;
        end

%% East
    case 'East'
        if strcmp(boundary.east, 'Dirichlet')
            stencil(index(i, j))     = 1;

        end

%% West
    case 'West'
        if strcmp(boundary.west, 'Dirichlet')
            stencil(index(i, j))     = 1;
        end
end

