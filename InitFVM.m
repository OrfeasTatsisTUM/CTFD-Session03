%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the parameters of the simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Geometry:
%
%    |---
%    |   ---
%    |      ----
%    |          |
% h1 |----------|  h2  <- symmetry axis
%    |          |
%    |      ----
%    |   ---
%    |---
%
%    |<--  l -->|


% Define dimension of the trapezoidal domain
% h2 <= h1 !

shape = 'linear';  % 'linear' or 'quadratic'

h1 = 4;
hm = 4;            % only necessary for quatratic option 
h2 = 2;
l = 5;

% Number of degrees of freedom (number of nodes per length)
dimX = 35;
dimY = 30;

% Shape of the Cooling Fin
% h2 <= h1 !

switch shape
    
    case 'linear'

        formfunction = @(xnorm) (1-xnorm)*h1/2 + xnorm*h2/2;

    case 'quadratic'

        c1 = h2+2*h1/2-2*hm;
        c2 = 2*hm - 3*h1/2 - h2/2;
        c3 = h1/2;

        formfunction = @(xnorm) c1*xnorm.^2 +c2*xnorm + c3;
        
    case 'crazy'

        d1 = 3;
        d2 = 4;
        
        formfunction = @(xnorm) (1-xnorm)*h1/2 + xnorm*h2/2+ (sin(2*pi*d1*xnorm)).*(1-(1-1/d2)*xnorm);
  
    % Add other cases of form function
end

%% Parameter for Conjugated Heat Transfer (For Session 04)
alpha = 2;
Tinf = 90;

%% Boundary conditions (Only Dirichlet applied in Session 03) 
% Type: 1) Dirichlet    2) Neumann    3) Robin
boundary.south = 'Neumann';
boundary.north = 'Dirichlet';
boundary.east  = 'Dirichlet';
boundary.west  = 'Dirichlet';

% Values for Dirichlet BC

TD.north = 10;
TD.south = 50;
TD.west  = 50;
TD.east  = 10;

%% Thermal Conductivity Parameters:
% Thermal conductivity Coefficient 0.05(ice) - 400(pure copper) [W/(m*K)]
% 1) homgenous      2) non_homogenous (region with different K)
% 3) random         4) linear (changing through x)
heat_conduc = 'homogenous';

% Define the Heat conductivity coefficient values
% for non_homogenous & linear cases it has been assumed that lamda changes on x axis

minlamda = 250;                     % minimum lamda value
deltalamda = 100;                   % lamda difference from side to side (x axis)

maxlamda = minlamda + deltalamda;   % maximum lamda value

switch heat_conduc

    case 'homogenous'
        lamda(1:dimY,1:dimX) = minlamda;

    case 'non_homogenous'
        lamda(1:dimY,1:round(dimX/2)) = minlamda;
        lamda(1:dimY,(round(dimX/2)+1):dimX) = maxlamda;

    case 'random'
        lamda = (deltalamda).*rand(dimY,dimX) + minlamda;

    case 'linear'
        lamda = zeros (dimY,dimX);
        dlamda = deltalamda/l;
        for i=1:dimX
            lamda(1:dimY,i) = minlamda + dlamda * (i-1)*l/(dimX-1);
        end
        
end
clear minlamda maxlamda deltalamda dlamda
