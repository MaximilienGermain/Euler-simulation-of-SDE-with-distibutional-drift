% Simulation of a realisation of a fractional Brownian motion with Hurst 
% index H on a Nx size grid
function [random,xgrid,B,M] = fbm(H,K,N,Nx,A,seed)

% Setting the seed to 1
rng(seed,'twister');

xgrid = linspace(A,K,Nx);
dx = 1/2^(N+1); % 2 times more precise than the grid for b^N

% Simulation of Nx - 1 independent gaussian variables
random = randn(1,Nx-1)';

% Construction of the correlation matrix
Gamma = zeros(Nx-1);
for i=1:Nx-1
    for j=1:i
        Gamma(i,j) = dx^(2*H)*(i^(2*H)+j^(2*H)-abs(i-j)^(2*H))/2;
        Gamma(j,i) = Gamma(i,j);
    end
end

M = chol(Gamma)'; % Here we must transpose M because of the definition chosen in Cholesky Matlab algorithm
B = [0 ; M*random];

end