% Simulation of a realisation of a fractional Brownian motion with Hurst 
% index H on a Nx size grid
function [xgrid,B] = fbm(H,N)

% Setting the seed to 1
%rng(10,'twister');
rng(1000,'twister');
Nx = 1 + (N+1)*2^(N+2);
xgrid = linspace(-N,N+1,Nx);
dx = 1/2^(N+1);

% Simulation of Nx independent gaussian variables
random = randn(1,Nx-1);

% Construction of the correlation matrix
Gamma = zeros(Nx-1);
for i=1:Nx-1
    for j=1:i
        Gamma(i,j) = dx^(2*H)*(i^(2*H)+j^(2*H)-abs(i-j)^(2*H))/2;
        Gamma(j,i) = Gamma(i,j);
    end
end

M = chol(Gamma);
% Here we must transpose M because of the definition chosen in Cholesky Matlab algorithm
B = [0  random*M];
% figure
% plot(xgrid,B)
% xlabel('x')
% xlim([min(xgrid) max(xgrid)])
% ylabel('B_x^H')
% chn = ['Sample path of a fractional Brownian motion B^H_x (Nx = ',num2str(Nx),' ; H = ',num2str(H),')'];
% title(chn)

end