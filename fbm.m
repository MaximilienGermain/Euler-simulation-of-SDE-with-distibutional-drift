function B = fbm(H,Nx,xmax)

% Setting the seed to 1
rng(10,'twister');
x = xmax/Nx*(0:Nx);
dx = xmax/Nx;

% Simulation of Nx independent gaussian variables
random = randn(1,Nx);

% Construction of the correlation matrix
Gamma = zeros(Nx);
for i=1:Nx
    for j=1:i
        Gamma(i,j) = dx^(2*H)*(i^(2*H)+j^(2*H)-abs(i-j)^(2*H))/2;
        Gamma(j,i) = Gamma(i,j);
    end
end

M = chol(Gamma);
% Here we must transpose M because of the definition chosen in Cholesky Matlab algorithm
B = [0 ; M'*random'];
plot(x,B)
xlabel('x')
ylabel('B_x^H')
chn = ['One realisation of a fractional Brownian motion (Nx = ',num2str(Nx),' ; H = ',num2str(H),')'];
title(chn)

end