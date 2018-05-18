K = 4;
N = 6;
Nx = 1+K*2^(N+2);
sqB = linspace(-K,K,Nx);
%B = 0.5*sqB.^2;
B = cos(sqB);
x = linspace(-K-1,K+1,1000);

% Construction of the coefficients matrix
Mu = computeMu(B,N,0,K);
value = b(Mu,K,x);
figure
plot(x,value)
grid on
grid minor
xlabel('$x$','Interpreter','latex')
ylabel('$Id^N(x)$','Interpreter','latex')
chn = ['Approximation of the derivative of B by Haar wavelets (N = ',num2str(N),')'];
title(chn,'Interpreter','latex')