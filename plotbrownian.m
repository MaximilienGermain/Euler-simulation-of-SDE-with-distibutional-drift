startN = 2; 
N = 3; 
Kmax = 4;
Nx = 1+Kmax*2^(startN+2) % 2 times more precise than the grid for b^N
H = 0.85;
[xgrid,B,~] = createfBm(H,Kmax,N,startN,Nx,-Kmax,1000);
figure
plot(xgrid,B)
grid on
grid minor
xlabel('$x$','Interpreter','latex')
ylabel('$B^H_x$','Interpreter','latex')
title(['Sample path of a fractional brownian motion with Hurst index ',num2str(H)])