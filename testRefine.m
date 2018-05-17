N = 2;
K = 3;
H = 0.95;
Nx = 1+K*2^(N+2);
[random,xgrid,B,M] = fbm(H,K,Nx,N);
[newgrid,newB,newM] = refinefbm(random,xgrid,H,N,M);

figure
plot(xgrid,B)
xlabel('x')
xlim([min(xgrid) max(xgrid)])
ylabel('B_x^H')
chn = ['Sample path of a fractional Brownian motion B^H_x (Nx = ',num2str(Nx),' ; H = ',num2str(H),')'];
title(chn)

figure
plot(newgrid,newB)
xlabel('x')
xlim([min(xgrid) max(xgrid)])
ylabel('B_x^H')
chn = ['Sample path of a fractional Brownian motion B^H_x (Nx = ',num2str(2*Nx-1),' ; H = ',num2str(H),')'];
title(chn)