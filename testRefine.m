N = 3;
K = 1;
H = 0.95;
[random,xgrid,B,M] = fbm(H,K,N);
[newrandom,newgrid,newB,newM] = refinefbm(random,xgrid,H,N,M);

figure
plot(xgrid,B)
grid on
grid minor
xlabel('x')
xlim([min(xgrid) max(xgrid)])
ylabel('B_x^H')
chn = ['Sample path of a fractional Brownian motion B^H_x (N = ',num2str(N),' ; H = ',num2str(H),')'];
title(chn)

figure
plot(newgrid,newB)
grid on
grid minor
xlabel('x')
xlim([min(newgrid) max(newgrid)])
ylabel('B_x^H')
chn = ['Sample path of a fractional Brownian motion B^H_x (N = ',num2str(N+1),' ; H = ',num2str(H),')'];
title(chn)