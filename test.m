X0 = 1;
startNT = 8;
NT = 14; % Time grid precision 2^(-NT)
N = 5; % Space grid precision 2^(-N)
T = 2;
%K = ceil(exp(2)+T*3);
K = 5;
H = 0.95;
%H = 0.55;
%H = 0.985;

%[random,xgrid,B,M] = fbm(H,K,N);
[xgrid,B,M] = createfBm(H,K,N);

%X = eulerMethod(X0,startNT,NT,N,T,H,B,xgrid,1,K);
X = eulerMethod(X0,startNT,NT,N,T,H,B,xgrid,0,K);