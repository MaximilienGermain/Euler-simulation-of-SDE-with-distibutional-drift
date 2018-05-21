% To fix omega, we must fix startNT, start, and Kmax
X0 = 1;
startNT = 1;
NT = 6; % Time grid precision 2^(-N)
N = 4; % Space grid precision 2^(-N)
start = 1; 
T = 1.2;
%K = ceil(exp(2)+T*3);

%K = start+5;
Kmax = start+4;
%K = 5; % We cannot take K = N yet without changing omega
H = 0.95;
%H = 0.95;
%H = 0.985;

Nx = 1+K*2^(start+2); % 2 times more precise than the grid for b^N
[xgrid,B,M] = createfBm(H,K,N,start,Nx,-K,1000);

%X = eulerMethod(X0,startNT,NT,N,T,H,B,xgrid,1,K);
for NT=6:10
    X = eulerMethod(X0,startNT,NT,N,T,H,B,xgrid,0,Kmax,0);
end