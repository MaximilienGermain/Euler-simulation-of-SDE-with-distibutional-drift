% To fix omega, we must fix seed and startNT. To fix b we must fix startN, and Kmax
X0 = 1;
startNT = 10;
startN = 3; 
NT = 10; % Time grid precision 2^(-N)
N = 5; % Space grid precision 2^(-N)
T = 1.2;
Kmax = startN+9;
Nx = 1+Kmax*2^(startN+2); % 2 times more precise than the grid for b^N
H = 0.85;
graphHaar = 0;
control = 0;
testId = 0;
frames = [];
seed = 2;
PlotActive = 1;

[xgrid,B,M] = createfBm(H,Kmax,N,startN,Nx,-Kmax,1000);
Mu = computeMu(B,N,testId,Kmax);

% graphHaar = 1;
% % control = 0;
% [X,~,frame,haar,control,usual] = eulerMethod(X0,startNT,NT,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,seed,PlotActive);
% imshow(['Haar H = ',num2str(H),' ; N = ',num2str(N),'.png'])

% Test with b=id
% testId = 1;
% graphHaar = 1;
% Kmax = 5;
% N = 4;
% B = linspace(-Kmax,Kmax,1+Kmax*2^(N+2));
% Mu = computeMu(B.^2/2,N,testId,Kmax);
% [X,~,frame,haar,control,usual] = eulerMethod(X0,startNT,NT,N,T,H,B.^2/2,Mu,B,testId,Kmax,graphHaar,control,seed,PlotActive);
% imshow(['Haar H = ',num2str(H),' ; N = ',num2str(N),'.png'])
% imshow(['Identity test n = ',num2str(NT),' ; N = ',num2str(N),'.png'])

% Convergence in n
% for NT=2:12
%     startNT = 2;
%     %control = 1;
%     %graphHaar = 1;
%     [X,~,frame,haar,control,usual] = eulerMethod(X0,startNT,NT,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,seed,PlotActive);
%     frames = [frames frame];
% end
% video('Convergence time grid',frames)
% close all

% Convergence in N
for N=startN:startN+7
    [xgrid,B,M] = createfBm(H,Kmax,N,startN,Nx,-Kmax,1000);
    Mu = computeMu(B,N,testId,Kmax);
    control = 1;
    [X,~,frame,haar,control,usual] = eulerMethod(X0,startNT,NT,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,seed,PlotActive);
    frames = [frames frame];
end
video('Convergence x grid',frames)
close all

%  Monte-Carlo N
% M = 2;
% expectations = [];
% var = [];
% NT = 10;
% startNT = 10;
% Nmax = 8;
% Kmax = Nmax;
% minN = 2;
% startN = minN;
% Nx = 1+Kmax*2^(startN+2);
% maxN = Nmax-4;
% [xref,Bref,M] = createfBm(H,Kmax,Nmax,startN,Nx,-Kmax,1000);
% Muref = computeMu(Bref,Nmax,testId,Kmax);
% tic
% for N=minN:maxN
%     values = [];    
%     for o=1:M
%         o
%         [xgrid,B,M] = createfBm(H,Kmax,N,startN,Nx,-Kmax,1000);
%         Mu = computeMu(B,N,testId,Kmax);
%         [X,Y,~,~,~,~] = eulerMethod(X0,startNT,NT,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,2^(ceil(log2(M))+N)+o,PlotActive);
%         [Xref,Yref,~,~,~,~] = eulerMethod(X0,startNT,NT,Nmax,T,H,Bref,Muref,xref,testId,Kmax,graphHaar,control,2^(ceil(log2(M))+N)+o,PlotActive);       
%         X
%         Y
%         Xref
%         Yref
%         values = [values 0.5*(max((X-Xref).^2)+max((Y-Yref).^2))];
%         close all
%     end     
%     meanvalues = mean(values);
%     expectations = [expectations meanvalues];
%     var = [var M/(M+1)*(mean(values.^2)-meanvalues.^2)]; % Ckeck for bias
% end
% toc
% figure
% Ns = minN:maxN;
% loglog(2.^Ns,expectations,'o')
% grid on 
% grid minor
% hold on
% xlabel('$N$','Interpreter','latex')
% ylabel('$E[{\sup}|X^{N,n}_t-X^N_t|^2]$','Interpreter','latex')
% a = loglog(2.^Ns,expectations+1.96*sqrt(var)/sqrt(M),'--b');
% loglog(2.^Ns,expectations-1.96*sqrt(var)/sqrt(M),'--b')
% legend([a],'95$\%$ confidence interval','Interpreter','latex')
% title('Error when $N$ varies','Interpreter','latex')

%  Monte-Carlo n
% M = 4;
% expectations = [];
% var = [];
% NTmax = 12;
% minN = 4;
% maxN = NTmax-2;
% PlotActive = 0;
% tic
% for NT=minN:maxN
%     values = [];
%     diff = NTmax - NT;
%     for o=1:M 
%         startNT = NT;
%         [X,Y,~,~,~,~] = eulerMethod(X0,startNT,NT,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,2^(ceil(log2(M))+NT)+o,PlotActive);
%         [Xref,Yref,~,~,~,~] = eulerMethod(X0,startNT,NTmax,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,2^(ceil(log2(M))+NT)+o,PlotActive);       
%         values = [values 0.5*(max((X-Xref([1 1+2^diff*(1:length(X)-1)])).^2)+max((Y-Yref([1 1+2^diff*(1:length(X)-1)])).^2))]; 
%         o
%     end     
%     meanvalues = mean(values);
%     expectations = [expectations meanvalues];
%     var = [var M/(M+1)*(mean(values.^2)-meanvalues.^2)];
% end
% toc
% figure
% Ns = minN:maxN;
% loglog(2.^Ns,expectations,'or')
% grid on 
% grid minor
% hold on
% xlabel('$n$','Interpreter','latex')
% ylabel('$E[{\sup}|X^{N,n}_t-X^{N,n_0}_t|^2]$','Interpreter','latex')
% xlim([0.8*min(2.^Ns) 1.2*max(2.^Ns)])
% a = loglog(2.^Ns,expectations+1.96*sqrt(var)/sqrt(M),'--b');
% loglog(2.^Ns,expectations-1.96*sqrt(var)/sqrt(M),'--b')
% legend(a,'95% confidence interval')
% title('Error when $n$ varies','Interpreter','latex') 