% To fix omega, we must fix seed and startNT. To fix b we must fix startN, and Kmax
X0 = 1;
startNT = 10;
startN = 3; 
NT = 10; % Time grid precision 2^(-N)
N = 5; % Space grid precision 2^(-N) Limited by 8 ????
T = 1.2;
Kmax = startN+7;
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
rng('shuffle');

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
for NT=2:13
    startNT = 2;
    PlotActive = 1;
    %control = 1;
    %graphHaar = 1;
    [X,~,frame,haar,control,usual] = eulerMethod(X0,startNT,NT,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,seed,PlotActive);
    frames = [frames frame];
end
video('Convergence time grid',frames)
close all

% Convergence in N
% for N=startN:startN+5
%     [xgrid,B,M] = createfBm(H,Kmax,N,startN,Nx,-Kmax,1000);
%     Mu = computeMu(B,N,testId,Kmax);
%     control = 1;
%     graphHaar = 1;
%     PlotActive = 1;
%     [X,~,frame,haar,control,usual] = eulerMethod(X0,startNT,NT,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,seed,PlotActive);
%     frames = [frames frame];
% end
% video('Convergence x grid',frames)
% close all

%  Monte-Carlo N
% MC = 1000;
% expectations = [];
% var = [];
% NT = 10;
% Nmax = 9;
% Kmax = Nmax;
% minN = 3;
% startN = minN;
% Nx = 1+Kmax*2^(startN+2);   
% maxN = Nmax-1;
% [xref,Bref,M] = createfBm(H,Kmax,Nmax,startN,Nx,-Kmax,1000);
% Muref = computeMu(Bref,Nmax,testId,Kmax);
% PlotActive = 0;
% tic
% for N=minN:maxN
%     values = [];    
%     for o=1:MC  
%         seed = randi(10^8);
%         [xgrid,B,M] = createfBm(H,Kmax,N,startN,Nx,-Kmax,1000);
%         Mu = computeMu(B,N,testId,Kmax);
%         [X,Y,~,~,~,~] = eulerMethod(X0,startNT,NT,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,seed,PlotActive);
%         [Xref,Yref,~,~,~,~] = eulerMethod(X0,startNT,NT,Nmax,T,H,Bref,Muref,xref,testId,Kmax,graphHaar,control,seed,PlotActive);       
%         values = [values 0.5*abs(max(X(end)-Xref(end))+max(Y(end)-Yref(end)))];
%         o
%     end     
%     meanvalues = mean(values);
%     expectations = [expectations meanvalues];
%     var = [var MC/(MC+1)*(mean(values.^2)-meanvalues.^2)];
% end
% toc
% figure
% Ns = minN:maxN;
% plot(log(2.^Ns),log(expectations),'o')
% grid on 
% grid minor
% hold on
% xlim([min(log(2.^Ns)) max(log(2.^Ns))])
% %ylim([min(log(expectations-1.96*sqrt(var)/sqrt(MC))) max(log(expectations+1.96*sqrt(var)/sqrt(MC)))])
% xlabel('$\log(N)$','Interpreter','latex')
% ylabel('$\log|E[X^{N}_t-X_t]|$','Interpreter','latex')
% a = plot(log(2.^Ns),log(expectations+1.96*sqrt(var)/sqrt(MC)),'--b');
% plot(log(2.^Ns),log(expectations-1.96*sqrt(var)/sqrt(MC)),'--b')
% legend([a],'95% confidence interval','Interpreter','latex')
% title('Error when $N$ varies estimated with Monte-Carlo method','Interpreter','latex')
% [beta0,beta1] = linearRegression(log(expectations)',log(2.^Ns)');
% order = - beta1

%  Monte-Carlo n
% M = 400;
% X0 = 0;
% NTmax = 12;
% minN = 5;
% maxN = NTmax-2;
% expectations = zeros(1,maxN-minN+1);
% var = zeros(1,maxN-minN+1);
% PlotActive = 0;
% stop = 0;
% tic
% for NT=minN:maxN
%     values = zeros(1,M);
%     diff = NTmax - NT;
%     for o=1:M 
%         startNT = NT;
%         seed = randi(10^8);
%         [X,Y,~,~,~,~] = eulerMethod(X0,startNT,NT,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,seed,PlotActive);
%         [Xref,Yref,~,~,~,~] = eulerMethod(X0,startNT,NTmax,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,seed,PlotActive);       
%         values(o) = 0.5*(max((X-Xref([1 1+2^diff*(1:length(X)-1)])).^2)+max((Y-Yref([1 1+2^diff*(1:length(X)-1)])).^2)); 
%         o
%         %if o>=2 && values(o) > 1.5
% %             [X,Y,~,~,~,~] = eulerMethod(X0,startNT,NT,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,seed,1);
% %             [Xref,Yref,~,~,~,~] = eulerMethod(X0,startNT,NTmax,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,seed,1);
% %             stop = 1;
% %             break
%         %end
%     end 
% %     if stop == 1
% %         break
% %     end
%     meanvalues = mean(values);
%     expectations(NT - minN + 1) = meanvalues;
%     var(NT - minN + 1) = M/(M+1)*(mean(values.^2)-meanvalues.^2);
% end
% toc
% figure
% Ns = minN:maxN;
% plot(log(2.^Ns),log(expectations),'or')
% grid on 
% grid minor
% hold on
% xlabel('$\log(n)$','Interpreter','latex')
% ylabel('$\log(E[{\sup}|X^{N,n}_t-X^{N,n_0}_t|^2])$','Interpreter','latex')
% xlim([min(log(2.^Ns)) max(log(2.^Ns))])
% ylim([min(log(expectations-1.96*sqrt(var)/sqrt(M))) max(log(expectations+1.96*sqrt(var)/sqrt(M)))])
% a = plot(log(2.^Ns),log(expectations+1.96*sqrt(var)/sqrt(M)),'--b');
% plot(log(2.^Ns),log(expectations-1.96*sqrt(var)/sqrt(M)),'--b')
% cont = linspace(min(log(2.^Ns)),max(log(2.^Ns)),1000);
% plot(cont,1.1226 - 0.8521*cont)
% legend(a,'95% confidence interval')
% title('Error when $n$ varies estimated with Monte-Carlo method for 400 paths','Interpreter','latex') 
% [beta0,beta1] = linearRegression(log(expectations)',log(2.^Ns)');
% order = - beta1