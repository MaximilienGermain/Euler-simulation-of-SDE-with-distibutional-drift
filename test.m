% To fix omega, we must fix seed and startNT. To fix b we must fix startN, and Kmax
X0 = 1;
startNT = 1;
startN = 1; 
NT = 6; % Time grid precision 2^(-N)
N = 6; % Space grid precision 2^(-N)
T = 1.2;
Kmax = startN+9;
Nx = 1+Kmax*2^(startN+2); % 2 times more precise than the grid for b^N
H = 0.85;
graphHaar = 0;
control = 0;
testId = 0;
frames = [];
seed = 2;

[xgrid,B,M] = createfBm(H,Kmax,N,startN,Nx,-Kmax,1000);

% graphHaar = 1;
% control = 0;
% [X,frame,haar,control,usual] = eulerMethod(X0,startNT,NT,N,T,H,B,xgrid,testId,Kmax,graphHaar,control,seed);
% imshow(['Haar H = ',num2str(H),' ; N = ',num2str(N),'.png'])

% Test with b=id
% testId = 1;
% graphHaar = 1;
% B = linspace(-Kmax,Kmax,1+Kmax*2^(N+2));
% [X,frame,haar,control,usual] = eulerMethod(X0,startNT,NT,N,T,H,B.^2/2,B,testId,Kmax,graphHaar,control,seed);
%imshow(['Haar H = ',num2str(H),' ; N = ',num2str(N),'.png'])
%imshow(['Identity test n = ',num2str(NT),' ; N = ',num2str(N),'.png'])

% Convergence in n
% for NT=2:12
%     startNT = 2;
%     %control = 1;
%     %graphHaar = 1;
%     [X,frame,haar,control,usual] = eulerMethod(X0,startNT,NT,N,T,H,B,xgrid,testId,Kmax,graphHaar,control,seed);
%     frames = [frames frame];
% end
% video('Convergence time grid',frames)
% close all

% Convergence in N
% for N=startN:startN+6
%     [xgrid,B,M] = createfBm(H,Kmax,N,startN,Nx,-Kmax,1000);
%     NT = 8;
%     control = 1;
%     [X,frame,haar,control,usual] = eulerMethod(X0,startNT,NT,N,T,H,B,xgrid,testId,Kmax,graphHaar,control,seed);
%     frames = [frames frame];
% end
% video('Convergence x grid',frames)
% close all

%  Monte-Carlo N
% M = 100;
% expectations = [];
% var = [];
% NT = 10;
% startNT = 10;
% Nmax = 11;
% minN = 5;
% maxN = NTmax-2;
% for N=minN:maxN
%     values = [];
%     diff = N - N;
%     for o=1:M % Check the seed
%         [X,frame,haar,control,usual] = eulerMethod(X0,startNT,NT,N,T,H,B,xgrid,testId,Kmax,graphHaar,control,2^(ceil(log2(M))+N)+o);
%         [Xref,frame,haar,control,usual] = eulerMethod(X0,startNT,NT,Nmax,T,H,B,xgrid,testId,Kmax,graphHaar,control,2^(ceil(log2(M))+N)+o);       
%         values = [values max((X-Xref).^2)]; 
%         close all
%     end     
%     meanvalues = mean(values);
%     expectations = [expectations meanvalues];
%     var = [var M/(M+1)*(mean(values.^2)-meanvalues.^2)]; % Ckeck for bias
% end
% 
% figure
% plot(minN:maxN,expectations,'o')
% grid on 
% grid minor
% hold on
% xlabel('$N$','Interpreter','latex')
% ylabel('$\mathbb{E}[\underset{0\leq t\leq T}{\sup}|X^{N,n}_t-X^N_t|^2]$','Interpreter','latex')
% a = plot(minN:maxN,expectations+1.96*sqrt(var)/sqrt(M),'--b');
% plot(minN:maxN,expectations-1.96*sqrt(var)/sqrt(M),'--b')
% legend([a],'95$\%$ confidence interval','Interpreter','latex')
% title('Error when $N$ varies','Interpreter','latex')

%  Monte-Carlo n
M = 100;
expectations = [];
var = [];
NTmax = 11;
startNT = 5;
minN = 5;
maxN = NTmax-2;
for NT=minN:maxN
    values = [];
    diff = NTmax - NT;
    for o=1:M % Check the seed
        2^(ceil(log2(M))+NT)+o;
        [X,frame,haar,control,usual] = eulerMethod(X0,startNT,NT,N,T,H,B,xgrid,testId,Kmax,graphHaar,control,2^(ceil(log2(M))+NT)+o);
        [Xref,frame,haar,control,usual] = eulerMethod(X0,startNT,NTmax,N,T,H,B,xgrid,testId,Kmax,graphHaar,control,2^(ceil(log2(M))+NT)+o);       
        values = [values max((X-Xref([1 1+2^diff*(1:length(X)-1)])).^2)]; 
        close all
    end     
    meanvalues = mean(values);
    expectations = [expectations meanvalues];
    var = [var M/(M+1)*(mean(values.^2)-meanvalues.^2)];
end

figure
plot(minN:maxN,expectations,'o')
grid on 
grid minor
hold on
xlabel('$n$','Interpreter','latex')
ylabel('$\mathbb{E}[\underset{0\leq t\leq T}{\sup}|X^{N,n}_t-X^N_t|^2]$','Interpreter','latex')
a = plot(minN:maxN,expectations+1.96*sqrt(var)/sqrt(M),'--b');
plot(minN:maxN,expectations-1.96*sqrt(var)/sqrt(M),'--b')
legend([a],'95$\%$ confidence interval','Interpreter','latex')
title('Error when $n$ varies','Interpreter','latex')  