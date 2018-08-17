function [expectations,var] = monteCarlo_N(X0,T,H,Nmax,Kmax,graphHaar,control,testId,minN,maxN,PlotActive,MC,NT,startNT,startN)

rng('shuffle');
Nx = 1+Kmax*2^(startN+2); 
[xref,Bref,~] = createfBm(H,Kmax,Nmax,startN,Nx,-Kmax,1000);
[Muref,Mureferr] = computeMu(Bref,Nmax,testId,Kmax,H);
expectations = [];
var = [];
f = waitbar(0,'Please wait...');
err = [];

for N=minN:maxN
    values = []; 
    [xgrid,B,~] = createfBm(H,Kmax,N,startN,Nx,-Kmax,1000);
    [Mu,~] = computeMu(B,N,testId,Kmax,H);
    %err = [err errorb(Muref,N,Nmax)];
    err = [err errorb(Mureferr,N,Nmax)];
    
    for o=1:MC  
        seed = randi(10^8);
        [X,Y,~,~,~,~] = eulerMethod(X0,startNT,NT,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,seed,PlotActive);
        waitbar(((N-minN)+(o-0.5)/MC)/(maxN-minN+1),f,['Computing Monte-Carlo simulation... (',num2str(floor(((N-minN)+(o-0.5)/MC)/(maxN-minN+1)*100)),' %)'])
        [Xref,Yref,~,~,~,~] = eulerMethod(X0,startNT,NT,Nmax,T,H,Bref,Muref,xref,testId,Kmax,graphHaar,control,seed,PlotActive);       
        values = [values 0.5*(abs(X(end)-Xref(end))+abs(Y(end)-Yref(end)))];
        %values = [values abs(X(end)-Xref(end))];
        %fprintf('%3.1f %%\n',100*((N-minN+1)+o/M)/(maxN-minN+1))
        waitbar(((N-minN)+o/MC)/(maxN-minN+1),f,['Computing Monte-Carlo simulation... (',num2str(floor(((N-minN)+o/MC)/(maxN-minN+1)*100)),' %)'])
    end     
    meanvalues = mean(values);
    expectations = [expectations meanvalues];
    var = [var MC/(MC+1)*(mean(values.^2)-meanvalues^2)];
end
figure
Ns = minN:maxN;
plot(log(Ns),log(expectations),'o')
grid on 
grid minor
hold on
xlim([min(log(Ns)) max(log(Ns))])
%ylim([min(log(expectations-1.96*sqrt(var)/sqrt(MC))) max(log(expectations+1.96*sqrt(var)/sqrt(MC)))])
xlabel('$\log(N)$','Interpreter','latex')
ylabel('$\log|E[X^{N,n_0}_T-X^{N_0,n_0}_T]|$','Interpreter','latex')
a = plot(log(Ns),log(expectations+1.96*sqrt(var)/sqrt(MC)),'--b');
plot(log(Ns),log(expectations-1.96*sqrt(var)/sqrt(MC)),'--b')
legend(a,{'$95\%$ confidence interval'},'Interpreter','latex')
title(['Error when $N$ varies estimated with Monte-Carlo method with ',num2str(MC),' paths'],'Interpreter','latex')
[beta0,beta1] = linearRegression(log(Ns)',log(expectations)');
order = - beta1
dispgrid = linspace(log(minN),log(maxN),1000);
plot(dispgrid,beta0+beta1*dispgrid)
delete(f)
name = ['MC = ',num2str(MC),'.eps'];
set(gcf,'PaperPositionMode','auto')
print(name,'-depsc','-tiff')

figure
plot(log(err),log(expectations),'o')
grid on 
grid minor
hold on
xlim([min(log(err)) max(log(err))])
%ylim([min(log(expectations-1.96*sqrt(var)/sqrt(MC))) max(log(expectations+1.96*sqrt(var)/sqrt(MC)))])
xlabel('$\log||b-b^N||$','Interpreter','latex')
ylabel('$\log|E[X^{N}_T-X^{N_0}_T]|$','Interpreter','latex')
a = plot(log(err),log(expectations+1.96*sqrt(var)/sqrt(MC)),'--b');
plot(log(err),log(expectations-1.96*sqrt(var)/sqrt(MC)),'--b')
legend(a,{'$95\%$ confidence interval'},'Interpreter','latex')
title(['Error when $N$ varies estimated with Monte-Carlo method with ',num2str(MC),' paths'],'Interpreter','latex')
[beta4,beta5] = linearRegression(log(err)',log(expectations)');
order = beta5
dispgrid = linspace(log(min(err)),log(max(err)),1000);
plot(dispgrid,beta4+beta5*dispgrid)
name = ['err MC = ',num2str(MC),'.eps'];
set(gcf,'PaperPositionMode','auto')
print(name,'-depsc','-tiff')

% figure
% plot(log(Ns),log(err))
% xlabel('$\log(N)$','Interpreter','latex')
% ylabel('$\log||b-b^N||$','Interpreter','latex')
% title(['Error in $H^s_2$ when $N$ varies'],'Interpreter','latex')
% [beta2,beta3] = linearRegression(log(Ns)',log(err)');
% dispgrid = linspace(log(minN),log(maxN),1000);
% plot(dispgrid,beta2+beta3*dispgrid)
% orderc = - beta3
% err
% name = ['error b ',num2str(MC),'.eps'];
% set(gcf,'PaperPositionMode','auto')
% print(name,'-depsc','-tiff')

end