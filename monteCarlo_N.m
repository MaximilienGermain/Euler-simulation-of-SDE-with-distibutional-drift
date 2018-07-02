function [expectations,var] = monteCarlo_N(X0,T,H,Nmax,Kmax,graphHaar,control,testId,minN,maxN,PlotActive,MC,NT,startNT,startN)

Nx = 1+Kmax*2^(startN+2); 
[xref,Bref,~] = createfBm(H,Kmax,Nmax,startN,Nx,-Kmax,1000);
Muref = computeMu(Bref,Nmax,testId,Kmax);
expectations = [];
var = [];
f = waitbar(0,'Please wait...');

for N=minN:maxN
    values = [];    
    for o=1:MC  
        seed = randi(10^8);
        [xgrid,B,~] = createfBm(H,Kmax,N,startN,Nx,-Kmax,1000);
        Mu = computeMu(B,N,testId,Kmax);
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
plot(log(2.^Ns),log(expectations),'o')
grid on 
grid minor
hold on
xlim([min(log(2.^Ns)) max(log(2.^Ns))])
%ylim([min(log(expectations-1.96*sqrt(var)/sqrt(MC))) max(log(expectations+1.96*sqrt(var)/sqrt(MC)))])
xlabel('$\log(N)$','Interpreter','latex')
ylabel('$\log|E[X^{N}_t-X^{N_0}_t]|$','Interpreter','latex')
a = plot(log(2.^Ns),log(expectations+1.96*sqrt(var)/sqrt(MC)),'--b');
plot(log(2.^Ns),log(expectations-1.96*sqrt(var)/sqrt(MC)),'--b')
legend(a,'95% confidence interval','Interpreter','latex')
title(['Error when $N$ varies estimated with Monte-Carlo method with ',num2str(MC),' paths'],'Interpreter','latex')
[beta0,beta1] = linearRegression(log(2.^Ns)',log(expectations)');
order = - beta1
dispgrid = linspace(log(2^minN),log(2^maxN),1000);
plot(dispgrid,beta0+beta1*dispgrid)
delete(f)

end