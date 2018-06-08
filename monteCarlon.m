function [expectations,var] = monteCarlon(X0,xgrid,B,N,T,Mu,H,Kmax,graphHaar,control,testId,minN,maxN,PlotActive,M,NTmax)

expectations = zeros(1,maxN-minN+1);
var = zeros(1,maxN-minN+1);
f = waitbar(0,'Please wait...');

for NT=minN:maxN
    values = zeros(1,M);
    diff = NTmax - NT;
    
    % M independent realizations
    for o=1:M 
        startNT = NT;
        seed = randi(10^8);
        [X,Y,~,~,~,~] = eulerMethod(X0,startNT,NT,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,seed,PlotActive);
        %fprintf('%3.1f %%\n',100*((NT-minN)+(o-1/2)/M)/(maxN-minN+1))
        waitbar(((NT-minN)+o/M)/(maxN-minN+1),f,['Computing Monte-Carlo simulation... (',num2str(floor(((NT-minN)+(o-0.5)/M)/(maxN-minN+1)*100)),' %)'])
        [Xref,Yref,~,~,~,~] = eulerMethod(X0,startNT,NTmax,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,seed,PlotActive);       
        values(o) = 0.5*(max((X-Xref([1 1+2^diff*(1:length(X)-1)])).^2)+max((Y-Yref([1 1+2^diff*(1:length(X)-1)])).^2)); 
        %fprintf('%3.1f %%\n',100*((NT-minN)+o/M)/(maxN-minN+1))
        waitbar(((NT-minN)+o/M)/(maxN-minN+1),f,['Computing Monte-Carlo simulation... (',num2str(floor(((NT-minN)+o/M)/(maxN-minN+1)*100)),' %)'])
    end
    
    % Expectation and unbiased variance
    meanvalues = mean(values);
    expectations(NT - minN + 1) = meanvalues;
    var(NT - minN + 1) = M/(M+1)*(mean(values.^2)-meanvalues.^2);
end
figure

% Display and linear regression
Ns = minN:maxN;
plot(log(2.^Ns),log(expectations),'or')
grid on 
grid minor
hold on
xlabel('$\log(n)$','Interpreter','latex')
ylabel('$\log(E[{\sup}|X^{N,n}_t-X^{N,n_0}_t|^2])$','Interpreter','latex')
xlim([min(log(2.^Ns)) max(log(2.^Ns))])
%ylim([min(log(expectations-1.96*sqrt(var)/sqrt(M))) max(log(expectations+1.96*sqrt(var)/sqrt(M)))])
a = plot(log(2.^Ns),log(expectations+1.96*sqrt(var)/sqrt(M)),'--b');
plot(log(2.^Ns),log(expectations-1.96*sqrt(var)/sqrt(M)),'--b')
legend(a,'95% confidence interval')
title(['Error when $n$ varies estimated with Monte-Carlo method for ',num2str(M),' paths'],'Interpreter','latex') 
[beta0,beta1] = linearRegression(log(2.^Ns)',log(expectations)');
order = - beta1
cont = linspace(min(log(2.^Ns)),max(log(2.^Ns)),1000);
plot(cont,beta0+beta1*cont)
delete(f)

end