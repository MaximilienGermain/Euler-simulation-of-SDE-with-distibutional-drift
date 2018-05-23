function plotUsualSolution(t,X,NT,N,T,Xeuler)

figure
modifiedEuler = plot(t,X,'r') ;
grid on
grid minor
hold on
xlim([0 T])
euler = plot(t,Xeuler,'b');
xlabel('$t$','Interpreter','latex')
ylabel('$X_t$','Interpreter','latex')
chn = ['Approximation of a sample path of the SDE solution ($n =\ $',num2str(NT),' ; $N =\ $',num2str(N),')'];
title(chn,'Interpreter','latex')
legend([modifiedEuler,euler],'Approximation with Haar wavelets','Usual Euler scheme')
    
end