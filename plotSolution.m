function frame = plotSolution(t,X,NT,N,T)

figure
modifiedEuler = plot(t,X,'b');
grid on
grid minor
xlim([0 T])
ylim([-3 2])
xlabel('$t$','Interpreter','latex')
ylabel('$X_t$','Interpreter','latex')
chn = ['Approximation of a sample path of the SDE solution ($n =\ $',num2str(NT),' ; $N =\ $',num2str(N),')'];
title(chn,'Interpreter','latex')
%name = ['n = ',num2str(NT),' ; N = ',num2str(N),'.png'];
%saveas(gcf,name);
frame = getframe(gcf);
    
end