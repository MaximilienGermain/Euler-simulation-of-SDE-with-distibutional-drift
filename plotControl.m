function frame = plotControl(t,X,NT,N,T,W)

figure
n = length(t);
modifiedEuler = plot(t,X-W(1:n)','b') ; %%%%% Test if there isn't saturation
grid on
grid minor
xlim([0 T])
xlabel('$t$','Interpreter','latex')
ylabel('$X_t-W_t$','Interpreter','latex')
chn = ['Approximation of the drift part of the SDE solution ($n =\ $',num2str(NT),' ; $N =\ $',num2str(N),')'];
title(chn,'Interpreter','latex')
%name = ['drift part n = ',num2str(NT),' ; N = ',num2str(N),'.png'];
%saveas(gcf,name);
frame = getframe(gcf);
    
end