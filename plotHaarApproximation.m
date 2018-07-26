function frame = plotHaarApproximation(K,Mu,N,B,xgrid,H)

figure 
dispgrid=linspace(-K,K,max(2*K*2^(N+2),2*K*30));  
approx = b(Mu,K,dispgrid);
max(abs(approx))
[AX,H1,H2] = plotyy(dispgrid,approx,dispgrid,continuous(B,xgrid,dispgrid));
%H1.LineStyle = '--';
ylabel(AX(1),'Approximation of $\frac{\mathrm{d}}{\mathrm{d} x}B^H_x$ by Haar wavelets','Interpreter','latex')
ylabel(AX(2),'$B^H_x$','Interpreter','latex')
set(AX(1),'ycolor','r')
set(AX(2),'ycolor','b')
set(H1,'Color','r')
set(H2,'Color','b')
set(AX(2),'Xgrid','on','Ygrid','on','YMinorGrid','on')
xlabel('$x$','Interpreter','latex')
xlim(AX(2), [min(xgrid) max(xgrid)])
xlim(AX(1), [min(xgrid) max(xgrid)])
ylim(AX(1), [min(approx)-0.2*abs(min(approx)) max(approx)+0.2*abs(min(approx))])
ylim(AX(2), [min(B)-0.2*abs(min(B)) max(B)+0.2*abs(min(B))])
chn = ['Sample path of $B^H_x$ and approximation of $\frac{\mathrm{d}}{\mathrm{d} x}B^H_x$ ($H = ',num2str(H),'\ ;\ N = ',num2str(N),'$)'];
title(chn,'Interpreter','latex')
name = ['Haar H = ',num2str(H),' ; N = ',num2str(N),'.eps'];
set(gcf,'PaperPositionMode','auto')
print(name,'-depsc','-tiff')
frame = getframe(gcf);
    
end