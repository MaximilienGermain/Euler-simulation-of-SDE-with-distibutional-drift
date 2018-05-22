function [xgrid,B,M] = createfBm(H,K,N,start,Nx,A,seed)

[random,xgrid,B,M] = fbm(H,K,start,Nx,A,seed);

% figure
% plot(xgrid,B)
% grid on
% grid minor
% xlabel('$x$','Interpreter','latex')
% xlim([min(xgrid) max(xgrid)])
% ylabel('$B_x^H$','Interpreter','latex')
% chn = ['Sample path of a fractional Brownian motion $B^H_x$ (N = ',num2str(start),' ; H = ',num2str(H),')'];
% title(chn,'Interpreter','latex')

Ninit = 1+K*2^(start+2);

for p=1:N-start
    [random,xgrid,B,M] = refinefbm(random,xgrid,H,start+p,M,Nx,p,seed); 
%     figure
%     plot(xgrid,B)
%     grid on
%     grid minor
%     xlabel('$x$','Interpreter','latex')
%     xlim([min(xgrid) max(xgrid)])
%     ylabel('$B_x^H$','Interpreter','latex')
%     chn = ['Sample path of a fractional Brownian motion $B^H_x$ (N = ',num2str(start+p),' ; H = ',num2str(H),')'];
%     title(chn,'Interpreter','latex')
end 

end