function [xgrid,B,M] = createfBm(H,K,N) % Crashes if N >= start + 2

start = 4;
[random,xgrid,B,M] = fbm(H,K,start);
% figure
% plot(xgrid,B)
% grid on
% grid minor
% xlabel('x')
% xlim([min(xgrid) max(xgrid)])
% ylabel('B_x^H')
% chn = ['Sample path of a fractional Brownian motion B^H_x (N = ',num2str(start),' ; H = ',num2str(H),')'];
% title(chn)
for p=1:N-start
    [random,xgrid,B,M] = refinefbm(random,xgrid,H,start,M);
    min(eig(M))
    max(eig(M))
%     figure
%     plot(xgrid,B)
%     grid on
%     grid minor
%     xlabel('x')
%     xlim([min(xgrid) max(xgrid)])
%     ylabel('B_x^H')
%     chn = ['Sample path of a fractional Brownian motion B^H_x (N = ',num2str(start+p),' ; H = ',num2str(H),')'];
%     title(chn)
%     end    
end