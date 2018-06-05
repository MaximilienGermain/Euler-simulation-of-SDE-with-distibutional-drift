delta = 0.000000000000000000001;
b = 1;
c = 1;
epsilon = (1:1000)/1000;
plot(epsilon,delta + b*epsilon.^(-0.1) + c*epsilon + delta./epsilon)
hold on
%plot(epsilon,delta + c*epsilon + (b+ delta)./epsilon)