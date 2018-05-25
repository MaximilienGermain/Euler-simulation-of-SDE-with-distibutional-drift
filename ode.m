n=100000;
eps = logspace(-80,-1);
t = linspace(0,1,n+1);
dt = 1/n;
alpha = 0.01;
values = [];
x0 = 0.01;

for epsilon=eps
    
    %x0 = epsilon;
    x=zeros(n+1,1);
    x(1) = x0;

    for i=1:n
        x(i+1) = x(i) + (x(i) + x(i)^alpha)*dt;
    end
    
    values = [values max(x)];
end
plot(t,x)
% min(values)
% semilogx(eps,values)

[x,y] = meshgrid(0:0.1:2,0:0.1:2);
u = ones(size(y));
v = x + x.^alpha;
figure
quiver(x,y,u,v)
