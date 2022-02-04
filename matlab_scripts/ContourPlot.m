clear all


%Plug in the ranges for x and t
x = linspace(-5,5);
t = linspace(0,2);

[X,T] = meshgrid(x,t);

%Enter the function u(x,t) with capital X and T
U = exp(-(x-2*t).^2);

contourf(X,T,U,10)

xlabel('Position x','FontSize',14);
ylabel('Time t','FontSize',14);
