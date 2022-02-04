
clear all


%Plug in the ranges for x and t 
%Leave the step size = 0.1 for a clear image
[X,T] = meshgrid(-8:0.1:8, 1:0.1:10);  

%Enter the function u(x,t) with capital X and T
U = (1/2)*sin(X-T) + (1/2)*sin(X+T);

surf(X,T,U)

xlabel('Position x','FontSize',14);
ylabel('Time t','FontSize',14);
zlabel('u(x,t)');