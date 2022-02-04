
clear all


% For the moving plot, you'll have to enter the time range (in t)
% and the function *twice* !


i = 1;

%Set the domain (in x)
%Leave the step size = 0.01 for a clear movie 
x = -8:0.01:8;

%Set the time range in t (leave the step size = 0.1)
for t = 0:0.1:10
    
    %Enter the function u(x,t) with lowercase x and t (only change the RHS)
    U(i,:) = 0.5*cos(x-t)+0.5*cos(x+t);
    
    i = i+1;
end

y1 = min(min(U));
y2 = max(max(U));

%Set the time range in t again (leave the step size = 0.1)
for t = 0:0.1:10 
    
  %Enter the function u(x,t) with lowercase x and t (only change the RHS)
  u = 0.5*cos(x-t)+0.5*cos(x+t);
    

  plot(x, u);
  title( sprintf('t = %.1f', t) );
  xlabel('Position x','FontSize',14);
  
  xlim([x(1) x(end)])
  ylim([y1-0.5 y2+0.5])
  pause( 0.1 );
end

