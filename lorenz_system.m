%Lorenz system from the Steve Brunton video
%link: https://www.youtube.com/watch?v=EnsB1wP3LFM

Beta = [10; 28; 8/3]; %chaotic values
x0 = [0; 1; 20]; %initial condition
dt = 0.001; %resolution
timeSpan = dt:dt:50; %time vector

options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));
[t,x] = ode45(@(t,x)lorenz(t,x,Beta),timeSpan,x0, options);

plot3(x(:,1), x(:,2), x(:,3));