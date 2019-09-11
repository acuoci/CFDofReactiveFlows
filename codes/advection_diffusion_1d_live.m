%% Advection-diffusion equation in 1D
%% 
% %% 1. Case description
% We want to numerically solve with the Finite Difference (FDM) method the *advection-diffusion 
% equation in 1D* on a domain with length $L$:
% 
% $$\frac{\partial f}{\partial t}+u\frac{\partial f}{\partial x}=\Gamma \frac{\partial^2 
% f}{\partial x^2 }$$
% 
% with *periodic boundary conditions* $f\left(x=0,\forall t\right)=f\left(x=L,\forall 
% t\right)$ and initial condition $f\left(t=0,\forall x\right)=A\;\mathrm{sin}\left(2\pi 
% k\;x\right)$, where $A$ (amplitude) and $k$ (wave-number) are user-defined parameteres.
% 
% The chosen initial and boundary conditions correspond to the following analytical 
% solution: 
% 
% $$a\left(x,t\right)=A\;\mathrm{exp}\left(-4\pi^2 k^2 \Gamma \;t\right)\;\mathrm{sin}\left(2\pi 
% k\left(x-u\;t\right)\right)$$
% 
% In particular, we are going to use the *explicit (or forward) Euler method* 
% in time, while for spatial discretization of advective and diffusive terms we 
% are going to use *centered, 2nd order discretizations*. We assume constant, 
% uniform velocity $u$ and diffusion coefficient $\Gamma$. The discretized equation 
% becomes:
% 
% $$\frac{f_i^{n+1} -f_i^n }{\Delta t}+u\frac{f_{i+1}^n -f_{i-1}^n }{2h}=\Gamma 
% \frac{f_{i+1}^n -2f_i^n +f_{i-1}^n }{h^2 }$$
% 
% where $n$ is the time level and $i$ the grid point.
% 
% We want also to monitor the error between the analytical $a\left(x,t\right)$ 
% and the numerical solution $f\left(x,t\right)$:
% 
% $$E={\left\|a\left(x,t\right)-f\left(x,t\right)\right\|}_2$$
% 
% where ${\left\|\mathit{\mathbf{v}}\right\|}_2$ indicates the Euclidean norm 
% of vector $\mathit{\mathbf{v}}$. 
% 
% We also monitor the squared areas below the analytical and numerical curves, 
% to better understand the differences between the analytical and numerical curves:
% 
% $$a_{\mathrm{int}}^2 \left(t\right)=\int_0^L a^2 \left(x,t\right)\mathrm{dx}$$
% 
% $$f_{\mathrm{int}}^2 \left(t\right)=\int_0^L f^2 \left(x,t\right)\mathrm{dx}$$
% 
% From the theory, we know that we can have stable solutions only if the following 
% constraints on the Courant number $Co$ and the Diffusion number $Di$ are satisfied:
% 
% $$\begin{array}{l}\mathrm{Co}=\frac{u\Delta t}{h}<1\\\mathrm{Di}=\frac{\Gamma 
% \Delta t}{h^2 }<\frac{1}{2}\end{array}$$
% 
% This means that the time step we can consider must satisfy the following expression:
% 
% $$\Delta t<\mathrm{min}\left(1\frac{h}{u},\frac{1}{2}\frac{h^2 }{\Gamma }\right)$$
%% 2. Implementation
% Numerical setup
% We start with the definition of the input data: $n$ (number of grid points), 
% $n_{step}$ (number of time steps), $\Delta t$ (time step), $u$ (velocity), $L$ 
% (length), $\Gamma \;$(diffusion equation):

np=21;              % number of grid points                
nstep=100;          % number of time steps
L=2.0;              % domain length [m]
dt=0.05;            % time step [s]
u=1;                % velocity [m/s]
D=0.05;             % diffusion coefficient [m2/s]
%% 
% We also define the 2 parameters for the initial condition:

A=0.5;              % amplitude of initial solution
k=1;                % wave number [1/m]
%% 
% We calculate now the spatial step $h$:

h=L/(np-1);         % grid step [m]
%% 
% In order to ensure a good computational efficiency, it is convenient to pre-allocate 
% memory for the relevant vectors $y$ (new solution at time $n+1$),$f$(old solution 
% at time $n$) and $a$ (the analytical solution):

fo=zeros(np,1);     % temporary numerical solution
f=zeros(np,1);      % current numerical solution
a=zeros(np,1);      % exact solution
%% 
% Let's define the initial solution:

for i=1:np
	f(i)=A*sin(2*pi*k*h*(i-1)); 
end
%% 
% We can now plot the initial solution:

plot(0:h:L, f);
%% 
% It is convenient to check if the stability constraints are satisfied:

Co = u*dt/h;                        % Courant number
Di = D*dt/h^2;                      % Diffusion number
dt_max = min(1*h/u, 0.5*h*h/D);     % Maximum allowed time step
fprintf('Co=%f, Di=%f, dt=%f, dt(max)=%f\n', Co, Di, dt, dt_max);
% b. Advancing the solution in time
% We can now proceed with the numerical solution, advancing in time, starting 
% from $t=0$. Let's apply the forward Euler method over all the internal points:
% 
% $$f_i^{n+1} =f_i^n -\frac{u\Delta t}{2h}\left(f_{i+1}^n -f_{i-1}^n \right)+\frac{\Gamma 
% \Delta t}{h^2 }\left(f_{i+1}^n -2f_i^n +f_{i-1}^n \right)$$
% 
% The corresponding code is reported in the following:
%%
% 
%   fo=f;      
%   for i=2:np-1 
%   	f(i) = fo(i)-(u*dt/2/h)*(fo(i+1)-fo(i-1))+...      % advection
%   		   D*(dt/h^2)*(fo(i+1)-2*fo(i)+fo(i-1));   % diffusion
%   end 
%
%% 
% Let's update the boundary conditions. In particular, we apply the same discretization 
% formula above on the last point of the grid, by considering that the point on 
% the right side is point $2$, because of periodicity:
% 
% $$f_N^{n+1} =f_N^n -\frac{u\Delta t}{2h}\left(f_2^n -f_{N-1}^n \right)+\frac{\Gamma 
% \Delta t}{h^2 }\left(f_2^n -2f_N^n +f_{N-1}^n \right)$$
% 
% The corresponding code is reported in the following:
%%
% 
%   f(np) = fo(np)-(u*dt/2/h)*(fo(2)-fo(np-1)) + D*(dt/h^2)*(fo(2)-2*fo(np)+fo(np-1)); 
%   f(1)  = f(np);
%
%% 
% Let's update the analytical solution:
%%
% 
%   for i=1:np 
%   	a(i) = A*exp(-4*pi*pi*k*k*D*t)*sin(2*pi*k*(h*(i-1)-u*t)); 
%   end  
%
%% 
% In order to compare the analytical and numerical solutions it is convenient 
% to define some integral quantities, for example the squared area below the curves 
% (which is a measure of the energy of the system):
%%
% 
%   a2_int = 0.;
%   f2_int = 0.;
%   for i=1:np-1
%        a2_int = a2_int + h/2*(a(i)^2+a(i+1)^2);
%        f2_int = f2_int + h/2*(f(i)^2+f(i+1)^2);
%   end 
%
%% 
% We also want to compare the error between the numerical and analytical solution 
% through the norm-2 of their difference:
%%
% 
%   E = 0;
%   for i=1:np 
%        E = E + (f(i)-a(i))^2;
%   end
%   E = h*sqrt(E); 
%
%% 
% We are now ready to move to the next time level:
%%
% 
%   t = t+dt;
%
% c. Graphical output
% Before starting the calculation, we prepare and open a video stream to register 
% the evolution of the solution in time:
%%
% 
%   video_name = 'advection_diffusion_1d.mp4';
%   videompg4 = VideoWriter(video_name, 'MPEG-4');
%   open(videompg4);
%
%% 
% At each time step we can update the plots showing the current numerical solution 
% (compared with the analytical profile):
%%
% 
%   message = sprintf('time=%d\na^2(int)=%d\ny^2(int)=%d', t, a2_int, f2_int);
%   hold off; plot([0:h:L],f,'linewidth',2); axis([0 L -1, 1]); % plot num. 
%   hold on; plot([0:h:L],a,'r','linewidth',2);                 % plot exact
%   hold on; legend('numerical', 'exact');                      % legend
%   xlabel('spatial coordinate [m]');
%   ylabel('solution');    
%   time = annotation('textbox',[0.15 0.8 0.1 0.1],'String',message,'EdgeColor','none');
%   frame = getframe(gcf);
%   writeVideo(v,frame);
%   delete(time);
%
%% 
% At the end of the calculations we have to close the video stream:
%%
% 
%   close(videompg4);
%
%% 3. Bringing it all together
% Let's bring it all together.

% Cleaning the MATLAB environment
%-------------------------------------------------------------------------%
close all;
clear variables;

% User-defined data
%-------------------------------------------------------------------------%
np=21;              % number of grid points                
nstep=100;          % number of time steps
L=2.0;              % domain length [m]
dt=0.05;            % time step [s]
u=1;                % velocity [m/s]
D=0.05;             % diffusion coefficient [m2/s]
A=0.5;              % amplitude of initial solution
k=1;                % wave number [1/m]

% Pre-processing of user-defined data
%-------------------------------------------------------------------------%
% Grid step calculation
h=L/(np-1);         % grid step [m]

% Memory allocation
fo=zeros(np,1);     % temporary numerical solution
f=zeros(np,1);      % current numerical solution
a=zeros(np,1);      % exact solution

% Initial solution
for i=1:np
	f(i)=A*sin(2*pi*k*h*(i-1)); 
end

% Check the stability conditions on time step
Co = u*dt/h;                        % Courant number
Di = D*dt/h^2;                      % Diffusion number
dt_max = min(1*h/u, 0.5*h*h/D);     % Maximum allowed time step
fprintf('Co=%f, Di=%f, dt=%f, dt(max)=%f\n', Co, Di, dt, dt_max);

% Video setup
%-------------------------------------------------------------------------%
video_name = 'advection_diffusion_1d.mp4';
videompg4 = VideoWriter(video_name, 'MPEG-4');
open(videompg4);

% Advancing in time
%-------------------------------------------------------------------------%
t = 0.;
for m=1:nstep
    
    % Update the analytical solution
    for i=1:np 
		a(i) = A*exp(-4*pi*pi*k*k*D*t)*sin(2*pi*k*(h*(i-1)-u*t)); 
    end  
    
    % Squared areas below the analytical and numerical solutions
    a2_int = 0.;
    f2_int = 0.;
    for i=1:np-1
         a2_int = a2_int + h/2*(a(i)^2+a(i+1)^2);
         f2_int = f2_int + h/2*(f(i)^2+f(i+1)^2);
    end  
    
    % Graphical output
    message = sprintf('time=%d\na^2(int)=%d\ny^2(int)=%d', t, a2_int, f2_int);
    hold off; plot(0:h:L,f,'linewidth',2); axis([0 L -1, 1]); % plot num. 
    hold on; plot(0:h:L,a,'r','linewidth',2);                 % plot exact
    hold on; legend('numerical', 'exact');                    % legend
    xlabel('spatial coordinate [m]');
    ylabel('solution');    
    time = annotation('textbox',[0.15 0.8 0.1 0.1],'String',message,'EdgeColor','none');
    frame = getframe(gcf);
    writeVideo(videompg4,frame);
    delete(time);
    
    % Forward Euler method
    fo=f;   
    for i=2:np-1 
		f(i) = fo(i)-(u*dt/2/h)*(fo(i+1)-fo(i-1))+...      % advection
			   D*(dt/h^2)*(fo(i+1)-2*fo(i)+fo(i-1));   % diffusion
    end 
    
    % Periodic boundary condition
    f(np) = fo(np)-(u*dt/2/h)*(fo(2)-fo(np-1))+...
            D*(dt/h^2)*(fo(2)-2*fo(np)+fo(np-1)); 
    f(1)  = f(np);
    
    % Update the error between numerical and analytical solution
    E = 0;
    for i=1:np 
        E = E + (f(i)-a(i))^2;
    end
    E = h*sqrt(E); 
    
    % New time step
    t=t+dt; 
    
    % Print the current time (every 25 steps)
    if (mod(m,25)==1), fprintf('time=%d E=%e\n', t, E); end
end

% Closing the video stream
close(videompg4);
% Play the video plotting the solution

implay(video_name);
%% 4. Exercises
%% 
% # Demonstrate that the time discretization error decreases linearly with the 
% time step adopting the same procedure used for studying the spatial accuracy.
% # Analyze the behavior of the solution by considering higher Courant and/or 
% Diffusion numbers. What happens when they are larger than 1?