%%
% Jordan Bayles
% URISC/Spring 2012 Research Project
% Model based on diffusion equation with fading memory, expressed as
% \[ v_{t} = \kappa H_A G_{\infty} v_{zz} + \kappa H_A \dot{Q} \]
% \[ Q_{t} = G_d/\tau v_{zz} - 1/\tau Q \]
% subbing Q_t into the first equation we get
% \begin{align}
% v_{t} &= \kappa H_A G_{\infty} v_{zz} + \kappa H_A (G_d/\tau v_{zz} - 1/\tau Q) \\
%       &= \kappa H_A (G_{\infty}+G_d/\tau) v_{zz} - \kappa H_A/\tau Q
% \end{align}
% Discretizing with \vec{v}=v(\vec(x)) and \vec(Q)=Q(\vec(x)), and Av\approx v_{zz}
% we can rewrite as a system ($h=\Delta z$)
% \dot{\begin{bmatrix}\vec{v}\\\vec{Q}\end{bmatrix}} = 
% \begin{bmatrix}\kappa H_A (G_{\infty}+G_d/\tau)A/h^2 &  - \kappa H_A/\tau I \\
% G_d/\tau A/h^2 &  - 1/\tau I\end{bmatrix}

%% Physical parameters
h = 1.28*10^(-3);                    % m, may need conversion
kap = 1.72*10^(-5);                    % m^4/(Ns), may need conversion
ha = 1.29*10^(6);                    % Pa, may need conversion

z0 = 0;                                % start of computational domain
zf = h;                                % end of computational domain
%T=5e3*h^2/kap/ha                    % end time
%T=3600
T = 100;
M = 40;                                % number of nodes (minus boundary)
dz = (zf-z0)/(M+1);                    % spatial step
c = 1;                                % From Hader paper
ginfty = 1;                            % From our notes
gd = c;                                % From our notes


%% Characterization of tau, the relaxation time
tau_1 = .01;
tau_2 = 200;
%tau = tau_1 + (tau_2 - tau_1).*rand(1, 1)
tau = 200;

% Time-related variables are placed here.
dt = .25;
N = ceil(T/dt)                        % number of time steps
theta = 0.5;

% t0, the time at which we switch from "ramp" to "hold" for the boundary condition,
% is defined in Haider, Schutgart based on the condition u0 = v0t0 and u0 = h/10.
% Discrepancy between Stuebner/Haider and Schutgart/Haider... Stubner paper says
% that v0 = h/10, specifies t0 = 10 seconds, tf = 3600 seconds.
v0 = h/10;
t0 = 10;

z = z0:dz:zf;                        % spatial nodes
t = dt:dt:T;                        % time steps

A = (diag(ones(M,1),-1)-2*diag(ones(M+1,1))+diag(ones(M,1),1));
% Identity matrix of size similar to A (double check syntax)
I = eye(M+1, M+1);

% Based off of this A, we edit some values to get correct boundaries.
A(M+1, M) = 2;

%% Initialize solution vectors (including boundary)
vold = zeros(M+1,1);
uold = vold;
pold = vold;
u = vold;
p = vold;
v = vold;


% Create a matrix for storing the values of u over time.
uxt = zeros(M+1, N);

% Implicit additions
Am=I-kap*ha*dt*(ginfty+gd/tau)*theta/dz/dz*A...
  +kap*ha*gd*(dt*theta/tau/dz)^2/(1+dt*theta/tau)*A; %A_minus

Ap=I+kap*ha*dt*(ginfty+gd/tau)*(1-theta)/dz/dz*A...
  -kap*ha*gd*(dt*(1-theta)/tau/dz)^2/(1+dt*theta/tau)*A; %A_plus

AforP=(kap*ha*dt*(1-theta)/tau...
	  +kap*ha*dt*theta/tau*(1-dt*(1-theta)/tau)/(1+dt*theta/tau))/dz/dz*A;


%% Before beginning the actual process, check if time step is OK
reply = input(sprintf('N = %d. Do you want to continue? Y/N [Y]: ', N), 's');
if isempty(reply)
    reply = 'Y'; % Default to yes
end

if reply == 'Y'

p0=0; %initialize
e1=0*v;e1(1)=1;
v0 = h/10;
sv0=dt*kap*ha*(ginfty+gd/tau)/dz/dz*v0*e1; %scaled v0 at n+theta

for n = 1:N,

  if (n+theta)*dt > t0
    v0 = 0;
    sv0 = 0*e1;   
  end
  p0=gd*v0+(p0-gd*v0)*exp(-dt/tau);
  sp0=dt*kap*ha/tau/dz/dz*p0*e1; %scaled p0 at n+theta


  v=Am\(Ap*vold-AforP*pold+sv0-sp0);
  p=(gd*dt/tau*(theta*v+(1-theta)*vold)+(1-dt*(1-theta)/tau)*pold)/(1+dt*theta/tau);
  u = uold+dt/2*(v+vold);

  % Here we finish the step by moving current values to old ones.
  pold = p;
  vold = v;
  uold = u;
  
  if(mod(n,floor(N/5))==0)

    subplot(1,2,1), plot(z,[v0;v].','-o')
    title(sprintf('V (or du/dt), with n=%d',n))
    %av = [z0 zf 1.27777 1.281111];
    %axis(av)
    subplot(1,2,2), plot(z,[0;u].','-o')
    title(sprintf('Displacement, with n=%d', n))
    %at = [0 N 0 N/20];
    %axis(at);

    pause
  end

    uxt(:, n) = u;

end
end

% Uxt is now completed, so we can graph to compare with Haider/Schugart
figure()
%hold on
%for i = 1:(M+1)
%	  plot(uxt(i, :), '-o')
%end
%hold off
	   plot(t,uxt(round(M/2),:),'-o')
title(sprintf('Plot of normalized displacements over time'))
xlabel(sprintf('Value of N: %d, T: %d', N, T))