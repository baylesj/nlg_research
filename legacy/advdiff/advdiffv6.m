%%
% Jordan Bayles
% URISC/Spring 2012 Research Project
% Model based on advection-diffusion equation, expressed as
% \[ v_{t} - \kappa H_A G_{\infty} v_{zz} = \kappa H_A \dot{Q} \]

%% Physical parameters
h = 1.28*10^(-3);					% m, may need conversion
kap = 1.72*10^(-5);					% m^4/(Ns), may need conversion
ha = 1.29*10^(6);					% Pa, may need conversion

z0 = 0;								% start of computational domain
zf = h;								% end of computational domain
%T=5e3*h^2/kap/ha					% end time
%T=3600;
N=100000;
M = 4;								% number of nodes (minus boundary)
dz = (zf-z0)/(M+1);					% spatial step

c = 1;								% From Hader paper
ginfty = 1;							% From our notes
gd = c;								% From our notes

D = kap*ha*ginfty;					% diffusion coefficient

%% Characterization of tau, the relaxation time
% In the case of tau uniform, equal chances from tau_1 to tau_2
% Stuebner Haider gives cases tau_1=.01,.1,1.0 and tau_2=10,100,200
% From earlier, we use tau in the definition of G(t), and thus tau_1, tau_2
% should be used as the limits of integration of \int g(t, tau) dF(tau).
% However, later we say tau Pdot + P = epsilon*E, where tau (~F) is the uniform
% distribution, and is given by r xi + m, xi ~ U[-1, 1].
tau_1 = .01;
tau_2 = 200;
%tau = tau_1 + (tau_2 - tau_1).*rand(1, 1)
tau = 100

% Time-related variables are placed here.
sf = 0.9;							% safety factor (Mingham)
% The *proper* theoretical time step is below (really slow...)
dt = sf*(dz*dz/(2*D))			% time step (how formulated? cfl condition)
%dt=tau/10
%dt = T/100;
%N = ceil(T/dt)						% number of time steps
T = dt*N							% rescale end time to be multiple of dt

% t0, the time at which we switch from "ramp" to "hold" for the boundary condition,
% is defined in Haider, Schutgart based on the condition u0 = v0t0 and u0 = h/10.
% Discrepancy between Stuebner/Haider and Schutgart/Haider... Stubner paper says
% that v0 = h/10, specifies t0 = 10 seconds, tf = 3600 seconds. Due to clarity of
% the paper, I will go with this definition for now. TODO check with Gibson.
v0 = h/10;
t0 = 10;

z = z0:dz:zf;						% spatial nodes
t = dt:dt:T;						% time steps

gam = dt/(dz*dz);					% constant for diffusion term
A = (diag(ones(M+1,1),-1)-2*diag(ones(M+2,1))+diag(ones(M+1,1),1));

% Based off of this A, we edit some values to get correct boundaries.
A(M+2, M+1) = 2;
A(1,2) = 0;							% Clear the 1 stored here from earlier gen.
A(1,1) = 0;
%% Initialize solution vectors (including boundary)
vold = zeros(M+2,1);
uold = vold;
u = vold;
qn = vold;
v = vold;

% The initial value for the ramp and hold function is set here.
vold(1)=v0;

%% Before beginning the actual process, check if time step is OK
reply = input(sprintf('N = %d. Do you want to continue? Y/N [Y]: ', N), 's');
if isempty(reply)
    reply = 'Y'; % Default to yes
end

if reply == 'Y'
for n = 1:N,

  % Pseudo-implementation of Kronecker delta, not available in MATLAB < 2012.
  % We want this evaluated at each running time. Note that function handles
  % are *not* supported in standard matrices, hence we set it here.
  if abs(n*dt-t0) < dt/2
    A(1,1) = -1/gam/D;
  end
%  A(1,1) = -(abs((n*dt-t0)) < dt/2)/gam;

  vzz = gam*A*vold; %this is actually dt*vzz

  % Update Qdot based upon the values of vzz (second order derivative of v, new).
  qdot = gd*vzz/tau-qn*(dt/tau); %this is actually dt*qdot
  %qdot1=qdot(1) % should be zero?

  % Update the value of Q_{n+1} based on Qdot.
  qnp1 = qn + qdot; %since qdot is actually dt*qdot
  
  % Finally, update v based on old copies of qdot, qnp1.
  % From 5/11 notes, we use:
  % (vnp1 - vn)/dt = D(v_{i-1} - 2v_i + v_{i+1})/dz^2 = DA \vec{v}/dz^2
  % v = vold+dt*D*A*v/dz^2+dt*kap*ha*qdot;
  v = vold+D*vzz+kap*ha*qdot; %since qdot is actually dt*qdot
                                  % and vzz is actually dt*vzz

  % First order approximation to u, the time integral of v.
  u = uold+dt/2*(v+vold);
  % Here we finish the step by moving current values to old ones.
  qn = qnp1;
  vold = v;
  uold = u;

  %sprintf('Value of A(1,1) = %d', A(1,1))
  %sprintf('Value of V(1,1) = %d', v(1,1))
  %sprintf('Value of n = %d', n)
  %sprintf('Value of dt = %d', dt)
  %sprintf('Value of t0 = %d', t0)
  
  if(mod(n,floor(N/10))==0)
    subplot(1,2,1), plot(z,v,'-o')
    title(sprintf('V (or du/dt), with n=%d',n))
    %az = [0 h 0 h/20];
    %axis(az)
    subplot(1,2,2), plot(z,u,'-o')
    title(sprintf('Displacement, with n=%d', n))
    %at = [0 N 0 N/20];
    %axis(at);

disp('hi')
pause
  end
  
end
end
