%%
% Jordan Bayles
% URISC/Spring 2012 Research Project
% Model based on advection-diffusion equation, expressed as
% \[ v_{t} - \kappa H_A G_{\infty} v_{zz} = \kappa H_A \dot{Q} \]

%% Physical parameters
h = 1.28*10^(-3);                    % m, may need conversion
kap = 1.72*10^(-5);                    % m^4/(Ns), may need conversion
ha = 1.29*10^(6);                    % Pa, may need conversion

z0 = 0;                                % start of computational domain
zf = h;                                % end of computational domain
%T=5e3*h^2/kap/ha                    % end time
pow = 5;
N = 10^pow
M = 16;                                % number of nodes (minus boundary)
dz = (zf-z0)/(M+1);                    % spatial step

c = 1;                                % From Hader paper
ginfty = 1;                            % From our notes
gd = c;                                % From our notes

% TODO find out what epsilon should be. Do we really want du/dz?
epsilon = 1;

D = kap*ha*ginfty;                    % diffusion coefficient

% Addition of Crank-Nicolson to increase stability
% In the case that theta = 0, dt for stability is leq dz^2/2D
theta = 0.5;

%% Characterization of tau, the relaxation time
tau_1 = .01;
tau_2 = 200;
%tau = tau_1 + (tau_2 - tau_1).*rand(1, 1)
tau = 100

% Time-related variables are placed here.
sf = 0.9;                            % safety factor (Mingham)
dt = sf*(dz*dz/(2*D))            % time step (how formulated? cfl condition)
% From notes, dt = tau / 10. Different tau? unstable for tau = 100.
%dt = tau/10;
%N = ceil(T/dt)                        % number of time steps
T = dt*N                            % rescale end time to be multiple of dt

% t0, the time at which we switch from "ramp" to "hold" for the boundary condition,
% is defined in Haider, Schutgart based on the condition u0 = v0t0 and u0 = h/10.
% Discrepancy between Stuebner/Haider and Schutgart/Haider... Stubner paper says
% that v0 = h/10, specifies t0 = 10 seconds, tf = 3600 seconds. Due to clarity of
% the paper, I will go with this definition for now. TODO check with Gibson.
v0 = h/10;
t0 = 10;

z = z0:dz:zf;                        % spatial nodes
t = dt:dt:T;                        % time steps

gam = dt/(dz*dz);                    % constant for diffusion term
A = (diag(ones(M+1,1),-1)-2*diag(ones(M+2,1))+diag(ones(M+1,1),1));
% Identity matrix of size similar to A (double check syntax)
I = eye(M+2, M+2);

% Based off of this A, we edit some values to get correct boundaries.
A(M+2, M+1) = 2;
A(1,2) = 0;                            % Clear the 1 stored here from earlier gen.
A(1,1) = 0;

%% Initialize solution vectors (including boundary)
vold = zeros(M+2,1);
uold = vold;
u = vold;
qn = vold;
v = vold;

% Create a matrix for storing the values of u over time.
res = 10;
index = 1;
uxt = zeros(M+2, N);

% The initial value for the ramp and hold function is set here.
vold(1)=v0;

%% Before beginning the actual process, check if time step is OK
reply = input(sprintf('N = %d. Do you want to continue? Y/N [Y]: ', N), 's');
if isempty(reply)
    reply = 'Y'; % Default to yes
end

%% Polynomial Chaos modifications, addl. variables.
% TODO: Check these values!
r = 99.9;
m = 100;
p = ?
ai = zeros(p, 1);
ainp1 = ai;
aM = zeros(p, M);
B = (r*M + m*I)

if reply == 'Y'
for n = 1:N,
  % Pseudo-implementation of Kronecker delta, not available in MATLAB < 2012.
  % We want this evaluated at each running time. Note that function handles
  % are *not* supported in standard matrices, hence we set it here.
  if abs(n*dt-t0) < dt/2
    A(1,1) = -1/gam/D;
  else % Not sure if this case is needed, I am thinking so
    A(1,1) = 0;
  end

  vzz = gam*A*vold; %this is actually dt*vzz
  
  % Polynomial Chaos goes here
  f = D*p;
  for iz = 1:M
    ai = aM(:, iz);
    f(1) = vzz(iz);
    ainp1= I*ai - dt*B\ai + epsilon*dt*B\f;
    aM(:, iz) = ainp1;
  end
  qn = qnp1;
  qnp1 = aM(1,:);

  %qnp1 = (1-dt/tau)*qn + epsilon/tau*vzz;
  % Update qdot (really dt*Qdot at the n+1/2 time step)
  qdot = (qnp1-qn);

  % Finally, update v based on old copies of qdot, qnp1.
  % Note that x = A\b is the same as x = inv(A)*b, except faster and more accurate.
  v = (I-theta*D*dt*A)\(vold + (1-theta)*D*vzz + qdot);
  % First order approximation to u, the time integral of v.
  u = uold+dt/2*(v+vold);
  % Here we finish the step by moving current values to old ones.
  qn = qnp1;
  vold = v;
  uold = u;
   
  if(1 == 2)
  %if(mod(n,floor(N/res))==0)
    %unxt(:, index) = u/h;
    %index = index + 1;
    
    subplot(1,2,1), plot(z,v,'-o')
    title(sprintf('V (or du/dt), with n=%d',n))
    %av = [z0 zf 1.27777 1.281111];
    %axis(av)
    subplot(1,2,2), plot(z,u,'-o')
    title(sprintf('Displacement, with n=%d', n))
    %at = [0 N 0 N/20];
    %axis(at);

    pause
  end

    uxt(:, index) = u/h;
    index = index + 1;
end
end

% Uxt is now completed, so we can graph to compare with Haider/Schugart
figure()
hold on
for i = 1:(M+2)
  plot(uxt(i, :), '-o')
end
hold off
title(sprintf('Plot of normalized displacements over time'))
xlabel(sprintf('Value of N: %d, T: %d', N, T))
uxtdiff = zeros(M+2, res-1);

