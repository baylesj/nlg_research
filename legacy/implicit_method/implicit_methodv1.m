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
T=3600

M = 16;                                % number of nodes (minus boundary)
dz = (zf-z0)/(M+1);                    % spatial step
c = 1;                                % From Hader paper
ginfty = 1;                            % From our notes
gd = c;                                % From our notes

D = kap*ha*ginfty;                    % diffusion coefficient

%% Characterization of tau, the relaxation time
tau_1 = .01;
tau_2 = 200;
%tau = tau_1 + (tau_2 - tau_1).*rand(1, 1)
tau = 200;

% Time-related variables are placed here.
dt = .25;
N = ceil(T/dt)                        % number of time steps
theta = 0.5;
index = 1;

% t0, the time at which we switch from "ramp" to "hold" for the boundary condition,
% is defined in Haider, Schutgart based on the condition u0 = v0t0 and u0 = h/10.
% Discrepancy between Stuebner/Haider and Schutgart/Haider... Stubner paper says
% that v0 = h/10, specifies t0 = 10 seconds, tf = 3600 seconds.
v0 = h/10;
t0 = 10;

z = z0:dz:zf;                        % spatial nodes
t = dt:dt:T;                        % time steps

gam = dt/(dz*dz);                    % constant for diffusion term

A = (diag(ones(M+1,1),-1)-2*diag(ones(M+2,1))+diag(ones(M+1,1),1));
% Identity matrix of size similar to A (double check syntax)
I = eye(M+1, M+1);

% Based off of this A, we edit some values to get correct boundaries.
A(M+2, M+1) = 2;

% TODO clean this up
A = A(2:end, 2:end);

%% Initialize solution vectors (including boundary)
vold = zeros(M+1,1);
uold = vold;
qold = vold;
u = vold;
q = vold;
v = vold;

w = [q; v];
wold = [qold; vold];

% Create a matrix for storing the values of u over time.
uxt = zeros(M+1, N);
% The initial value for the ramp and hold function is set here.
vold(1)=v0;

% Implicit additions
Aimp = [ (1+theta/tau)*I, -dt*gd*theta/dz/dz/tau*A;
          dt/tau*theta*I, I-dt*gam*(D+gd/tau/dt)*A ];

Bimp = [ (1-dt/tau*(1-theta))*I, gd*dt/tau/dz/dz*(1-theta)*A;
         -dt/tau*(1-theta)*I, I+dt*gam*(D+gd/tau/dt)*(1-theta)*A ];

% Store the LU factorization of Aimp for later
[Limp, Uimp] = lu(Aimp);

%% Before beginning the actual process, check if time step is OK
reply = input(sprintf('N = %d. Do you want to continue? Y/N [Y]: ', N), 's');
if isempty(reply)
    reply = 'Y'; % Default to yes
end

if reply == 'Y'

zhalf = (z0+dz/2):dz:(zf-dz/2)
for n = 1:N,
  % Pseudo-implementation of Kronecker delta, not available in MATLAB < 2012.
  if n*dt < t0
    v0 = h/10;
  else 
    v0 = 0;
  end

  w = Uimp\(Limp\Bimp*wold);
  q = w(1:M+1);
  v = w(M+2:end);
  u = uold+dt/2*(v-vold);
  % Here we finish the step by moving current values to old ones.
  qold = q;
  vold = v;
  uold = u;
  wold = w;
  
  if(mod(n,floor(N/2))==0)

    subplot(1,2,1), plot(zhalf,v.','-o')
    title(sprintf('V (or du/dt), with n=%d',n))
    %av = [z0 zf 1.27777 1.281111];
    %axis(av)
    subplot(1,2,2), plot(zhalf,u.','-o')
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
for i = 1:(M+1)
  plot(uxt(i, :), '-o')
end
hold off
title(sprintf('Plot of normalized displacements over time'))
xlabel(sprintf('Value of N: %d, T: %d', N, T))
