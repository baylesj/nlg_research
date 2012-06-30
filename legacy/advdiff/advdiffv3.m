%%
% Jordan Bayles
% URISC/Spring 2012 Research Project
% Model based on advection-diffusion equation, expressed as
% \[ v_{t} - \kappa H_A G_{\infty} v_{zz} = \kappa H_A \dot{Q} \]

%% Physical parameters
h=1.28*10^(-3);          			% m, may need conversion
kap=1.72*10^(-5); 					% m^4/(Ns), may need conversion
ha=1.29*10^(6);						% Pa, may need conversion

z0=0;								% start of computational domain
zf=h;								% end of computational domain
%T=5e3*h^2/kap/ha					% end time
T = 1e-7
M=64;								% number of nodes (minus boundary)
dz=(zf-z0)/(M+1);					% spatial step

c=1;								% From Hader paper
ginfty=1;							% From our notes
gd=c;								% From our notes

v=0;								% advection velocity
D=kap*ha*ginfty						% diffusion coefficient

%% Characterization of tau, the relaxation time
% In the case of tau uniform, equal chances from tau_1 to tau_2
% Stuebner Haider gives cases tau_1=.01,.1,1.0 and tau_2=10,100,200
% From earlier, we use tau in the definition of G(t), and thus tau_1, tau_2
% should be used as the limits of integration of \int g(t, tau) dF(tau).
% However, later we say tau Pdot + P = epsilon*E, where tau (~F) is the uniform
% distribution, and is given by r xi + m, xi ~ U[-1, 1], this interpretation
% is given above, but probably incorrect still.
tau_1 = .01
tau_2 = 200
tau = tau_1 + (tau_2 - tau_1).*rand(1, 1)

%% Time considerations
sf=0.9;								% safety factor (Mingham)
dt=sf*(dz*dz/(v*dz+2*D))			% time step (how formulated?)
N=ceil(T/dt)						% number of time steps
T=dt*N								% rescale end time to be multiple of dt

z=z0:dz:zf;							% spatial nodes
t=dt:dt:T;							% time steps

gam=dt/(dz*dz);					% constant for diffusion term
%A=gam*(diag(ones(M-1,1),-1)-2*diag(ones(M,1))+diag(ones(M-1,1),1));
A=gam*(diag(ones(M+1,1),-1)-2*diag(ones(M+2,1))+diag(ones(M+1,1),1));

%% Initialize solution vectors (including boundary)
vold=zeros(1,M+2).';

qn=vold;
v=vold;

%% Before beginning the actual process, check if time step is OK
reply = input(sprintf('N = %d. Do you want to continue? Y/N [Y]: ', N), 's');
if isempty(reply)
    reply = 'Y'; % Default to yes
end

if reply == 'Y'
for n=1:N,
  vzz=A*v; 
  %% Initial value analysis
  % Considering that qnp1 defaults to zero (qn is initialized to zero)
  % Vzz is set to zero as well, then qdot equals zero, vold is zero
  % Thus v is zero, and nothing changes...

  %update qdot:
  qdot=gd/tau*vzz-qn/tau;
  qnp1=qn + dt*qdot;
  %qnp1=qn+dt*gd/tau*vzz-dt/tau*qn;
  %qdot=gd/tau*vzz-qnp1/tau;
  v=vold+D*vzz+dt*kap*ha*qdot;
  
  for j=[],%2:M+1

  v(j)= vold(j) ...
	+ D*gam*(vold(j+1)-2*vold(j)+vold(j-1)) ...
	+ dt*kap*ha*qdot(j);%f(j,n);
  end

  qn=qnp1;
  vold=v;
  

  if(mod(n,floor(N/20))==0)
	n
	v
	qdot
	plot(z,v,'-o')
	title(sprintf('n=%d',n))
	az = [0 h 0 h/20];
	axis(az)
	pause
  end
  
end
end
