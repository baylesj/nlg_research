%%
% Jordan Bayles
% URISC/Spring 2012 Research Project
% Model based on advection-diffusion equation, expressed as
% \[ v_{t} - \kappa H_A G_{\infty} v_{zz} = \kappa H_A \dot{Q} \]

%% Question for Gibson:
% I did some research trying to get a feel for advection/diffusion
% models after our last discussion, but I am still a little unsure how
% everything fits together.
% I think that in our model we currently have no advection; just diffusion, and a
% forcing term of the form $\kappa * H_A * \dot{Q}$ however I am not
% *completely* sure if this is the right approach. WIth that in mind, this
% may not be a completely correct interpretation

%% Physical parameters
h=1.28*10^(-3);          			% m, may need conversion
kap=1.72*10^(-5); 					% m^4/(Ns), may need conversion
ha=1.29*10^(6);						% Pa, may need conversion

z0=0;								% start of computational domain
zf=h;								% end of computational domain
T=85.1;								% end time
M=64;								% number of nodes (minus boundary)
dz=(zf-z0)/(M+1);					% spatial step

c=1;								% From Hader paper
ginfty=1;							% From our notes
gd=c;								% From our notes
qdot=0;								% Assumption for now, will be changed

v=0;								% advection velocity
D=kap*ha*ginfty						% diffusion coefficient

% Define external force in space-time domain
f=kap*ha*qdot;
% Need characterization of qdot, for now try
f= z.*t.*0
%% Time considerations
sf=0.9;								% safety factor (Mingham)
%dt=sf*(dz*dz/(v*dz+2*D))    		% time step (how formulated?)
dt = 0.1;
N=ceil(T/dt)						% number of time steps
T=dt*N								% rescale end time to be multiple of dt

z=z0:dz:zf;							% spatial nodes
t=dt:dt:T;							% time steps

nu=v*dt/dz;							% cfl number
gam=D*dt/(dz*dz);					% constant for diffusion term

%f=(abs(z-0.5-dz/4)<dz/2)'*((t-pi<dt/2).*sin(t));

%% Initialize solution vectors (including boundary)
uold=zeros(1,M+2);
u=uold;
pause

for n=1:N,
  for j=2:M+1
    u(j)= uold(j)-nu*(uold(j)-uold(j-1)) ...
	+ gam*(uold(j+1)-2*uold(j)+uold(j-1)) ...
	+ dt*f(j,n);
  end
  
  uold=u

  if(mod(n,floor(N/20))==0)
    n
    plot(z,u,'-o')
    title(sprintf('n=%d',n))
    az = [0 h 0 0.04];
	axis(az)
    pause
  end
  
end
