L=20;
x_min=-L/2;
x_max=L/2;
J=998+1;
x=x_min:L/J:x_max;
dx=L/J;

%%%initial wavefunction
x_0=0;
k_0=5;
sigmax=0.4;
psi0=initialwavefunction(x,x_0,k_0,sigmax);

dt=2/(2*pi/sigmax+k_0)^2;
dt0=0.001;
dt1=0.0005;
total_time=2;

T1=ones(J+1,1);
T=-spdiags([-T1,2*T1,-T1],-1:1,J+1,J+1)/(2*dx^2);