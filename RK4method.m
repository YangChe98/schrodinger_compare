function psi_rk4=RK4method(psi0,dt0,x,J,time_step,dx)
T1=ones(J+1,1);
psi_rk4=zeros(J+1,time_step+1);
psi_rk4(:,1)=psi0.';
T=spdiags([-T1,2*T1,-T1],-1:1,J+1,J+1)/(2*dx^2);
t_now=0;
% v=zeros(size(x));
% V_now1=diag(v,0);
% V_now=sparse(V_now1);
% V_now=functionv(x,t_now);
sparsei = speye(J+1);
for i=1:time_step
    t1=t_now;
    t2=t_now+dt0/2;
    t3=t_now+dt0/2;
    t4=t_now+dt0;
    t_now=t_now+dt0;
   
    V1=functionv(x,t1);
    V2=functionv(x,t2);
    V3=functionv(x,t3);
    V4=functionv(x,t4);
    k1=dt0*(-1i*(T+V1))*psi_rk4(:,i);
    k2=dt0*(-1i*(T+V2))*(psi_rk4(:,i)+k1/2);
    k3=dt0*(-1i*(T+V3))*(psi_rk4(:,i)+k2/2);
    k4=dt0*(-1i*(T+V4))*(psi_rk4(:,i)+k3);

 
    psi_rk4(:,i+1)=psi_rk4(:,i)+(k1 + 2*k2 + 2*k3 + k4) / 6;
end
end