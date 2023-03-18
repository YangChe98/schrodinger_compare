function psi_cn=CNmethod(psi0,dt0,x,J,time_step,dx)
T1=ones(J+1,1);
psi_cn=zeros(J+1,time_step+1);
psi_cn(:,1)=psi0.';
T=spdiags([-T1,2*T1,-T1],-1:1,J+1,J+1)/(2*dx^2);
t_now=0;
% v=zeros(size(x));
% V_now1=diag(v,0);
% V_now=sparse(V_now1);
% V_now=functionv(x,t_now);
sparsei = speye(J+1);
for i=1:time_step
    t_forward=t_now;
    t_now=t_now+dt0;
    V_now=functionv(x,t_now);
    V_forward=functionv(x,t_forward);
    U1=sparsei -1i*dt0*(T+V_forward)/2;
    U2=sparsei +1i*dt0*(T+V_now)/2;
    psi_now=U2\(U1*psi_cn(:,i));
    psi_cn(:,i+1)=psi_now;
end
end