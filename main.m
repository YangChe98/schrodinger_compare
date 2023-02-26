L=20;
x_min=-L/2;
x_max=L/2;
J=998+1;
x=x_min:L/J:x_max;
dx=L/J;
time_step=2000;

%%%initial wavefunction
x_0=0;
k_0=5;
sigmax=0.4;
psi0=initialwavefunction(x,x_0,k_0,sigmax);

dt=2/(2*pi/sigmax+k_0)^2;
dt0=0.001;
dt1=0.0005;
total_time=2;
% ;
% 
% T1=ones(J+1,1);
% psi_cn=zeros(J+1,time_step+1);
% psi_cn(:,1)=psi0.';
% T=-spdiags([-T1,2*T1,-T1],-1:1,J+1,J+1)/(2*dx^2);
% t_now=0;
% % v=zeros(size(x));
% % V_now1=diag(v,0);
% % V_now=sparse(V_now1);
% % V_now=functionv(x,t_now);
% sparsei = speye(J+1);
% for i=1:time_step
%     t_forward=t_now;
%     t_now=t_now+dt0;
%     V_now=functionv(x,t_now);
%     V_forward=functionv(x,t_forward);
%     U1=sparsei +1i*dt0*(T+V_forward)/2;
%     U2=sparsei -1i*dt0*(T+V_now)/2;
%     psi_now=U2\(U1*psi_cn(:,i));
%     psi_cn(:,i+1)=psi_now;
% end
% psi_cn=CNmethod(psi0,dt0,x,J,time_step,dx);
time_step1=4000;
psi_rk4=RK4method(psi0,dt1,x,J,time_step1,dx);
% printgif("cn_method.gif",psi_cn,time_step,x,1e-5)
% 
% fig = figure;
% 
% for i=1:(time_step+1)
%     plot(x,abs(psi_cn(:,i)))
%     xlim([-10,10])
%     ylim([0,1.2])
%     drawnow
%     frame = getframe(fig);
%     im{i} = frame2im(frame);
% end
% 
% filename = "testgif.gif"; % Specify the output file name
% for idx = 1:(time_step+1)
%     [A,map] = rgb2ind(im{idx},256);
%     if idx == 1
%         imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",dt0);
%     else
%         imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",dt0);
%     end
% end