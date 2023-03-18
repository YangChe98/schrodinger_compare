function dy=schrodingereq(t,y,x,J)
T1=ones(J+1,1);
dx=0.02;
T=spdiags([-T1,2*T1,-T1],-1:1,J+1,J+1)/(2*dx^2);
V=functionv(x,t);
H=T+V;
dy=-1i*H*y;
end