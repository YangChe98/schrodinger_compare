function V_now=functionv(x,t)
%v=zeros(size(x));
l=10;
v=pi^2*(1-cos(pi*x/l))/(2*l^2)+(sin(t))^2*pi*(sin(pi*x/l))/l;
V_now1=diag(v,0);
V_now=sparse(V_now1);
end