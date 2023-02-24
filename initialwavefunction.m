function psi0=initialwavefunction(x,x_0,k,sigma)
g=sqrt(1/(sqrt(pi)*sigma))*exp(-(x-x_0).^2/(2*sigma^2));
psi0=exp(1i*k*(x-x_0)).*g;
end