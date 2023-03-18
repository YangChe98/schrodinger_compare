% Define the time-dependent Hamiltonian
% function H = H(t)
%     H = [0, exp(-1i*t); exp(1i*t), 0];
% end

% Define the initial wave function
psi0 = [1; 0];

% Define the time interval and time step
t0 = 0;
tf = 10;
dt = 0.01;
Nt = floor((tf-t0)/dt);

% Define the evolution operator for H(t) over a time step dt
% function U = U(dt)
%     H1 = H(dt/2);
%     H2 = H(dt);
%     H3 = H(dt/2);
%     U = expm(-1i*H2*dt)*expm(-1i*H1*dt)*expm(-1i*H3*dt);
% end

% Initialize arrays to store the wave function at each time step
psi = zeros(2, Nt+1);
psi(:,1) = psi0;

% Time evolution using the split-operator method
for i = 1:Nt
    psi(:,i+1) = U(dt)*psi(:,i);
end

% Plot the results
t = linspace(t0,tf,Nt+1);
figure;
plot(t, abs(psi(1,:)).^2, 'b-', 'LineWidth', 2);
hold on;
plot(t, abs(psi(2,:)).^2, 'r--', 'LineWidth', 2);
xlabel('Time');
ylabel('Probability');
legend('|psi_1|^2', '|psi_2|^2');

