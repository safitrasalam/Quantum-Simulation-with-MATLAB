% This is a program to simulate the Rydberg atom with a triangle configuration.
% Define variables:
N = 100000;
n = 3;
t = linspace(0,30e-6,N);
dt = t(2)-t(1);
hbar = 1.05457182e-34;
Omega = 2*pi*1e+6;
C_6 = 2*pi*516*1e9;
d13 = 6.6;
d12 = 6.6;

% Matrices and constants that will be use
sx = [0,1;1,0];
sz = [1,0;0,-1];
gamma = 2*pi*30*1e3;
F = sqrt(gamma/2)*sz;
I = [1,0;0,1];
P = [1,0;0,0];

% Van der-Waals interactions between Rydberg atoms
U_matrix = zeros(8,8);
U12 = C_6/d12^6;
U23 = U12;
U13 = C_6/d13^6;
U_matrix(7,7) = U12;
U_matrix(4,4) = U23;
U_matrix(6,6) = U13;
U_matrix(8,8) = U13 + U23 + U12;

% Hamiltonian of the system
H = hbar*Omega*(kron(kron(sx,P),P) + kron(kron(P,sx),P) + kron(kron(P,P),sx))/2 + hbar*U_matrix;

% Density matrix and its initial value
rho = zeros(size(H,1),size(H,2)); % rho_t 2^N_atom, 2^N_atom 
rho(1,1) = 1;

% Lindbladian for each atom
L1 = kron(kron(F,I),I);
L2 = kron(kron(I,F),I);
L3 = kron(kron(I,I),F); 
L0 = L1 + L2 + L3;

% Runge-Kutta 4th order method to solve the Lindblad-Master equation
for i=1:N
    % {L'*L, rho}
    A1 = L1'*L1*rho + rho*L1'*L1;
    A2 = L2'*L2*rho + rho*L2'*L2;
    A3 = L3'*L3*rho + rho*L3'*L3;


    % Individual and collective dephasing
    L_ind = L1*rho*L1' + L2*rho*L2' + L3*rho*L3' - 1/2*(A1 + A2 + A3);
    L_c = L0*rho*L0' - 1/2*(L0'*L0*rho + rho*L0'*L0);

    % Density matrix equation
    rho_dot = @(t,rho) (-sqrt(-1)/hbar)*(H*rho-rho*H) + L_ind + L_c;

    k_1 = rho_dot(t(i),rho);
    k_2 = rho_dot(t(i)+0.5*dt,rho+0.5*dt*k_1);
    k_3 = rho_dot((t(i)+0.5*dt),(rho+0.5*dt*k_2));
    k_4 = rho_dot((t(i)+dt),(rho+k_3*dt));
    rho = rho + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;
    prob_0(i)=rho(1,1); % Probabilitas atom kembali ke keadaan awal
end

% plot probability of an atom return to ground state (P0) vs time (t)
figure
plot(t,prob_0)
xlim([0 5e-6])
xlabel('Time (second)')
ylabel('P_0(t)')

% Fast Fourier transform of probability of an atom return to the ground state,
spectrum_dynamics=fftshift(fft(prob_0));
fs=1/dt;
w = (-N/2:N/2-1)*(fs/N);

% Plot FFT(P0) vs frequency (f)
figure
plot(w/1e6,abs(spectrum_dynamics))
xlabel('Frequency (MHz)')
ylabel('FFT(P_0(t))')
xlim([0 5])
ylim([0 4e3])
