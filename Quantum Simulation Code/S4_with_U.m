clear
% variables
N = 100000;
n = 3;
t = linspace(0,30e-6,N);
dt = t(2)-t(1);
hbar = 1.05457182e-34;
Omega = 2*pi*1.0e+6;
C_6 = 2*pi*516e9;
lattice_const=6.6;
d12 = lattice_const;
d13 = lattice_const;
d14 = lattice_const;
d23 = lattice_const*sqrt(3);
d24 = lattice_const*sqrt(3);
d34 = lattice_const*sqrt(3);

% Matrices and constant
sx = [0,1;1,0];
sz = [1,0;0,-1];
gamma = 2*pi*30*1e3;
F = sqrt(gamma/2)*sz;
I = [1,0;0,1];
P = [1,0;0,0];

% Van der-Waals interaction 
U_matrix = zeros(16,16);
U12 = C_6/d12^6;
U13 = C_6/d13^6;
U14 = C_6/d14^6;
U23 = C_6/d23^6;
U24 = C_6/d24^6;
U34 = C_6/d34^6;

U_matrix(4,4) = U34;
U_matrix(6,6) = U24;
U_matrix(7,7) = U23;
U_matrix(8,8) = U23 + U24 + U34;
U_matrix(10,10) = U14;
U_matrix(11,11) = U13;
U_matrix(12,12) = U13 + U14 + U34;
U_matrix(13,13) = U12;
U_matrix(14,14) = U12 + U14 + U24;
U_matrix(15,15) = U12 + U13 + U23;
U_matrix(16,16) = U12 + U13 + U14 + U23 + U24 + U34;

%Hamiltonian
H = hbar*Omega*(kron(kron(kron(sx,P),P),P) + kron(kron(kron(P,sx),I),I) + kron(kron(kron(P,I),sx),I) + kron(kron(kron(P,I),I),sx))/2 + hbar*U_matrix;

% rho
rho = zeros(size(H,1),size(H,2)); % rho_t 2^N_atom, 2^N_atom 
rho(1,1) = 1;

% Linbladian

L1 = kron(kron(kron(F,I),I),I);
L2 = kron(kron(kron(I,F),I),I);
L3 = kron(kron(kron(I,I),F),I);
L4 = kron(kron(kron(I,I),I),F);
L0 = L1 + L2 + L3 + L4;

%the equation
% rho_dot = @(t,rho) (-sqrt(-1)/hbar)*(H.*rho-rho.*H);   
dt=t(2)-t(1);

for i=1:N
    % {L'*L, rho}
    A1 = L1'*L1*rho + rho*L1'*L1;
    A2 = L2'*L2*rho + rho*L2'*L2;
    A3 = L3'*L3*rho + rho*L3'*L3;
    A4 = L4'*L4*rho + rho*L4'*L4;

    % Linblad Ind and C 
    L_ind = L1*rho*L1' + L2*rho*L2' + L3*rho*L3' - 1/2*(A1 + A2 + A3 + A4);
    L_c = L0*rho*L0' - 1/2*(L0'*L0*rho + rho*L0'*L0);
    % Runge-Kutta
    rho_dot = @(t,rho) (-sqrt(-1)/hbar)*(H*rho-rho*H) + L_ind + L_c;

    k_1 = rho_dot(t(i),rho);
    k_2 = rho_dot(t(i)+0.5*dt,rho+0.5*dt*k_1);
    k_3 = rho_dot((t(i)+0.5*dt),(rho+0.5*dt*k_2));
    k_4 = rho_dot((t(i)+dt),(rho+k_3*dt));
    rho = rho + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;
    prob_0(i)=rho(1,1);
end

%plot

%for i=1:N
 %    rho=rho+rho_dot(t,rho)*dt;
 %    prob_0(i)=rho(1,1);
%end

%plot
% plot(t,rho)
figure
plot(t,prob_0)
xlim([0 5e-6])
xlabel('Time (second)')
ylabel('P_0(t)')

% fft
spectrum_dynamics=fftshift(fft(prob_0));
fs=1/dt;
w = (-N/2:N/2-1)*(fs/N);
figure
plot(w/1e6,abs(spectrum_dynamics))
xlabel('Frequency (MHz)')
ylabel('FFT(P_0(t))')
xlim([0 5])
ylim([0 1e4])
