# Quantum Simulation
During my internship in Korea Advanced Institute of Science and Technology (KAIST), I joined the Quantum Computing Laboratory under the supervision of Prof. Jaewook Ahn and Mr. Andrew Byun. In this project, I simulate the Rydberg atom quantum computing by using the PXP model. It is a specific quantum model that describes the interaction between neighboring spins. The purpose of this work is to predict the output and the wavelength needed to excite the atoms to the Rydberg state. Here are the steps of the simulation:

1. Define the variables,
2. Determine the Hamiltonian from the PXP model,
3. Calculate the Eigenenergies and Eigenstate,
4. Determine the density matrix from Lindbladian-Master equation,
5. Solve the equation using 4th order Runge-Kutta method,
6. Plot the density matrix against time,
7. Apply the fast Fourier transform to density matrix,
8. Plot the FFT of density matrix against frequency.

Here are the results for triangle:

![Results for triangle](https://github.com/safitrasalam/Quantum-Simulation-with-MATLAB/blob/main/Quantum%20Simulation%20Code/Triangle.png)

Right triangle:

![Results for right triangle](https://github.com/safitrasalam/Quantum-Simulation-with-MATLAB/blob/main/Quantum%20Simulation%20Code/Right%20triangle.png)

Linear: 

![Results for linear](https://github.com/safitrasalam/Quantum-Simulation-with-MATLAB/blob/main/Quantum%20Simulation%20Code/Linear.png)

Star (S4):

![Results for S4](https://github.com/safitrasalam/Quantum-Simulation-with-MATLAB/blob/main/Quantum%20Simulation%20Code/S4.png)

Complete (K4):

![Results for K4](https://github.com/safitrasalam/Quantum-Simulation-with-MATLAB/blob/main/Quantum%20Simulation%20Code/K4.png)

Cycle (C4):

![Results for C4](https://github.com/safitrasalam/Quantum-Simulation-with-MATLAB/blob/main/Quantum%20Simulation%20Code/C4.png)

Diamond (K4e):

![Results for K4e](https://github.com/safitrasalam/Quantum-Simulation-with-MATLAB/blob/main/Quantum%20Simulation%20Code/K4e.png)
