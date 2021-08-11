# FC-ElectromangeticShowerSimulation

This is the final project for the computational physics course of MEFT.
The objective of the project was to simulate an electromagnetic particle shower initiated by a high energy photon, using monte-carlo methods.

Description of the files:

- Material- structure that stores constants relative to the material where the particles propagate.
- Particle- structure that stores the energy, linear momentum, initial and final position, name and vivo (if the particle has reached the final position, stopped propagating) and has the constructor, destructor, print, amd methods to get ans set the values mentioned along with the mass of the particle.
- Propagator - structure made to propagate particles P, in a medium Al. It public methods to propagate photons, positrons and electrons. It has private methods for pair production (photon creates an electron and a positron), loss of energy by bremsstrahlung ratiation and inelastic colisions (for electron and positrons), electron-positron annihilation with the creation of photons.
- Funcoes - methods to draw histograms, graphs and store values in files.
- Derivator - sub class of Func1D, inplements numerical derivation
- EqSolver - implements gauss elimination method
- FCmatrix - structure for matrix 
- FCmatrixFull - sub class of FCmatrix where all the values in the matrix are stored. It has medothds to get a column, row, maximum of a column, or row, the calculation of the determinant and a methods to swap rows, along with the =,+,-,*,[] methods.
- Func1D - structure for functions
- Integrator - sub class of Func1D. implements numerical integration 
- projeto.C - main of the project
- makefile
- cFCgraphics - made by the course prof. 
