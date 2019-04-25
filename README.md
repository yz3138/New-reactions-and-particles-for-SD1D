# New-reactions-and-particles-for-SD1D
Some new reactions and new species of particles are included into SD1D module through these files.

By Y.L.Zhou
24/04/2019
This branch is created to add new reactions and particles in SD1D.

New particle species: 
1. hydrogen molecules (H2), 
2. charged hydrogen molecules (H2+).

New reactions: 
1. Dissociation ionization rate coefficient (DI); T_min = 3.0eV.

   e- + H2 -> e- + e- + H2+ -> H+ + H(1s) + 2e-
2. Molecules dissociated by collision with atoms (Diss_H_H2); T_min = 3.0eV.
   H(1s) + H2(v) -> H(1s) + 2H(1s)
3. Molecule ionized by collision with electrons; T_min = 3.0eV.
   e + H2 -> 2e + H2+
4. H2+ ionized by collision with electrons; T_min = 3.0eV
   e + H2+ -> 2e + 2H+

The data used in these files are from:
[1] R. K. Janev, D. Reiter, Collision Processes in Low-Temperature Hydrogen Plasmas, Zentralbibliothek Jülich (2003)
[2] R. K. Janev, W. D. Langer, K. Evans Jr, D.E. Post Jr, Elementary Processes In Hydrogen-Helium Plasmas, Springer-Verlag Berlin Heidelberg (1987)
