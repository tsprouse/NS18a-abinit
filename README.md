# NS18a-abinit
Ab initio code for final project in NucStruc-SPR2018

Run (or import) abinit - can set mj,parity, nmax truncation / subspaces, and identify nucleus by Z and A for basis enumeration, hamiltonian generation, and solution. Gets the ground state 'right.'
End of script shows results for He4 and H2 binding energies. In principle, *could* solve any nucleus, but for the wrong (right?) nucleus / computer combination, it could take a few million years.

Probably maxes out around A=5 before basis enumeration would take prohibitively long to perform. For higher values of A, would need a "smarter" basis enumeration routine - but that wasn't so much of the point for this exercise...

tbme.py just solves the deuteron problem using A.S. TBME's provided from Morten (N2LO for hbar-omega = 10MeV?), in contrast to the full-blown second-quantization approach used in abinit.py. We wrote this first before moving on to the abinit.py script for arbitrary nuclear systems.

--Trevor and Patrick



