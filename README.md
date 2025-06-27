The two UMATs define general finite strain  viscoelasticity constitutive equations with damage. 

One version is standard, the second is defined for HYBRID elements.

This is a general version of Reese and Govindjee (1998) model applied for any type of (I_1,I_2)-invariants strain energy for the deviatoric part and any viscoelastic hydrostatic energy U(J,Je). 

The constitutive equations as the damage parameters are detailed in [2.	F. Gouhier, J. Diani, A. vandenbroucke, 2024. A finite strain viscoelastic model with damage and tension-compression asymmetry considerations for solid propellants. Mechanics of Materials, 199, 105152, https://doi.org/10.1016/j.mechmat.2024.105152.]

The HowtoUse_UMAT.pdf file shows several validation cases and explains how to modify the UMAT according to the strain energy of interest.

Several input files are provided for testing on a 1 element cube or on structure rectangular bar for torsion (DMAtorsion). Note that for the latest test, DMAtorsionS.inp and DMAtrosionR.inp allows comparing the result between Simo's model directly implemented in Abaqus and Govindjee and Reese's model implemented here.
