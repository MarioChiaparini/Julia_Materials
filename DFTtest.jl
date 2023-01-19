using DFTK 
using Plots 
using Unitful
using UnitfulAtomic
# Silicon Ground State Energy - LDA basis set 
print("Classical example for LDA Ground State Calculation - Silicon Crystal")
#parameters os silicon atom 
a = 5.431u"angstrom" # silicon lattice constant 
lattice = a / 2 * [[0 1 1.]; # Silicon lattice vectors
                   [1 0 1.]; # Specified columns
                   [1 1 0.];];
#load HGH pseudopotentials for the Silicon atom 
Silicon = ElementPsp(:Si, psp=load_psp("hgh/lda/Si-q4"))
# Specify the atoms position 
atoms = [Silicon, Silicon]
position = [ones(3)/8, -ones(3)/8]
# Select the type of model and it basis set  
type = model_LDA(lattice, atoms, position)
kgrid = [4 4 4] # k-points grid (Regular Monkhorst-Pack grid)
Kinetic_cutoff = 7 # Kninetic energy cut off 
basis = PlaneWaveBasis(model; Kinetic_cutoff; kgrid)
# Run the SCF procedure to take the ground state 
scf = self_consistent_field(basis, tol=1e-4)