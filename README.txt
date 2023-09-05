Updated Graphical User Interface on Atom/Molecule Descriptor Profile

The original graphical user interface (Guidemol) was set available in github framework:

https://github.com/jairesdesousa/guidemol/tree/main

This graphical user interface consists on a Python computer program based on the RDKit software (RDKit: Open-source cheminformatics. https://www.rdkit.org). 
The original version processes molecular structures and calculate molecular descriptors with a graphical user interface using the tkinter package. 
It can calculate descriptors already implemented in RDKit as well as grid representations of 3D molecular structures using the electrostatic potential or voxels.
The complete description was released in ChemRxiv:

https://chemrxiv.org/engage/chemrxiv/article-details/6448fd24e4bbbe4bbf475860

The new version of this open source graphical user interface (Atom-Molecule-Descriptor-Profile.py) involves the generation of atomic descriptors of a generical chemical:

1- Atomic number.
2- Atomic contribution to Molecular refractivity.
3- LogP contribution of atoms.
4- MMFF Charges.
5- Gasteiger charges.
6- Mulliken-Jaffe electronegativities:(https://www.knowledgedoor.com/2/elements_handbook/mulliken-jaffe_electronegativity.html).
7- Jazzy set of atomic properties (https://jazzy.readthedocs.io/en/latest/cookbook.html) :
	a. atomic number (z)
	b. formal charge (q)
	c. partial charge (eeq)
	d. atomic-charge dependent dynamic atomic polarizabilities (alp)
	e. number of lone pairs (num_lp)
	f. Acceptor strength (sda)
	g. C-H donor strength (sdc)
	h. X-H donor strength when X is a non-carbon atom (sdx)
	i. hybridisation (hyb)

Author: Gonçalo V. S. M. Carrera

This work was supported by the Associate Laboratory for Green Chemistry – LAQV which is financed by National Funds from FCT/MCTES (UIDB/50006/2020 and UIDP/50006/2020).
