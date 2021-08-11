*NanoFilm* is an easy-to-use package whose main functionality is to allow the user 
to position and manipulate one or more molecules on substrates, thus generating structures 
that will be used as input coordinates for first-principles calculations or atomistic simulations.

In this very first version, the possibility of creating a molecular film object, i.e., Lattice objects, 
has not yet been implemented. The documentation is still limited to comments in the code. 
A "how-to" of the package's objects and methods can be seen in the Python script given as an example below.

As you surely know, citations are the currency of science. Therefore, 
if you find this package useful for your own research on adsorption of molecules on 
surfaces, please consider citing my articles about this topic:

- Adsorption of metal-phthalocyanine molecules onto the Si(111) surface passivated by delta-doping: ab initio calculations, R. G. A. Veiga, R. H. Miwa, and A. B. McLean. Phys. Rev. B, 93, 115301 (2016).
- Self-assembly of NiTPP on Cu(111): a transition from disordered 1D wires to 2D chiral domains, S. Fatayer, R. G. A. Veiga, M. J. Prieto, E. Perim, R. Landers, R. H. Miwa, and A. de Siervo. Phys. Chem. Chem. Phys., 17, 18344 (2015).
- Structural, electronic, and magnetic properties of pristine and oxygen-adsorbed graphene nanoribbons, R. H. Miwa, R. G. A. Veiga, and G. P. Srivastava. Appl. Surf. Sci., 256, 5776 (2010).
- Quenching of local magnetic moment in oxygen adsorbed graphene nanoribbons, R. G. A. Veiga, R. H. Miwa, and G. P. Srivastava. J. Chem. Phys, 128, 201101 (2008).
- Ab initio study of TCNQ-doped carbon nanotubes, R. G. A. Veiga and R. H. Miwa. Phys. Rev. B, 73, 245422 (2006).

This helps me, among other things, to get funding for my research, which includes working on other codes like this.

If you have any question, just drop an e-mail to *roberto.veiga@ufabc.edu.br*.

# Basic Usage

---

## Adsorbing molecules onto graphene


```python
from nanofilm import AtomType,Atom,Molecule,Slab,Film   # NanoFilm package
from copy import deepcopy
```

#### Sets the AtomType objects


```python
hydrogen=AtomType("H","hydrogen",1,1.007)
carbon=AtomType("C","carbon",12,12.0107)
```

#### Creates first H2 molecule from scratch


```python
mol1=Molecule("H2")                 # H2 molecule, initially empty (no atoms)
h1=Atom(hydrogen,[-0.44,0.0,0.0])   # First H atom
h2=Atom(hydrogen,[0.44,0.0,0.0])    # Second H atom
mol1.add_atom(h1)                   # Adds the first H atom to the molecule
mol1.add_atom(h2)                   # Adds the second H atom to the molecule
mol1.write_xyz()
```

#### Creates the substrate (Slab object) taking the coordinates from a .xyz file


```python
slab=Slab("graphene","graphene.xyz")  # Slab object (in this case, graphene, with coordinates taken from an XYZ file)
slab.replicate_slab(4,3)              # Replicates the slab in the XY plane
slab.add_ads_site([3.69,5.68128])     # First adsorption site
slab.add_ads_site([4.92,4.97114])     # Second adsorption site
slab.write_xyz()
```

#### Adsorbs the H2 molecule on the first adsorption site (creating first a Film object)


```python
film=Film()
film.mindistance=0.75                     # Minimum separation between the molecule and the substrate
film.set_substrate(slab)                  # Sets graphene as the substrate
film.add_molecule(mol1,0,vert_sep=0.75)   # Places the H2 molecule over the first adsorption site
film.write_xyz("film1.xyz")
```

#### Rotates the H2 molecule and enforces a minimum distance from the substrate to prevent overlap


```python
film.molecules[0].rotate_molecule(90,45,0)   # Rotates the H2 molecule
film.enforce_mindist()                       # Required to enforce that the molecule will not overlap with the substrate
film.write_xyz("film2.xyz")
```

#### Makes a (deep) copy of the H2 molecule, placing it above the second adsorption site


```python
mol2=deepcopy(mol1)
film.add_molecule(mol2,1,vert_sep=1.0)
film.molecules[1].rotate_molecule(15,30,45)
film.enforce_mindist()
film.write_xyz("film3.xyz")
```

#### Removes both H2 molecules and adsorbs a benzene molecule


```python
film.remove_molecule(0)                               # Removes the first H2 molecule
film.remove_molecule(0)                               # Second H2 molecule became the first (and unique), then index zero
mol3=Molecule("Benzene","benzene.xyz")                # Benzene molecule (coordinates taken from an .xyz file)
film.mindistance=2.0                                  # Redefines the minimum molecule-substrate separation
film.add_molecule(mol3,pos=[2.5,3.6],vert_sep=2.0)    # Benzene molecule is placed over the substrate at an arbitrary location
film.write_xyz("film4.xyz")
```

#### Rotates the benzene molecule


```python
film.molecules[0].rotate_molecule(30,30,30)
film.write_xyz("film5.xyz")
film.enforce_mindist()
film.write_xyz("film6.xyz")
```

#### Makes a copy of the benzene molecule, places it over the second adsorption site and displaces it


```python
mol4=deepcopy(mol3)
film.add_molecule(mol4,1,vert_sep=2.0)
film.write_xyz("film7.xyz")
film.molecules[1].displace_molecule([1.0,1.5,-5.0])
film.write_xyz("film8.xyz")
film.enforce_mindist()
film.write_xyz("film9.xyz")
```


```python

```
