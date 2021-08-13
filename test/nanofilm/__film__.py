from nanofilm.__molecule__ import Molecule
from nanofilm.__lattice__ import Lattice
from nanofilm.__slab__ import Slab
from math import inf
from nanofilm.__exception__ import BasicException

class Film:
    def __init__(self):
        '''
        __init__() -> Class constructor.

        Returns
        -------
        None.
        '''
        self.slab=None
        self.molecules=[]
        self.mindistance=1.0
    
    def set_substrate(self,slab):
        '''
        set_substrate(slab) -> sets the substrate.

        Parameters
        ----------
        slab : Slab object
            Substrate.

        Returns
        -------
        None.
        '''
        if isinstance(slab,Slab):
            self.slab=slab
        else:
            raise FilmError("'slab' must be a Slab object!")
            
            return
        
    def add_molecule(self,molecule,ads_site="",ads_site_index=-1,pos=[],
                     vert_sep=1.0):
        '''
        add_molecule(molecule,ads_site,ads_site_index,pos,vert_sep) -> adsorbs 
            a molecule onto the previously defined substrate.

        Parameters
        ----------
        molecule : Molecule object
            Adsorbed molecule.
        ads_site : string
            Label of the adsorption site. The default is "".
        ads_site_index : int, optional
            Index of an adsorption site on the substrate. The default is -1,
            which means that no adsorption site will be used.
        pos : Python list, optional
            Coordinates in the XY plane where the center of mass of the 
            molecule will be placed over. The default is [].
        vert_sep : float, optional
            Nearest vertical distance separating the molecule from the 
            substrate. The default is 1.0 (Angstrom).

        Returns
        -------
        None.
        '''       
        if self.slab is None:
            raise FilmError("No substrate has been defined!")
            
            return
        elif not isinstance(molecule,Molecule):
            raise FilmError("'molecule' must be a Molecule object!")
            
            return
        
        ads_site=str(ads_site)
        xnew=[0.0,0.0,0.0]
        
        if not isinstance(ads_site_index,int):
            raise FilmError("'ads_site_index' must be an integer!")
            
            return
        
        if (not isinstance(vert_sep,(int,float))) or vert_sep<=0.0:
            raise FilmError("'vert_sep' must be a number greater than zero!")
            
            return
        elif vert_sep<self.mindistance:
            vert_sep=self.mindistance
        
        if len(ads_site)>0 and (ads_site in self.slab.ads_sites.keys()) and\
            ads_site_index>=0 and ads_site_index<len(self.slab.ads_sites[ads_site]):
            xnew[0]=self.slab.ads_sites[ads_site][ads_site_index][0]
            xnew[1]=self.slab.ads_sites[ads_site][ads_site_index][1]
        elif len(pos)==2:
            xnew[0]=pos[0]
            xnew[1]=pos[1]
        else:
            raise FilmError("Either 'ads_site' and 'ads_site_index' or 'pos' must be properly set!")
            
            return
                   
        self.molecules.append(molecule)
        
        zmin=self.slab.top+vert_sep
        dz=zmin-molecule.minz
        xnew[2]=molecule.com[2]+dz
                
        molecule.move_molecule_to(xnew)
        
    def remove_molecule(self,molid):
        '''
        remove_molecule(molid) -> removes an adsorbed molecule.

        Parameters
        ----------
        molid : int
            ID of the adsorbed molecule to be removed.

        Returns
        -------
        None.
        '''
        if isinstance(molid,int) and molid>=0 and molid<Molecule.__currmolid__:
            for i in range(len(self.molecules)):
                if self.molecules[i].id==molid:
                    self.molecules.pop(i)
                    
                    break
        else:
            raise FilmError("'molid' must be an integer greater than or equal to zero!")
            
            return
        
    def enforce_mindist(self,moltype=None,mindist=None):
        '''
        enforce_mindist(mindist) -> if necessary, displaces vertically a molecule
            such that the nearest separation of the molecule with respect to 
            the top of the substrate is the prescribed minimum distance.

        Parameters
        ----------
        moltype : string, optional
            Label identifying the molecule. The default is None.
        mindist : float, optional
            Minimum distance separating the molecule from the substrate. The 
            default is None.

        Returns
        -------
        None.
        '''
        if mindist is None or mindist<self.mindistance:
            mindist=self.mindistance
        else:
            if not isinstance(mindist,(float,int)):
                raise FilmError("'mindist' must be a number greater than zero!")
                
                return
        
        for molecule in self.molecules:
            if moltype is None or molecule.moltype==moltype:
                if molecule.minz-self.slab.top<mindist:
                    dispz=mindist-(molecule.minz-self.slab.top)
                
                    molecule.displace_molecule([0.0,0.0,dispz])
                
    def write_xyz(self,filename="film.xyz",vacuum=10.0):
        '''
        write_xyz(filename,vacuum) -> saves the atomic coordinates of the 
            adsorbed system into an XYZ file.

        Parameters
        ----------
        filename : string, optional
            Name of the XYZ file. The default is "film.xyz".
        vacuum : float, optional
            Height of the vacuum layer to be added to the top of the film. 
            The default is 10.0 (Angstroms).

        Returns
        -------
        None.
        '''
        maxz=-inf
        natoms=len(self.slab.atoms)
        
        for molecule in self.molecules:
            natoms+=len(molecule.atoms)
            
            if molecule.maxz>maxz:
                maxz=molecule.maxz
                
        with open(filename,"w") as f:
            f.write("%d\n" % natoms)
            f.write('Lattice="%.6f %.6f 0.0 %.6f %.6f 0.0 0.0 0.0 %.6f" Properties=species:S:1:pos:R:3\n' 
                    % (self.slab.lattice_vectors[0][0],
                       self.slab.lattice_vectors[0][1],
                       self.slab.lattice_vectors[1][0],
                       self.slab.lattice_vectors[1][1],
                       maxz+vacuum))
            
            for atom in self.slab.atoms:
                f.write("%s %.6f %.6f %.6f\n" % (atom.symbol,atom.x[0],
                                                 atom.x[1],atom.x[2]))
                
            for molecule in self.molecules:
                for atom in molecule.atoms:
                    f.write("%s %.6f %.6f %.6f\n" % (atom.symbol,atom.x[0],
                            atom.x[1],atom.x[2]))
                           
class FilmError(BasicException):
    pass