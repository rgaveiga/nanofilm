from nanofilm.__atom__ import AtomType,Atom
from math import inf
from copy import deepcopy
from nanofilm.__exception__ import BasicException

class Slab:
    def __init__(self,slabname,filename="",atomlist=[],
                 lattice_vectors=[[1.0,0.0],[0.0,1.0]],ads_sites={}):
        '''
        __init__(slabname,filename,atomlist,lattice_vectors,ads_sites) -> class 
            constructor.

        Parameters
        ----------
        slabname : string
            A label allowing you to identify the slab, for instance, "graphene".
        filename : string, optional
            Name of an XYZ file containing the slab structure. The default 
            is "".
        atomlist : Python list, optional
            List of Atom objects to be added to the slab. The default is [].
        lattice_vectors : Python list
            List of slab lattice vectors in the XY plane. The default is 
            [[1.0,0.0],[0.0,1.0]].
        ads_sites : Python dictionary, optional
            Dictionary of adsorption sites on the surface. The default is {}.

        Returns
        -------
        None.
        '''
        if isinstance(slabname,str) and len(slabname)>0:
            self.slabname=slabname
        else:
            raise SlabError("'slabname' must be a string!")
            
            self.__del__()
            
        self.atoms=[]       # List of Atom objects that constitute the slab
        self.ads_sites={}   # Dictionary of adsorption sites
        self.top=-inf       # Maximum Z coordinate of the slab
        self.bottom=inf     # Miniumum Z coordinate of the slab
        
        if len(lattice_vectors)==2:
            if isinstance(lattice_vectors,list):
                if isinstance(lattice_vectors[0],list) and isinstance(lattice_vectors[1],list):
                    if len(lattice_vectors[0])==len(lattice_vectors[1])==2:
                        self.lattice_vectors=lattice_vectors.copy()
            else:
                raise SlabError("'lattice_vectors' must be a list with two elements!")
                
                self.__del__()
        
        if len(filename)>0:
            if isinstance(filename,str):
                self.read_xyz(filename)
            else:
                raise SlabError("'filename' must be a valid file name!")
            
                self.__del__()
            
        if len(atomlist)>0:
            if isinstance(atomlist,list):
                for atom in atomlist:
                    self.add_atom(atom)
            else:
                raise SlabError("'atomlist' must be a non-empyt list!")
            
                self.__del__()
        
        if bool(ads_sites):
            if isinstance(ads_sites,dict):
                for key,site in ads_sites.items():
                    if not key in self.ads_sites.keys():
                        self.ads_sites[key]=[]
                    
                    if isinstance(site,list) and len(site)==2:
                        if isinstance(site[0],(float,int)) and isinstance(site[1],(float,int)):
                            self.ads_sites[key].append(site)
                        else:                            
                            raise SlabError("'pos' must contain numeric coordinates!")
                
                            self.__del__()
            else:
                raise SlabError("'ads_sites' must be a non-empty dictionary!")
                
                self.__del__()
                
    def add_atom(self,atom):
        '''
        add_atom(atom) -> adds an atom to the slab.

        Parameters
        ----------
        atom : Atom object
            New atom added to the slab.

        Returns
        -------
        None.
        '''
        if isinstance(atom,Atom):
            self.atoms.append(atom)
            
            if atom.x[2]>self.top:
                self.top=atom.x[2]
            elif atom.x[2]<self.bottom:
                self.bottom=atom.x[2]
        else:
            raise SlabError("'atom' must be an Atom object!")

            return
        
    def remove_atom(self,atomid):
        '''
        remove_atom(atomid) -> removes an atom from the slab.

        Parameters
        ----------
        index : int
            ID of the atom to be removed.

        Returns
        -------
        None.
        '''
        if isinstance(atomid,int) and atomid>=0 and atomid<Atom.__curratomid__:
            for i in range(len(self.atoms)):
                if self.atoms[i].id==atomid:
                    self.atoms.pop(i)
                    
                    break
        else:
            raise SlabError("'atomid' must be an integer greater than or equal to zero!")
            
            return
        
        self.__get_slab_z__()
        
    def add_ads_site(self,ads_site,pos):
        '''
        add_ads_site(ads_site,pos) -> adds an adsorption site on the slab.

        Parameters
        ----------
        ads_site : string
            Label of the adsorption site.
        pos : Python list
            Coordinates of the adsorption site on the XY plane.

        Returns
        -------
        None.
        '''
        ads_site=str(ads_site)
        
        if not ads_site in self.ads_sites.keys():
            self.ads_sites[ads_site]=[]
        
        if isinstance(pos,list) and len(pos)==2:
            if isinstance(pos[0],(int,float)) and isinstance(pos[1],(int,float)):
                self.ads_sites[ads_site].append(pos)
            else:
                if len(self.ads_sites[ads_site])==0:
                    self.ads_sites.pop(ads_site)
                
                raise SlabError("'pos' must contain numeric coordinates!")
                
                return
        else:
            raise SlabError("'pos' must be a list with two elements!")
                
            return
            
    def remove_ads_site(self,ads_site,index):
        '''
        remove_ads_site(ads_site,index) -> removes an adsorption site.

        Parameters
        ----------
        ads_site : string
            Label of the adsorption site.
        index : int
            Index of the adsorption site.

        Returns
        -------
        None.
        '''
        ads_site=str(ads_site)
        
        if not ads_site in self.ads_sites.keys():
            raise SlabError("'ads_site' is not a valid label of adsorption sites!")
            
            return
        
        if isinstance(index,int) and index>=0 and index<len(self.ads_sites[ads_site]):
            self.ads_sites[ads_site].pop(index)
            
            if len(self.ads_sites[ads_site])==0:
                self.ads_sites.pop(ads_site)
        else:
            raise SlabError("'index' must be an integer greater than or equal to zero!")
            
            return
        
    def replicate_slab(self,n,m,replicate_ads_sites=False):
        '''
        replicate_slab(n,m,replicate_ads_sites) -> replicates the slab in the 
            XY plane.

        Parameters
        ----------
        n : int
            Number of replications along the first lattice vector.
        m : int
            Number of replications along the second lattice vector.
        replicate_ads_sites : logical, optional
            Whether also replicate the adsorption sites or not. The default is 
            False.

        Returns
        -------
        None.
        '''
        if not (isinstance(n,int) and isinstance(m,int) and n>=0 and m>=0):
            raise SlabError("'n' and 'm' must be positive integers!")
            
            return
        elif n==0 and m==0:
            raise SlabError("Either 'm' or 'n' must be greater than zero!")
            
            return
        
        ucell=deepcopy(self.atoms)
        
        if replicate_ads_sites and bool(self.ads_sites):
            ucell_ads_sites=deepcopy(self.ads_sites)
        else:
            ucell_ads_sites={}
        
        for i in range(n):
            for j in range(m):
                if i==0 and j==0:
                    continue
                
                for atom in ucell:
                    newatom=deepcopy(atom)
                    newatom.x[0]+=i*self.lattice_vectors[0][0]+j*self.lattice_vectors[1][0]
                    newatom.x[1]+=i*self.lattice_vectors[0][1]+j*self.lattice_vectors[1][1]

                    self.add_atom(newatom)
                    
                if bool(ucell_ads_sites):                    
                    for key,site in ucell_ads_sites.items():
                        newsite=site.copy()
                        newsite[0]+=i*self.lattice_vectors[0][0]+j*self.lattice_vectors[1][0]
                        newsite[1]+=i*self.lattice_vectors[0][1]+j*self.lattice_vectors[1][1]
                        
                        self.add_ads_site(key,newsite)
                    
        self.lattice_vectors[0][0]*=n
        self.lattice_vectors[0][1]*=n
        self.lattice_vectors[1][0]*=m
        self.lattice_vectors[1][1]*=m
                                
    def read_xyz(self,filename):
        '''
        read_xyz(filename) -> reads the structure of the slab from a file
            in XYZ format.

        Parameters
        ----------
        filename : string
            Name of the XYZ file.

        Returns
        -------
        None.
        '''
        with open(filename,"r") as f:
            lines=f.readlines()
            
        count=0
        extxyz=-1

        for line in lines:
            if count==0:
                count=1
                
                continue
            elif count==1:
                extxyz=line.find("Lattice")
                
                if extxyz>=0:
                    line=line.replace('Lattice','')
                    line=line.replace('=','')
                    line=line.replace('"','')
                
            l=line.split()
            
            if count==1 and extxyz>=0:
                self.lattice_vectors[0][0]=float(l[0])
                self.lattice_vectors[0][1]=float(l[1])
                self.lattice_vectors[1][0]=float(l[3])
                self.lattice_vectors[1][1]=float(l[4])
                count=2
                
                continue
        
            if len(l)>=4:
                for atomtype in AtomType.__AtomTypes__:
                    if l[0]==atomtype.symbol:
                        self.atoms.append(Atom(atomtype,x=[float(l[1]),
                                                           float(l[2]),
                                                           float(l[3])]))

                        break

        self.__get_slab_z__()
    
    def write_xyz(self,filename="slab.xyz"):
        '''
        write_xyz(filename) -> saves the atomic coordinates of the slab into
            an XYZ file.

        Parameters
        ----------
        filename : string, optional
            Name of the XYZ file. The default is "slab.xyz".

        Returns
        -------
        None.
        '''
        dimz=(self.top-self.bottom)+20.0
        
        with open(filename,"w") as f:
            f.write("%d\n" % (len(self.atoms)))
            f.write('Lattice="%.6f %.6f 0.0 %.6f %.6f 0.0 0.0 0.0 %.6f" Properties=species:S:1:pos:R:3\n' 
                    % (self.lattice_vectors[0][0],self.lattice_vectors[0][1],
                       self.lattice_vectors[1][0],self.lattice_vectors[1][1],
                       dimz))
            
            for atom in self.atoms:
                f.write("%s %.6f %.6f %.6f\n" % (atom.symbol,atom.x[0],
                                                 atom.x[1],atom.x[2]))
        
    def __get_slab_z__(self):
        '''
        __get_slab_z__() -> sets the minimum and maximum Z coordinates of 
            the slab.

        Returns
        -------
        None.
        '''        
        for atom in self.atoms:
            if atom.x[2]>self.top:
                self.top=atom.x[2]
            elif atom.x[2]<self.bottom:
                self.bottom=atom.x[2]
                
    def __del__(self):
        '''
        __del__() -> class destructor

        Returns
        -------
        None.
        '''
        pass
    
class SlabError(BasicException):
    pass