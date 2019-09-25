"""
Copyright 2019 Ruud van Deursen, Firmenich SA.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

from rdkit.Chem import MolToSmiles,MolFromSmiles
from rdkit.Chem import RWMol,Atom,Bond,BondType
import g6_reader

class LGI:

    def __init__(self):
        self.d2atno = [5,9,8,7,6,15,16]
        self.replacements = [("B","@"),("F","a"),("O","b"),("N","c"),("C","d"),("P","e"),("S","f")]
        self.allowed = set([0,1,2,3,4,5,6])

    def Write(self,degrees,edges,canonical=True):
        if set(degrees).issubset(self.allowed):
            # Define the molecule
            cp = RWMol()
            _ = [cp.AddAtom(Atom(self.d2atno[D])) for D in degrees]
            _ = [cp.AddBond(f,t,BondType.SINGLE) for f,t in edges]

            # Export as canonical SMILES or a random SMILES
            if canonical:
                out = MolToSmiles(cp,canonical=True)
            else:
                out = MolToSmiles(cp,canonical=False,doRandom=True)

            # Carry out replacements
            for src,dst in self.replacements:
                out = out.replace(src,dst)
            return out.upper()
        else:
            return None
    
class SmiToLgi(LGI):
    """ 
    Class translates SMILES to graphs. This translator is based on
    organic chemistry and can produce graphs with vertex degrees
    in the range [0,6].
    """
    
    def __init__(self):
        """
        Constructor of Smi
        """
        super(SmiToLgi,self).__init__()
    
    def Translate(self,smi,canonical=True):
        """
        Method translates a SMILES-string to a undirected
        graph G(V,E) with featureless vertices and unweighted
        edges, e.g. the graph equivalent of a saturated hydrocarbon.
        Input:
        smi
        """
        # Make a copy of the molecule to address the degrees
        mol = MolFromSmiles(smi)
        degrees = [atom.GetDegree() for atom in mol.GetAtoms()]
        edges = [(bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()) for bond in mol.GetBonds()]
        return self.Write(degrees,edges,canonical=canonical)
            
class G6ToLgi(LGI):
    """ 
    Class translates SMILES to graphs. This translator is based on
    organic chemistry and can produce graphs with vertex degrees
    in the range [0,6].
    """
    
    def __init__(self):
        super(G6ToLgi,self).__init__()
        
    def Translate(self,g6,canonical=True):
        G = g6_reader.decodegeng(g6)
        V = G.number_of_nodes()
        degrees = [G.degree(idx) for idx in range(V)]
        return self.Write(degrees,G.edges(),canonical)
            
"""
Statically constructed instance of the class.
"""
smi_to_lgi = SmiToLgi()

"""
Static instance to G6ToLGI.
"""
g6_to_lgi = G6ToLgi()
