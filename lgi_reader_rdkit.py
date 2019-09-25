"""
Copyright 2019 Ruud van Deursen, Firmenich SA.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

from rdkit.Chem import MolFromSmiles,MolToInchiKey
import networkx as nx

def IsValid(G,degree):
    """
    Method checks if a graph is valid, i.e. the real
    degree of the vertices matches the character
    in the lgi-string.
    G      -- Graph.
    degree -- Expected degree
    Return:
    True if all degrees match the character in the string.
    """
    for idx in range(G.number_of_nodes()):
        expt = degree[idx]
        act = G.degree(idx)
        if act != expt:
            return False
    return True

class LGIReader:
    """ 
    Class LGIReader defines a class to read LGI-graph strings.
    This class can currently only read lgi-strings with a degree of 1-6.
    """
    
    def __init__(self):
        """ Constructor of LGIReader """
        super(LGIReader,self).__init__()
        self.replacements = [("B","O"),("C","N"),("D","C"),("E","P"),("F","S"),("A","F")]
        self.degree = dict([(char,idx+1) for idx,char in enumerate(["A","B","C","D","E","F"])])
        self.known = [src for src,dst in self.replacements]
        
    def Read(self,lgi):
        """
        Method Read imports an lgi to Graph.
        """
        try:
            # Extract the degree
            D = [self.degree[lgi[idx]] for idx in range(len(lgi)) if lgi[idx] in self.known]

            # Translate to smiles and import using RDKit
            smi = "%s"%(lgi)            
            for src,dst in self.replacements:
                smi = smi.replace(src,dst)
            mol = MolFromSmiles(smi)

            # Define the graph
            G = nx.Graph()
            for bond in mol.GetBonds():
                f,t = bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
                G.add_edge(f,t)            

            # Done
            if IsValid(G,D):
                return G
            else: 
                return None
        except:
            return None
        
"""
Static instance of the class.
"""
lgireader = LGIReader()
