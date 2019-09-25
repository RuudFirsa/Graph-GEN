"""
Copyright 2019 Ruud van Deursen, Firmenich SA.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

from lgi_translator import smi_to_lgi,g6_to_lgi
import tqdm
import concurrent.futures
import multiprocessing
import sys

def Run(lines):
    """
    Method runs the translation of graphs from g6 to lgi.
    Input:
    lines -- List with lines in g6-format
    Return:
    List with lines in lgi-format.
    """
    num_in = len(lines)
    num_processes = multiprocessing.cpu_count()
    with concurrent.futures.ProcessPoolExecutor(num_processes) as pool:
        translated = list(tqdm.tqdm(
            pool.map(g6_to_lgi.Translate, lines, chunksize=16), total=num_in)) 
    return list(filter(lambda x: x is not None,translated))

if __name__ == "__main__":
    # Define
    fin = sys.argv[1]
    fout = fin.replace(".g6",".lgi")
    with open(fin,"r") as f:
        lines = [line.strip().split("\t")[0] for line in f.readlines()]
        
    # Translate
    lgi_lines = Run(lines)
    
    # Write
    with open(fout,"w") as f:
        _ = [f.write("%s\n"%(lgi)) for lgi in lgi_lines]
    print("Output written to %s - %s lines"%(fout,len(lgi_lines)))
