To do
============

## rotifer/devel/alpha/sequence.py
- [x] Substitute the __init__ function to one that will use a load df function  
- [x] Create a load function (similar to the actual __init__)  
- [x] Use the load function to imporove the realing function. Ex: return load(realing)  
- [ ] Create function to_hmm and to_hhm (The hhm could use information of secondary structure. see hhpred wiki)
- [ ] Pairwise comparasion of all x all, network, find the most connected (different edge values e cutoff), create community and order by it. 
- [ ] __thinking about it__ Remove freq table from the __init__. Just call it when it needed.
- [ ] Create function to run hhsearch with the aligment. To decided wich PDB should be add to the sequence aligment.   
- [ ] Create the function to add fetch and add sequence to the aligment.
