To do
============

## rotifer/devel/alpha/sequence.py
- [x] Substitute the __init__ function to one that will use a load df function  
- [x] Create a load function (similar to the actual __init__)  
- [x] Use the load function to imporove the realing function. Ex: return load(realing)  
- [x] Create function to_hmm and to_hhm (The hhm could use information of secondary structure. see hhpred wiki)
- [x] Pairwise comparasion of all x all, network, find the most connected (different edge values e cutoff), create community and order by it. 
- [x] Create function to run hhsearch with the aligment. To decided wich PDB should be add to the sequence aligment.   
- [x] Create the function to add fetch and add sequence to the aligment.
- [ ] ADD the two HHSEARCH line resutls (the one with |, +, and the one with consensus ) to the annotation function, to use it blast result.
- [ ] ADD  option to the to_file to be able to export it in tab format with all annotation, and selecting the delimiter ( use np.savetxt, to be able to use two charcater as delimiter).
- [ ] Fix the sort function to automatically detect a list or string in the by options.
- [ ] Sort aln by specific position (more than one is also possible)
- [ ] ADD to fig, use the ggnicastro/bin/aln2fig2 as a start

