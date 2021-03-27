for s in conseq pfetch trimal2msa acc2operon rconfig rmsd_comparer uniprot2gi rnexplorer
do
	eval "$(register-python-argcomplete $s)"
done
