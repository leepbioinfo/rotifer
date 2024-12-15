-- Rotifer's shell configuration

-- Local variables
local myversion = "0.2"
local myfn  = myFileName()
local mymfn = myModuleFullName()
local mydir = realpath(myfn:gsub(mymfn .. ".lua",""):gsub("/etc/lmod/modules/",""))
local datadir = mydir .. "/share/rotifer/data"

-- Environment variables
prepend_path("PYTHONPATH", mydir .. "/lib")
prepend_path("PERL5LIB", mydir .. "/perl/lib")
if (os.getenv("DATABASES") == nil) then
	setenv("DATABASES",datadir)
else
	datadir = os.getenv("DATABASES")
end
if (os.getenv("ROTIFER_DATA") == nil) then
	setenv("ROTIFER_DATA",datadir)
end

-- PATH
prepend_path("PATH", mydir .. "/bin")

-- ----------
-- Functions
-- ----------

-- acc2pfam
set_shell_function("acc2pfam"," \
    SEQ=$(efetch -db protein -format fasta -id $1);\
    ( echo \"$SEQ\" | hmmscan --cpu 4 /databases/pfam/Pfam - | hmmer2table -c model=version;\
    echo \"$SEQ\" | phobius 2> /dev/null | phobius2table -e 0.0101 ) | domain2architecture -e 0.0101 | architecture2table | padtable\
","")

-- acc2profiledb
set_shell_function("acc2profiledb"," \
    ( efetch -db protein -format fasta -id $1 | rpsblast -db $ROTIFER_DATA/rpsdb/allprofiles | blast2table -c profiledb -s;\
    efetch -db protein -format fasta -id $1 | phobius 2> /dev/null | phobius2table -e 0.0101 ) | domain2architecture -e 0.0101 | architecture2table | padtable\
","")

-- aln_order_by_tree
set_shell_function("aln_order_by_tree"," \
    aln=$1;\
    tree=$2;\
    tmp=$(mktemp -p /tmp aln_order_by_tree.XXXXXX);\
    treeutil -of leaves $tree > $tmp;\
    aln2seqrows -r '\\t' $aln | tjoin -r '\\n' -i1 1 -i2 1 -1 1 -f1 '$F[1] ne \"ID\"' -f '$F[0]=\">$F[0]\";1' -p $tmp -;\
    rm -f $tmp\
","")

-- distinct_dom_per_prot
set_shell_function("distinct_dom_per_prot"," \
    arch=$1;\
    shift;\
    architecture2table $arch | tgroup -g 0 -a 1=count:distinct -i 0..1 \"$@\"\
","")

-- dom_count
set_shell_function("dom_count"," \
    arch=$1;\
    shift;\
    architecture2table $arch | tgroup -f '$F[2]=$F[1];1' -g 0 -g 1 -a 2=count -i 0..2 \"$@\"\
","")

-- dom_per_prot
set_shell_function("dom_per_prot"," \
    arch=$1;\
    shift;\
    architecture2table $arch | tgroup -g 0 -a 1=count -i 0..1 \"$@\"\
","")

-- tdesc
set_shell_function("tdesc"," \
    sep=$2;\
    if [ \"$sep\" == \"\" ]; then\
        sep='\\t';\
    fi;\
    head -n1 $1 | ttranspose -s $sep | tfilter -f 'unshift(@F,$H{c}++);1'\
","")

-- Whatis
whatis("Module: " .. mymfn)
whatis("Filename: " .. myfn)
whatis("Directory: " .. mydir)
whatis("Version: " .. myversion)
whatis("Category: data analysis ")
whatis("Description: Rotifer's shell setup.")
whatis("URL: https://github.com/leepbioinfo/rotifer ")

-- Help message
local helpMsg = [[
ROTIFER
=======

Rapid Open-source Tools and Infrastructure For data Exploration and Research
----------------------------------------------------------------------------

ROTIFER is a multi-language collection of high-level programming libraries
for the development of data analysis pipelines, mostly targeting problems in
comparative genomics and the computational analysis of biological sequences.

It also provides a collection of easy to use command line tools based on this framework.

Current installation directory: ]] .. mydir

