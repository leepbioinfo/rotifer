def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def count_series(series, normalize=False, cut_off=False, count='domain'):
    ''' 
    Function to flatten a pd.Series.value_counts and count the number of domains recognized in the whole proteins set of proteins
    Mind that by defaul, it count domain not architecture. If one wants to count architecture, it should change the count option to 'architecture'.
    The normalize options is to write the results as frequency.
    The cut_off options is to print results only above a given threashold. If the normalize function is True, you should give the cutoff value as pct. 
    Usage:
    in: count_series(df.pfam, normalize=True, cut_off=10)
    out:'T6SS_HCP(15.15%), SP(11.36%)'
    in: df.groupby('c80e3').agg({'pfam':count_series})
                                                 pfam
    c80e3                                            
    AAN66278.1                                       
    ABO22067.1                           SP(1), TM(1)
    ABU79335.1  Gal_mutarotas_2(1), Glyco_hydro_31(1)
    ABV15410.1                           SP(1), TM(1)
    ABV15414.1                                  TM(1)

    '''
    l = []
    if count == 'architecture':
        s = series.dropna().value_counts(normalize=normalize)
    else:
        s = series.str.split('+').explode().dropna().value_counts(normalize=normalize)

    if cut_off:
        cut_off = cut_off/100
        s = s.where(lambda x: x>= cut_off).dropna()
    for y,z in s.iteritems():
        if normalize:
            l.append(f'{y}({100 * z:.2f}%)')
        else:
            l.append(f'{y}({z})')
    return ', '.join(l)

 def fetch_seq(seqs):
     import time
     import rotifer.core.functions as cf
     from rotifer.devel.alpha.sequence import sequence
     from rotifer.devel.alpha.gian_func import chunks
     if isinstance(seqs, list):
         if len(seqs) > 200:
             seqs = chunks(seqs, 200)
     else:
         seqs = list(seqs)
     seq_string = ''
     for x in seqs:
         f = Entrez.efetch(db = 'protein', rettype = 'fasta', retmode = 'text', id = ','.join(x), api_key = cf.loadAPI()).read()
         seq_string = seq_string + f
         time.sleep(1)
     return sequence(seq_string)

