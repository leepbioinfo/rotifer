def arch_2_svg(df, length=False, font_size=8):
    import yaml
    import numpy as np
    import seaborn as sns
    import pandas as pd
    import pygraphviz as pgv
    from pygraphviz.agraph import re
    import rotifer.interval.utils as riu
    from rotifer.devel.alpha.sequence import complete

    #loading the domain dict to select the domains that shoul be kept and rename it (if nescssary)
    #with open('../network/laks_domains.yaml', 'r') as f:
    #    domain_dict = yaml.load(f, Loader=yaml.SafeLoader)

    # As reference Nature's standard figure sizes are 89 mm (3.50 inches) wide (single column) and 183 mm (7.20472 inches) wide (double column) 
    # If a protein with 5000 aminoacids would fit a whole line of doublw colum figure, each aa would then use 0.001440944 inches.
    # Therefor we could use this number as conversion factor to draw the domains shapes.
    # Example : a domain of 300 aa would have a shape of 300 *0.001440944 = 0.4322832 inches
    #Therefore to use calculate the number of characters that fit a shape of given size we have to use the followin forumule:
    # [size of shape]/(font_point * [character point numeber])
    # Using the full_space_character one can easyly check 

    # Size in inches of each point in font size, to calculate how much fit n the shapes:
    font_point = float(0.013837)
    full_space_character = '━'
    aa_scale = float(0.001440944) # This scale would be enough for a protein of 5000 aa on a double column figure (seems better to scalonate *4)

    colors = sns.color_palette("pastel").as_hex() +sns.color_palette("deep").as_hex()  + sns.color_palette('muted').as_hex() + sns.color_palette('colorblind')
    shapes = ['box',
              'octagon',
              'circle',
              'hexagon',
              'ellipse',
              'diamond',
              'doubleoctagon',
              'tripleoctagon']
    color_code = pd.DataFrame(shapes, columns=['shapes']).join(pd.DataFrame(colors, columns=['colors']), how='cross')

    a = df.copy() 
    if add_length and 'length' not in df.columns:
        import rotifer.devel.beta.sequence as rdbs
        a['length'] = a.ID.map(rdbs.sequence(a.ID.drop_duplicates().tolist()).df.set_index('id').length.to_dict())


    a = complete.complete(a, values={'domain':'unk'})
    d = complete.complete(a, values={'domain':'unk'})

    a = a.reset_index()
    #a.rename({'index':'pid_order'}, axis=1, inplace=True)
    a['pid_order'] = a.ID.map({ k:v for (k,v) in zip(a.ID,range(0,len(a)))})
    a = a.reset_index(drop=True).reset_index()

    #Selecting the colors
    b = a.domain.value_counts().reset_index().join(color_code, how='left')
    cc = b.set_index('index').colors.to_dict()
    sc = b.set_index('index').shapes.to_dict()

    #Preparing to create different lines

    #m = a.groupby('pid')['index'].apply(lambda x: min(x))
    #li = list(a.groupby('pid')['index'].apply(lambda x: list(x)))
    #for x in range(len(li)):
    #    li[x].insert(0,a.iloc[li[x]].pid.unique()[0])

    #Creating the graph.
    A = pgv.AGraph()
    # adding nodes  some edges

    li =[]
    gb = a.groupby('pid_order')
    for x  in gb:
        l=[]
        for y,z in x[1].iterrows():
            A.add_node(z['index'],  label=z.domain, shape=sc[z['domain']], style='round,filled',  color='lightgray', fillcolor=cc[z['domain']])
            l.append(z['index'])    
        li.append(l)
        li.append(x[1].ID.drop_duplicates().tolist())
        A.add_node(x[1].ID.unique()[0],  label=x[1].ID.unique()[0], shape='none')
        
    for x in range(len(li)):
        A.add_subgraph(li[x], rank="same")

    # v = list(m)
    v = [li[x][0] for x in range(len(li))]
    [A.add_edge(v.pop(0), v[0], penwidth = 0) for x in range(len(v)-1)]
    A.graph_attr.update(nodesep= 0.02)
    A.graph_attr.update(ranksep= 0.05)
    A.draw("subgraph7.svg", prog="dot")

    d = d.query('domain !="unk"')
    d['size']  = d.end - d.start
    d['esc'] = d['size'] * (aa_scale*4)
    #d.domain = d.domain.replace(domain_dict['Domain'])
    d['dl'] = d.domain.str.len()
    d['f_space']= round(d.esc/(font_point*font_size))
    d = d[d.f_space >= 1]
    d.f_space = d.f_space.astype(int)
    d.loc[d.domain =='unk', 'domain'] = d.loc[d.domain =='unk', 'f_space'].apply(lambda x: x * full_space_character)
    d['pid_order'] = d.ID.map(d.ID.drop_duplicates().reset_index(drop=True).reset_index().set_index('ID')['index'].to_dict())
    d = d.reset_index(drop=True).reset_index()
    #if size_median:
    #     d['esc'] = d.groupby('domain').esc.transform('median')

    gb = d.groupby('pid_order')
    li =[]
    A = pgv.AGraph()
    for x in gb:
        l =[]
        for y,z in x[1].iterrows():
            if z.domain.__contains__(full_space_character):
                A.add_node(z['index'],  label=z.domain, width=z['esc'],fixedsize='true', shape='none',height='0.4', fontsize=font_size)
                l.append(z['index'])
            else:
                A.add_node(z['index'],  label=z.domain,style='rounded,filled', width=z['esc'],fixedsize='true',height='0.4', color='lightgray', fontsize=font_size)
                l.append(z['index'])
        li.append(l)
        li.append(x[1].ID.drop_duplicates().tolist())
        A.add_node(x[1].ID.unique()[0],  label=x[1].ID.unique()[0], shape='none')
        
    for x in range(len(li)):
        A.add_subgraph(li[x], rank="same")
        
    v = [li[x][0] for x in range(len(li))]
    [A.add_edge(v.pop(0), v[0], penwidth = 0) for x in range(len(v)-1)]
    A.graph_attr.update(nodesep= 0.02)
    A.graph_attr.update(ranksep= 0.05)
    A.draw("subgraph10.svg", prog="dot")

def operon_fig(df, domain_dict=False):
    ''' saves a svg figure of neighbourhood dataframe
    '''

    import pygraphviz as pgv
    from pygraphviz.agraph import re
    import yaml
    import numpy as np
    import pandas as pd
    fontsize = 8    
    height=0.4

    def split_domain(df, column='pfam',fill ='unk', domain_dict=domain_dict, strand=False, after=10, before=10, remove_tm=False):
        '''
        It will use the query as anchor to delect the given number of after and before genes to
        create a dataframe of domains.
        One can use a dcitionary of domain to select the ones that would be kept and rename it
        '''
        def reverse_domains(pid):
            def fix_direction(df):
                df = df.copy()
                q = df.query('query ==1').query('strand ==-1').block_id
                df.strand = np.where(df.block_id.isin(q), df.strand * -1, df.strand)
                return df
            pid = fix_direction(pid)
            strand = pid.strand.unique()[0]
            if strand == -1:
                pid.domp = pid.domp.iloc[::-1].to_list()
            return pid.reset_index(drop=True)

        domdf = df.reset_index(drop=True).select_neighbors(strand=strand, after=after, before=before).copy().reset_index(drop=True).query('type == "CDS"')
        domdf['dom'] = domdf[column].str.split('+')
        domdf = domdf.explode('dom')
        domdf['domp'] = pd.Series(np.where(domdf.pid == domdf.pid.shift(1), 1,0))
        domdf.domp = domdf.groupby(['pid','block_id'])['domp'].cumsum()

        if fill:
            domdf.dom = domdf.dom.fillna(fill)

        # removing tm sig
        if remove_tm:
            domdf['tmc'] = np.where(domdf.dom.isin(["TM","SIG"]),0,1)
            only_tm_sig = list(domdf.groupby('pid').tmc.sum().where(lambda x : x==0).dropna().index)
            domdf.dom = np.where(domdf.pid.isin(only_tm_sig), fill, domdf.dom)
            domdf = domdf.query('dom no in ["TM", "SIG"]')

        domdf = domdf.groupby('pid').apply(
                reverse_domains
                ).reset_index(
                        drop=True
                        ).sort_values([
            'nucleotide',
            'block_id',
            'start',
            'end',
            'domp'])

        if not domain_dict:
            return domdf
            
        #domdf = domdf[domdf.dom.isin(list(domain_dict.keys()))]
        domdf.dom = domdf.dom.replace(domain_dict)


        return domdf

    te = df.query('type =="CDS"').copy()
    te = split_domain(te)
    te['shape'] = 'rectangle'
    # Removing duplicates sequencial domain in same protein
    te = te.loc[~((te.pid.shift() == te.pid) & (te.dom.shift() == te.dom))].reset_index(drop=True)

    l = te.drop_duplicates(['block_id','pid'], keep='first').query('strand ==-1').index
    r = te.drop_duplicates(['block_id','pid'], keep='last').query('strand ==1').index
    te.loc[r,'shape'] = 'rarrow'
    te.loc[l,'shape'] = 'larrow'
    te['height'] = np.where(te['shape'].isin(['larrow', 'rarrow']),  height, None)
    te = te.reset_index(drop=True).reset_index()

    m = te.groupby('block_id')['index'].apply(lambda x: min(x))
    li = list(te.groupby('block_id')['index'].apply(lambda x: list(x)))

    A = pgv.AGraph()
    # add some edges
    for x,y in te.iterrows():
        A.add_node(y['index'],  label=y.dom, shape=y['shape'],style='filled',height=y['height'],  color='lightgray', fillcolor='lightgray', fontsize=fontsize)
        
    for x in range(len(li)):
        A.add_subgraph(li[x], rank="same")
        
    v = list(m)
    [A.add_edge(v.pop(0), v[0], penwidth = 0) for x in range(len(v)-1)]
    A.graph_attr.update(nodesep= 0.02)
    A.draw("subgraph3.svg", prog="dot")
    return print(A.string())

