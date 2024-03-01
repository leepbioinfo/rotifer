import pygraphviz as pgv
from pygraphviz.agraph import re
import yaml
import numpy as np
import seaborn as sns
import yaml
import pandas as pd



#loading the domain dict to select the domains that shoul be kept and rename it (if nescssary)
with open('../network/laks_domains.yaml', 'r') as f:
    domain_dict = yaml.load(f, Loader=yaml.SafeLoader)
    

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


# Load architecture files
with open('./selected_arch_operon.yaml', 'r') as f:
    
    ttt = yaml.load(f, Loader=yaml.SafeLoader)
    
a = pd.DataFrame()
for x in ttt['Architectures'].keys():
    a1 = pd.DataFrame()
    for y in ttt['Architectures'][x].keys():
        a2 =  pd.DataFrame(ttt['Architectures'][x][y]).T
        a2['y'] = y
        a2['x'] = x
        a1 = pd.concat([a1,a2])    
    a = pd.concat([a, a1])
    
a.reset_index(inplace=True)
a.columns =['pid', 'assembly', 'arch', 't1','t2']

a['domain'] = a.arch.str.split('+')
a = a.reset_index()
a.rename({'index':'pid_order'}, axis=1, inplace=True)
a = a.explode('domain')
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
        A.add_node(z['index'],  label=z.domain, shape=sc[z['domain']],style='rounded,filled',  color='lightgray', fillcolor=cc[z['domain']])
        l.append(z['index'])    
    li.append(l)
    li.append(x[1].pid.drop_duplicates().tolist())
    A.add_node(x[1].pid.unique()[0],  label=x[1].pid.unique()[0], shape='none')
    
for x in range(len(li)):
    A.add_subgraph(li[x], rank="same")
    
# v = list(m)
v = [li[x][0] for x in range(len(li))]
[A.add_edge(v.pop(0), v[0], penwidth = 0) for x in range(len(v)-1)]
A.graph_attr.update(nodesep= 0.02)
A.graph_attr.update(ranksep= 0.05)
A.draw("subgraph7.svg", prog="dot")


import rotifer.interval.utils as riu
import rotifer.devel.alpha.sequence.complete as complete


    
font_size= 8
d = complete(a)
if not unk:
   d = d.query('domain !="unk"')
d['size']  = d.end - d.start
d['esc'] = d['size'] * (aa_scale*4)
d.domain = d.domain.replace(domain_dict['Domain'])
d['dl'] = d.domain.str.len()
d['f_space']= round(d.esc/(font_point*font_size))
d = d[d.f_space >= 1]
d.f_space = d.f_space.astype(int)
d.loc[d.domain =='unk', 'domain'] = d.loc[d.domain =='unk', 'f_space'].apply(lambda x: x * full_space_character)
d['pid_order'] = d.pid.map(d.pid.drop_duplicates().reset_index(drop=True).reset_index().set_index('pid')['index'].to_dict())
d = d.reset_index(drop=True).reset_index()
if size_median:
     d['esc'] = d.groupby('domain').esc.transform('median')



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
    li.append(x[1].pid.drop_duplicates().tolist())
    A.add_node(x[1].pid.unique()[0],  label=x[1].pid.unique()[0], shape='none')
    
for x in range(len(li)):
    A.add_subgraph(li[x], rank="same")
    
v = [li[x][0] for x in range(len(li))]
[A.add_edge(v.pop(0), v[0], penwidth = 0) for x in range(len(v)-1)]
A.graph_attr.update(nodesep= 0.02)
A.graph_attr.update(ranksep= 0.05)
A.draw("subgraph10.svg", prog="dot")
