# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 09:09:48 2017

@author: mmoskon
"""

#import plotly.plotly as py
import plotly
from plotly.graph_objs import Line, Scatter, Marker, Data, Layout, Figure, Annotation, Annotations, Font, XAxis, YAxis
import networkx as nx
import make_pairs
#import model_data

def draw_plotly(model, met_names, pairs, fluxes):        
    if not pairs:
        return
    
    G = nx.DiGraph()
    G.add_edges_from(pairs)
    
    if fluxes and type(list(fluxes.values())[0]) == tuple:           
        fl = {}
        for i in fluxes:
            fl[i] = fluxes[i][1]
        producing = set(make_pairs.direct_producers(model, met_names, fl))
        consuming = set(make_pairs.direct_consumers(model, met_names, fl))            
    else:
        producing = set(make_pairs.direct_producers(model, met_names, fluxes))
        consuming = set(make_pairs.direct_consumers(model, met_names, fluxes))
    
    val_map = {'reaction': '#666699',
               'metabolite': '#666699',
               'producing': '#006600',
               'both': '#663300',
               'consuming': '#ff0000',
               'observed': "#0066cc",
               'path': '#ff66ff'}
    
    """
    val_map = {'reaction': 'b',
               'metabolite': 'c',
               'producing': 'g',
               'both': 'y',
               'consuming': 'r',
               'observed': 'm',
               'path': 'm'}
    """
    
    node_colors = []
    for node in G.nodes():          
        if node in met_names:
            color = val_map['observed']
        elif node in producing & consuming:
            color = val_map['both']
        elif node in producing:
            color = val_map['producing']
        elif node in consuming:
            color = val_map['consuming']
        elif node in model.reactions:
            color = val_map['reaction'] 
        else:
            color = val_map['metabolite']             
        node_colors.append(color)            
    
    node_types = ['observed metabolites' if node in met_names else 'producing/consuming reactions' if node in producing & consuming else 'producing reactions' if node in producing else 'consuming reactions' if node in consuming else 'reactions' if node in model.reactions else 'metabolites' for node in G.nodes()]
    
    #node_ids = {node:i for i, node in enumerate(G.nodes())}
    
    node_labels = {}
    for node in G.nodes():                             
        if node in model.reactions:
            f = fluxes[node]
            if type(f) == tuple:
                node_labels[node] = node + " :" + "{:.3f}\n{:.3f}".format(round(f[0],3),round(f[1],3))
            else:
                node_labels[node] = node + " :" + "{:.3f}".format(round(fluxes[node],3))
        else:
            node_labels[node] = node
           
    node_sizes = [20 if a in model.reactions else 10 for a in G.nodes()]
            
    #pos = nx.nx_pydot.graphviz_layout(G, prog="dot")
    pos = nx.fruchterman_reingold_layout(G) 
    
    N = G.nodes()
    E = G.edges()
    
    Xv=[pos[k][0] for k in pos]
    Yv=[pos[k][1] for k in pos]
    labels = [node_labels[l] for l in node_labels]
    #full_labels = [model.reactions[a] if a in model.reactions else model.metabolites[a] for a in G.nodes()]
    
    full_labels = []
    for i, node in enumerate(G.nodes()):                             
        if node in model.reactions:
            f = fluxes[node]
            if type(f) == tuple:
                full_labels.append(model.reactions[node] + " :" + "{:.3f}\n{:.3f}".format(round(f[0],3),round(f[1],3)))
            else:
                full_labels.append(model.reactions[node] +  " :" + "{:.3f}".format(round(fluxes[node],3)))
        else:
            full_labels.append(model.metabolites[node])
    
    
    
    
    Xed=[]
    Yed=[]
    for edge in E:
        Xed+=[pos[edge[0]][0],pos[edge[1]][0], None]
        Yed+=[pos[edge[0]][1],pos[edge[1]][1], None] 
    
    edge_trace = Scatter(x=Xed,
                           y=Yed,
                           mode='lines+markers',
                           line=Line(color='rgb(210,210,210)', width=1),
                           hoverinfo='none'
                           ) 
    
    traces = [edge_trace]
    """
    node_trace = Scatter(x=Xv[1:3],
                       y=Yv[1:3],
                       mode='markers',
                       name='dot',
                       marker=Marker(symbol='dot',
                                     size=node_sizes[1:3], 
                                     color= node_colors[1:3],
                                     line=Line(color='rgb(50,50,50)', width=0.5)
                                     ),
                       text=labels[1:3],
                       hoverinfo='text'
                       )
       
    
    traces.append(node_trace)
    
       
    node_trace = Scatter(x=Xv[3:],
                       y=Yv[3:],
                       mode='markers',
                       name='dot',
                       marker=Marker(symbol='dot',
                                     size=node_sizes[3:], 
                                     color= node_colors[3:],
                                     line=Line(color='rgb(50,50,50)', width=0.5)
                                     ),
                       text=labels[3:],
                       hoverinfo='text'
                       )
       
    
    traces.append(node_trace)
    """
    
    N_len = len(N)
    
    for i in range(N_len):
        traces.append(Scatter(x=[Xv[i]],
                       y=[Yv[i]],
                       mode='markers',
                       name=node_types[i],
                       marker=Marker(symbol='dot',
                                     size=[node_sizes[i]], 
                                     color= [node_colors[i]],                                     
                                     line=Line(color='rgb(50,50,50)', width=0.5)
                                     ),
                       text=[full_labels[i]],
                       hoverinfo='text'
                       ))
                   


    
            
    annot=""#"This networkx.Graph has the Fruchterman-Reingold layout<br>Code:"+"<a href='http://nbviewer.ipython.org/gist/empet/07ea33b2e4e0b84193bd'> [2]</a>"
    
    annotations=Annotations([       
           Annotation(
           showarrow=False, 
            text=labels[i],  
            xref='x1',     
            yref='y1',     
            x=Xv[i],#+1,  
            y=Yv[i],#+1,  
            xanchor='left',   
            yanchor='bottom',  
            font=Font(
            size=14
            )     
            )
            for i in range(N_len)
        ] +[Annotation(
           showarrow=False, 
            text=labels[0],  
            xref='x1',     
            yref='y1',     
            x=Xv[0],#+1,  
            y=Yv[0],#+1,  
            xanchor='left',   
            yanchor='bottom',  
            font=Font(
            size=14
            )     
            )]
            )
    #print(annotations)
    axis=dict(showline=False, # hide axis line, grid, ticklabels and  title
          zeroline=False,
          showgrid=False,
          showticklabels=False,
          title='' 
          )
    
    
    layout = Layout(
            title='',
            titlefont=dict(size=16),
            showlegend=False,
            #showlegend=True,
            hovermode='closest',
            margin=dict(b=20,l=5,r=5,t=40),          
            annotations=annotations,
            xaxis=XAxis(axis),
            yaxis=YAxis(axis))
    
    data1=Data(traces)#Data([edge_trace, node_trace])
    fig1=Figure(data=data1, layout=layout)
    fig1['layout']['annotations'][0]['text']=annot
    #py.iplot(fig1, filename='Coautorship-network-nx')
    plotly.offline.plot(fig1)  