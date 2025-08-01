from operator import pos
from dash import Dash, dcc, html, Input, Output, ctx, callback, State
import dash_cytoscape as cyto
import networkx as nx
import numpy as np
import pandas as pd
from rotifer.devel.alpha import gian_func as gf

# Enable SVG export
cyto.load_extra_layouts()

def network_annotation(networkdf,
                     node_color=None,
                     node_shape='ellipse',
                     second_shape={},  
                     second_shape_list = [],
                     node_outline=None,
                     node_outline_color={},
                     node_outline_list=[],
                     edge_color=True,
                     unidirected_graph=True):

    import networkx as nx
    import pandas as pd

    networkdf['edge_total'] = networkdf.groupby(['source', 'target'])['edge_count'].transform('sum')

    if unidirected_graph:
        G = nx.from_pandas_edgelist(networkdf, edge_attr=['edge_type', 'edge_count', 'edge_total'])
    else:
        G = nx.from_pandas_edgelist(networkdf, edge_attr=['edge_type', 'edge_count', 'edge_total'], create_using=nx.MultiDiGraph())


    if node_color:
        nx.set_node_attributes(G, node_color, name="node_color") 
    else:
        node_color={item: "#D3D3D3" for item in G.nodes()}
        nx.set_node_attributes(G, node_color, name="node_color") 

    #Getting the size of nodes
    # In the future we can add more ways of getting the node size
    node_degree = dict(G.degree())
    nx.set_node_attributes(G, node_degree, name="degree") 
    

    #Getting the shape of nodes
    #Creating a dictionary to shape nodes
    if second_shape == {}:
        second_shape={item: "hexagon" for item in second_shape_list}
    elif isinstance(second_shape, str):
        second_shape={item: second_shape for item in second_shape_list}
    else:
        pass

    node_shape_dict  = {node: second_shape[node] if node in second_shape_list else node_shape for node in G.nodes()}
    nx.set_node_attributes(G, node_shape_dict, name="node_shape")


    # Getting the outline of nodes:
    if node_outline_color == {}:
        node_outline_color={item: "#000000" for item in node_outline_list}
    elif isinstance(node_outline_color, str):
        node_outline_color={item: node_outline_color for item in node_outline_list}
    else:
        pass
    if node_outline:
        node_outline_dict = {node: node_outline_color[node] if node in node_outline_list else node_color[node] for node in G.nodes()}
        nx.set_node_attributes(G, node_outline_dict, name="node_outline") 
    
    if edge_color:
        edge_colors=[]
        edge_alpha=[]
        edge_styles=[]
    #Getting the color of edges():
        if isinstance(G, nx.MultiDiGraph):
            edges_obj = G.edges(keys=True) 
            for source,target, key in edges_obj:
                if node_color[source] == node_color[target]:
                    edge_colors.append(node_color[source])
                    edge_styles.append('solid')
                    edge_alpha.append(0.7)
                else:
                    edge_colors.append('gray')
                    edge_styles.append('dotted')
                    edge_alpha.append(0.2)    
                   
            dict_edge_colors = { (u, v, key): color for (u, v, key), color in zip(edges_obj, edge_colors)}
            dict_edge_style = { (u, v, key): style for (u, v, key), style in zip(edges_obj, edge_styles)}
            dict_edge_alpha = { (u, v, key): alpha for (u, v, key), alpha in zip(edges_obj, edge_alpha)}

        else:    
            edges_obj = G.edges() 
            for source,target in edges_obj:
                if node_color[source] == node_color[target]:
                    edge_colors.append(node_color[source])
                    edge_styles.append('solid')
                    edge_alpha.append(0.7)
                else:
                    edge_colors.append('gray')
                    edge_styles.append('dotted')
                    edge_alpha.append(0.2)    
            dict_edge_colors = {k: v for k, v in zip(G.edges(), edge_colors)}
            dict_edge_style = {k: v for k, v in zip(G.edges(), edge_styles)}
            dict_edge_alpha = {k: v for k, v in zip(G.edges(), edge_alpha)}

        nx.set_edge_attributes(G, dict_edge_style, 'line_style')
        nx.set_edge_attributes(G, dict_edge_colors, 'edge_color')
        nx.set_edge_attributes(G, dict_edge_alpha, 'opacity')        

    return G

def calculate_positions(G, iterations):
    # Recalculate positions as needed
    pos1 = nx.spring_layout(G, iterations=iterations, seed=1)
    pos2 = nx.kamada_kawai_layout(G, pos=pos1)
    pos3 = nx.kamada_kawai_layout(G)
    pos4 = nx.circular_layout(G)
    return {"Spring Layout": pos1,
            "Kamada-Kawai from Spring layout": pos2,
            "Kamada-Kawai": pos3,
            "Circular": pos4}

def convert_to_cytoscape(G, positions, layout_key, scaling_factor=900, x_scale=1, y_scale=1):
    cy = nx.cytoscape_data(G)
    # Modify node positions and style
    for node in cy["elements"]["nodes"]:
        node_id = node["data"]["id"]
        if node_id in positions[layout_key]:
            p = positions[layout_key][node_id]
            node["position"] = {
                "x": int(p[1] * scaling_factor * x_scale),
                "y": int(p[0] * scaling_factor * y_scale),
            }
    return cy["elements"]["nodes"] + cy["elements"]["edges"]


def network_dash(G,pos_dict, xx, function):
    from dash import Dash, dcc, html, Input, Output, ctx, callback, State
    import dash
    import dash_cytoscape as cyto
    import networkx as nx
    from networkx.readwrite import json_graph
    import numpy as np
    import pandas as pd
    from rotifer.devel.alpha import net_functions
    from dash.exceptions import PreventUpdate


    # Serializing G to save it od dash  as dcc.store, that way we could manipulate in any way that we need
    G_serialized = json_graph.node_link_data(G.copy())

    app = dash.Dash(__name__)
    default_stylesheet = [
        {
            'selector': 'node',
            'style': {
                'content': 'data(name)',
                'background-color': 'data(node_color)',
                'border-color': 'data(node_outline)',
                "border-width": 3,
                'shape': 'data(node_shape)',
                'width': 'data(size)',
                'height': 'data(size)',
                'text-valign': 'bottom',
                'text-halign': 'center',
                'font-size': '12px', 
                'color': 'black'
            }
        },
        {
            'selector': 'edge',
            'style': {
                'line-color': 'data(edge_color)',
                'line-style': 'data(line_style)',
                'opacity': 'data(opacity)'
            }
        }
    ]
    
    
    app.layout = html.Div([
        dcc.Dropdown(
            id='layout-selector',
            options=[{'label': key, 'value': key} for key, value in pos_dict.items()],
            value=list(pos_dict.keys())[0],
            clearable=False,
            style={'width': '30%', 'margin-right': '10p'}
        ),
        html.Label("Spring Layout iterations:", style={'margin-right': '5px'}),
        dcc.Input(
                id='iterations-input',
                type='number',
                placeholder='Iterations',
                value=10,
                min=1,
                step=1,
                style={'width': '30px'}
            ),
        html.Label("Minimmun connection:", style={'margin-right': '5px'}),
        dcc.Input(
                id='connections-input',
                type='number',
                placeholder='min_connections',
                value=0,
                min=0,
                step=1,
                style={'width': '30px'}
            ),
        html.Label("x scale:", style={'margin-right': '5px'}),
        dcc.Input(
                id='x_scale',
                type='number',
                placeholder='x_scale',
                value=1,
                min=1,
                step=0.1,
                style={'width': '30px'}
            ),
        html.Label("y scale:", style={'margin-right': '5px'}),
        dcc.Input(
                id='y_scale',
                type='number',
                placeholder='y scale',
                value=1,
                min=1,
                step=0.1,
                style={'width': '30px'}
            ),
        html.Button("Update layouts", id="up_lay"),
        dcc.Input(id='new-group-key', type='text', placeholder='Group name'),
        dcc.Input(id='new-group-color', type='text', placeholder='Color hex'),
        html.Button('Add group', id='add-group-button', n_clicks=0),
    html.Div([
        html.Div([
            # Source option with radio and input
            html.Div([
                dcc.RadioItems(
                    id="select-mode",
                    options=[
                        {"label": "Source", "value": "source"},
                        {"label": "Target", "value": "target"}
                    ],
                    value="source",
                    labelStyle={'display': 'inline-block', 'margin-right': '15px'}
                )
            ], style={"margin-bottom": "10px"}),
            html.Div([
                html.Label("Source Node", htmlFor="select_source", style={"margin-right": "5px"}),
                dcc.Input(id="select_source", type="text", placeholder="Source Node", style={"width": "150px"})
            ], style={"display": "inline-block", "margin-right": "30px"}),
            html.Div([
                html.Label("Target Node", htmlFor="select_target", style={"margin-right": "5px"}),
                dcc.Input(
                    id="select_target",
                    type="text",
                    placeholder="Target Node",
                    style={"width": "150px", "margin-right": "10px"}
                    ),
                html.Button("Add Edge", id="add-edge-button", n_clicks=0, style={"margin-right": "5px"}),
                html.Button("Remove Edge", id="remove-edge-button", n_clicks=0)
                ], style={"display": "inline-block"})
        ])
    ]),

    html.Div(id="node-input-box"),


        html.Div([
            dcc.Dropdown(
                id="colorscheme",
                options=[],   # preenchido no callback
                value=None,   # será definido no callback também
                placeholder="Select color scheme",
                style={"flex": "0 0 10%", "margin-right": "10px"}
            ),
            dcc.Dropdown(
                id='group_selection',    
                options=[],
                value=[],
                multi=True,
                placeholder="Select groups",
                style={"flex": "1"}
            )
        ], style={"display": "flex", "align-items": "center", "margin-bottom": "10px"}),

        cyto.Cytoscape(
            id='cytoscape',
            layout={'name': 'preset'},
            style={'height': '85vh', 'width': '100%'},
            stylesheet=default_stylesheet,
            elements=[]  # Initially empty
        ),
        html.Label("Font Size:"),
        dcc.Input(
            id='font-size-input',
            type='number',
            value=12,  # Default font size
            min=0,
            max=24,
            step=0.5,
            style={'width': '30px'}
        ),
        html.Label("Edge width:"),
        dcc.Input(
            id='edge-size',
            type='number',
            value=1.5,  # Default font size
            min=0,
            max=10,
            step=0.1,
            style={'width': '30px'}
        ),
        html.Label("Outline width:"),
        dcc.Input(
            id='outline-size-input',
            type='number',
            value=1,  # Default font size
            min=0,
            max=20,
            step=0.1,
            style={'width': '30px'}
        ),
        html.Label("Node size:", style={'margin-right': '5px'}),
        dcc.Input(
            id='Node_scale',
            type='number',
            value=3,  # Default font size
            min=0,
            max=50,
            step=0.1,
            style={'width': '30px'}
        ),
        html.Label("Edge distances:",style={'margin-right': '5px'}),
        dcc.Input(
            id='scale-size-input',
            type='range',
            value=900,  # Default font size
            min=100,
            max=3000,
            step=100,
            style={'width': '90px'}
        ),
        dcc.Input(
            id='Saved_position_text',
            value='Pos_name',  # Default font size
            type='text',
            style={'width': '70px'}
        ),
        html.Button("Save Position", id="save-button", n_clicks=0),
        html.Button("Print Stored Data", id="print-button"),  # Button to trigger action
        html.Label("Save network as:"),
        html.Button("JPG", id="btn-get-jpg" ),
        html.Button("PNG", id="btn-get-png"),
        html.Button("SVG", id="btn-get-svg"),
        #### Store  to do the double clik
        dcc.Store(id="edge-color-cycle", data=[]),
        dcc.Store(id="last-edge-click", data={"id": None, "time": 0}),
        dcc.Store(id="last-node-click", data={"id": None, "time": 0}),
        dcc.Store(id="node-color-cycle", data=[]),
        ### Adding stroe for groups
        dcc.Store(id="xx-store", data=xx),
        dcc.Store(id="function-store", data=function),
        ### Saving position
        dcc.Store(id="pos-dict-store", data=pos_dict),
        #### Saving the networkX G object 
        dcc.Store(id="graph-store", data=G_serialized),

        ### Dummy store to update add and remove edges on fly
        dcc.Store(id="edge-modification-trigger", data=0),




        dcc.Store(id="saved-pos-store", data={}, clear_data=True, storage_type='memory'),
        html.Div(id="output"),  # Output for display
        
    ])
   


    @app.callback(
        Output('cytoscape', 'stylesheet'),
        [Input('font-size-input', 'value'),  # Font size input
         Input('outline-size-input', 'value'),
         Input('edge-size', 'value')]       # Edge size input
    )
    def update_stylesheet(font_size,outline_size, edge_size):
        return [
            {
                'selector': 'node',
                'style': {
                    'content': 'data(name)',
                    'background-color': 'data(node_color)',
                    'border-color': 'data(node_outline)',
                    'shape': 'data(node_shape)',
                    "border-width": f'{outline_size}px',
                    'width': 'data(size)',
                    'height': 'data(size)',
                    'text-valign': 'bottom',
                    'text-halign': 'center',
                    'font-size': f'{font_size}px',  # Dynamic font size
                    'color': 'black'
                }
            },
            {
                'selector': 'edge',
                'style': {
                    'line-color': 'data(edge_color)',
                    'line-style': 'data(line_style)',
                    'opacity': 'data(opacity)',
                    'width': f'{edge_size}px'  # Dynamic edge width
                }
            }
        ]
    """    
    @app.callback(
        Output('cytoscape', 'elements'),
        [Input('layout-selector', 'value'),
        Input('group_selection', 'value'),
        Input('connections-input', 'value'),
        Input('Node_scale', 'value'),
        Input('x_scale', 'value'),
        Input('y_scale', 'value'),
        Input('scale-size-input', 'value')]
    )
    def update_graph(selected_option, groups,min_connections,node_scale,x_scale,y_scale, scale_value):
        scale_value = int(scale_value)
       # y_scale = int(scale_value)
       # x_scale = int(scale_value)
        # Dynamically generate elements based on selected_option
        #### Node filter based on groups:
        filter_nodes = list({key: value for key, value in function.items() if value in groups}.keys())
        G_filtered = G.subgraph(filter_nodes)
        ### Edge filtered based in total Edge count
        filter_edges = [(u, v) for u, v, attrs in G_filtered.edges(data=True) if attrs.get('edge_total', 0) >= min_connections]
        subgraph = G_filtered.edge_subgraph(filter_edges).copy()
        subgraph.remove_nodes_from(list(nx.isolates(subgraph)))
        for edge in subgraph.edges(data=True):
            if function[edge[0]] == function[edge[1]]:
                edge[2]['line_style']='solid'
                edge[2]['opacity']=1
                edge[2]['edge_color']= xx[function[edge[0]]]
            else:
                edge[2]['line_style']='dotted'
                edge[2]['opacity']=0.3
                edge[2]['edge_color']= 'gray'

        ## Scale the node size:
        for node, attrs in subgraph.nodes(data=True):
            attrs['size'] = np.sqrt(attrs['degree']) * node_scale
            attrs['node_color'] = xx[function[node]]
            
        element = net_functions.convert_to_cytoscape(subgraph, pos_dict, selected_option,x_scale=x_scale, y_scale=y_scale, scaling_factor=scale_value)
        return element
    """ 

    @callback(
        Output("select_source", "value"),
        Output("select_target", "value"),
        Input("cytoscape", "tapNodeData"),
        State("select-mode", "value"),
    )
    def update_input_boxes(node_data,toggle):
        if node_data and "id" in node_data:
            node_id = node_data["id"]
            if  toggle == "source":
                return node_id, dash.no_update
            elif toggle == "target":
                return dash.no_update, node_id
        return dash.no_update, dash.no_update


    # Dicionários para armazenar estilos personalizados
    edge_style_overrides = {}
    node_style_overrides = {}

    @callback(
        Output("cytoscape", "elements"),
        Output("last-edge-click", "data"),
        Output("edge-color-cycle", "data"),
        Output("last-node-click", "data"),
        Output("node-color-cycle", "data"),
        Output("function-store", "data"),
        Input("cytoscape", "tapEdgeData"),
        Input("cytoscape", "tapNodeData"),
        Input("layout-selector", "value"),
        Input("group_selection", "value"),
        Input("connections-input", "value"),
        Input("Node_scale", "value"),
        Input("x_scale", "value"),
        Input("y_scale", "value"),
        Input("scale-size-input", "value"),
        Input("pos-dict-store", "data"),
        Input("edge-modification-trigger", "data"),
        State("cytoscape", "elements"),
        State("last-edge-click", "data"),
        State("edge-color-cycle", "data"),
        State("last-node-click", "data"),
        State("node-color-cycle", "data"),
        State('xx-store', 'data'),
        State("function-store", "data"),
        State("pos-dict-store", "data"),
        State("graph-store", "data"),
    )
    def update_elements_or_cycle(
        edge_data, node_data,
        selected_option, groups, min_connections, node_scale, x_scale, y_scale, scale_value,pos_dict_input,dummy_call,
        elements, last_click, color_cycle, last_node_click, node_color_cycle, func_color_dict, function_data, pos_dict_data, serial_graph
    ):

        from dash import callback, Input, Output, State, ctx, no_update
        import time
        import numpy as np
        triggered_id = ctx.triggered_id

        # Duplo clique em aresta
        if triggered_id == "cytoscape" and ctx.triggered[0]["prop_id"] == "cytoscape.tapEdgeData":
            if not edge_data or "id" not in edge_data:
                return no_update, last_click, color_cycle, last_node_click, node_color_cycle, function_data

            edge_id = f"{edge_data['source']}|{edge_data['target']}"
            now = time.time()
            double_click = (last_click["id"] == edge_id and (now - last_click["time"]) < 0.3)

            if not double_click:
                return no_update, {"id": edge_id, "time": now}, color_cycle, last_node_click, node_color_cycle, function_data

            # Captura estilos únicos de arestas visíveis
            if not color_cycle:
                seen = set()
                color_cycle = []
                for el in elements:
                    if "source" not in el.get("data", {}):
                        continue
                    color = el["data"].get("edge_color", "#ccc")
                    style = el["data"].get("line_style", "solid")
                    opacity = round(float(el["data"].get("opacity", 1)), 2)
                    key = (color, style, opacity)
                    if key not in seen:
                        seen.add(key)
                        color_cycle.append(key)

            # Aplica estilo seguinte com contador e salva override
            new_elements = []
            for el in elements:
                data = el.get("data", {})
                if f"{data.get('source')}|{data.get('target')}" == edge_id:
                    cycle_index = (data.get("cycle_index", -1) + 1) % len(color_cycle)
                    data["cycle_index"] = cycle_index
                    next_style = color_cycle[cycle_index]
                    data["edge_color"], data["line_style"], data["opacity"] = next_style

                    el["style"] = el.get("style", {})
                    el["style"].update({
                        "line-color": next_style[0],
                        "target-arrow-color": next_style[0],
                        "line-style": next_style[1],
                        "opacity": next_style[2],
                    })

                    edge_style_overrides[edge_id] = next_style

                new_elements.append(el)

            return new_elements, {"id": edge_id, "time": now}, color_cycle, last_node_click, node_color_cycle, function_data

        # Duplo clique em nó
        if triggered_id == "cytoscape" and ctx.triggered[0]["prop_id"] == "cytoscape.tapNodeData":
            if not node_data or "id" not in node_data:
                return no_update, last_click, color_cycle, last_node_click, node_color_cycle, function_data

            node_id = node_data["id"]
            now = time.time()
            double_click = (last_node_click["id"] == node_id and (now - last_node_click["time"]) < 0.3)

            if not double_click:
                return no_update, last_click, color_cycle, {"id": node_id, "time": now}, node_color_cycle, function_data
            """
            if not node_color_cycle:
                node_color_cycle = list(func_color_dict.values())
            """
            node_color_cycle = list(func_color_dict.values())

            new_elements = []
            for el in elements:
                data = el.get("data", {})
                if data.get("id") == node_id:
                    current_color = data.get("node_color", "#ccc")
                    try:
                        next_color = node_color_cycle[(node_color_cycle.index(current_color) + 1) % len(node_color_cycle)]
                    except ValueError:
                        next_color = node_color_cycle[0]

                    data["node_color"] = next_color
                    el["style"] = el.get("style", {})
                    el["style"]["background-color"] = next_color
                    node_style_overrides[node_id] = next_color
                    """
                    for func_name, color in func_color_dict.items():
                        if color == next_color:
                            function_data[node_id] = func_name
                            break
                    """        
                    # Procura o nome da função correspondente à cor
                    matched_group = next((func_name for func_name, color in func_color_dict.items() if color == next_color), None)

                    if matched_group is not None:
                        function_data[node_id] = matched_group
                    else:
                        # Cor não está no dicionário original, tenta buscar no xx-store atualizado
                        for func_name, color in func_color_dict.items():
                            if color == next_color:
                                function_data[node_id] = func_name
                                break

                new_elements.append(el)

            return new_elements, last_click, color_cycle, {"id": node_id, "time": now}, node_color_cycle, function_data

        # Atualização normal do grafo
        scale_value = int(scale_value)
        filter_nodes = list({key: value for key, value in function_data.items() if value in groups}.keys())
        stored_G = json_graph.node_link_graph(serial_graph)
        G_filtered = stored_G.subgraph(filter_nodes)

        filter_edges = [(u, v) for u, v, attrs in G_filtered.edges(data=True)
                        if attrs.get('edge_total', 0) >= min_connections]
        subgraph = G_filtered.edge_subgraph(filter_edges).copy()
        subgraph.remove_nodes_from(list(nx.isolates(subgraph)))

        for edge in subgraph.edges(data=True):
            uid = f"{edge[0]}|{edge[1]}"
            if uid in edge_style_overrides:
                color, style, opacity = edge_style_overrides[uid]
                edge[2]['line_style'] = style
                edge[2]['opacity'] = opacity
                edge[2]['edge_color'] = color
            elif function_data[edge[0]] == function_data[edge[1]]:
                edge[2]['line_style'] = 'solid'
                edge[2]['opacity'] = 1
                edge_group = function_data.get(edge[0], '')
                edge[2]['edge_color'] = func_color_dict.get(edge_group, '#D3D3D3')
            else:
                edge[2]['line_style'] = 'dotted'
                edge[2]['opacity'] = 0.3
                edge[2]['edge_color'] = 'gray'

        for node, attrs in subgraph.nodes(data=True):
            attrs['size'] = np.sqrt(attrs['degree']) * node_scale
            if node in node_style_overrides:
                attrs['node_color'] = node_style_overrides[node]
            else:
                attrs['node_color'] = func_color_dict[function_data[node]]

        new_elements = convert_to_cytoscape(
            subgraph, pos_dict_data, selected_option,
            x_scale=x_scale, y_scale=y_scale, scaling_factor=scale_value
        )
        return new_elements, last_click, color_cycle, last_node_click, node_color_cycle, function_data

    @app.callback(
        [Output("saved-pos-store", "data"),
        Output("layout-selector", "options"),
        Output("pos-dict-store","data"),
        Output("layout-selector", "value")],
        [Input("save-button", "n_clicks"),
        State("Saved_position_text", "value"),
        State('x_scale', 'value'),
        State('y_scale', 'value'),
        State('scale-size-input', 'value'),
        State("cytoscape", "elements"),
        State("pos-dict-store", "data"),],
        prevent_initial_call=True,
    )
    def save_positions(n_clicks, text,x_scale,y_scale, scale, n_elements, current_pos_dict):
        # Extract positions from the Cytoscape graph
        positions = {
            node["data"]["id"]: node.get("position", {})
            for node in n_elements if "position" in node
        }
        positions = {key: np.array([value['y']/scale/y_scale, value['x']/scale/x_scale]) for key, value in positions.items()}
        #global p
        current_pos_dict[text] = positions
        to_selector = [{'label': key, 'value': key} for key, value in current_pos_dict.items()]
        return positions, to_selector, current_pos_dict, text
    
    
    @app.callback(
        Output("output", "children"),
        Input("print-button", "n_clicks"),
        State("saved-pos-store", "data"),
        prevent_initial_call=True,  # Prevent execution on page load
    )
    def print_stored_data(n_clicks, toprint):
        
        print("Stored data:", toprint)  # Print to console
        return f"Printed data: {toprint}"  # Show feedback in the app
    
    
    @app.callback(
        Output("cytoscape", "generateImage"),
        [
            Input("btn-get-jpg", "n_clicks"),
            Input("btn-get-png", "n_clicks"),
            Input("btn-get-svg", "n_clicks"),
        ],
        prevent_initial_call=True
    )
    def get_image(get_jpg_clicks, get_png_clicks, get_svg_clicks):
        triggered_id = ctx.triggered_id
        file_type_map = {
            'btn-get-jpg': 'jpg',
            'btn-get-png': 'png',
            'btn-get-svg': 'svg'
        }
        ftype = file_type_map.get(triggered_id, 'jpg')
        action = 'download' if ftype == 'svg' else 'both'
        return {'type': ftype, 'action': action}
  
    ##### Adding dinamically groups to dropdown menu:
    @app.callback(
            Output("group_selection", "options"),
            Output("group_selection", "value"),
            Input("xx-store", "data"),
            Input("function-store", "data")
            )
    def update_group_dropdown(xx_data, function_data):
        all_keys = set(xx_data.keys()) | set(function_data.values())
        options = [{"label": k, "value": k} for k in sorted(all_keys)]
        return options, list(sorted(all_keys))

    @app.callback(
        Output('xx-store', 'data'),
        Input('add-group-button', 'n_clicks'),
        State('new-group-key', 'value'),
        State('new-group-color', 'value'),
        State('xx-store', 'data'),
        prevent_initial_call=True
    )
    def add_group(n_clicks, new_key, new_color, current_xx):
        if not new_key or not new_color:
            return dash.no_update
        current_xx[new_key] = new_color
        return current_xx



    @app.callback(
        Output("layout-selector", "options", allow_duplicate=True),
        Output("pos-dict-store","data", allow_duplicate=True),
        Input("up_lay", "n_clicks"),
        State('group_selection', 'value'),
        State('iterations-input', 'value'),
        State('connections-input', 'value'),
        State("pos-dict-store", "data"),
        State("graph-store", "data"),
        State("function-store", "data"),
        prevent_initial_call=True,  # Prevent execution on page load
    )
    def update_layouts(n_clicks, groups,Spring_lay_iteractions, min_connections, pos_dict, serial_graph, saved_func_dict):
        stored_G = json_graph.node_link_graph(serial_graph)
        filter_nodes = list({key: value for key, value in saved_func_dict.items() if value in groups}.keys())
        G_filtered = stored_G.subgraph(filter_nodes)
        ### Edge filtered based in total Edge count
        filter_edges = [(u, v) for u, v, attrs in G_filtered.edges(data=True) if attrs.get('edge_total', 0) >= min_connections]
        subgraph = G_filtered.edge_subgraph(filter_edges).copy()
        subgraph.remove_nodes_from(list(nx.isolates(subgraph)))
        newp = net_functions.calculate_positions(subgraph, Spring_lay_iteractions)
        pos_dict['updated Spring'] = newp['Spring Layout']
        pos_dict['updated Kamada from spring'] = newp['Kamada-Kawai from Spring layout']
        pos_dict['updated Kamada'] = newp['Kamada-Kawai']
        pos_dict['updated Circular'] = newp['Circular']
        to_selector = [{'label': key, 'value': key} for key, value in pos_dict.items()]
        return to_selector, pos_dict


    @callback(
            Output("graph-store", "data", allow_duplicate=True),
            Output("edge-modification-trigger", "data", allow_duplicate=True),
            Output("pos-dict-store", "data", allow_duplicate=True),
            Output("layout-selector", "options", allow_duplicate=True),
            Output("layout-selector", "value", allow_duplicate=True),
            Input("add-edge-button", "n_clicks"),
            Input("remove-edge-button", "n_clicks"),
            State("select_source", "value"),
            State("select_target", "value"),
            State("graph-store", "data"),
            State("edge-modification-trigger", "data"),
            State("cytoscape", "elements"),
            State("pos-dict-store", "data"),
            prevent_initial_call=True
            )
    def modify_edge_and_store_current(
            add_clicks, remove_clicks, source, target,
            graph_data, trigger_count, elements, pos_dict
            ):
        from dash.exceptions import PreventUpdate
        from dash import ctx
        import numpy as np

        if not source or not target:
            raise PreventUpdate

        G = json_graph.node_link_graph(graph_data)

        if ctx.triggered_id == "add-edge-button":
            G.add_edge(source, target, edge_type="manual", edge_count=1, edge_total=1)
        elif ctx.triggered_id == "remove-edge-button":
            try:
                G.remove_edge(source, target)
            except:
                pass

        # Save current Cytoscape positions as "Current" layout
        positions = {
                el["data"]["id"]: [el["position"]["y"] / 900, el["position"]["x"] / 900]
                for el in elements if "position" in el
                }
        pos_dict["Current"] = {k: np.array(v) for k, v in positions.items()}
        options = [{'label': key, 'value': key} for key in pos_dict.keys()]

        return json_graph.node_link_data(G), trigger_count + 1, pos_dict, options, "Current"



    app.run_server(jupyter_mode="external",debug=True)        
    if __name__ == '__main__':
        app.run_server(jupyter_mode="external",debug=True)        
 

def get_network_community(G,
                 community_to_color='leidein', weight='edge_total', resolution_parameter=0.05):
    import community
    import seaborn as sns
    import leidenalg as la
    import igraph as ig
    import matplotlib.pyplot as plt
    from matplotlib import cm as cm
    import numpy as np
    import networkx as nx

    #Community detection by Louvain
    partition = community.best_partition(G,weight=weight)
    #Df to easily map community to node:
    c = pd.DataFrame.from_dict(partition,orient='index').reset_index().rename(
       {'index': 'cluster', 0: 'Louvain'}, axis=1)
    #Leiden partition:
    part = la.find_partition(ig.Graph.from_networkx(G), la.ModularityVertexPartition, weights=weight)
    c['leidein'] = pd.Series(part.membership)
    #Leiden partition playing with resolution parameter:
    part2 = la.find_partition(ig.Graph.from_networkx(G), la.CPMVertexPartition, weights=weight, resolution_parameter = resolution_parameter)
    c['leidein_2'] = pd.Series(part2.membership)


    c.Louvain  = 'Louvain_' + c.Louvain.astype(str)
    c.leidein  = 'Leiden_' + c.leidein.astype(str)
    c.leidein_2 = 'Leiden_' + c.leidein_2.astype(str)

    community_dict = {"Louvain" : [c.set_index('cluster').Louvain.to_dict(),pd.Series(sns.color_palette('pastel',c.Louvain.nunique()).as_hex(), index=c.Louvain.unique()).to_dict()],
                      "Leidein":  [c.set_index('cluster').leidein.to_dict(),pd.Series(sns.color_palette('pastel',c.leidein.nunique()).as_hex(), index=c.leidein.unique()).to_dict()],
                      "Leidein_2": [c.set_index('cluster').leidein_2.to_dict(),pd.Series(sns.color_palette('pastel',c.leidein_2.nunique()).as_hex(), index=c.leidein_2.unique()).to_dict()]}

    return community_dict

def clique_jaccard_community(G,
                             min_jaccard_index=0.4,
                             resolution_parameter=0.5,
                             min_clique_size=10):
    import net_functions
    cliques = list(nx.find_cliques(G))
    cliques = pd.Series(cliques).rename('cliques').to_frame()
    cliques['clique_size'] = cliques.cliques.apply(lambda x: len(x))
    cliques = cliques.query('clique_size >= @min_clique_size')
    cliques['t']  = 1
    cliques.reset_index(inplace=True)
    cliques.rename({'index':'clique_n'}, axis=1, inplace=True)
    cm = cliques.merge(cliques, on='t')
    cm['uni'] = cm.apply(lambda x: set(x.cliques_x).union(set(x.cliques_y)), axis=1)
    cm['inter'] = cm.apply(lambda x: set(x.cliques_x).intersection(set(x.cliques_y)), axis=1)
    cm['Jaccard_index'] = cm.apply(lambda x: len(x.inter)/len(x.uni), axis=1)
    cm.rename({"clique_n_x" :'source', 'clique_n_y' :'target'}, axis=1, inplace=True)
    cm = cm[cm.Jaccard_index >= min_jaccard_index]
    G = nx.from_pandas_edgelist(cm, edge_attr='Jaccard_index')
    z = net_functions.get_network_community(G, weight='Jaccard_index', resolution_parameter=resolution_parameter)
    cm['source_leidein_resolution'] = cm.source.replace(z['Leidein_2'][0])
    cm['target_leidein_resolution'] = cm.target.replace(z['Leidein_2'][0])
    cm['source_leidein'] = cm.source.replace(z['Leidein'][0])
    cm['target_leidein'] = cm.target.replace(z['Leidein'][0])
    cm['source_louvain'] = cm.source.replace(z['Louvain'][0])
    cm['target_louvain'] = cm.target.replace(z['Louvain'][0])
    # result = {"Louvain" : cm.explode(
    return cm

def subgraph_stats(net,
                   subnet,
                   theme_dict,
                   theme_to_stats,
                   alternative = "two-sided",
                   return_df = False):
    from  scipy.stats import fisher_exact
    import pandas as pd
    to_test = pd.Series(list(subnet.nodes),name="node",dtype='object').to_frame()
    to_test['theme'] = to_test.node.map(theme_dict)
    total = pd.Series(list(net.nodes),name="node").to_frame()
    total['theme'] = total.node.map(theme_dict)
    sub_total = to_test.node.nunique()
    sub_theme = to_test.query('theme == @theme_to_stats').node.nunique()
    ttt = to_test.node.tolist()
    total = total[~total.node.isin(ttt)]
    net_total = total.node.nunique()
    net_theme = total.query('theme == @theme_to_stats').node.nunique()
    to_fisher = [[sub_theme, sub_total - sub_theme], [net_theme, net_total - net_theme]] 
    if return_df:
        print(fisher_exact(to_fisher, alternative=alternative))
        return pd.DataFrame(to_fisher, columns=['Theme', 'Total'], index=['Subnetwork', 'net - subnet'] )
    return fisher_exact(to_fisher,alternative=alternative)

def get_relations(domdf,
                  dom='full',
                  min_connections = 5,
                  self_relation = False,
                  to_yaml=False,
                  iteration=1,
                  filter_rename_yaml  = False):
    from rotifer.devel.alpha import gian_func as gf
    """ 
    domdf to network, it can create dictionary to create a YAML file or a dataframe.
    dom: list of domains to be used as query in the network, if full all domdf will be used
    min_connections: Minimum number of relationship (Operon or fusion) anode should have to be presente in the network.
    yaml: If decided to transform the network in dictionary to latter be exported as YAML file.
    Iteration: If a list of domains was supplied, the number of iteration in the network to collect other nodes.
    """
    domdf.dom = domdf.dom.str.strip()
    if filter_rename_yaml:
        import yaml
        with open(filter_rename_yaml, "r") as file:
            data = yaml.safe_load(file)
            data = pd.DataFrame(data).T
            data.loc[data.Display_name =="", 'Display_name'] = data.loc[data.Display_name ==""].index.tolist()
            rename_dict = data.Display_name.to_dict()
            domdf.dom = domdf.dom.replace(rename_dict)
            to_display = data.query('Include ==1').Display_name.unique().tolist()
            domdf = domdf.query('dom in @to_display')

    if isinstance(dom, str):
        dom = [dom]
    if dom[0] =='full':
        block_id = domdf.block_id.tolist()
    else:
        for x in range(iteration):
            block_id = domdf.query('dom in @dom').block_id.tolist()
            dom = domdf.query('block_id in @block_id').dom.unique().tolist()
    to_net = domdf.query('block_id in @block_id')[['pid', 'block_id', 'dom', 'strand']]
    to_net['block_pos'] = 1
    to_net['block_pos'] = to_net.groupby(['block_id'])['block_pos'].cumsum()
    to_net = to_net.merge(to_net, on ='block_id')
    to_net.loc[((to_net.strand_x == 1 ) & (to_net.block_pos_x < to_net.block_pos_y) & (to_net.pid_x != to_net.pid_y)), 'edge_type'] = 'Downstream'
    to_net.loc[((to_net.strand_x == 1 ) & (to_net.block_pos_x > to_net.block_pos_y) & (to_net.pid_x != to_net.pid_y)), 'edge_type'] = 'Upstream'
    to_net.loc[((to_net.strand_x == 1 ) & (to_net.block_pos_x > to_net.block_pos_y) & (to_net.pid_x == to_net.pid_y)), 'edge_type'] = 'Nterminal'
    to_net.loc[((to_net.strand_x == 1 ) & (to_net.block_pos_x < to_net.block_pos_y) & (to_net.pid_x == to_net.pid_y)), 'edge_type'] = 'Cterminal'
    to_net.loc[((to_net.strand_x == -1 ) & (to_net.block_pos_x < to_net.block_pos_y) & (to_net.pid_x != to_net.pid_y)), 'edge_type'] = 'Downstream'
    to_net.loc[((to_net.strand_x == -1 ) & (to_net.block_pos_x > to_net.block_pos_y) & (to_net.pid_x != to_net.pid_y)), 'edge_type'] = 'Upstream'
    to_net.loc[((to_net.strand_x == -1 ) & (to_net.block_pos_x > to_net.block_pos_y) & (to_net.pid_x == to_net.pid_y)), 'edge_type'] = 'Cterminal'
    to_net.loc[((to_net.strand_x == -1 ) & (to_net.block_pos_x < to_net.block_pos_y) & (to_net.pid_x == to_net.pid_y)), 'edge_type'] = 'Nterminal'

    to_net= to_net.query('edge_type.notna()')
    to_net = to_net[['block_id','dom_x','dom_y','edge_type']]
    to_net.columns = ['block_id','source','target','edge_type']
    if to_yaml:
        to_net['edge_count'] = to_net.groupby(['source', 'target']).edge_type.transform('count')
        to_net = to_net.query('edge_count >= @ min_connections').drop(['edge_count'], axis=1)
        yamldf = pd.pivot_table(to_net,
                            index='source',
                            columns='edge_type',
                            values = ['target', 'block_id'],
                            aggfunc={
                                'target':lambda x : gf.count_series(x, as_list=True),
                                'block_id': 'nunique'}
                            )
        yamldf = yamldf.loc[:,'block_id'].sum(axis=1).rename('total').to_frame().join(yamldf.loc[:,'target']).sort_values('total')
        yamldf.fillna('', inplace=True)
        yamldf.Nterminal = yamldf.Nterminal.apply(lambda x: x.split(' ') if isinstance(x, str) else x)
        yamldf.Cterminal = yamldf.Cterminal.apply(lambda x: x.split(' ') if isinstance(x, str) else x)
        yamldf.Downstream = yamldf.Downstream.apply(lambda x: x.split(' ') if isinstance(x, str) else x)
        yamldf.Upstream = yamldf.Upstream.apply(lambda x: x.split(' ') if isinstance(x, str) else x)
        yamldf['L'] = yamldf[['Cterminal', 'Downstream', 'Nterminal', 'Upstream']].apply(
                lambda row: [item for sublist in row for item in sublist if item != ''],
                axis=1
                )
        yamldf.L = yamldf.L.apply(len)
        yamldf = yamldf.query('L >0')
        yamldf = yamldf.drop(['L'], axis=1)

        return(yamldf.to_dict(orient='index'))
    else:
        to_net['edge_total'] = to_net.groupby(['source', 'target']).edge_type.transform('count')
        to_net = to_net.groupby(['source', 'target', 'edge_type']).agg(
                {'block_id' :'count',
                 'edge_total' : 'first'}).reset_index().sort_values('block_id')
        to_net = to_net.rename({'block_id':'edge_count'},axis=1)        
        to_net = to_net.query('edge_total >= @min_connections').drop(['edge_total'], axis=1)
        if self_relation:
            pass
        else:
            to_net = to_net[to_net.source != to_net.target]
        return(to_net)
