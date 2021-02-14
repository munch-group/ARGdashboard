import plotly.graph_objects as go

import json
from textwrap import dedent as d

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State

import dash_bootstrap_components as dbc

#from app import app

# external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
# external_stylesheets = [dbc.themes.GRID]
external_stylesheets = [dbc.themes.BOOTSTRAP]

from jupyter_dash import JupyterDash
app = JupyterDash(__name__, external_stylesheets=external_stylesheets)

# app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

server = app.server
app.config.suppress_callback_exceptions = True

import json

###########################################e

import pandas as pd
import networkx as nx

import arg
from arg import Coalescent, Recombination, interval_sum, get_breakpoints, get_child_lineages, rescale_positions, marginal_arg, traverse_marginal, marginal_trees



# G = nx.random_geometric_graph(20, 0.125)

# edge_x = []
# edge_y = []
# for edge in G.edges():
#     x0, y0 = G.nodes[edge[0]]['pos']
#     x1, y1 = G.nodes[edge[1]]['pos']
#     edge_x.append(x0)
#     edge_x.append(x1)
#     edge_x.append(None)
#     edge_y.append(y0)
#     edge_y.append(y1)
#     edge_y.append(None)


# node_x = []
# node_y = []
# for node in G.nodes():
#     x, y = G.nodes[node]['pos']
#     node_x.append(x)
#     node_y.append(y)


# node_adjacencies = []
# node_text = []
# for node, adjacencies in enumerate(G.adjacency()):
#     node_adjacencies.append(len(adjacencies[1]))
#     node_text.append('# of connections: '+str(len(adjacencies[1])))

# # node_trace.marker.color = node_adjacencies
# # node_trace.text = node_text

###########################################


layout = html.Div(
    [
        # dbc.Row(
        #     [
        #         dbc.Col(
        #             [
        #                 dbc.Container(
        #                     [
        #                         html.H3("Ancestral recombination graph"),
        #                     ], fluid=True, style={'padding-left': 20,
        #                                         'padding-top': 10,
        #                                         'padding-bottom': 0,
        #                                         }               
        #                 ),
        #             ], width=8
        #         ),
        #         dbc.Col(
        #             [
        #                 dbc.Container(
        #                     [
        #                     dcc.Link(
        #                        html.H3("Dashboards"),
        #                        href='/',
        #                        style={'color': 'lightgrey'}
        #                     ),
        #                     ], fluid=True, style={'text-align': 'right',
        #                                         'padding-right': 40,
        #                                         'padding-top': 10,
        #                                         'padding-bottom': 0,
        #                                         'color': 'lightgrey',
        #                                         }               
        #                 ),
        #             ], width=4
        #         )
        #     ], justify='center', no_gutters=True,
        # ),

        # Hidden div inside the app that stores the intermediate value
        html.Div(id='intermediate-value', style={'display': 'none'}),

        # row for arg and marginal trees
        dbc.Row(
            [
                # column for arg 
                dbc.Col(
                    [
                        # arg
                        dbc.Container(
                            [
                                dbc.Container(
                                    [

                                        dbc.Row(
                                            [
                                                dbc.Col(
                                                    [
                                                        html.B("Simulation:"),
                                                        # "Simulation:",
                                                        dcc.Dropdown(
                                                            id='sim-dropdown',
                                                            options=[
                                                                {'label': "ARG", 'value': 'arg'},
                                                                # {'label': "SMC", 'value': 'smc'},
                                                                {'label': "SMC'", 'value': 'smcprime'}
                                                            ],
                                                            value='arg', searchable=False, clearable=False
                                                        ),
                                                    ], width=2
                                                ),                                                        
                                                dbc.Col(
                                                    [
                                                        html.B("Nr samples:"),
                                                        dcc.Dropdown(
                                                            id='samples-dropdown',
                                                            options=[
                                                                {'label': "3", 'value': 3},
                                                                {'label': "4", 'value': 4},
                                                                {'label': "5", 'value': 5}
                                                            ],
                                                            value=5, searchable=False, clearable=False,
                                                            style={
                                                                    # 'height': '20px', 
                                                                    # 'width': '80px', 
                                                                    'font-size': "0.85rem",
                                                                    # 'min-height': '1px',
                                                                    },
                                                        ),
                                                    ], width=2
                                                ),  
                                                dbc.Col(
                                                    [
                                                        html.B("Length:"),
                                                        dcc.Dropdown(
                                                            id='seqlen-dropdown',
                                                            options=[
                                                                {'label': "1kb", 'value': 1e+3},
                                                                {'label': "2kb", 'value': 2e+3},
                                                                {'label': "4kb", 'value': 4e+3}
                                                            ],
                                                            value=2e+3, searchable=False, clearable=False,
                                                            style={
                                                                    # 'height': '20px', 
                                                                    # 'width': '80px', 
                                                                    'font-size': "0.85rem",
                                                                    # 'min-height': '1px',
                                                                    },
                                                        ),
                                                    ], width=2
                                                ),                                                                                                        
                                                dbc.Col(
                                                    [ 
                                                        # html.Button('New simulation', id='new-arg-button')
                                                        dbc.Button('New', id='new-arg-button', 
                                                            color="primary", #size="sm", #outline=True,
                                                            style={'height': 35, 'font-size': "0.85rem"},
                                                            className="mr-1"
                                                            )
                                                    ], width=3
                                                ),    
                                                dbc.Col(
                                                    [
                                                        html.Div(id='arg-header'),

                                                        # dcc.Markdown(d("""
                                                        # **Ancestral recombination graph:**   
                                                        # Nodes are colored by amount of ancestral sequence.
                                                        # """), ),                    
                                                    ], width=3
                                                ),                                                                                                                                                    
                                                # dbc.Col(
                                                #     [ 
                                                #         dcc.Loading(
                                                #             id="loading-1",
                                                #             type="default",
                                                #             children=html.Div(id="arg-figure"),
                                                #         ),  
                                                #     ], width=1
                                                # )
                                                    
                                            ], justify="between", align="end", #no_gutters=True, 
                                        ),


                                        dcc.Graph(id='arg-figure',
                                                clear_on_unhover=True,
                                                figure={'layout': {
                                                            'height': 570,
                                                            # 'margin': {'l': 0, 'b': 0, 't': 0, 'r': 0},
                                                                }
                                                            },
                                                    ),

  

                                    ], className='pretty_container', fluid=True,
                                ),
                            ], style={'padding': 20}
                        ),
                    ], width=8, 
                ),

                # column for marginal trees
                dbc.Col(
                    [
                        dbc.Container(
                            [
                                dbc.Container(
                                    [
                                        dcc.Markdown(d("""
                                        **Marginal tree(s):** Hover over an ARG node.
                                        """), ),                    
                                        dcc.Graph(id='marginal-tree',
                                                    figure={'layout': {
                                                            # 'title': 'Marginal tree',
                                                            'height': 250,
                                                            # 'margin': {'l': 10, 'b': 0, 't': 0, 'r': 0},
                                                                }
                                                            },),
                                    ], className='pretty_container'
                                ),
                            ], style={'padding': 20, 'padding-left': 0, 'padding-bottom': 0}
                        ),
                        dbc.Container(
                            [
                                dbc.Container(
                                    [
                                        dcc.Markdown(d("""
                                        **Ancestral sequences:** Hover over an ARG node.
                                        """), ),                    
                                        dcc.Graph(id='ancestral-sequence',
                                                figure={'layout': {
                                                    'height': 250,
                                                        }
                                                        },),
                                    ], className='pretty_container'
                                ),            
                            ], style={'padding': 20, 'padding-left': 0}
                        ),            
                    ], width=4, 
                ),        
            ], no_gutters=True
        ),

        dbc.Row(
            [
                dbc.Col(
                    [
                        dbc.Container(
                            [
                                dbc.Container(
                                    [
                                        dbc.Container(
                                            [
                                                dcc.Markdown(d("""
                                                **Coalesce and recombination events:** 
                                                Slide to see progression of events.
                                                """)),
                                            ], 
                                        ),
                                        dbc.Container(
                                            [
                                                dcc.Slider(
                                                    id='event-slider',
                                                    min=0, max=40, value=0, 
                                                    marks={str(i): str(i) for i in range(0, 40)}, 
                                                    step=None,
                                                    ),
                                            ], style={'padding-bottom': 20}
                                        )
                                    ], className='pretty_container'
                                )
                            ], style={'padding': 20, 'padding-top': 0},
                        )
                    ], width=6
                ),  
                dbc.Col(
                    [
                        dbc.Container(
                            [
                                dbc.Container(
                                    [
                                        dbc.Container(
                                            [
                                                dcc.Markdown(d("""
                                                    **Recombination points:** 
                                                    Slide to see graph for only part of the sequence.
                                                """)),
                                            ]
                                        ),
                                        dbc.Container(
                                            [
                                                dcc.RangeSlider(
                                                    id='seq-slider',
                                                    min=0,
                                                    max=1000,
                                                    value=[0, 1000],
                                                    # step=None,
                                                    marks={0: '0', 1000: '1'},
                                                    pushable=30,
                                                )
                                            ], style={'padding-bottom': 20}
                                        ),
                                    ], className='pretty_container', 
                                ),
                            ], style={'padding': 20, 'padding-left': 0, 'padding-top': 0}
                        ),
                    ], width=6, align='start',
                ),
            ], no_gutters=True, 
        ),
    ], style={'padding': 20}
)


def arg_figure_data(nodes):

    traces = []

    edge_x = []
    edge_y = [] 
    # for lineage in get_parent_lineages(nodes, root=False):
    for lineage in get_child_lineages(nodes):
        # start
        edge_x.append(lineage.down.xpos)
        edge_y.append(lineage.down.height)
        # end
        edge_x.append(lineage.up.xpos)
        edge_y.append(lineage.up.height)
        # gap
        edge_x.append(None)
        edge_y.append(None)

    traces.append(dict(
        x=edge_x,
        y=edge_y,
        mode='lines',
        opacity=1,
        hoverinfo = 'skip',
        line={
            'color': 'grey',
        },
        name=''
    ))

    node_x = []
    node_y = []    
    node_text = []
    node_color = []
    for node in nodes:
        node_x.append(node.xpos)
        node_y.append(node.height)
        prop_ancestral = 1
        if type(node) is Coalescent:
            prop_ancestral = interval_sum(node.parent.intervals)
        elif type(node) is Recombination:
            prop_ancestral = interval_sum(node.child.intervals)
        node_text.append(f"Fraction ancestral: {round(prop_ancestral, 2)}<br>Event: {type(node).__name__}")

        node_color.append(prop_ancestral)

    traces.append(dict(
        x=node_x,
        y=node_y,
        text=node_text,
        # range_color=[0, 1],
        # cmin=0,
        # cmax=1,
        mode='markers',
        opacity=1,
        hoverinfo ='text',
        marker={
            'size': 10,
            'color': node_color,
            'cmin': 0,
            'cmax': 1,
            'line': {'width': 0.7, 'color': 'white'},
            # 'colorscale': 'Viridis',
            'colorscale': 'Rainbow',
            'colorbar': {'title': 'Fraction<br>ancestral<br>sequence',
                        'titleside': 'top',
                        'thickness': 15,
                        'len': 0.5,
                        # 'tickmode': 'array',
                        'tickvals': [0, 0.5, 1],
                        # 'ticktext': ['0', '1'],
                        'ticks': 'outside',
                        },
        },
        name=''
    ))

    return dict(data=traces,
                layout=dict(xaxis=dict(fixedrange=True, 
                                       range=[-0.1, 1.1], #title='Samples',
                                       showgrid=False, showline=False, 
                                       zeroline=False, showticklabels=False
                                       ),
                            yaxis=dict(fixedrange=True, 
                                       range=[-0.1, 1.1], #title='Time',
                                       showgrid=False, showline=False, 
                                       zeroline=False, showticklabels=False
                                       ),
                            hovermode='closest',
                            range_color=[0,1],
                            margin= {'l': 50, 'b': 20, 't': 20, 'r': 20},
                            transition = {'duration': 0},
                            showlegend=False
                            )
                )


def tree_figure_data(node_lists):

    traces = []

    edge_x = []
    edge_y = [] 
    node_x = []
    node_y = []    
    node_color = []

    for i, nodes in enumerate(node_lists):

        # for lineage in get_parent_lineages(nodes, root=False):
        for lineage in get_child_lineages(nodes):
            # start
            edge_x.append(lineage.down.xpos)
            edge_y.append(lineage.down.height)
            # end
            edge_x.append(lineage.up.xpos)
            edge_y.append(lineage.up.height)
            # gap
            edge_x.append(None)
            edge_y.append(None)

        for node in nodes:
            node_x.append(node.xpos)
            node_y.append(node.height)

            node_color.append(i/len(node_lists))

    traces.append(dict(
        x=edge_x,
        y=edge_y,
        mode='lines',
        opacity=1,
        hoverinfo = 'skip',
        line={
            'color': 'grey',
        },
        name=''
    ))

    traces.append(dict(
        x=node_x,
        y=node_y,
        mode='markers',
        opacity=1,
        hoverinfo ='text',
        marker={
            'size': 7,
            'color': node_color,
            'cmin': 0,
            'cmax': 1,
            'colorscale': 'Rainbow',
            'line': {'width': 0.3, 'color': 'white'},
        },
        name=''
    ))

    return dict(data=traces,
                layout=dict(xaxis=dict(fixedrange=True, 
                                       range=[-0.05, 1.05], #title='Samples',
                                       showgrid=False, showline=False, 
                                       zeroline=False, showticklabels=False
                                       ),
                            yaxis=dict(fixedrange=True, 
                                       range=[-0.1, 1.1], #title='Time',
                                       showgrid=False, showline=False, 
                                       zeroline=False, showticklabels=False
                                       ),
                            hovermode='closest',
                            range_color=[0,1],
                            margin= {'l': 50, 'b': 20, 't': 20, 'r': 20},
                            transition = {'duration': 0},
                            showlegend=False
                            )
                )


@app.callback(
    Output('arg-header', 'children'),
    [Input('new-arg-button', 'n_clicks')])
def update_header(n_clicks):

    # if n_clicks is None:
    #     return dcc.Markdown(d("""
    #             **Click the button**   
    #             to show a simulation.
    #             """))
    # else:
    #     return dcc.Markdown(d("""
    #                 **Simulation #{}:**   
    #                 Node color show proportion of ancestral material.
    #                 """.format(n_clicks)))

    if n_clicks is None:
        n_sim = 1
    else:
        n_sim = n_clicks + 1

    return dcc.Markdown(d("""
                **Simulation #{}:**   
                """.format(n_sim)))

@app.callback(Output('intermediate-value', 'children'), 
    [Input('new-arg-button', 'n_clicks'),
     Input('sim-dropdown', 'value'),
     Input('samples-dropdown', 'value'),
     Input('seqlen-dropdown', 'value')])
def new_data(n_clicks, sim, samples, length):

    nodes = arg.get_arg_nodes(L=length, n=samples, simulation=sim)
    rescale_positions(nodes)
    json_str = arg.arg2json(nodes)
    return json_str


    # arg = list(range(value))#[1,2,3]

    # return json.dumps(arg)
#     return cleaned_df.to_json(date_format='iso', orient='split')


@app.callback(
    [Output(component_id='event-slider', component_property='min'),
     Output(component_id='event-slider', component_property='max'),
     Output(component_id='event-slider', component_property='step'),
     Output(component_id='event-slider', component_property='value')],
    [Input('intermediate-value', 'children')])    
def update_event_slider(jsonified_data):
    if jsonified_data:
        nodes = arg.json2arg(jsonified_data)
    else:
        nodes = []

    nr_leaves = len([n for n in nodes if type(n) is arg.Leaf])
    nr_events = len(nodes)-nr_leaves
    return 0, nr_events, 1, nr_events

@app.callback(
    [Output(component_id='seq-slider', component_property='min'),
     Output(component_id='seq-slider', component_property='max'),
     Output(component_id='seq-slider', component_property='value'),
     Output(component_id='seq-slider', component_property='marks')],
    [Input('intermediate-value', 'children')])    
def update_seq_slider(jsonified_data):
    if jsonified_data:
        nodes = arg.json2arg(jsonified_data)
    else:
        nodes = []

    breakpoints = get_breakpoints(nodes)
    # print(breakpoints)
    marks = dict((b*1000, str(i+1)) for i, b in enumerate(breakpoints))
    # print(marks)
    # marks={0: '0', 500:'0.5', 1000: '1'},

    return 0, 1000, [0, 1000], marks


@app.callback(
    Output('arg-figure', 'figure'),
    [Input('intermediate-value', 'children'),
     Input('event-slider', 'value'),
     Input('seq-slider', 'value')])
def update_arg_figure(jsonified_data, event, interval):

    if jsonified_data:
        nodes = arg.json2arg(jsonified_data)

        interval = [i/1000 for i in interval]

        # Get marginal arg for interval
        marg_arg_nodes = marginal_arg(nodes, interval)
        # print(interval)
        # get only subset of events
        nr_leaves = len([n for n in nodes if type(n) is arg.Leaf])
        new_nodes = marg_arg_nodes[:nr_leaves+event]
    else:
        new_nodes = []

    return arg_figure_data(new_nodes)


@app.callback(
    Output('marginal-tree', 'figure'),
    [Input('intermediate-value', 'children'),
     Input('arg-figure', 'hoverData')])
def update_marg_tree_figure(jsonified_data, hover):

    # EVEN COOLER: THE ANCESTRAL SEQUENCES PANEL COULD SHOW EACH INTERVAL (between brakpoints)
    # OF THE FOCUS NODE AS COLORED BARS REPRESENTING ANCESTRAL SEQUENCE (and non-ancestral as gray)
    # AND THEN THE TREE PANEL COULD SHOW THE CORRESPONDING MARGINAL TREES IN THE SAME COLORS.
    # FOR COALESCENCE THAT SHOULD SHOW INTERVALS FOR THE PARENT LINEAGE
    # FOR RECOMBINATION IT SHOULD LSHOW INTERVALS FOR THE CHILD LINEAGE

    marg_tree_list = []
    if hover and jsonified_data:
        nodes = arg.json2arg(jsonified_data)
        focus_node_idx = hover['points'][0]['pointIndex']
        focus_node = nodes[focus_node_idx]

        if type(focus_node) is Recombination:
            intervals = focus_node.child.intervals
        else:
            intervals = focus_node.parent.intervals

        for interval in intervals:
            # get marginal arg under focus node
            new_nodes = traverse_marginal(focus_node, interval)

            new_nodes = list(new_nodes)
            new_nodes.sort(key=lambda x: x.height)
            # new_nodes = sorted(new_nodes, key=lambda x: int(repr(x)))

            marg_trees = marginal_trees(new_nodes)
            marg_tree_list.extend(marg_trees)

        nr_cols = len(marg_tree_list)

        # TODO: make rows and cols, and make sure maginal trees has same aspect ratio, and that the same nodes are placed at same height in each tree.
        space = 0.5
        
    
#         heights = set()
#         for tree in marg_tree_list:
#             if tree:
#                 heights.add(tree[-1].height)
#         same_height = len(heights) <= 1
        
        for i in range(nr_cols):
            tree = marg_tree_list[i]
#             if tree and same_height:
#                 rescale_positions(tree)
            if not tree:
                print('Empty tree...')
            for node in tree:
                # TODO: add some space between trees
                node.xpos = node.xpos/(nr_cols+(nr_cols-1)*space) + i/nr_cols
            marg_tree_list[i] = tree

    # TODO: Keep "dangling root" branch here

    if marg_tree_list:
        return(tree_figure_data(marg_tree_list))
        # return(arg_figure_data([x for sublist in marg_tree_list for x in sublist]))
    else:
        return(tree_figure_data([]))

import plotly.colors

def get_continuous_color(colorscale, intermed):
    """
    Plotly continuous colorscales assign colors to the range [0, 1]. This function computes the intermediate
    color for any value in that range.

    Plotly doesn't make the colorscales directly accessible in a common format.
    Some are ready to use:
    
        colorscale = plotly.colors.PLOTLY_SCALES["Greens"]

    Others are just swatches that need to be constructed into a colorscale:

        viridis_colors, scale = plotly.colors.convert_colors_to_same_type(plotly.colors.sequential.Viridis)
        colorscale = plotly.colors.make_colorscale(viridis_colors, scale=scale)

    :param colorscale: A plotly continuous colorscale defined with RGB string colors.
    :param intermed: value in the range [0, 1]
    :return: color in rgb string format
    :rtype: str
    """
    if len(colorscale) < 1:
        raise ValueError("colorscale must have at least one color")

    if intermed <= 0 or len(colorscale) == 1:
        return colorscale[0][1]
    if intermed >= 1:
        return colorscale[-1][1]

    for cutoff, color in colorscale:
        if intermed > cutoff:
            low_cutoff, low_color = cutoff, color
        else:
            high_cutoff, high_color = cutoff, color
            break

    # noinspection PyUnboundLocalVariable
    return plotly.colors.find_intermediate_color(
        lowcolor=low_color, highcolor=high_color,
        intermed=((intermed - low_cutoff) / (high_cutoff - low_cutoff)),
        colortype="rgb")

@app.callback(
    Output('ancestral-sequence', 'figure'),
    [Input('intermediate-value', 'children'),
     Input('arg-figure', 'hoverData')])    
def update_ancestral_seq_figure(jsonified_data, hover):

    figure_data = dict(
                    layout=dict(xaxis=dict(fixedrange=True,
                                           range=[-0.01, 1.01], #title='Samples',
                                           showgrid=False, showline=False, 
                                           zeroline=False, showticklabels=False
                                           ),
                               yaxis=dict(fixedrange=True,
                                           range=[0, 1], #title='Time',
                                           showgrid=False, showline=False, 
                                           zeroline=False, showticklabels=False
                                           ),
                               hovermode='closest',
                               margin= {'l': 0, 'b': 0, 't': 20, 'r': 0},
                               transition = {'duration': 0},
                               showlegend=False,
                               shapes=[],
                               )
                    )
    
    if hover and jsonified_data:
        nodes = arg.json2arg(jsonified_data)
        focus_node_idx = hover['points'][0]['pointIndex']
        focus_node = nodes[focus_node_idx]

        if type(focus_node) is Recombination:
            intervals = focus_node.child.intervals
        else:
            intervals = focus_node.parent.intervals
        
        segments = list()
        for interval in intervals:        

            new_nodes = traverse_marginal(focus_node, interval)
            new_nodes = list(new_nodes)
            new_nodes.sort(key=lambda x: x.height)

            breakpoints = get_breakpoints(new_nodes)
            starts = sorted(set([interval[0]] + breakpoints))
            ends = sorted(set([interval[1]] + breakpoints))

            segments.extend(list(zip(starts, ends)))

        colors, _ = plotly.colors.convert_colors_to_same_type(plotly.colors.sequential.Rainbow)
        colorscale = plotly.colors.make_colorscale(colors)

        shape = dict(type='rect', xref='x', yref='y', fillcolor=None, line= {'width': 1},
                    x0=0/3, y0=0.5, x1=1/3, y1=0.6)
        figure_data['layout']['shapes'].append(shape)
            
        for i, segment in enumerate(segments):
        
            # Add the shapes you need...
            color=get_continuous_color(colorscale, intermed=i/len(segments))

            offset = 0
            shape = dict(type='rect', xref='x', yref='y', fillcolor= color, line= {'width': 1},
                        x0=segment[0]/3, y0=0.5, x1=segment[1]/3, y1=0.6)
            figure_data['layout']['shapes'].append(shape)

    return figure_data
    
# @app.callback(
#     Output('ancestral-sequence', 'figure'),
#     [Input('event-slider', 'value'),
#      Input('seq-slider', 'value'),
#      Input('arg-figure', 'hoverData')])
# def update_ancestral_seq_figure(node, interval, hover):
#     figure_data = dict(
#                         layout=dict(xaxis=dict(fixedrange=True,
#                                                range=[0, 1], #title='Samples',
#                                                showgrid=False, showline=False, 
#                                                zeroline=False, showticklabels=False
#                                                ),
#                                    yaxis=dict(fixedrange=True,
#                                                range=[0, 1], #title='Time',
#                                                showgrid=False, showline=False, 
#                                                zeroline=False, showticklabels=False
#                                                ),
#                                    hovermode='closest',
#                                    margin= {'l': 0, 'b': 0, 't': 20, 'r': 0},
#                                    transition = {'duration': 0},
#                                    showlegend=False,
#                                    shapes=[],
#                                    )
#                         )
#     if hover is None:
#         return figure_data
#     else:
#         for i in range(1):
            
#             # Add the shapes you need...
#             color='blue'
            
#             shape = dict(type='rect', xref='x', yref='y', fillcolor= color, line= {'width': 0},
#                         x0=0.2, y0=0.2, x1=0.6, y1=0.6)
#             figure_data['layout']['shapes'].append(shape)
#         return figure_data


# @app.callback(
#     Output('hover-data', 'children'),
#     [Input('arg-figure', 'hoverData')])
# def display_hover_data(hoverData):
#     return json.dumps(hoverData, indent=2)


# @app.callback(
#     Output('click-data', 'children'),
#     [Input('arg-figure', 'clickData')])
# def display_click_data(clickData):
#     return json.dumps(clickData, indent=2)

app.layout = layout

if __name__ == '__main__':
    app.run_server(debug=True)