from dash import html, dcc, callback, Output, Input, State, dash_table
from algorithms import *
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import dash_bio as bio
import plotly.figure_factory as ff
from copy import deepcopy

tversky_parametrs = html.Div([
                        dbc.Label("Weight a", className="ms-2 mb-2"),
                        dbc.Input(id="a-input", type="text", placeholder="Value between 0 and 1"),
                        dbc.Label("Weight b", className="mt-2 ms-2 mb-2"),
                        dbc.Input(id="b-input", type="text", placeholder="Value between 0 and 1")
                    ], id="tversky-parameters-content")

similarity_options = dbc.Select(options = [
                    {'label': 'Tanimoto', 'value': 'Tanimoto'},
                    {'label': 'Dice', 'value': 'Dice'},
                    {'label': 'Cosine', 'value': 'Cosine'},
                    {'label': 'Sokal', 'value': 'Sokal'},
                    {'label': 'Russel', 'value': 'Russel'},
                    {'label': 'Kulczynski', 'value': 'Kulczynski'},
                    {'label': 'McConnaughey', 'value': 'McConnaughey'},
                    {'label': 'Tversky', 'value': 'Tversky'}
                ], 
                id="similarity-coefficient",
                className="mb-2")

rdkit = html.Div([
            dbc.Alert("Fingerprint parameters", className="mb-3", color="info"),
            html.Div([
                dbc.Label("Min Path Length", className="mt-2 ms-2 mb-2"),
                dcc.Slider(id="min-path-slider", min=1, max=9, step=1, value=1, className="mb-3"),
                dbc.Label("Max Path Length", className="ms-2 mb-2"),
                dcc.Slider(id="max-path-slider", min=2, max=10, step=1, value=2,
                           marks={i: str(i) for i in range(2,11)}, className="mb-3"),
                dbc.Label("Fingerprint size", className="ms-2 mb-2"),
                dcc.Slider(id="fps-slider-rdkit", min=1024, max=4096, step=None,
                           marks = {1024:"1024", 2048:"2048", 4096:"4096"},
                           value=2048, included=False, className="mb-3"),
            ], 
            className="mb-3")
        ], 
        hidden=True,
        id="rdkit-parameters")

morgan = html.Div([
            dbc.Alert("Fingerprint parameters", className="mb-3", color="info"),
            html.Div([
                dbc.Label("Radius Length", className="mt-2 ms-2 mb-2"),
                dcc.Slider(id="radius-slider", min=1, max=10, step=1, value=1, className="mb-3"),
                dbc.Label("Fingerprint size", className="ms-2 mb-2"),
                dcc.Slider(id="fps-slider-morgan", min=1024, max=4096, step=None,
                           marks = {1024:"1024", 2048:"2048", 4096:"4096"},
                           value=2048, included=False, className="mb-3"),
            ], 
            className="mb-3")
        ], 
        hidden=True,
        id="morgan-parameters")

atompair = html.Div([
            dbc.Alert("Fingerprint parameters", className="mb-3", color="info"),
            html.Div([
                dbc.Label("Fingerprint size", className="ms-2 mb-2"),
                dcc.Slider(id="fps-slider-atompairs", min=1024, max=4096, step=None,
                           marks = {1024:"1024", 2048:"2048", 4096:"4096"},
                           value=2048, included=False, className="mb-3"),
            ], 
            className="mb-3")
        ], 
        hidden=True,
        id="atompairs-parameters")

maccs = html.Div([
            dbc.Alert("Fingerprint parameters", className="mb-3", color="info"),
            html.Div([], 
            className="mb-3")
        ],
        hidden=True, 
        id="maccs-parameters")

layout = dbc.Container([
    dcc.Store(id='fingerprint-store-data'),
    dcc.Store(id='data-frame-data'),
    dbc.Row([
        dbc.Navbar(
            html.H2('Molecular Similarity Visualiser', className="ms-3 mt-2"), 
            color="primary",
            className="mb-3")
    ]),

    dbc.Row([
        dbc.Col([
            dbc.Card([
                html.H4("Welcome to the Dashboard",
                        className="text-left mb-3"),
                dbc.Textarea(id="textarea-input",
                             className="mb-3", 
                             placeholder="SMILES strings separated by comma"),
                html.Div([
                    dbc.Label("Select fingerprint type"),
                    dbc.Select(options = [
                        {'label': 'RDKit', 'value': 'RDKit'},
                        {'label': 'AtomPairs', 'value': 'AtomPairs'},
                        {'label': 'Morgan', 'value': 'Morgan'},
                        {'label': 'MACCS Keys', 'value': 'MACCS Keys'}
                    ],
                    id="fingerprint-type",
                    className="mb-3")]),
                dbc.Collapse(
                    children=html.Div([
                        rdkit, 
                        morgan,
                        atompair,
                        maccs
                    ], id="fingerprint-parameters-content"),
                    is_open=False,
                    id="fingerprint-parameters-collapse",
                ),
                dbc.Button('Submit', id='submit-button', color='primary', className="mt-3 mb-3", n_clicks=0),
                html.Div(id="generation-alert"),
                html.Div(id="validation-alert")
            ], body = True)
        ], width = 4),
        dbc.Col([
            html.Div(id="visuals-column", style={"overflow": "hidden"})
        ])
    ]),
], fluid=True)



@callback(
    Output('fingerprint-parameters-collapse', 'is_open'),
    Output('rdkit-parameters', 'hidden'),
    Output('morgan-parameters', 'hidden'),
    Output('atompairs-parameters', 'hidden'),
    Output('maccs-parameters', 'hidden'),
    Output('rdkit-parameters', 'children'),
    Output('morgan-parameters', 'children'),
    Output('atompairs-parameters', 'children'),
    Output('maccs-parameters', 'children'),
    Input('fingerprint-type', 'value'),
    prevent_initial_call = True
)
def update_fingerprint_parameters(fingerprint_type: str):
    hide_rdkit = hide_morgan = hide_maccs = hide_atompairs = True
    rdkit_children = None
    morgan_children = None
    atompairs_children = None
    maccs_children = None

    children = html.Div([
        dbc.Label("Select similarity coefficient", className="ms-2 mb-2"),
        similarity_options,
        dbc.Collapse(
            id="tversky-parameters-collapse",
            is_open=False,
            children=tversky_parametrs
        ),
    ])

    if fingerprint_type == "RDKit":
        hide_rdkit = False
        rdkit_children = html.Div([
            dbc.Alert("Fingerprint parameters", className="mb-3", color="info"),
            html.Div([
                dbc.Label("Min Path Length", className="mt-2 ms-2 mb-2"),
                dcc.Slider(id="min-path-slider", min=1, max=9, step=1, value=1, className="mb-3"),
                dbc.Label("Max Path Length", className="ms-2 mb-2"),
                dcc.Slider(id="max-path-slider", min=2, max=10, step=1, value=2,
                           marks={i: str(i) for i in range(2,11)}, className="mb-3"),
                dbc.Label("Fingerprint size", className="ms-2 mb-2"),
                dcc.Slider(id="fps-slider-rdkit", min=1024, max=4096, step=None,
                           marks = {1024:"1024", 2048:"2048", 4096:"4096"},
                           value=2048, included=False, className="mb-3"),
                children
            ], 
            className="mb-3"),
        ])

    elif fingerprint_type == "Morgan":
        hide_morgan = False
        morgan_children = html.Div([
            dbc.Alert("Fingerprint parameters", className="mb-3", color="info"),
            html.Div([
                dbc.Label("Radius Length", className="mt-2 ms-2 mb-2"),
                dcc.Slider(id="radius-slider", min=1, max=10, step=1, value=1, className="mb-3"),
                dbc.Label("Fingerprint size", className="ms-2 mb-2"),
                dcc.Slider(id="fps-slider-morgan", min=1024, max=4096, step=None,
                           marks = {1024:"1024", 2048:"2048", 4096:"4096"},
                           value=2048, included=False, className="mb-3"),
                children
            ], 
            className="mb-3")
        ])

    elif fingerprint_type == "AtomPairs":
        hide_atompairs = False
        atompairs_children = html.Div([
            dbc.Alert("Fingerprint parameters", className="mb-3", color="info"),
            html.Div([
                dbc.Label("Fingerprint size", className="ms-2 mb-2"),
                dcc.Slider(id="fps-slider-atompairs", min=1024, max=4096, step=None,
                           marks = {1024:"1024", 2048:"2048", 4096:"4096"},
                           value=2048, included=False, className="mb-3"),
                children
            ], 
            className="mb-3"),
        ])


    elif fingerprint_type == "MACCS Keys":
        hide_maccs = False
        maccs_children = html.Div([
            dbc.Alert("Fingerprint parameters", className="mb-3", color="info"),
            html.Div([
                children
            ], 
            className="mb-3")
        ])
    
    return True, hide_rdkit, hide_morgan, hide_atompairs, hide_maccs, rdkit_children, morgan_children, atompairs_children, maccs_children


@callback(
    Output("max-path-slider", "min"),
    Output("max-path-slider", "marks"),
    Output("max-path-slider", "value"),
    Input("min-path-slider", "value"),
    State("max-path-slider", "value")
)
def update_max_path_slider(min_val: int, current_max_val: int):
    new_min = min_val + 1
    new_max = 10

    marks = {i: str(i) for i in range(new_min, new_max + 1)}

    if current_max_val < new_min:
        current_max_val = new_min

    return new_min, marks, current_max_val


@callback(
    Output('tversky-parameters-collapse', 'is_open'),
    Input('similarity-coefficient', 'value'),
    suppress_callback_exceptions=True
)
def show_tversky_parameters(similarity_coefficient: str):
    if similarity_coefficient == "Tversky":
        return True
    else:
        return False


@callback(
    Output('validation-alert', 'children'),
    Output('visuals-column', 'children'),
    Output('fingerprint-store-data', 'data'),
    Output('data-frame-data', 'data'),
    Output('generation-alert', 'children'),
    Input('submit-button', 'n_clicks'),
    State('fingerprint-type', 'value'),
    State('similarity-coefficient', 'value'),
    State('textarea-input', 'value'),
    State('min-path-slider', 'value'),
    State('max-path-slider', 'value'),
    State('fps-slider-rdkit', 'value'),
    State('fps-slider-atompairs', 'value'),
    State('radius-slider', 'value'),
    State('fps-slider-morgan', 'value'),
    State('a-input', 'value'),
    State('b-input', 'value'),
    prevent_intial_call=True
)
def submit_form(n_clicks: int, fingerprint_type:str=None, similarity_coefficient:str=None, text_value:str=None, min_path:int=None, 
                max_path:int=None, fps_rdkit:int = None, fps_atompairs:int = None, radius:int = None, fps_morgan:int = None, weight_a:str=None, weight_b:str=None):
    
    data = {
        'min_path': min_path,
        'max_path': max_path,
        'fps_rdkit': fps_rdkit,
        'fps_atompairs': fps_atompairs,
        'radius': radius, 
        'fps_morgan': fps_morgan,
        'a': weight_a,
        'b': weight_b
    }

    if n_clicks > 0:
        if not fingerprint_type:
            return dbc.Alert("You must choose a fingerprint type!", className="mb-3", color="warning"), None, None, None, None
        if not text_value or text_value.strip() == "":
            return dbc.Alert("You must input at least two SMILES strings!", className="mb-3", color="warning"), None, None, None, None
        if not similarity_coefficient:
            return dbc.Alert("You must choose a similarity coefficient!", className="mb-3", color="warning"), None, None, None, None
        

        if weight_a is not None and weight_b is not None:
            if float(weight_a) < 0 or float(weight_a) > 1 or float(weight_b) < 0 or float(weight_b) > 1:
                return dbc.Alert("Tversky parameters a and b must be between 0 and 1", className="mb-3", color="warning"), None, None, None, None

        if check_textarea_input(text_value) == False:
            return dbc.Alert("There exists an invalid SMILES string!", className="mb-3", color="warning"), None, None, None, None
        elif len(check_textarea_input(text_value)) < 2:
            return dbc.Alert("There must be at least 2 SMILES strings!", className="mb-3", color="warning"), None, None, None, None
        else:
            data_frame_generator = DataFrameGenerator(smiles=check_textarea_input(text_value),
                                                    generation_strategy=fingerprint_type,
                                                    similarity_strategy=similarity_coefficient,
                                                    data=data)
            
            df = data_frame_generator.get_data_frame()
            fingerprint_dict = data_frame_generator.get_fingerprint_indices()

            head_length = 0
            if len(df.columns) > 10:
                head_length = 10
            else:
                head_length = len(df.columns)

            
            
            masked_df = df.where(np.triu(np.ones(df.shape), k=1).astype(bool)).head(head_length)
            similarity_series = masked_df.stack()
            top_similar = similarity_series.sort_values(ascending=False)
            top_similar_df = top_similar.reset_index()
            top_similar_df.columns = ['Molecule 1', 'Molecule 2', 'Similarity']

            swapped_df = top_similar_df.rename(columns={
                                                        'Molecule 1': 'Molecule 2',
                                                        'Molecule 2': 'Molecule 1'
                                                    }
                                                )[['Molecule 1', 'Molecule 2', 'Similarity']]
            
            full_symetric_df = pd.concat([top_similar_df, swapped_df], ignore_index=True)
            full_symetric_df = full_symetric_df.sort_values(by='Similarity', ascending=False)


            heatmap = go.Figure(
                data=go.Heatmap(
                    z=df.values,
                    x=df.columns,
                    y=df.index,
                    colorscale='RdBu',
                    colorbar=dict(title=f'{similarity_coefficient} Similarity'),
                    text=df.values,
                    texttemplate="%{text: .2f}",
                    hovertemplate="Similarity: %{z: .2f}<br>y: %{y}<br>x: %{x}<extra></extra>"
                )
            )

            smiles_list=df.index.astype(str).tolist()
            dendrogram = ff.create_dendrogram(
                df.values,
                labels=smiles_list,
                orientation='left')
            dendrogram.update_traces(hoverinfo='x+y')
            dendrogram.update_layout(
                yaxis=dict(showticklabels=False),
                width=900,
                height=600
            )

            if len(df.columns) <= 6 and all(map(lambda col: len(col) <= 35, df.columns)):
                heatmap.update_layout(
                    xaxis_nticks=len(df.columns),
                    yaxis_nticks=len(df.index)
                )
            else:
                heatmap.update_layout(
                    xaxis=dict(showticklabels=False),
                    yaxis=dict(showticklabels=False),
                    hoverlabel=dict(
                        font=dict(
                            size=12
                        )
                    )
                )

            if len(df.columns) <= 6 and all(map(lambda col: len(col) <= 35, df.columns)):
                clustergram = bio.Clustergram(
                            data=df,
                            column_labels=list(df.columns.values),
                            row_labels=list(df.index),
                            hidden_labels=['row'],
                            display_ratio=[0.1, 0.75]
                        )
            else:
                clustergram = bio.Clustergram(
                            data=df,
                            column_labels=list(df.columns.values),
                            row_labels=list(df.index),
                            hidden_labels=['row', 'col'],
                            display_ratio=[0.1, 0.75],
                            height=600,
                            width=900
                        )
            
            rdkit_legend = html.Div([
                                html.Span(style={
                                    "display": "inline-block",
                                    "width": "12px",
                                    "height": "12px",
                                    "borderRadius": "50%",
                                    "backgroundColor": "#e3e334",
                                    "marginRight": "8px"
                                }),
                                "Aromatic atoms"
                            ], className="mb-3 d-none", id="fingerprint-legend")
            
            morgan_legend = html.Div([
                html.Div([
                    html.Span(style={
                        "display": "inline-block",
                        "width": "12px",
                        "height": "12px",
                        "borderRadius": "50%",
                        "backgroundColor": "#e3e334",
                        "marginRight": "6px"
                    }),
                    html.Span("Aromatic atoms", style={"marginRight": "20px"}),

                    html.Span(style={
                        "display": "inline-block",
                        "width": "12px",
                        "height": "12px",
                        "borderRadius": "50%",
                        "backgroundColor": "#9c9ce3",
                        "marginRight": "6px"
                    }),
                    html.Span("Central atom in the environment", style={"marginRight": "20px"}),

                    html.Span(style={
                        "display": "inline-block",
                        "width": "12px",
                        "height": "12px",
                        "borderRadius": "50%",
                        "backgroundColor": "#cfcfcf",
                        "marginRight": "6px"
                    }),
                    html.Span("Aliphatic atoms")
                ], style={"display": "flex", "justifyContent": "center", "alignItems": "center"})
            ], className="mb-3 d-none", id="fingerprint-legend")
            
            fingerprint_bits = dbc.Card([
                dbc.CardBody([
                    dbc.Row([
                        dbc.Col([
                            dbc.Label("Select molecule", className="mt-2 ms-2 mb-2"),
                                dbc.Select(
                                    options=[{"label": key, "value": key} for key in fingerprint_dict.keys()],
                                    id="smiles-select"
                                )],
                                    width=6  
                                    ),
                        dbc.Col([
                            dbc.Label("Select fingerprint bit", className="mt-2 ms-2 mb-2"),
                                dbc.Select(
                                    options=None,
                                    value=None,
                                    id="bit-select"
                                )],
                                    width=6
                                    ),
                                ], className="mb-3")]),
                dbc.CardImg(
                    src=None,
                    bottom=True,
                    id="card-img-bit",
                    style={"width": "40%", "height": "auto"} ,
                    className="mx-auto d-block mb-1"
                ),
                html.Div(
                    morgan_legend if fingerprint_type == "Morgan" else rdkit_legend,
                    className="text-center"
                    )
                ])
            
            similarity_map_legend = html.Div([
                html.Div([
                    html.Span(style={
                        "display": "inline-block",
                        "width": "12px",
                        "height": "12px",
                        "borderRadius": "50%",
                        "backgroundColor": "#6d8964",
                        "marginRight": "6px"
                    }),
                    html.Span("Positive difference (removing bits decreases similarity)", style={"marginRight": "20px"}),

                    html.Span(style={
                        "display": "inline-block",
                        "width": "12px",
                        "height": "12px",
                        "borderRadius": "50%",
                        "backgroundColor": "#ad447e",
                        "marginRight": "6px"
                    }),
                    html.Span("Negative difference (removing bits increases similarity)", style={"marginRight": "20px"}),

                ], style={"display": "flex", "justifyContent": "center", "alignItems": "center"})
            ], className="mb-3 d-none", id="similarity-map-legend")
            
            similarity_map = dbc.Card([
                dbc.CardBody([
                    dbc.Row([
                        dbc.Col([
                            dbc.Label("Select reference molecule", className="mt-2 ms-2 mb-2"),
                                dbc.Select(
                                    options=[{"label": molecule, "value": molecule} for molecule in df.columns],
                                    id="molecule-select-1"
                                )],
                                    width=5  
                                    ),
                        dbc.Col([
                            dbc.Label("Select probe molecule", className="mt-2 ms-2 mb-2"),
                                dbc.Select(
                                    options=[{"label": molecule, "value": molecule} for molecule in df.columns],
                                    id="molecule-select-2"
                                )],
                                    width=5
                                    ),
                        dbc.Col([
                            dbc.Label("Select something!", style={"visibility": "hidden"}),
                            dbc.Button('Generate', id='generate-button', color='primary', className="ms-2 mt-1", n_clicks=0)
                        ])], className="mb-3")
                    ]),
                dbc.CardImg(
                    src=None,
                    bottom=True,
                    id="card-img-similarity",
                    style={"width": "65%", "height": "auto"} ,
                    className="mx-auto d-block mb-1"
                ),
                html.Div(
                    similarity_map_legend,
                    className="text-center"
                )
            ])

            

            table = dash_table.DataTable(
                id='similarity-table',
                columns=[{"name": col, "id": col, "type": "text"} for col in full_symetric_df.columns],
                data=full_symetric_df.to_dict('records'),
                page_size=10,
                sort_action='native',
                filter_action='native',
                cell_selectable=True,
                active_cell=None,
                tooltip_header={
                    "Molecule 1": "Surround SMILES in quotation marks (\"\")",
                    "Molecule 2": "Surround SMILES in quotation marks (\"\")",
                    "Similarity": "Use numbers, >, =>, <, <=, = etc. e.g. >0.4"
                }, 
                tooltip_delay=0,
                tooltip_duration=None,
                style_table={
                    'overflowX': 'auto',
                    'maxWidth': '100%'
                },
                style_cell={
                    'whiteSpace': 'normal',  
                    'height': 'auto',            
                    'textAlign': 'left',
                    'minWidth': '50px',
                    'width': '200px',
                    'maxWidth': '150px',
                    'overflow': 'hidden',
                    'textOverflow': 'ellipsis'
                },
                style_cell_conditional=[
                    {
                        'if': {'column_id': 'Similarity'},
                        'width': '40px',
                        'minWidth': '20px',
                        'maxWidth': '40px',
                        'textAlign': 'center'
                    }
                ],
                style_data={
                    'whiteSpace': 'normal',
                    'height': 'auto',
                }
            )

            molecule_image = dbc.Card([
                dbc.CardBody([
                    dbc.Row([
                        dbc.Col([
                            dbc.Label("Select molecule", className="mt-2 ms-2 mb-2"),
                                dbc.Select(
                                    options=[{"label": key, "value": key} for key in fingerprint_dict.keys()],
                                    id="molecule-select"
                                )])
                            ], className="mb-3")]),
                dbc.CardImg(
                    src=None,
                    bottom=True,
                    id="card-img-molecule",
                    style={"width": "40%", "height": "auto"} ,
                    className="mx-auto d-block mb-1"
                )])
            
            show_bits_tab = fingerprint_type in ["RDKit", "Morgan"]
            show_similarity_map_tab = fingerprint_type in ["Morgan", "AtomPairs"] and similarity_coefficient != "Tversky"
            tabs = dbc.Tabs(
                [
                    dbc.Tab(
                        dbc.Card(
                            dcc.Graph(id='heatmap', figure=heatmap, style={'width': '98%', 'height': 'auto'}),
                            className="mb-3 mr-3"
                            ),
                        label='Heatmap'),
                    dbc.Tab(
                        dbc.Card(
                            dcc.Graph(id='dendrogram', figure=dendrogram, style={'width': '98%', 'height': 'auto'}),
                            className="mb-3 mr-3"
                        ),
                        label="Dendrogram"),
                    dbc.Tab(
                        dbc.Card(
                            dcc.Graph(id="clustergram", figure=clustergram, style={'width': '98%', 'height': 'auto'}),
                            className="mb-3 mr-3"
                            ),
                        label="Clustergram"),
                    dbc.Tab(
                        dbc.Card(
                            dbc.CardBody(table),
                        ),
                        label="Table",
                        className="mb-3 mr-3"),
                    dbc.Tab(
                        molecule_image,
                            label="Molecule Image", 
                            id="molecule-image-tab",
                            tab_id="molecule-image",
                            className="mb-3 mr-3"),
                    dbc.Tab(
                        fingerprint_bits,
                            label="Fingerprint Bits", 
                            id="fingerprint-bits-tab",
                            tab_style={"display": "none"} if not show_bits_tab else {},
                            tab_id="fingerprint-bits",
                            className="mb-3 mr-3"),
                    dbc.Tab(
                        similarity_map,
                            label="Similarity Map",
                            id="similarity-map-tab",
                            tab_style={"display": "none"} if not show_similarity_map_tab else {},
                            tab_id="similarity-map-tab-id",
                            className="mb-3 mr-3")
                ]
            )

            generation_alert = dbc.Alert("NOTE: Changing parameters on the dashboard doesn't change the visualisations until you press the SUBMIT button.", className="mb-3", color="info")

            return None, tabs, fingerprint_dict, list(df.columns), generation_alert


@callback(
    Output('card-img-molecule', 'src'),
    Input('molecule-select', 'value'),
    State('fingerprint-type', 'value'),
    State('similarity-coefficient', 'value'),
    State('textarea-input', 'value'),
    State('min-path-slider', 'value'),
    State('max-path-slider', 'value'),
    State('fps-slider-rdkit', 'value'),
    State('fps-slider-atompairs', 'value'),
    State('radius-slider', 'value'),
    State('fps-slider-morgan', 'value'),
    State('a-input', 'value'),
    State('b-input', 'value'),
    prevent_intial_call=True
)
def get_molecule_image(smiles:str, fingerprint_type: str, similarity_coefficient: str, text_value: str, min_path:int=None, 
                max_path:int=None, fps_rdkit:int = None, fps_atompairs:int = None, radius:int = None, fps_morgan:int = None, weight_a:str=None, weight_b:str=None):
    
    data = {
        'min_path': min_path,
        'max_path': max_path,
        'fps_rdkit': fps_rdkit,
        'fps_atompairs': fps_atompairs,
        'radius': radius, 
        'fps_morgan': fps_morgan,
        'a': weight_a,
        'b': weight_b
    }

    data_frame_generator = DataFrameGenerator(
        smiles=check_textarea_input(text_value),
        generation_strategy=fingerprint_type,
        similarity_strategy=similarity_coefficient,
        data=data)
            
    return data_frame_generator.get_molecule_image(smiles)

@callback(
    Output('bit-select', 'options'),
    Output('bit-select', 'value'),
    Output('card-img-bit', 'src', allow_duplicate=True),
    Output('fingerprint-legend', 'className', allow_duplicate=True),
    Input('smiles-select', 'value'),
    State('fingerprint-store-data', 'data'),
    prevent_initial_call='initial_duplicate'
)
def get_fingerprint_bit_select(value:str, data:dict):
    fingerprint_indices = data[value]
    return [{"label": bit, "value": bit} for bit in fingerprint_indices], None, None, "d-none"

@callback(
    Output('card-img-bit', 'src'),
    Output('fingerprint-legend', 'className'),
    Input('bit-select', 'value'),
    State('smiles-select', 'value'),
    State('fingerprint-type', 'value'),
    State('similarity-coefficient', 'value'),
    State('textarea-input', 'value'),
    State('min-path-slider', 'value'),
    State('max-path-slider', 'value'),
    State('fps-slider-rdkit', 'value'),
    State('fps-slider-atompairs', 'value'),
    State('radius-slider', 'value'),
    State('fps-slider-morgan', 'value'),
    State('a-input', 'value'),
    State('b-input', 'value'),
    prevent_intial_call=True
)
def get_fingerprint_image(bit_value: int, smiles:str, fingerprint_type: str, similarity_coefficient: str, text_value: str, min_path:int=None, 
                max_path:int=None, fps_rdkit:int = None, fps_atompairs:int = None, radius:int = None, fps_morgan:int = None, weight_a:str=None, weight_b:str=None):
    
    data = {
        'min_path': min_path,
        'max_path': max_path,
        'fps_rdkit': fps_rdkit,
        'fps_atompairs': fps_atompairs,
        'radius': radius, 
        'fps_morgan': fps_morgan,
        'a': weight_a,
        'b': weight_b
    }

    data_frame_generator = DataFrameGenerator(
        smiles=check_textarea_input(text_value),
        generation_strategy=fingerprint_type,
        similarity_strategy=similarity_coefficient,
        data=data)
            
    return data_frame_generator.get_fingerprint_bit_image(smiles, bit_value), "ms-3 mb-3 mr-3"
            
@callback(
    Output("molecule-select-1", "options"),
    Output("molecule-select-2", "options"),
    Input("molecule-select-1", "value"),
    Input("molecule-select-2", "value"),
    State('data-frame-data', 'data')  
)
def update_molecule_select_options(selected1, selected2, molecules):
    options1 = [{"label": mol, "value": mol, "disabled": (mol == selected2)} for mol in molecules]
    options2 = [{"label": mol, "value": mol, "disabled": (mol == selected1)} for mol in molecules]

    return options1, options2

@callback(
    Output('card-img-similarity', 'src'),
    Output('similarity-map-legend', 'className'),
    Input('generate-button', 'n_clicks'),
    State('molecule-select-1', 'value'),
    State('molecule-select-2', 'value'),
    State('fingerprint-type', 'value'),
    State('similarity-coefficient', 'value'),
    State('textarea-input', 'value'),
    State('min-path-slider', 'value'),
    State('max-path-slider', 'value'),
    State('fps-slider-rdkit', 'value'),
    State('fps-slider-atompairs', 'value'),
    State('radius-slider', 'value'),
    State('fps-slider-morgan', 'value'),
    State('a-input', 'value'),
    State('b-input', 'value'),
    prevent_intial_call=True
)
def get_similarity_map_image(n_clicks: int, smiles1:str, smiles2:str, fingerprint_type: str, similarity_coefficient: str, text_value: str, min_path:int=None, 
                max_path:int=None, fps_rdkit:int = None, fps_atompairs:int = None, radius:int = None, fps_morgan:int = None, weight_a:str=None, weight_b:str=None):
    
    data = {
        'min_path': min_path,
        'max_path': max_path,
        'fps_rdkit': fps_rdkit,
        'fps_atompairs': fps_atompairs,
        'radius': radius, 
        'fps_morgan': fps_morgan,
        'a': weight_a,
        'b': weight_b
    }

    data_frame_generator = DataFrameGenerator(
        smiles=check_textarea_input(text_value),
        generation_strategy=fingerprint_type,
        similarity_strategy=similarity_coefficient,
        data=data)
            
    return data_frame_generator.get_similarity_map(smiles1, smiles2), "mb-3"
    