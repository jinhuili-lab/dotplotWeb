import base64  # Add this import statement
from dash import Dash, dcc, html, Input, Output, ctx, callback
import dash
from dash.dependencies import Input, Output, State
import dash_daq as daq
import plotly.express as px
import pandas as pd
from textwrap import wrap
from flask_caching import Cache
import io
# import dash_bio as dashbio
import urllib.request as urlreq
# from dash.dependencies import Input, Output
from dash import dash_table

# Initialize Dash app and Flask server
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
server = app.server

# Initialize Flask-Caching
cache = Cache(server, config={'CACHE_TYPE': 'simple'})

# Function to visualize blastn results
@cache.memoize()
def visualize_blastn_results(file_contents, color,line_weight,marker_size,getDirection):
    file_like_object = io.BytesIO(file_contents)
    df = pd.read_csv(file_like_object, comment='#', sep="\t", header=None)
    file_like_object.seek(0)
    title = [line.strip() for  line in file_like_object.readlines() if line.decode().startswith("#")]
    print(type(title))#.strip("# Query:")
    title = "<br>".join(wrap(title[1].decode().strip("# Query: "), width=60))
    df['line_num'] = df.reset_index().index + 1
    df.columns = ["F" + str(i) for i in range(16)]
    if getDirection == "All":
        df = df
    elif getDirection == "Revise":
        df = getRevise(df)
    else:
        df = getForward(df)
    dic = {}
    dic2 = {}
    num = 0
    for index, row in df.iterrows():
        x1 = row["F6"]
        x2 = row["F7"]
        y1 = row["F8"]
        y2 = row["F9"]
        seq_q = row["F12"]
        line_number = row[15]

        # Store data in dictionary
        dic[num] = [x1, y1, seq_q, line_number]
        # dic2[line_number] = [str(x1) + "_" + str(y1) + "_" + str(x2) + "_" + str(y2), seq_q, line_number]

        num += 1
        dic[num] = [x2, y2, seq_q, line_number]
        num += 1

    df1 = pd.DataFrame(dic).T
    df1.columns = ["Pos" + str(i) for i in range(2)] + ["seq_q", "No"]
    df1 = df1.sort_values(by="No")

    # df2 = pd.DataFrame(dic2).T
    # df2.columns = ["F" + str(i) for i in range(3)] #["label","seq","No."]
    # df2 = df2.sort_values(by="F0")
    # fasta_list=""
    # for i in range(len(df2)):
    #     fasta_list+=(">" + df2["F0"][i+1] + '\n' + df2["F1"][i+1]+ '\n' )
    # fasta_list = urlreq.urlopen(
    #     'https://git.io/alignment_viewer_p53.fasta'
    # ).read().decode('utf-8')
    # align = dashbio.AlignmentChart(
    #         id='alignment-viewer',
    #         data=fasta_list
    #     )
    # align = "1"
    fig = px.line(df1, x="Pos0", y="Pos1", title=title, #line_dash_sequence =["1000px"],
                  color="seq_q", width=800, height=800, markers=True, symbol="No",hover_data={"seq_q": True})
    fig.update_layout(
        showlegend=False,
        # size = 10,
        plot_bgcolor="#fff",
        uniformtext_minsize=1,
        uniformtext_mode='hide',
        xaxis=dict(zeroline=False),
        yaxis=dict(zeroline=False),
        paper_bgcolor="white",
        title_xanchor="center",
        title_xref="paper",
        title_x=0.5,
        title_yanchor='top'
    )
    fig.update_traces(
        # ModuleNotFoundError
        line=dict(width=line_weight),
        marker_size = marker_size, 
        marker_line = dict(width=1, color=color['hex']), 
        # selector = dict(mode = "markers")
        )
    # fig.update_geos(framecolor="black", framewidth=5)
    fig.update_xaxes(showline=True, linewidth=3, linecolor='black', title=None)
    fig.update_yaxes(showline=True, linewidth=3, linecolor='black', title=None)
    return fig #,align
styles = {
    'output': {
        'overflow-y': 'scroll',
        'overflow-wrap': 'break-word',
        'height': 'calc(100% - 25px)',
        'border': 'thin lightgrey solid'
    },
    'tab': {'height': 'calc(98vh - 115px)'}
}
def getForward(df):
    result = df[(df.iloc[:, 6] <= df.iloc[:, 7]) & (df.iloc[:, 8] <= df.iloc[:, 9])]
    return result
def getRevise(df):
    result = df[(df.iloc[:, 6] <= df.iloc[:, 7]) & (df.iloc[:, 8] > df.iloc[:, 9])]
    return result
# Define Dash layout
app.layout = html.Div([
    html.H1(["Blastn2Dotplot v1.0"]),
    html.Hr(),
    html.Div(className="three columns", style={'color':'blue'}, children= [
        html.Div([
            daq.ColorPicker(id = 'color', label = "shadow color ",value = dict(hex="#000000")),
            html.H4("line size"),
            dcc.Slider( 
                0, 
                20, 
                step=0.1, 
                value=6, 
                id='line_size' 
                ),
            html.H4("marker size"),
            dcc.Slider( 
                0, 
                20, 
                step=0.1, 
                value=13, 
                id='marker_size' 
                ),
            html.H4("Sequence Direction"),
            dcc.Dropdown(['All', 'Forward', 'Revise'], 'All', id='getDirection'),
            dcc.Upload(
                id='upload-data',
                children=html.Div([
                    'Drag and Drop or ',
                    html.A('Select Files')
                ]),
                style={
                    'width': '90%',
                    'height': '60px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '10px',
                },
                # Allow multiple files to be uploaded
                multiple=False
            ),
            
        ]),
    ]),
    html.Div(className="eight columns", children= [
        # html.Hr(),
        html.Div(id='output-data-upload'),
        html.Hr(),
        html.Div(id='output-container'),
    ]),
    
])


@app.callback(
    Output('output-container', 'children'),
    [Input('blastn-plot', 'hoverData')])
def display_hover_data(hoverData):
    if hoverData is not None:
        # point_index = hoverData['points'][0]['pointIndex']
        # # 获取该点的文本信息
        # text_info = hoverData['points'][0]['text']
        # return f"You hovered over point {point_index + 1}, text: {text_info}"
        x = hoverData['points'][0]['x']
        y = hoverData['points'][0]['y']
        seq = hoverData['points'][0]['customdata']
        dic= [{ "x":x,   "y":y,   "seq": seq,      }]

        # 获取该点的文本信息
        # text_info = hoverData['points'][0]['text']
        # html_table
        # dtable = dash_table.DataTable(dic, [{"name": i, "id": i} for i in dic.keys()])
        return  f"You hovered over: {dic}"#x: {x}, y: {y},{hoverData}
        # return f"You hovered over: {hoverData}"
    else:
        return "Hover over the graph to see data"


# Callback to update output
@app.callback(Output('output-data-upload', 'children'),
              Input('upload-data', 'contents'),
              Input('color', 'value'),
              Input('line_size', 'value'),
              Input('marker_size', 'value'),
              Input('getDirection', 'value'),
              prevent_initial_call=True,
              )
def update_output(contents, color,line_size, marker_size,getDirection):#
    if contents is not None:
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        figs = visualize_blastn_results(decoded, color,line_size,marker_size,getDirection)
        pl = dcc.Graph(
            id='blastn-plot',
            figure=figs#[0],
        )
        return pl
if __name__ == '__main__':
    app.run_server(debug=True)