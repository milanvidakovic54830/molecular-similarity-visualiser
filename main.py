from dash import Dash
import dash_bootstrap_components as dbc
from layout import layout

external_stylesheets = [dbc.themes.LUMEN]
app = Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = layout

if __name__ == '__main__':
    app.run()