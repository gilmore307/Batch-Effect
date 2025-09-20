# ===============================
# File: _3_welcome.py
# ===============================
from dash import html
import dash_bootstrap_components as dbc
from _1_components import build_navbar


def welcome_layout(active_path: str):
    return html.Div(
        [
            build_navbar(active_path),
            dbc.Container(
                [
                    html.H2("Welcome"),
                    html.P("Placeholder for introduction content."),
                ],
                fluid=True,
            ),
        ]
    )
