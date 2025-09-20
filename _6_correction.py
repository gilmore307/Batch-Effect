# ===============================
# File: _6_correction.py
# ===============================
from typing import Sequence
from dash import html
from dash.dependencies import Input, Output, State
import dash
import dash_bootstrap_components as dbc

from _1_components import build_navbar
from _2_utils import get_session_dir, run_methods, SUPPORTED_METHODS, DEFAULT_METHODS


def correction_layout(active_path: str):
    checklist_options = [
        {"label": display, "value": code} for code, display in SUPPORTED_METHODS
    ]
    return html.Div(
        [
            build_navbar(active_path),
            dbc.Container(
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                html.H2("Run batch correction"),
                                html.P("Select the correction methods to run, then start the pipeline."),
                                html.H4("Select correction methods"),
                                dcc.Checklist(
                                    id="method-selection",
                                    options=checklist_options,
                                    value=list(DEFAULT_METHODS),
                                    inputStyle={"margin-right": "0.5rem"},
                                    labelStyle={"display": "block", "margin-bottom": "0.4rem"},
                                ),
                                dbc.Button("Run correction", id="run-correction", color="primary", className="mb-3", disabled=True),
                            ],
                            xs=12,
                            md=8,
                            lg=9,
                        ),
                    ],
                    align="start",
                ),
                fluid=True,
            ),
        ]
    )


def register_correction_callbacks(app):
    @app.callback(Output("selected-methods", "data"), Input("method-selection", "value"), prevent_initial_call=True)
    def sync_method_selection(selected):
        return selected or []

    @app.callback(
        Output("correction-complete", "data"),
        Input("run-correction", "n_clicks"),
        State("selected-methods", "data"),
        State("session-id", "data"),
        prevent_initial_call=True,
    )
    def perform_correction(n_clicks: int, methods: Sequence[str], session_id: str):
        if not n_clicks:
            raise dash.exceptions.PreventUpdate
        if not session_id:
            return False
        session_dir = get_session_dir(session_id)
        if not (session_dir / "raw.csv").exists() or not (session_dir / "metadata.csv").exists():
            return False
        methods = methods or []
        if not methods:
            return False

        success, log = run_methods(session_dir, methods)
        return bool(success)
