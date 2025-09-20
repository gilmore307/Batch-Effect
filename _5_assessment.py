# ===============================
# File: _5_assessment.py
# ===============================
from dash import html, dcc
from dash.dependencies import Input, Output, State
import dash
import dash_bootstrap_components as dbc

from _1_components import build_navbar
from _2_utils import (
    get_session_dir,
    run_r_scripts,
    PRE_SCRIPTS,
    POST_SCRIPTS,
    PRE_FIGURES,
    POST_FIGURES,
    render_assessment_tabs,
)


def assessment_layout(active_path: str, stage: str):
    figure_specs = PRE_FIGURES if stage == "pre" else POST_FIGURES
    header = "Pre-correction assessment" if stage == "pre" else "Post-correction assessment"
    button_id = "run-pre-assessment" if stage == "pre" else "run-post-assessment"
    status_id = "pre-assessment-log" if stage == "pre" else "post-assessment-log"
    gallery_id = "pre-assessment-gallery" if stage == "pre" else "post-assessment-gallery"

    return html.Div(
        [
            build_navbar(active_path),
            dbc.Container(
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                html.H2(header),
                                dbc.Button(
                                    "Run assessment",
                                    id=button_id,
                                    color="primary",
                                    className="mb-3",
                                    disabled=True,
                                ),
                                html.Div(
                                    "Run to populate the tabs on the right.",
                                    className="text-muted",
                                ),
                                html.Div(style={"height": "70vh"}),
                            ],
                            xs=12,
                            md=8,
                            lg=8,
                        ),
                        dbc.Col(
                            [
                                html.H6("Assessment Results", className="mb-2"),
                                html.Div(
                                    id=gallery_id,
                                    children=html.Div("Results will appear here once generated."),
                                    style={"position": "sticky", "top": "1rem", "maxHeight": "80vh", "overflowY": "auto"},
                                ),
                            ],
                            xs=12,
                            md=4,
                            lg=4,
                        ),
                    ],
                    align="start",
                ),
                fluid=True,
            ),
        ]
    )


def register_pre_post_callbacks(app):
    @app.callback(
        Output("pre-assessment-gallery", "children"),
        Output("pre-complete", "data"),
        Input("run-pre-assessment", "n_clicks"),
        State("session-id", "data"),
        prevent_initial_call=True,
    )
    def perform_pre_assessment(n_clicks: int, session_id: str):
        if not n_clicks:
            raise dash.exceptions.PreventUpdate
        session_dir = get_session_dir(session_id)
        if not (session_dir / "raw.csv").exists() or not (session_dir / "metadata.csv").exists():
            message = html.Div("Upload both raw.csv and metadata.csv before running the assessment.")
            return message, False

        success, scripts_log = run_r_scripts(PRE_SCRIPTS, session_dir)
        tabs = render_assessment_tabs(session_dir, PRE_FIGURES)
        return tabs, bool(success)

    @app.callback(
        Output("post-assessment-gallery", "children"),
        Output("post-complete", "data"),
        Input("run-post-assessment", "n_clicks"),
        State("session-id", "data"),
        prevent_initial_call=True,
    )
    def perform_post_assessment(n_clicks: int, session_id: str):
        if not n_clicks:
            raise dash.exceptions.PreventUpdate
        session_dir = get_session_dir(session_id)
        if not (session_dir / "raw.csv").exists() or not (session_dir / "metadata.csv").exists():
            message = html.Div("Upload both raw.csv and metadata.csv before running the assessment.")
            return message, False

        success, scripts_log = run_r_scripts(POST_SCRIPTS, session_dir)
        tabs = render_assessment_tabs(session_dir, POST_FIGURES)
        return tabs, bool(success)
