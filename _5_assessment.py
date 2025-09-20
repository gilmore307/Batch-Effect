# ===============================
# File: _5_assessment.py
# ===============================
from dash import html
from dash.dependencies import Input, Output, State
import dash
import dash_bootstrap_components as dbc

from _1_components import build_navbar, logger_box
from _2_utils import get_session_dir, run_r_scripts, render_figures, PRE_SCRIPTS, POST_SCRIPTS, PRE_FIGURES, POST_FIGURES


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
                        logger_box(status_id, "Assessment log"),
                        dbc.Col(
                            [
                                html.H2(header),
                                dbc.Button("Run assessment", id=button_id, color="primary", className="mb-3", disabled=True),
                                html.Div(id=gallery_id, children=html.Div("Results will appear here once generated.")),
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


def register_pre_post_callbacks(app):
    @app.callback(
        Output("pre-assessment-log", "children"),
        Output("pre-assessment-log", "is_open"),
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
            message = "Upload both raw.csv and metadata.csv before running the assessment."
            return message, True, html.Div(), False

        success, scripts_log = run_r_scripts(PRE_SCRIPTS, session_dir)
        gallery = render_figures(session_dir, PRE_FIGURES)
        return scripts_log, True, gallery, bool(success)

    @app.callback(
        Output("post-assessment-log", "children"),
        Output("post-assessment-log", "is_open"),
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
            message = "Upload both raw.csv and metadata.csv before running the assessment."
            return message, True, html.Div(), False

        success, scripts_log = run_r_scripts(POST_SCRIPTS, session_dir)
        gallery = render_figures(session_dir, POST_FIGURES)
        return scripts_log, True, gallery, bool(success)