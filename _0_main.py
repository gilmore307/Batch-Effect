# ===============================
# File: _0_main.py
# ===============================
from pathlib import Path
import uuid
import shutil

import dash
from dash import Dash, dcc, html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc

# Local modules
from _1_components import build_navbar, NAV_LINKS, NAV_ID_MAP
from _2_utils import cleanup_old_sessions, get_session_dir
from _3_welcome import welcome_layout
from _4_upload import upload_layout, register_upload_callbacks
from _5_assessment import assessment_layout, register_pre_post_callbacks
from _6_correction import correction_layout, register_correction_callbacks
from _7_ranking import ranking_layout, register_ranking_callbacks


app: Dash = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.FLATLY],
    suppress_callback_exceptions=True,
)
server = app.server


# ---- Layout factory ----
PAGE_FACTORY = {
    "/": welcome_layout,
    "/upload": upload_layout,
    "/pre": lambda active: assessment_layout(active, stage="pre"),
    "/correction": correction_layout,
    "/post": lambda active: assessment_layout(active, stage="post"),
    "/ranking": ranking_layout,
}


def serve_layout() -> html.Div:
    session_id = str(uuid.uuid4())
    cleanup_old_sessions()
    return html.Div(
        [
            dcc.Location(id="page-url", refresh=False),

            # Session + state stores
            dcc.Store(id="session-id", storage_type="session", data=session_id),
            dcc.Store(id="selected-methods", storage_type="session", data=[]),
            dcc.Store(id="upload-complete", storage_type="session", data=False),
            dcc.Store(id="pre-complete", storage_type="session", data=False),
            dcc.Store(id="correction-complete", storage_type="session", data=False),
            dcc.Store(id="post-complete", storage_type="session", data=False),

            # Download
            dcc.Download(id="download-results"),

            # Page mount
            html.Div(id="page-content"),
        ]
    )


app.layout = serve_layout


# ---- Page routing ----
@app.callback(Output("page-content", "children"), Input("page-url", "pathname"))
def render_page(pathname: str):
    layout_factory = PAGE_FACTORY.get(pathname, welcome_layout)
    return layout_factory(pathname)


# ---- Enable/disable navigation buttons based on stage completion ----
@app.callback(
    Output(NAV_ID_MAP["/"], "disabled"),          # Welcome (never disabled)
    Output(NAV_ID_MAP["/upload"], "disabled"),    # Upload
    Output(NAV_ID_MAP["/pre"], "disabled"),       # Pre
    Output(NAV_ID_MAP["/correction"], "disabled"),# Correction
    Output(NAV_ID_MAP["/post"], "disabled"),      # Post
    Output(NAV_ID_MAP["/ranking"], "disabled"),   # Ranking
    Input("upload-complete", "data"),
    Input("pre-complete", "data"),
    Input("correction-complete", "data"),
    Input("post-complete", "data"),
)
def gate_nav_buttons(upload_done, pre_done, correction_done, post_done):
    upload_done = bool(upload_done)
    pre_done = bool(pre_done)
    correction_done = bool(correction_done)
    post_done = bool(post_done)
    return (
        False,                 # Welcome
        False,                 # Upload always enabled
        not upload_done,       # Pre requires upload
        not pre_done,          # Correction requires pre
        not correction_done,   # Post requires correction
        not post_done,         # Ranking requires post
    )


# ---- Gate action buttons by upload state ----
@app.callback(
    Output("run-pre-assessment", "disabled"),
    Input("upload-complete", "data"),
)
def toggle_pre_assessment_button(upload_complete: bool) -> bool:
    return not bool(upload_complete)


@app.callback(
    Output("run-correction", "disabled"),
    Input("upload-complete", "data"),
)
def toggle_correction_button(upload_complete: bool) -> bool:
    return not bool(upload_complete)


# ---- Download results (ZIP session dir) ----
@app.callback(
    Output("download-results", "data"),
    Input("download-results-btn", "n_clicks"),
    State("session-id", "data"),
    prevent_initial_call=True,
)
def download_results(n_clicks: int, session_id: str):
    if not n_clicks:
        raise dash.exceptions.PreventUpdate
    session_dir = get_session_dir(session_id)
    if not session_dir.exists():
        raise dash.exceptions.PreventUpdate
    archive_path = session_dir / "results"
    zip_path = shutil.make_archive(str(archive_path), "zip", session_dir)
    return dcc.send_file(zip_path)


# ---- Register page-specific callbacks ----
register_upload_callbacks(app)
register_pre_post_callbacks(app)
register_correction_callbacks(app)
register_ranking_callbacks(app)


if __name__ == "__main__":
    app.run_server(debug=True, host="127.0.0.1", port=8051)

