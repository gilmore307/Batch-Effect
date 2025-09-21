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
            dcc.Store(id="preprocess-complete", storage_type="session", data=False),
            dcc.Store(id="pre-complete", storage_type="session", data=False),
            dcc.Store(id="pre-started", storage_type="session", data=False),
            dcc.Store(id="correction-complete", storage_type="session", data=False),
            dcc.Store(id="post-complete", storage_type="session", data=False),
            dcc.Store(id="help-shown", storage_type="session", data=False),
            dcc.Store(id="runlog-path", storage_type="session", data=""),

            # Download
            dcc.Download(id="download-results"),

            # Page mount (preload Home so navbar buttons exist for callbacks)
            html.Div(id="page-content", children=welcome_layout("/")),

            # Background log poller
            dcc.Interval(id="runlog-interval", interval=800, n_intervals=0, disabled=True),

            # Confirm restart modal (shown when clicking Upload from other pages)
            dbc.Modal(
                [
                    dbc.ModalHeader(dbc.ModalTitle("Restart session?")),
                    dbc.ModalBody(
                        "Starting a new session clears your current uploads and results. Continue?"
                    ),
                    dbc.ModalFooter(
                        [
                            dbc.Button("Cancel", id="confirm-restart-no", color="secondary", className="me-2"),
                            dbc.Button("Restart", id="confirm-restart-yes", color="danger"),
                        ]
                    ),
                ],
                id="confirm-restart-modal",
                is_open=False,
                backdrop="static",
                keyboard=False,
                centered=True,
            ),

            # Help modal (placeholder)
            dbc.Modal(
                [
                    dbc.ModalHeader(dbc.ModalTitle("Help")),
                    dbc.ModalBody(
                        html.Div([
                            html.P("Help content coming soon."),
                            html.P("Resources:"),
                            html.Ul([
                                html.Li(html.A("GitHub repository", href="https://github.com/gilmore307/Batch-Effect", target="_blank")),
                                html.Li(html.A("Report an issue", href="https://github.com/gilmore307/Batch-Effect/issues", target="_blank")),
                            ])
                        ])
                    ),
                    dbc.ModalFooter(
                        dbc.Button("Close", id="help-close", color="secondary")
                    ),
                ],
                id="help-modal",
                is_open=False,
                centered=True,
            ),

            # Run log modal
            dbc.Modal(
                [
                    dbc.ModalHeader(dbc.ModalTitle("Run Log")),
                    dbc.ModalBody(
                        html.Pre(
                            id="runlog-content",
                            children="Log will appear here...",
                            style={
                                "whiteSpace": "pre-wrap",
                                "fontFamily": "monospace",
                                "maxHeight": "70vh",
                                "overflowY": "auto",
                                "marginBottom": 0,
                            },
                        )
                    ),
                    dbc.ModalFooter(
                        dbc.Button("Close", id="runlog-close", color="secondary")
                    ),
                ],
                id="runlog-modal",
                is_open=False,
                centered=True,
                size="xl",
                scrollable=True,
            ),
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
    Output(NAV_ID_MAP["/"], "disabled"),          # Home (never disabled)
    Output(NAV_ID_MAP["/upload"], "disabled"),    # Upload
    Output(NAV_ID_MAP["/pre"], "disabled"),       # Pre
    Output(NAV_ID_MAP["/correction"], "disabled"),# Correction
    Output(NAV_ID_MAP["/post"], "disabled"),      # Post
    Input("preprocess-complete", "data"),
    Input("pre-started", "data"),
    Input("pre-complete", "data"),
    Input("correction-complete", "data"),
    Input("post-complete", "data"),
)
def gate_nav_buttons(preprocess_done, pre_started, pre_done, correction_done, post_done):
    preprocess_done = bool(preprocess_done)
    pre_started = bool(pre_started)
    pre_done = bool(pre_done)
    correction_done = bool(correction_done)
    post_done = bool(post_done)
    return (
        False,                 # Welcome
        False,                 # Upload always enabled
        not preprocess_done,   # Pre requires preprocess
        not pre_started,       # Correction enabled after pre is started
        not correction_done,   # Post requires correction
    )


@app.callback(
    Output("run-correction", "disabled"),
    Input("upload-complete", "data"),
)
def toggle_correction_button(upload_complete: bool) -> bool:
    return not bool(upload_complete)


# ---- Intercept Upload nav to confirm restart ----
@app.callback(
    Output("confirm-restart-modal", "is_open"),
    Input(NAV_ID_MAP["/upload"], "n_clicks"),
    State("page-url", "pathname"),
    prevent_initial_call=True,
)
def open_restart_modal(n_clicks: int, pathname: str) -> bool:
    if not n_clicks:
        raise dash.exceptions.PreventUpdate
    # Only ask when not already on the upload page
    # and not when navigating from Home ("/") to Upload
    if pathname in ("/upload", "/"):
        raise dash.exceptions.PreventUpdate
    return True


@app.callback(
    Output("session-id", "data"),
    Output("selected-methods", "data", allow_duplicate=True),
    Output("upload-complete", "data", allow_duplicate=True),
    Output("preprocess-complete", "data", allow_duplicate=True),
    Output("pre-started", "data", allow_duplicate=True),
    Output("pre-complete", "data", allow_duplicate=True),
    Output("correction-complete", "data", allow_duplicate=True),
    Output("post-complete", "data", allow_duplicate=True),
    Output("page-url", "pathname"),
    Output("confirm-restart-modal", "is_open", allow_duplicate=True),
    Input("confirm-restart-yes", "n_clicks"),
    Input("confirm-restart-no", "n_clicks"),
    State("session-id", "data"),
    prevent_initial_call=True,
)
def handle_restart(confirm_yes: int, confirm_no: int, current_session: str):
    ctx = dash.callback_context
    if not ctx.triggered:
        raise dash.exceptions.PreventUpdate
    trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]
    if trigger_id == "confirm-restart-yes":
        # Delete previous session folder if present
        try:
            if current_session:
                old_dir = get_session_dir(current_session)
                if old_dir.exists():
                    shutil.rmtree(old_dir, ignore_errors=True)
        except Exception:
            # Ignore cleanup errors; proceed with new session
            pass
        new_session = str(uuid.uuid4())
        return (
            new_session,  # session-id
            [],           # selected-methods
            False,        # upload-complete
            False,        # preprocess-complete
            False,        # pre-started
            False,        # pre-complete
            False,        # correction-complete
            False,        # post-complete
            "/upload",   # navigate to upload
            False,        # close modal
        )
    else:
        # Cancel: keep everything, just close modal
        return (
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,  # preprocess-complete
            dash.no_update,  # pre-started
            dash.no_update,
            dash.no_update,
            dash.no_update,
            dash.no_update,
            False,
        )

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

# Enable Download button only after post-assessment finished
@app.callback(
    Output("download-results-btn", "disabled"),
    Input("post-complete", "data"),
)
def enable_download(post_done: bool) -> bool:
    return not bool(post_done)


@app.callback(
    Output("help-modal", "is_open"),
    Output("help-shown", "data"),
    Input("help-open", "n_clicks"),
    Input("help-close", "n_clicks"),
    State("help-shown", "data"),
)
def toggle_help_modal(open_clicks, close_clicks, help_shown):
    # Open once at session start (initial call), then only on user action.
    ctx = dash.callback_context
    if not ctx.triggered:
        if not help_shown:
            return True, True
        return dash.no_update, dash.no_update
    trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]
    if trigger_id == "help-open":
        return True, True
    if trigger_id == "help-close":
        return False, dash.no_update
    return dash.no_update, dash.no_update


# Highlight next-step nav button when user can proceed
@app.callback(
    Output(NAV_ID_MAP["/"], "color"), Output(NAV_ID_MAP["/"], "outline"),
    Output(NAV_ID_MAP["/upload"], "color"), Output(NAV_ID_MAP["/upload"], "outline"),
    Output(NAV_ID_MAP["/pre"], "color"), Output(NAV_ID_MAP["/pre"], "outline"),
    Output(NAV_ID_MAP["/correction"], "color"), Output(NAV_ID_MAP["/correction"], "outline"),
    Output(NAV_ID_MAP["/post"], "color"), Output(NAV_ID_MAP["/post"], "outline"),
    Input("page-url", "pathname"),
    Input("upload-complete", "data"),
    Input("pre-started", "data"),
    Input("pre-complete", "data"),
    Input("correction-complete", "data"),
    Input("post-complete", "data"),
)
def highlight_next(pathname, upload_done, pre_started, pre_done, correction_done, post_done):
    upload_done = bool(upload_done)
    pre_started = bool(pre_started)
    pre_done = bool(pre_done)
    correction_done = bool(correction_done)
    post_done = bool(post_done)

    # Determine next step user can take
    if not upload_done:
        next_step = "/upload"
    elif not pre_started:
        next_step = "/pre"
    elif not correction_done:
        next_step = "/correction"
    elif not post_done:
        next_step = "/post"
    else:
        next_step = None

    def style_for(path):
        if path == pathname:
            return ("light", False)
        if path == next_step:
            return ("warning", True)  # highlight next step
        return ("light", True)

    home = style_for("/")
    upload = style_for("/upload")
    pre = style_for("/pre")
    corr = style_for("/correction")
    post = style_for("/post")
    return (
        home[0], home[1],
        upload[0], upload[1],
        pre[0], pre[1],
        corr[0], corr[1],
        post[0], post[1],
    )


# ---- Run log modal handlers ----
@app.callback(
    Output("runlog-modal", "is_open", allow_duplicate=True),
    Output("runlog-interval", "disabled", allow_duplicate=True),
    Input("log-open", "n_clicks"),
    Input("runlog-close", "n_clicks"),
    State("runlog-modal", "is_open"),
    prevent_initial_call=True,
)
def toggle_runlog_modal(open_clicks, close_clicks, is_open):
    # Only user-driven open/close here; job callbacks will also open it explicitly.
    ctx = dash.callback_context
    if not ctx.triggered:
        raise dash.exceptions.PreventUpdate
    trig = ctx.triggered[0]["prop_id"].split(".")[0]
    if trig == "log-open":
        return True, False
    if trig == "runlog-close":
        return False, True
    return is_open, dash.no_update


@app.callback(
    Output("runlog-content", "children"),
    Input("runlog-interval", "n_intervals"),
    State("runlog-path", "data"),
)
def update_runlog_content(n, log_path):
    if not log_path:
        raise dash.exceptions.PreventUpdate
    try:
        p = Path(log_path)
        if not p.exists():
            return "Waiting for log file..."
        text = p.read_text(encoding="utf-8", errors="replace")
        # limit to last 10000 chars to avoid huge payloads
        return text[-10000:]
    except Exception:
        return "(Failed to read log)"


# (Removed init_runlog_on_actions to avoid referencing page-specific IDs
#  that are not present in the current layout.)

if __name__ == "__main__":
    app.run_server(debug=True, host="127.0.0.1", port=8051)

