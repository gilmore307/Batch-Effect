# ===============================
# File: _1_components.py
# ===============================
from dash import html
import dash_bootstrap_components as dbc

# Nav items (path, label)
NAV_LINKS = (
    ("/", "Home"),
    ("/upload", "Upload Files"),
    ("/pre", "Pre-correction Assessment"),
    ("/correction", "Batch Effect Correction"),
    ("/post", "Post-correction Assessment"),
)

# Stable IDs for navbar buttons
NAV_ID_MAP = {
    "/": "nav-welcome-btn",
    "/upload": "nav-upload-btn",
    "/pre": "nav-pre-btn",
    "/correction": "nav-correction-btn",
    "/post": "nav-post-btn",
}


def logger_box(log_id: str, title: str = "Progress log") -> dbc.Col:
    return dbc.Col(
        dbc.Card(
            [
                dbc.CardHeader(html.Strong(title)),
                dbc.CardBody(
                    dbc.Alert(
                        id=log_id,
                        color="secondary",
                        is_open=False,
                        style={
                            "whiteSpace": "pre-wrap",
                            "fontFamily": "monospace",
                            "maxHeight": "70vh",
                            "overflowY": "auto",
                            "marginBottom": 0,
                        },
                    )
                ),
            ],
            className="h-100",
        ),
        xs=12,
        md=4,
        lg=3,
        className="mb-4",
    )


def build_navbar(active_path: str) -> dbc.Navbar:
    # Row 1: Home + Help + external links + Download
    top_items = []

    # Home button
    is_home = active_path == "/"
    top_items.append(
        dbc.Button(
            "Home",
            id=NAV_ID_MAP["/"],
            href="/",
            color="light",
            outline=not is_home,
            className="me-2 mb-2 be-nav-btn",
            size="sm",
        )
    )

    # Help modal trigger
    top_items.append(
        dbc.Button(
            "Help",
            id="help-open",
            color="info",
            outline=True,
            className="me-2 mb-2",
            size="sm",
        )
    )

    # External links
    top_items.append(
        dbc.Button(
            "GitHub",
            color="secondary",
            outline=True,
            className="me-2 mb-2",
            size="sm",
            href="https://github.com/gilmore307/Batch-Effect",
            target="_blank",
        )
    )
    top_items.append(
        dbc.Button(
            "Report issue",
            color="secondary",
            outline=True,
            className="me-2 mb-2",
            size="sm",
            href="https://github.com/gilmore307/Batch-Effect/issues",
            target="_blank",
        )
    )

    # Spacer then Download button on the right
    top_items.append(html.Div(className="flex-grow-1"))
    top_items.append(
        dbc.Button(
            "Download results",
            id="download-results-btn",
            color="warning",
            className="ms-3 mb-2",
            size="sm",
            disabled=True,
        )
    )

    top_row = html.Div(top_items, className="d-flex align-items-center flex-wrap py-1")

    # Row 2: Pipeline navigation (exclude Home)
    pipeline_paths = [p for p, _ in NAV_LINKS if p != "/"]
    bottom_items = []
    last_path = pipeline_paths[-1]

    for path, label in [x for x in NAV_LINKS if x[0] != "/"]:
        is_active = path == active_path
        # Upload: direct href only when on Home; otherwise intercept for restart modal
        if path == "/upload":
            href = "/upload" if is_home else None
        else:
            href = path
        bottom_items.append(
            dbc.Button(
                label,
                id=NAV_ID_MAP[path],
                href=href,
                color="light",
                outline=not is_active,
                className="me-2 mb-2 be-nav-btn",
                size="sm",
                disabled=False,
            )
        )
        if path != last_path:
            bottom_items.append(html.Span("鈫?, className="mx-2 text-white fw-bold d-none d-md-inline"))

    bottom_row = html.Div(bottom_items, className="d-flex align-items-center flex-wrap py-1")

    return dbc.Navbar(
        dbc.Container([top_row, bottom_row], fluid=True, className="py-1"),
        color="primary",  # brighter background for better contrast
        dark=True,
        className="mb-4 be-navbar",
    )

