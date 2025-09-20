# ===============================
# File: _1_components.py
# ===============================
from dash import html
import dash_bootstrap_components as dbc

# Nav items (path, label)
NAV_LINKS = (
    ("/", "Welcome"),
    ("/upload", "Upload"),
    ("/pre", "Pre-assessment"),
    ("/correction", "Correction"),
    ("/post", "Post-assessment"),
    ("/ranking", "Ranking"),
)

# Stable IDs for navbar buttons
NAV_ID_MAP = {
    "/": "nav-welcome-btn",
    "/upload": "nav-upload-btn",
    "/pre": "nav-pre-btn",
    "/correction": "nav-correction-btn",
    "/post": "nav-post-btn",
    "/ranking": "nav-ranking-btn",
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
    items = []
    last_path = NAV_LINKS[-1][0]

    for path, label in NAV_LINKS:
        is_active = (path == active_path)
        # Intercept Upload clicks via callback (no href) to confirm restart
        href = None if path == "/upload" else path
        items.append(
            dbc.Button(
                label,
                id=NAV_ID_MAP[path],
                href=href,
                color="light" if is_active else "secondary",
                outline=not is_active,
                className="me-2 mb-2",
                size="sm",
                disabled=False,
            )
        )
        if path != last_path:
            items.append(html.Span("→", className="mx-1 text-white-50 fw-bold d-none d-md-inline"))

    items.append(
        dbc.Button(
            "Download results",
            id="download-results-btn",
            color="warning",
            className="ms-3 mb-2",
            size="sm",
        )
    )

    return dbc.Navbar(
        dbc.Container(items, fluid=True, className="py-2"),
        color="dark",  # ensures high contrast
        dark=True,
        className="mb-4",
    )
