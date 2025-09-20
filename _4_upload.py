# ===============================
# File: _4_upload.py
# ===============================
from typing import List
from dash import dcc, html
from dash.dependencies import Input, Output, State
import dash
import dash_bootstrap_components as dbc

from _1_components import build_navbar
from _2_utils import (
    get_session_dir,
    save_uploaded_file,
    human_size,
    run_preprocess,
)


def upload_layout(active_path: str):
    return html.Div(
        [
            build_navbar(active_path),
            dbc.Container(
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                html.H2("Upload data"),
                                html.P("Provide the raw count matrix and the matching metadata table."),
                                dbc.Row(
                                    [
                                        dbc.Col(
                                            dbc.Card(
                                                [
                                                    dbc.CardHeader(html.Strong("Count matrix (CSV)")),
                                                    dbc.CardBody(
                                                        dcc.Upload(
                                                            id="upload-matrix",
                                                            children=html.Div(["Drag & drop or click to upload count matrix"]),
                                                            multiple=False,
                                                            className="border border-secondary rounded p-5 text-center bg-light",
                                                            accept=".csv,text/csv",
                                                        )
                                                    ),
                                                    dbc.CardFooter(html.Div(id="matrix-file-info", className="text-muted", children="No file uploaded yet.")),
                                                ],
                                                className="h-100",
                                            ),
                                            md=6,
                                            className="mb-3",
                                        ),
                                        dbc.Col(
                                            dbc.Card(
                                                [
                                                    dbc.CardHeader(html.Strong("Metadata (CSV)")),
                                                    dbc.CardBody(
                                                        dcc.Upload(
                                                            id="upload-metadata",
                                                            children=html.Div(["Drag & drop or click to upload metadata"]),
                                                            multiple=False,
                                                            className="border border-secondary rounded p-5 text-center bg-light",
                                                            accept=".csv,text/csv",
                                                        )
                                                    ),
                                                    dbc.CardFooter(html.Div(id="metadata-file-info", className="text-muted", children="No file uploaded yet.")),
                                                ],
                                                className="h-100",
                                            ),
                                            md=6,
                                            className="mb-3",
                                        ),
                                    ]
                                ),
                                html.Hr(),
                                html.H4("Process uploads"),
                                html.P("After both files are uploaded, click Process to review metadata columns and run preprocessing."),
                                dbc.Button("Process", id="process-uploads", color="primary", className="mb-3", disabled=True),
                                html.Div(id="metadata-columns-display", className="mb-2"),
                                html.Div(id="column-mapping-container", className="mb-3"),
                                html.Div(id="process-result", className="mb-2"),
                                html.P("Re-uploading a file will replace the previous one."),
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


def register_upload_callbacks(app):
    @app.callback(
        Output("upload-complete", "data"),
        Output("matrix-file-info", "children"),
        Output("metadata-file-info", "children"),
        Input("upload-matrix", "contents"),
        Input("upload-metadata", "contents"),
        State("upload-matrix", "filename"),
        State("upload-metadata", "filename"),
        State("session-id", "data"),
        prevent_initial_call=True,
    )
    def handle_upload(matrix_contents, metadata_contents, matrix_name, metadata_name, session_id):
        if not session_id:
            return False, "No file uploaded yet.", "No file uploaded yet."

        session_dir = get_session_dir(session_id)
        saved_items: List[str] = []
        matrix_info = "No file uploaded yet."
        metadata_info = "No file uploaded yet."

        if matrix_contents:
            save_uploaded_file(matrix_contents, session_dir, "raw.csv")
            size = (session_dir / "raw.csv").stat().st_size
            matrix_info = dbc.Badge("Uploaded", color="success", className="me-2"), html.Span(f"Source: {matrix_name} → Saved as: raw.csv • {human_size(size)}")
            saved_items.append(f"Matrix saved as raw.csv (source: {matrix_name})")

        if metadata_contents:
            save_uploaded_file(metadata_contents, session_dir, "metadata.csv")
            size = (session_dir / "metadata.csv").stat().st_size
            metadata_info = dbc.Badge("Uploaded", color="success", className="me-2"), html.Span(f"Source: {metadata_name} → Saved as: metadata.csv • {human_size(size)}")
            saved_items.append(f"Metadata saved as metadata.csv (source: {metadata_name})")

        upload_complete = (session_dir / "raw.csv").exists() and (session_dir / "metadata.csv").exists()
        return upload_complete, matrix_info, metadata_info

    @app.callback(
        Output("process-uploads", "disabled"),
        Input("upload-complete", "data"),
    )
    def toggle_process_button(upload_complete: bool):
        return not bool(upload_complete)

    @app.callback(
        Output("metadata-columns-display", "children"),
        Output("column-mapping-container", "children"),
        Input("process-uploads", "n_clicks"),
        State("session-id", "data"),
        prevent_initial_call=True,
    )
    def show_metadata_columns(n_clicks: int, session_id: str):
        if not n_clicks:
            raise dash.exceptions.PreventUpdate
        if not session_id:
            return dash.no_update, dash.no_update
        session_dir = get_session_dir(session_id)
        meta_path = session_dir / "metadata.csv"
        if not meta_path.exists():
            return dash.no_update, dash.no_update

        # Read header from metadata
        try:
            import csv
            with meta_path.open("r", encoding="utf-8") as fh:
                reader = csv.reader(fh)
                header = next(reader)
        except Exception:
            return dash.no_update, dash.no_update

        col_names = [name.strip() for name in header if name and name.strip()]
        if not col_names:
            return dash.no_update, dash.no_update

        chips = html.Div([
            html.Span(name, className="badge bg-secondary me-1 mb-1") for name in col_names
        ])

        opts = [{"label": name, "value": name} for name in col_names]
        mapping_ui = dbc.Card([
            dbc.CardHeader(html.Strong("Map metadata columns")),
            dbc.CardBody([
                dbc.Row([
                    dbc.Col([
                        dbc.Label("Sample ID column"),
                        dcc.Dropdown(id="map-sample-id", options=opts, placeholder="Select sample_id column"),
                    ], md=4),
                    dbc.Col([
                        dbc.Label("Batch ID column"),
                        dcc.Dropdown(id="map-batch-id", options=opts, placeholder="Select batch_id column"),
                    ], md=4),
                    dbc.Col([
                        dbc.Label("Phenotype column"),
                        dcc.Dropdown(id="map-phenotype", options=opts, placeholder="Select phenotype column"),
                    ], md=4),
                ], className="gy-2"),
                dbc.Button("Apply mapping and preprocess", id="apply-mapping", color="success", className="mt-3"),
            ])
        ], className="mt-2")

        info = html.Div([
            html.H6("Columns in metadata.csv:"),
            chips,
        ])
        return info, mapping_ui

    @app.callback(
        Output("process-result", "children"),
        Input("apply-mapping", "n_clicks"),
        State("map-sample-id", "value"),
        State("map-batch-id", "value"),
        State("map-phenotype", "value"),
        State("session-id", "data"),
        prevent_initial_call=True,
    )
    def apply_mapping_and_preprocess(n_clicks: int, sample_col: str, batch_col: str, pheno_col: str, session_id: str):
        if not n_clicks:
            raise dash.exceptions.PreventUpdate
        if not session_id:
            return "Session not initialised."
        if not sample_col or not batch_col or not pheno_col:
            return "Please select all three columns before applying."

        session_dir = get_session_dir(session_id)
        meta_path = session_dir / "metadata.csv"
        if not meta_path.exists():
            return "metadata.csv not found in session."

        # Rename columns by rewriting CSV
        try:
            import csv
            with meta_path.open("r", encoding="utf-8", newline="") as fh:
                reader = csv.DictReader(fh)
                orig_fieldnames = reader.fieldnames or []
                rows = list(reader)
            # Build new header mapping
            new_fieldnames = []
            for name in orig_fieldnames:
                if name == sample_col:
                    new_fieldnames.append("sample_id")
                elif name == batch_col:
                    new_fieldnames.append("batch_id")
                elif name == pheno_col:
                    new_fieldnames.append("phenotype")
                else:
                    new_fieldnames.append(name)
            # Write back with new headers
            with meta_path.open("w", encoding="utf-8", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=new_fieldnames)
                writer.writeheader()
                for row in rows:
                    out_row = {}
                    for old, new in zip(orig_fieldnames, new_fieldnames):
                        out_row[new] = row.get(old)
                    writer.writerow(out_row)
        except Exception as exc:
            return f"Failed to rename metadata columns: {exc}"

        # Run preprocess.R in the session directory
        ok, log = run_preprocess(session_dir)
        prefix = "Preprocess succeeded." if ok else "Preprocess failed."
        return html.Div([
            html.Div(prefix),
            html.Pre(log, style={"whiteSpace": "pre-wrap"}),
        ])
