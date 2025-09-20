# ===============================
# File: _7_ranking.py
# ===============================
from dash import html
from dash.dependencies import Input, Output, State
import dash
import dash_bootstrap_components as dbc

from _1_components import build_navbar
from _2_utils import get_session_dir, aggregate_rankings


def ranking_layout(active_path: str):
    return html.Div(
        [
            build_navbar(active_path),
            dbc.Container(
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                html.H2("Method ranking"),
                                dbc.Button("Refresh ranking", id="refresh-ranking", color="primary", className="mb-3"),
                                html.Div(id="ranking-table"),
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


def register_ranking_callbacks(app):
    @app.callback(
        Output("ranking-table", "children"),
        Input("refresh-ranking", "n_clicks"),
        State("session-id", "data"),
        prevent_initial_call=True,
    )
    def update_ranking(n_clicks: int, session_id: str):
        if not n_clicks:
            raise dash.exceptions.PreventUpdate
        session_dir = get_session_dir(session_id)
        metrics, rows = aggregate_rankings(session_dir)
        if not metrics:
            return html.Div("No ranking files found. Run the post-correction assessment first.")

        header = ["Method", "Average rank", "Metrics", *metrics]
        table_header = html.Thead(html.Tr([html.Th(col) for col in header]))
        table_body = html.Tbody([html.Tr([html.Td(row.get(col, "")) for col in header]) for row in rows])
        table = dbc.Table([table_header, table_body], bordered=True, hover=True, responsive=True, className="text-nowrap")
        return table
