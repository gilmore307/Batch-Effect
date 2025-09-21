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
    render_group_tabset,
    build_ranking_tab,
    build_raw_assessments_tab,
    build_overall_div,
)


def assessment_layout(active_path: str, stage: str):
    header = "Pre-correction Assessment" if stage == "pre" else "Post-correction Assessment"
    gallery_id = "pre-assessment-gallery" if stage == "pre" else "post-assessment-gallery"

    # Define groups (key, title, script)
    pre_groups = [
        ("mosaic", "Mosaic plot", "Mosaic.R"),
        ("pca", "PCA", "pca.R"),
        ("pcoa", "PCoA", "pcoa.R"),
        ("nmds", "NMDS", "NMDS.R"),
        ("dissimilarity", "Dissimilarity heatmaps", "Dissimilarity_Heatmaps.R"),
        ("r2", "R^2", "R2.R"),
        ("prda", "pRDA", "pRDA.R"),
        ("pvca", "PVCA", "pvca.R"),
        ("alignment", "Alignment score", "Alignment_Score.R"),
        ("auc", "AUC", "AUC.R"),
    ]
    post_extra = [
        ("lisi", "LISI", "LISI.R"),
        ("ebm", "Entropy score", "Entropy_Score.R"),
        ("silhouette", "Silhouette score", "Silhouette.R"),
    ]
    groups = pre_groups + (post_extra if stage == "post" else [])

    # Build preset tabs with a Run button per group and placeholder until run
    tab_items = []
    for key, title, _ in groups:
        run_id = f"run-{stage}-{key}"
        content_id = f"{stage}-{key}-content"
        placeholder = html.Div("Click Run to generate results.")
        tab_items.append(
            dcc.Tab(
                label=title,
                value=f"tab-{key}",
                children=html.Div([
                    dbc.Button("Run", id=run_id, size="sm", color="primary", className="mb-2"),
                    dcc.Loading(html.Div(id=content_id, children=placeholder), type="default"),
                ]),
        )
        )

    # Overall tab at the end (auto-updates after any run)
    tab_items.append(
        dcc.Tab(
            label="Overall",
            value="tab-overall",
            children=html.Div(id=f"{stage}-overall-content", children=html.Div("No results yet.")),
        )
    )

    tabs = dcc.Tabs(children=tab_items, value=(tab_items[0].value if tab_items else None), vertical=True, className="be-results-tabs")

    return html.Div(
        [
            build_navbar(active_path),
            dbc.Container(
                [
                    html.H2(header),
                    html.Div(
                        "Click Run inside a tab to compute only that assessment.",
                        className="text-muted mb-3",
                    ),
                    tabs,
                ],
                fluid=True,
            ),
        ]
    )


def register_pre_post_callbacks(app):
    # Group definitions for per-tab runs
    pre_groups = [
        ("mosaic", "Mosaic plot", "Mosaic.R"),
        ("pca", "PCA", "pca.R"),
        ("pcoa", "PCoA", "pcoa.R"),
        ("nmds", "NMDS", "NMDS.R"),
        ("dissimilarity", "Dissimilarity heatmaps", "Dissimilarity_Heatmaps.R"),
        ("r2", "R^2", "R2.R"),
        ("prda", "pRDA", "pRDA.R"),
        ("pvca", "PVCA", "pvca.R"),
        ("alignment", "Alignment score", "Alignment_Score.R"),
        ("auc", "AUC", "AUC.R"),
    ]
    post_extra = [
        ("lisi", "LISI", "LISI.R"),
        ("ebm", "Entropy score", "Entropy_Score.R"),
        ("silhouette", "Silhouette score", "Silhouette.R"),
    ]

    def _register_group(stage: str, key: str, script_name: str):
        run_id = f"run-{stage}-{key}"
        content_id = f"{stage}-{key}-content"

        outputs = [Output(content_id, "children")]
        if stage == "pre":
            outputs.append(Output("pre-started", "data", allow_duplicate=True))
        if stage == "post":
            outputs.append(Output("post-complete", "data", allow_duplicate=True))
        outputs.extend([
            Output("runlog-path", "data", allow_duplicate=True),
            Output("runlog-modal", "is_open", allow_duplicate=True),
            Output("runlog-interval", "disabled", allow_duplicate=True),
            Output(f"{stage}-overall-content", "children", allow_duplicate=True),
        ])

        @app.callback(*outputs, Input(run_id, "n_clicks"), State("session-id", "data"), prevent_initial_call=True)
        def _run_one(n_clicks: int, session_id: str, _stage=stage, _key=key, _script=script_name):
            if not n_clicks:
                raise dash.exceptions.PreventUpdate
            if not session_id:
                if _stage == "pre":
                    return html.Div("Session not initialised."), True, dash.no_update, dash.no_update, dash.no_update, dash.no_update
                else:
                    return html.Div("Session not initialised."), True, dash.no_update, dash.no_update, dash.no_update, dash.no_update
            session_dir = get_session_dir(session_id)
            if not (session_dir / "raw.csv").exists() or not (session_dir / "metadata.csv").exists():
                if _stage == "pre":
                    return html.Div("Upload both raw.csv and metadata.csv first."), True, dash.no_update, dash.no_update, dash.no_update, dash.no_update
                else:
                    return html.Div("Upload both raw.csv and metadata.csv first."), True, dash.no_update, dash.no_update, dash.no_update, dash.no_update

            log_path = session_dir / f"{_stage}_{_key}.log"
            try:
                log_path.write_text("", encoding="utf-8")
            except Exception:
                pass
            success, _ = run_r_scripts((_script,), session_dir, log_path=log_path)
            content = render_group_tabset(session_dir, _stage, _key)
            overall = build_overall_div(session_dir, _stage)
            if _stage == "pre":
                return content, True, str(log_path), True, False, overall
            else:
                return content, True, str(log_path), True, False, overall

    # Register all group callbacks
    for key, _, script in pre_groups:
        _register_group("pre", key, script)
    for key, _, script in pre_groups + post_extra:
        _register_group("post", key, script)

