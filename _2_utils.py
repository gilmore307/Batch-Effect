# ===============================
# File: _2_utils.py
# ===============================
import base64
import csv
import os
import shutil
import subprocess
import textwrap
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple, Optional

from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc

# Paths & constants
BASE_DIR = Path(__file__).resolve().parent
OUTPUT_ROOT = BASE_DIR / "output"
PLOTS_DIR = BASE_DIR / "plots"
METHODS_SCRIPT = BASE_DIR / "methods.R"
PREPROCESS_SCRIPT = BASE_DIR / "preprocess.R"
CLEANUP_HOURS = 6


# ---- Dataclasses & config ----
@dataclass
class FigureSpec:
    label: str
    filename: str


PRE_FIGURES: Sequence[FigureSpec] = (
    FigureSpec("Mosaic plot", "mosaic_plot.png"),
    FigureSpec("PCA", "pca.png"),
    FigureSpec("PCoA (Aitchison)", "pcoa_aitchison.png"),
    FigureSpec("PCoA (Bray-Curtis)", "pcoa_braycurtis.png"),
    FigureSpec("NMDS (Aitchison)", "nmds_aitchison.png"),
    FigureSpec("NMDS (Bray-Curtis)", "nmds_braycurtis.png"),
    FigureSpec("Dissimilarity heatmaps (Aitchison)", "dissimilarity_heatmaps_aitchison.png"),
    FigureSpec("Dissimilarity heatmaps (Bray-Curtis)", "dissimilarity_heatmaps_braycurtis.png"),
    FigureSpec("R^2 (Aitchison)", "R2_aitchison.png"),
    FigureSpec("R^2 (Bray-Curtis)", "R2_braycurtis.png"),
    FigureSpec("pRDA (Aitchison)", "pRDA_aitchison.png"),
    FigureSpec("pRDA (Bray-Curtis)", "pRDA_braycurtis.png"),
    FigureSpec("PVCA", "PVCA.png"),
    FigureSpec("Alignment score", "alignment_score.png"),
    FigureSpec("AUC", "auroc.png"),
)

POST_EXTRA_FIGURES: Sequence[FigureSpec] = (
    FigureSpec("LISI", "LISI.png"),
    FigureSpec("Entropy score", "ebm.png"),
    FigureSpec("Silhouette score", "silhouette.png"),
)

POST_FIGURES: Sequence[FigureSpec] = PRE_FIGURES + POST_EXTRA_FIGURES


PRE_SCRIPTS: Sequence[str] = (
    "Mosaic.R",
    "pca.R",
    "pcoa.R",
    "NMDS.R",
    "Dissimilarity_Heatmaps.R",
    "R2.R",
    "pRDA.R",
    "pvca.R",
    "Alignment_Score.R",
    "AUC.R",
)

POST_SCRIPTS: Sequence[str] = PRE_SCRIPTS + (
    "LISI.R",
    "Entropy_Score.R",
    "Silhouette.R",
)


SUPPORTED_METHODS: Sequence[Tuple[str, str]] = (
    ("QN", "Quantile normalization"),
    ("BMC", "Batch mean centering"),
    ("limma", "Limma removeBatchEffect"),
    ("ConQuR", "ConQuR"),
    ("PLSDA", "PLSDA-batch"),
    ("ComBat", "ComBat"),
    ("FSQN", "FSQN"),
    ("MMUPHin", "MMUPHin"),
    ("RUV", "RUV-III-NB"),
    ("MetaDICT", "Meta-DICT"),
    ("SVD", "Surrogate variable decomposition"),
    ("PN", "Percentile normalization"),
    ("FAbatch", "FAbatch"),
    ("ComBatSeq", "ComBat-Seq"),
    ("DEBIAS", "DEBIAS"),
)

# Map method identifiers (as they may appear in CSVs) to formal display names
_FORMAL_MAP: Dict[str, str] = {code.lower(): display for code, display in SUPPORTED_METHODS}
# Add common lowercase variants used by R outputs
_FORMAL_ALIASES: Dict[str, str] = {
    "qn": _FORMAL_MAP["qn"],
    "bmc": _FORMAL_MAP["bmc"],
    "limma": _FORMAL_MAP["limma"],
    "conqur": _FORMAL_MAP["conqur"],
    "plsda": _FORMAL_MAP["plsda"],
    "plsdabatch": _FORMAL_MAP["plsda"],
    "combat": _FORMAL_MAP["combat"],
    "fsqn": _FORMAL_MAP["fsqn"],
    "mmuphin": _FORMAL_MAP["mmuphin"],
    "ruv": _FORMAL_MAP["ruv"],
    "ruv-iiinb": _FORMAL_MAP["ruv"],
    "ruviiinb": _FORMAL_MAP["ruv"],
    "metadict": _FORMAL_MAP["metadict"],
    "meta-dict": _FORMAL_MAP["metadict"],
    "svd": _FORMAL_MAP["svd"],
    "pn": _FORMAL_MAP["pn"],
    "percentilenormalization": _FORMAL_MAP["pn"],
    "fabatch": _FORMAL_MAP["fabatch"],
    "combatseq": _FORMAL_MAP["combatseq"],
    "combat-seq": _FORMAL_MAP["combatseq"],
    "debias": _FORMAL_MAP["debias"],
    "debias-m": _FORMAL_MAP["debias"],
}

def method_formal_name(name: str) -> str:
    if not name:
        return name
    key = (name or "").strip()
    # Keep special baseline label as-is
    if key.lower() == "before correction":
        return "Before correction"
    lk = key.replace(" ", "").replace("_", "").lower()
    return _FORMAL_ALIASES.get(lk, _FORMAL_MAP.get(lk, key))

DEFAULT_METHODS: Sequence[str] = ("ComBat", "FSQN")


# ---- General helpers ----

def human_size(num_bytes: int) -> str:
    units = ["B", "KB", "MB", "GB", "TB"]
    size = float(num_bytes)
    for u in units:
        if size < 1024 or u == units[-1]:
            return f"{size:.1f} {u}" if u != "B" else f"{int(size)} {u}"
        size /= 1024


def ensure_output_root() -> None:
    OUTPUT_ROOT.mkdir(exist_ok=True)


def cleanup_old_sessions(max_age_hours: int = CLEANUP_HOURS) -> None:
    if not OUTPUT_ROOT.exists():
        return
    cutoff = time.time() - max_age_hours * 3600
    for session_dir in OUTPUT_ROOT.iterdir():
        try:
            if not session_dir.is_dir():
                continue
            if session_dir.stat().st_mtime < cutoff:
                shutil.rmtree(session_dir, ignore_errors=True)
        except OSError:
            continue


def get_session_dir(session_id: str) -> Path:
    ensure_output_root()
    session_path = OUTPUT_ROOT / session_id
    session_path.mkdir(exist_ok=True, parents=True)
    return session_path


def decode_contents(contents: str) -> bytes:
    header, encoded = contents.split(",", 1)
    return base64.b64decode(encoded)


def save_uploaded_file(contents: str, directory: Path, dest_name: str) -> None:
    data = decode_contents(contents)
    (directory / dest_name).write_bytes(data)


def run_command(command: Sequence[str], cwd: Path) -> Tuple[bool, str]:
    def _s(val) -> str:
        if val is None:
            return ""
        try:
            return val.strip()
        except Exception:
            try:
                return str(val).strip()
            except Exception:
                return ""

    try:
        result = subprocess.run(
            command,
            cwd=cwd,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=True,
        )
        log = textwrap.dedent(
            f"""
            $ {' '.join(command)}
            {_s(result.stdout)}
            {_s(result.stderr)}
            """
        ).strip()
        return True, log
    except subprocess.CalledProcessError as exc:
        stdout = getattr(exc, "stdout", getattr(exc, "output", None))
        stderr = getattr(exc, "stderr", None)
        log = textwrap.dedent(
            f"""
            $ {' '.join(command)}
            {_s(stdout)}
            {_s(stderr)}
            """
        ).strip()
        return False, log


def run_command_streaming(command: Sequence[str], cwd: Path, log_path: Path) -> Tuple[bool, str]:
    """Run a command streaming stdout to a log file for real-time viewing.

    Returns (success, combined_log_text_at_end).
    """
    log_path.parent.mkdir(parents=True, exist_ok=True)
    combined: List[str] = []
    try:
        with log_path.open("a", encoding="utf-8", errors="replace") as logf:
            header = "$ " + " ".join(command) + "\n"
            logf.write(header)
            logf.flush()
            proc = subprocess.Popen(
                command,
                cwd=cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
            )
            assert proc.stdout is not None
            for line in proc.stdout:
                logf.write(line)
                logf.flush()
                combined.append(line)
            rc = proc.wait()
            success = (rc == 0)
            return success, header + "".join(combined)
    except Exception as exc:
        return False, f"$ {' '.join(command)}\n{exc}"


def run_r_scripts(script_names: Sequence[str], output_dir: Path, log_path: Optional[Path] = None) -> Tuple[bool, str]:
    logs: List[str] = []
    for script in script_names:
        script_path = PLOTS_DIR / script
        if not script_path.exists():
            logs.append(f"Warning: Script not found: {script}")
            continue
        cmd = ("Rscript", str(script_path), str(output_dir))
        if log_path is not None:
            success, log = run_command_streaming(cmd, cwd=BASE_DIR, log_path=log_path)
        else:
            success, log = run_command(cmd, cwd=BASE_DIR)
        logs.append(log)
        if not success:
            return False, "\n\n".join(logs)
    return True, "\n\n".join(logs)


def run_methods(session_dir: Path, methods: Iterable[str], log_path: Optional[Path] = None) -> Tuple[bool, str]:
    method_arg = ",".join(methods)
    command = ("Rscript", str(METHODS_SCRIPT), method_arg, str(session_dir))
    if log_path is not None:
        return run_command_streaming(command, cwd=BASE_DIR, log_path=log_path)
    return run_command(command, cwd=BASE_DIR)


def run_preprocess(session_dir: Path, log_path: Optional[Path] = None) -> Tuple[bool, str]:
    """Run preprocess.R in the project root for a given session directory.

    It writes outputs into the provided session directory and reads the
    matrix from session_dir/raw.csv by default.
    """
    if not PREPROCESS_SCRIPT.exists():
        return False, f"Script not found: {PREPROCESS_SCRIPT.name}"
    matrix_path = session_dir / "raw.csv"
    command = ("Rscript", str(PREPROCESS_SCRIPT), str(session_dir), str(matrix_path))
    if log_path is not None:
        return run_command_streaming(command, cwd=BASE_DIR, log_path=log_path)
    return run_command(command, cwd=BASE_DIR)


def render_figures(session_dir: Path, figures: Sequence[FigureSpec]):
    cards = []
    for spec in figures:
        file_path = session_dir / spec.filename
        if not file_path.exists():
            continue
        encoded = base64.b64encode(file_path.read_bytes()).decode("ascii")
        src = f"data:image/png;base64,{encoded}"
        cards.append(
            dbc.Col(
                dbc.Card(
                    [
                        dbc.CardImg(src=src, top=True),
                        dbc.CardBody(html.H6(spec.label, className="card-title")),
                    ],
                    className="h-100",
                ),
                xs=12,
                sm=6,
                lg=4,
                className="mb-4",
            )
        )
    if not cards:
        return html.Div("No figures available yet. Run the analysis to generate outputs.")
    return dbc.Row(cards, className="gy-2")


def _read_csv_rows(csv_path: Path) -> Tuple[List[str], List[List[str]]]:
    import csv as _csv
    with csv_path.open("r", encoding="utf-8", newline="") as fh:
        reader = _csv.reader(fh)
        rows = list(reader)
    if not rows:
        return [], []
    header = rows[0]
    data = rows[1:] if len(rows) > 1 else []
    return header, data


def _find_file_case_insensitive(directory: Path, target_name: str) -> Path | None:
    t = target_name.lower()
    for p in directory.iterdir():
        if p.is_file() and p.name.lower() == t:
            return p
    return None


def _candidate_csvs_for_image(filename: str) -> List[str]:
    stem = Path(filename).stem
    s = stem.lower()
    bases: List[str] = []
    # image families with shared tables
    if s.startswith("pcoa_"):
        bases = ["pcoa"]
    elif s.startswith("nmds_"):
        bases = ["nmds"]
    elif s.startswith("dissimilarity_") or s.startswith("dissimilarity-") or s.startswith("dissimilarity"):
        bases = ["dissimilarity"]
    elif s.startswith("r2_"):
        bases = ["r2"]
    elif s.startswith("prda_") or s.startswith("prda"):
        bases = ["pRDA", "prda"]
    elif s == "pvca":
        bases = ["PVCA", "pvca"]
    elif s == "alignment_score":
        bases = ["alignment_score"]
    elif s == "auroc":
        bases = ["auroc"]
    elif s == "lisi":
        bases = ["LISI", "lisi"]
    elif s == "ebm":
        bases = ["ebm"]
    elif s == "silhouette":
        bases = ["silhouette"]
    elif s == "pca":
        bases = ["pca"]
    else:
        bases = [stem]

    # Generate candidate filenames (prefer ranking first)
    candidates: List[str] = []
    for b in bases:
        # ranking summary across methods
        candidates.append(f"{b}_ranking.csv")
    for b in bases:
        # baseline-only assessment
        candidates.append(f"{b}_raw_assessment.csv")
    # Special LISI grid file
    for b in bases:
        if b.lower() == "lisi":
            candidates.append(f"{b}_k_grid.csv")
    return candidates


def render_assessment_tabs(session_dir: Path, figures: Sequence[FigureSpec], stage: str = "pre", extra_tabs: Sequence = ()):  # extra dcc.Tab items appended
    """Render a vertical tab set of assessment outputs.

    Group metrics by base (e.g., PCoA, NMDS, R2, pRDA, Dissimilarity heatmaps).
    For groups with multiple geometries (Aitchison/Bray-Curtis), show a single
    top-level tab with sub-tabs: Aitchison (CLR), Bray-Curtis (TSS), and a third
    sub-tab that is either Assessment (stage='pre') or Rank (stage='post').
    Single-geometry plots also include the third sub-tab accordingly.
    """

    # Consistent styles for top-level and inner tabs (keeps size stable)
    TOP_TAB_STYLE = {
        "borderBottom": "1px solid #d6d6d6",
        "padding": "6px",
        "fontWeight": "bold",
        "width": "250px",
        "height": "60px",
        "display": "flex",
        "alignItems": "center",
        "justifyContent": "center",
    }
    TOP_TAB_SELECTED_STYLE = {
        "borderTop": "1px solid #d6d6d6",
        "borderBottom": "1px solid #d6d6d6",
        "backgroundColor": "#f8f9fa",
        "color": "#0d6efd",
        "padding": "6px",
        "width": "250px",
        "height": "60px",
        "display": "flex",
        "alignItems": "center",
        "justifyContent": "center",
    }
    SUBTAB_STYLE_BASE = {
        "borderBottom": "1px solid #d6d6d6",
        "padding": "6px",
        "fontWeight": "bold",
        "height": "60px",
        "display": "inline-flex",
        "alignItems": "center",
        "justifyContent": "center",
        # Fixed width to keep constant size, no shrinking
        "width": "250px",
        "minWidth": "250px",
        "flex": "0 0 250px",
        "boxSizing": "border-box",
    }
    SUBTAB_SELECTED_STYLE_BASE = {
        "borderTop": "1px solid #d6d6d6",
        "borderBottom": "1px solid #d6d6d6",
        "backgroundColor": "#f8f9fa",
        "color": "#0d6efd",
        "padding": "6px",
        "height": "60px",
        "display": "inline-flex",
        "alignItems": "center",
        "justifyContent": "center",
        # Fixed width to keep constant size, no shrinking
        "width": "250px",
        "minWidth": "250px",
        "flex": "0 0 250px",
        "boxSizing": "border-box",
    }

    def content_for_image(filename: str) -> html.Div:
        img_path = session_dir / filename
        if not img_path.exists():
            return html.Div("Image not found.")
        encoded = base64.b64encode(img_path.read_bytes()).decode("ascii")
        src = f"data:image/png;base64,{encoded}"
        img = html.Img(src=src, style={"maxWidth": "100%", "height": "auto"})
        return html.Div([img])

    tabs = []
    first_value = None

    # Build groups
    groups = {}

    def add_group_item(key: str, title: str, geom: str, filename: str):
        g = groups.setdefault(key, {"title": title, "ait": None, "bray": None, "single": None})
        if geom == "ait":
            g["ait"] = filename
        elif geom == "bray":
            g["bray"] = filename
        else:
            g["single"] = filename

    for spec in figures:
        fn = spec.filename
        low = fn.lower()
        if low.startswith("pcoa_aitchison"):
            add_group_item("pcoa", "PCoA", "ait", fn)
        elif low.startswith("pcoa_braycurtis"):
            add_group_item("pcoa", "PCoA", "bray", fn)
        elif low.startswith("nmds_aitchison"):
            add_group_item("nmds", "NMDS", "ait", fn)
        elif low.startswith("nmds_braycurtis"):
            add_group_item("nmds", "NMDS", "bray", fn)
        elif low.startswith("dissimilarity_heatmaps_aitchison"):
            add_group_item("dissimilarity", "Dissimilarity heatmaps", "ait", fn)
        elif low.startswith("dissimilarity_heatmaps_braycurtis"):
            add_group_item("dissimilarity", "Dissimilarity heatmaps", "bray", fn)
        elif low.startswith("r2_aitchison"):
            add_group_item("r2", "R^2", "ait", fn)
        elif low.startswith("r2_braycurtis"):
            add_group_item("r2", "R^2", "bray", fn)
        elif low.startswith("prda_aitchison"):
            add_group_item("prda", "pRDA", "ait", fn)
        elif low.startswith("prda_braycurtis"):
            add_group_item("prda", "pRDA", "bray", fn)
        elif low == "pca.png":
            add_group_item("pca", "PCA", "ait", fn)
        else:
            # fallback: one tab per figure (with its own assessment sub-tab if any)
            key = f"single:{low}"
            add_group_item(key, spec.label, "single", fn)

    # Build tabs with sub-tabs
    for idx, (key, g) in enumerate(groups.items()):
        title = g["title"]
        sub_tabs = []
        # choose a representative filename for assessment lookup
        rep = g["ait"] or g["bray"] or g["single"]

        # Build sub-tab definitions first so we can size them evenly later
        sub_defs = []
        if g["ait"]:
            sub_defs.append(("Aitchison (CLR)", f"{key}-ait", content_for_image(g["ait"])) )
        if g["bray"]:
            sub_defs.append(("Bray-Curtis (TSS)", f"{key}-bray", content_for_image(g["bray"])) )

        # Third sub-tab: Assessment (pre) or Rank (post)
        third_label = "Assessment" if stage == "pre" else "Rank"
        third_content = None
        if rep:
            for cand in _candidate_csvs_for_image(rep):
                if stage == "pre" and cand.endswith("_raw_assessment.csv"):
                    found = _find_file_case_insensitive(session_dir, cand)
                    if found and found.exists():
                        header, data = _read_csv_rows(found)
                        if header:
                            thead = html.Thead(html.Tr([html.Th(h) for h in header]))
                            tbody = html.Tbody([html.Tr([html.Td(cell) for cell in row]) for row in data])
                            third_content = dbc.Table([thead, tbody], bordered=True, hover=True, size="sm")
                        break
                if stage == "post" and cand.endswith("_ranking.csv"):
                    found = _find_file_case_insensitive(session_dir, cand)
                    if found and found.exists():
                        header, data = _read_csv_rows(found)
                        if header:
                            thead = html.Thead(html.Tr([html.Th(h) for h in header]))
                            tbody = html.Tbody([html.Tr([html.Td(cell) for cell in row]) for row in data])
                            third_content = dbc.Table([thead, tbody], bordered=True, hover=True, size="sm")
                        break
        if third_content is None:
            third_content = html.Div("No table found.")
        sub_defs.append((third_label, f"{key}-third", html.Div(third_content)))

        # Create sub-tabs with fixed width; allow horizontal scrolling in container
        SUBTAB_STYLE = dict(SUBTAB_STYLE_BASE)
        SUBTAB_SELECTED_STYLE = dict(SUBTAB_SELECTED_STYLE_BASE)
        sub_tabs = [
            dcc.Tab(label=lbl, value=val, children=child, style=SUBTAB_STYLE, selected_style=SUBTAB_SELECTED_STYLE)
            for (lbl, val, child) in sub_defs
        ]

        inner_default = sub_tabs[0].value if sub_tabs else None
        inner = dcc.Tabs(
            className="be-subtabs",
            children=sub_tabs,
            value=inner_default,
            vertical=False,
            mobile_breakpoint=0,
            # Keep in one row and allow horizontal scroll
            style={
                "width": "100%",
                "display": "block",
                "overflowX": "auto",
                "whiteSpace": "nowrap",
            },
        )

        top_value = f"tab-{key}"
        tabs.append(dcc.Tab(label=title, value=top_value, children=html.Div(inner), style=TOP_TAB_STYLE, selected_style=TOP_TAB_SELECTED_STYLE))
        if first_value is None:
            first_value = top_value

    if not tabs and not extra_tabs:
        return html.Div("No figures available yet. Run the analysis to generate outputs.")

    # Append any extra tabs (e.g., raw assessments, ranking)
    for extra in extra_tabs:
        if first_value is None:
            first_value = getattr(extra, 'value', None) or 'extra-0'
        tabs.append(extra)

    return dcc.Tabs(
        children=tabs,
        value=first_value or (extra_tabs[0].value if extra_tabs else None),
        vertical=True,
        className="be-results-tabs",
    )


def render_group_tabset(session_dir: Path, stage: str, key: str):
    """Render only one group's inner tabset (Aitchison/Bray + Assessment/Rank).

    key examples: 'pcoa','nmds','dissimilarity','r2','prda','pca','mosaic',
                  'pvca','alignment','auc','lisi','ebm','silhouette'
    """
    figures = PRE_FIGURES if stage == "pre" else POST_FIGURES

    # Keep sub-tab size consistent using inline styles
    SUBTAB_STYLE_BASE = {
        "borderBottom": "1px solid #d6d6d6",
        "padding": "6px",
        "fontWeight": "bold",
        "height": "80px",
        "display": "inline-flex",
        "alignItems": "center",
        "justifyContent": "center",
        # Fixed width to keep constant size, no shrinking
        "width": "250px",
        "minWidth": "250px",
        "flex": "0 0 250px",
        "boxSizing": "border-box",
    }
    SUBTAB_SELECTED_STYLE_BASE = {
        "borderTop": "1px solid #d6d6d6",
        "borderBottom": "1px solid #d6d6d6",
        "backgroundColor": "#f8f9fa",
        "color": "#0d6efd",
        "padding": "6px",
        "height": "80px",
        "display": "inline-flex",
        "alignItems": "center",
        "justifyContent": "center",
        # Fixed width to keep constant size, no shrinking
        "width": "250px",
        "minWidth": "250px",
        "flex": "0 0 250px",
        "boxSizing": "border-box",
    }

    def content_for_image(filename: str) -> html.Div:
        img_path = session_dir / filename
        if not img_path.exists():
            return html.Div("Image not found.")
        encoded = base64.b64encode(img_path.read_bytes()).decode("ascii")
        src = f"data:image/png;base64,{encoded}"
        img = html.Img(src=src, style={"maxWidth": "100%", "height": "auto"})
        return html.Div([img])

    # Resolve filenames in this group
    g: Dict[str, Optional[str]] = {"ait": None, "bray": None, "single": None, "title": key}
    for spec in figures:
        low = spec.filename.lower()
        if key == "pcoa":
            if low.startswith("pcoa_aitchison"):
                g["ait"] = spec.filename
            elif low.startswith("pcoa_braycurtis"):
                g["bray"] = spec.filename
            g["title"] = "PCoA"
        elif key == "nmds":
            if low.startswith("nmds_aitchison"):
                g["ait"] = spec.filename
            elif low.startswith("nmds_braycurtis"):
                g["bray"] = spec.filename
            g["title"] = "NMDS"
        elif key == "dissimilarity":
            if low.startswith("dissimilarity_heatmaps_aitchison"):
                g["ait"] = spec.filename
            elif low.startswith("dissimilarity_heatmaps_braycurtis"):
                g["bray"] = spec.filename
            g["title"] = "Dissimilarity heatmaps"
        elif key == "r2":
            if low.startswith("r2_aitchison"):
                g["ait"] = spec.filename
            elif low.startswith("r2_braycurtis"):
                g["bray"] = spec.filename
            g["title"] = "R^2"
        elif key == "prda":
            if low.startswith("prda_aitchison"):
                g["ait"] = spec.filename
            elif low.startswith("prda_braycurtis"):
                g["bray"] = spec.filename
            g["title"] = "pRDA"
        elif key == "pca" and spec.filename.lower() == "pca.png":
            g["single"] = spec.filename; g["title"] = "PCA"
        elif key == "mosaic" and spec.filename.lower() == "mosaic_plot.png":
            g["single"] = spec.filename; g["title"] = "Mosaic plot"
        elif key == "pvca" and spec.filename.lower() == "pvca.png":
            g["single"] = spec.filename; g["title"] = "PVCA"
        elif key == "alignment" and spec.filename.lower() == "alignment_score.png":
            g["single"] = spec.filename; g["title"] = "Alignment score"
        elif key == "auc" and spec.filename.lower() == "auroc.png":
            g["single"] = spec.filename; g["title"] = "AUC"
        elif key == "lisi" and spec.filename.lower() == "lisi.png":
            g["single"] = spec.filename; g["title"] = "LISI"
        elif key == "ebm" and spec.filename.lower() == "ebm.png":
            g["single"] = spec.filename; g["title"] = "Entropy score"
        elif key == "silhouette" and spec.filename.lower() == "silhouette.png":
            g["single"] = spec.filename; g["title"] = "Silhouette score"

    sub_tabs = []
    sub_defs = []
    if g["ait"]:
        sub_defs.append(("Aitchison (CLR)", f"{key}-ait", content_for_image(g["ait"])) )
    if g["bray"]:
        sub_defs.append(("Bray-Curtis (TSS)", f"{key}-bray", content_for_image(g["bray"])) )
    if g["single"] and not (g["ait"] or g["bray"]):
        sub_defs.append(("Plot", f"{key}-plot", content_for_image(g["single"])) )

    rep = g["ait"] or g["bray"] or g["single"]
    third_label = "Assessment" if stage == "pre" else "Rank"
    third_content = None
    if rep:
        for cand in _candidate_csvs_for_image(rep):
            if stage == "pre" and cand.endswith("_raw_assessment.csv"):
                found = _find_file_case_insensitive(session_dir, cand)
                if found and found.exists():
                    header, data = _read_csv_rows(found)
                    if header:
                        thead = html.Thead(html.Tr([html.Th(h) for h in header]))
                        tbody = html.Tbody([html.Tr([html.Td(cell) for cell in row]) for row in data])
                        third_content = dbc.Table([thead, tbody], bordered=True, hover=True, size="sm")
                    break
            if stage == "post" and cand.endswith("_ranking.csv"):
                found = _find_file_case_insensitive(session_dir, cand)
                if found and found.exists():
                    header, data = _read_csv_rows(found)
                    if header:
                        thead = html.Thead(html.Tr([html.Th(h) for h in header]))
                        tbody = html.Tbody([html.Tr([html.Td(cell) for cell in row]) for row in data])
                        third_content = dbc.Table([thead, tbody], bordered=True, hover=True, size="sm")
                    break
    if third_content is not None:
        sub_defs.append((third_label, f"{key}-third", html.Div(third_content)))

    SUBTAB_STYLE = dict(SUBTAB_STYLE_BASE)
    SUBTAB_SELECTED_STYLE = dict(SUBTAB_SELECTED_STYLE_BASE)
    sub_tabs = [
        dcc.Tab(label=lbl, value=val, children=child, style=SUBTAB_STYLE, selected_style=SUBTAB_SELECTED_STYLE)
        for (lbl, val, child) in sub_defs
    ]

    return dcc.Tabs(
        children=sub_tabs,
        value=sub_tabs[0].value if sub_tabs else None,
        className="be-subtabs",
        style={
            "width": "100%",
            "display": "block",
            "overflowX": "auto",
            "whiteSpace": "nowrap",
        },
        vertical=False,
        mobile_breakpoint=0,
    )


def aggregate_rankings(session_dir: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    ranking_files = sorted(session_dir.glob("*_ranking.csv"))
    if not ranking_files:
        return [], []

    metric_names: List[str] = []
    metric_ranks: Dict[str, Dict[str, float]] = {}
    method_ranks: defaultdict[str, List[float]] = defaultdict(list)

    for csv_path in ranking_files:
        metric_name = csv_path.stem.replace("_ranking", "").replace("_", " ").title()
        metric_names.append(metric_name)
        with csv_path.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                continue
            field_map = {name.lower(): name for name in reader.fieldnames}
            method_key = field_map.get("method")
            rank_key = field_map.get("rank")
            if not method_key or not rank_key:
                continue
            for row in reader:
                method = (row.get(method_key, "") or "").strip()
                val = (row.get(rank_key, "") or "").strip()
                if not method or not val or val.lower() == "na":
                    continue
                try:
                    rank = float(val)
                except ValueError:
                    continue
                disp = method_formal_name(method)
                metric_ranks.setdefault(metric_name, {})[disp] = rank
                method_ranks[disp].append(rank)

    if not metric_ranks:
        return [], []

    unique_metrics = list(metric_ranks.keys())
    methods = sorted(method_ranks)
    rows: List[Dict[str, str]] = []
    for method in methods:
        avg_rank = sum(method_ranks[method]) / len(method_ranks[method])
        row = {
            "Method": method,
            "Average rank": f"{avg_rank:.2f}",
            "Metrics": str(len(method_ranks[method])),
        }
        for metric in unique_metrics:
            value = metric_ranks.get(metric, {}).get(method)
            if value is None:
                row[metric] = "-"
            else:
                # Show integer ranks (no values) in overall table
                try:
                    ival = int(round(value))
                    row[metric] = str(ival)
                except Exception:
                    row[metric] = str(value)
        rows.append(row)
    rows.sort(key=lambda item: float(item["Average rank"]))
    return unique_metrics, rows


def build_ranking_tab(session_dir: Path):
    metrics, rows = aggregate_rankings(session_dir)
    if not metrics:
        content = html.Div("No ranking files found. Run the post-correction assessment first.")
    else:
        header = ["Method", "Average rank", "Metrics", *metrics]
        table_header = html.Thead(html.Tr([html.Th(col) for col in header]))
        table_body = html.Tbody([html.Tr([html.Td(row.get(col, "")) for col in header]) for row in rows])
        content = dbc.Table([table_header, table_body], bordered=True, hover=True, responsive=True, className="text-nowrap")
    return dcc.Tab(label="Method Ranking", value="tab-ranking", children=html.Div(content))


def build_assessment_overview_table(session_dir: Path, stage: str):
    # Build a single integrated table: Test | Criteria | Test Result | Comment
    csv_files = list(session_dir.glob("*.csv"))
    candidates = [p for p in csv_files if p.name.lower().endswith("_raw_assessment.csv")]
    rows: List[Dict[str, str]] = []

    def pick_value(header: List[str], data_row: List[str]) -> str:
        hl = [h.lower() for h in header]
        pri = ["score", "auc", "silhouette", "entropy", "alignment", "r2", "mean_between"]
        for key in pri:
            if key in hl:
                idx = hl.index(key)
                return data_row[idx]
        # fallback to first numeric-looking cell
        for i, v in enumerate(data_row):
            try:
                float(v)
                return v
            except Exception:
                continue
        return "-"

    for path in sorted(candidates):
        test_name = path.stem.replace("_raw_assessment", "").replace("_", " ").title()
        header, data = _read_csv_rows(path)
        if not header or not data:
            continue
        # Prefer the row with Method == 'Before correction'
        method_idx = None
        try:
            method_idx = [h.lower() for h in header].index("method")
        except ValueError:
            method_idx = None
        row = None
        if method_idx is not None:
            for r in data:
                if (r[method_idx] or "").strip().lower() == "before correction":
                    row = r; break
        if row is None:
            row = data[0]
        # Criteria
        crit = "-"
        if "criteria" in [h.lower() for h in header]:
            crit = row[[h.lower() for h in header].index("criteria")]
        # Comment from Needs_Correction if available
        comment = "-"
        if "needs_correction" in [h.lower() for h in header]:
            val = row[[h.lower() for h in header].index("needs_correction")]
            comment = "Need correction" if str(val).strip().upper() in ("TRUE", "1") else "No correction"
        rows.append({
            "Test": test_name,
            "Criteria": crit,
            "Test Result": pick_value(header, row),
            "Comment": comment,
        })

    columns = [
        {"name": "Test", "id": "Test"},
        {"name": "Criteria", "id": "Criteria"},
        {"name": "Test Result", "id": "Test Result"},
        {"name": "Comment", "id": "Comment"},
    ]
    return dash_table.DataTable(
        data=rows,
        columns=columns,
        sort_action="native",
        filter_action="native",
        page_size=20,
        style_table={"overflowX": "auto"},
        style_cell={"padding": "6px", "textAlign": "left"},
    )


def build_overall_div(session_dir: Path, stage: str):
    if stage == "post":
        # Build an interactive ranking table using all *_ranking.csv files
        metrics, rows = aggregate_rankings(session_dir)
        if not metrics:
            body = html.Div("No ranking files found. Run post-assessment tabs to generate *_ranking.csv files.")
        else:
            # Convert to typed values for proper sorting
            typed_rows: List[Dict[str, object]] = []
            for idx, r in enumerate(rows, start=1):
                tr: Dict[str, object] = {"Method": r.get("Method", "")}
                # Add Rank column (1-based by Average rank order)
                try:
                    avg = float(r.get("Average rank", "nan"))
                except Exception:
                    avg = float("nan")
                tr["Rank"] = idx
                tr["Average rank"] = avg
                # Metrics count
                try:
                    tr["Metrics"] = int(r.get("Metrics", 0))
                except Exception:
                    tr["Metrics"] = r.get("Metrics", 0)
                for m in metrics:
                    val = r.get(m, "-")
                    if val == "-" or val is None:
                        tr[m] = None
                    else:
                        try:
                            tr[m] = int(val)
                        except Exception:
                            try:
                                tr[m] = float(val)
                            except Exception:
                                tr[m] = val
                typed_rows.append(tr)

            columns = (
                [{"name": "Rank", "id": "Rank", "type": "numeric"},
                 {"name": "Method", "id": "Method"},
                 {"name": "Average rank", "id": "Average rank", "type": "numeric"},
                 {"name": "Metrics", "id": "Metrics", "type": "numeric"}]
                + [{"name": m, "id": m, "type": "numeric"} for m in metrics]
            )

            body = dash_table.DataTable(
                data=typed_rows,
                columns=columns,
                sort_action="native",
                sort_by=[{"column_id": "Rank", "direction": "asc"}],
                filter_action="native",
                page_size=25,
                style_table={"overflowX": "auto"},
                style_cell={
                    "padding": "6px",
                    "textAlign": "left",
                    "minWidth": "80px",
                    "width": "120px",
                    "maxWidth": "260px",
                    "whiteSpace": "normal",
                },
                style_header={"fontWeight": "600"},
                fill_width=True,
            )

        return html.Div([
            html.H6("Overall Ranking"),
            body,
        ])
    # Pre: show overall summary of raw assessment values
    return html.Div([
        html.H6("Overall Summary"),
        build_assessment_overview_table(session_dir, stage),
    ])


def build_raw_assessments_tab(session_dir: Path):
    # Collect all *_raw_assessment.csv files (case-insensitive)
    csv_files = list(session_dir.glob("*.csv"))
    candidates = [p for p in csv_files if p.name.lower().endswith("_raw_assessment.csv")]
    if not candidates:
        body = html.Div("No raw assessment files found.")
        return dcc.Tab(label="Raw Assessments", value="tab-raw", children=html.Div(body))

    items = []
    for path in sorted(candidates):
        # Derive metric name from filename
        metric = path.stem.replace("_raw_assessment", "").replace("_", " ").title()
        header, data = _read_csv_rows(path)
        if not header:
            continue
        thead = html.Thead(html.Tr([html.Th(h) for h in header]))
        tbody = html.Tbody([html.Tr([html.Td(cell) for cell in row]) for row in data])
        table = dbc.Table([thead, tbody], bordered=True, hover=True, size="sm")
        items.append(
            dbc.AccordionItem(table, title=metric)
        )

    if not items:
        body = html.Div("No valid raw assessment tables found.")
    else:
        body = dbc.Accordion(items, always_open=False)

    return dcc.Tab(label="Raw Assessments", value="tab-raw", children=html.Div(body))

