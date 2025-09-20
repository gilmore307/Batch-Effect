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
from typing import Dict, Iterable, List, Sequence, Tuple

from dash import html, dcc
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


def run_r_scripts(script_names: Sequence[str], output_dir: Path) -> Tuple[bool, str]:
    logs: List[str] = []
    for script in script_names:
        script_path = PLOTS_DIR / script
        if not script_path.exists():
            logs.append(f"Warning: Script not found: {script}")
            continue
        success, log = run_command(("Rscript", str(script_path), str(output_dir)), cwd=BASE_DIR)
        logs.append(log)
        if not success:
            return False, "\n\n".join(logs)
    return True, "\n\n".join(logs)


def run_methods(session_dir: Path, methods: Iterable[str]) -> Tuple[bool, str]:
    method_arg = ",".join(methods)
    command = ("Rscript", str(METHODS_SCRIPT), method_arg, str(session_dir))
    return run_command(command, cwd=BASE_DIR)


def run_preprocess(session_dir: Path) -> Tuple[bool, str]:
    """Run preprocess.R in the project root for a given session directory.

    It writes outputs into the provided session directory and reads the
    matrix from session_dir/raw.csv by default.
    """
    if not PREPROCESS_SCRIPT.exists():
        return False, f"Script not found: {PREPROCESS_SCRIPT.name}"
    matrix_path = session_dir / "raw.csv"
    command = ("Rscript", str(PREPROCESS_SCRIPT), str(session_dir), str(matrix_path))
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


def render_assessment_tabs(session_dir: Path, figures: Sequence[FigureSpec]):
    """Render a tab set where each tab shows an image and an optional table.

    The table is looked up by taking the PNG filename stem and using a CSV
    with the same stem (e.g., pca.png -> pca.csv). If not found, only the
    image is shown.
    """
    tabs = []
    first_value = None
    for idx, spec in enumerate(figures):
        img_path = session_dir / spec.filename
        if not img_path.exists():
            continue
        # Image
        encoded = base64.b64encode(img_path.read_bytes()).decode("ascii")
        src = f"data:image/png;base64,{encoded}"
        img = html.Img(src=src, style={"maxWidth": "100%", "height": "auto"})

        # Table (optional)
        table_comp = None
        # Try to resolve the matching CSV based on naming from R scripts
        for cand in _candidate_csvs_for_image(spec.filename):
            found = _find_file_case_insensitive(session_dir, cand)
            if found and found.exists():
                header, data = _read_csv_rows(found)
                if header:
                    thead = html.Thead(html.Tr([html.Th(h) for h in header]))
                    tbody = html.Tbody([html.Tr([html.Td(cell) for cell in row]) for row in data])
                    table_comp = dbc.Table([thead, tbody], bordered=True, hover=True, size="sm", className="mt-2")
                break

        content_children = [img]
        if table_comp is not None:
            content_children.append(table_comp)

        tab = dcc.Tab(label=spec.label, value=f"tab-{idx}", children=html.Div(content_children))
        tabs.append(tab)
        if first_value is None:
            first_value = f"tab-{idx}"

    if not tabs:
        return html.Div("No figures available yet. Run the analysis to generate outputs.")

    return dcc.Tabs(children=tabs, value=first_value)


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
            row[metric] = f"{value:.2f}" if value is not None else "-"
        rows.append(row)
    rows.sort(key=lambda item: float(item["Average rank"]))
    return unique_metrics, rows
