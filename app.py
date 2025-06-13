import dash
from dash import html, Output, Input, State, dcc
import dash_bootstrap_components as dbc
import subprocess
import uuid
import os
import shutil
import redis
import base64

# ------------------------------
# Redis Setup
# ------------------------------
r = redis.Redis(host='localhost', port=6379, decode_responses=True)
SESSION_TTL = 60  # 5 minutes

# ------------------------------
# Dash App Setup
# ------------------------------
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

methods_options = [
    {"label": "ConQuR", "value": "ConQuR"},
    {"label": "ALRA", "value": "ALRA"},
    {"label": "PLSDAbatch", "value": "PLSDA"},
    {"label": "ComBat", "value": "ComBat"},
    {"label": "FSQN", "value": "FSQN"}
]

def main_layout():
    return dbc.Container([
        html.H2("Upload Input Data"),

        html.Div([
            html.Label("Upload Count Matrix (.csv):"),
            dcc.Upload(
                id="upload-matrix",
                children=html.Div(["Drag & Drop or Click to Upload Matrix"]),
                style={"width": "100%", "height": "60px", "lineHeight": "60px",
                       "borderWidth": "1px", "borderStyle": "dashed", "borderRadius": "5px",
                       "textAlign": "center", "margin-bottom": "10px"},
                multiple=False
            ),
            html.Label("Upload Metadata (.csv):"),
            dcc.Upload(
                id="upload-metadata",
                children=html.Div(["Drag & Drop or Click to Upload Metadata"]),
                style={"width": "100%", "height": "60px", "lineHeight": "60px",
                       "borderWidth": "1px", "borderStyle": "dashed", "borderRadius": "5px",
                       "textAlign": "center"},
                multiple=False
            ),
            html.Div(id="upload-status", style={"marginTop": "10px", "color": "green"})
        ]),

        html.H2("Select Normalization Methods"),

        html.Div([
            dcc.Checklist(
                id="method-checklist",
                options=methods_options,
                value=[],
                inputStyle={"margin-right": "8px", "margin-left": "16px"}
            )
        ]),

        html.Br(),
        dbc.Button("Run Selected Methods", id="run-button", color="primary", n_clicks=0),
        html.Br(),

        dcc.Textarea(id="output-area", style={"width": "100%", "height": 300})
    ])

def visualization_layout(uid):
    output_dir = os.path.join("output", uid)
    image_files = [
        ("Alignment Score", "alignment_score.png"),
        ("AUC Heatmap", "auc_heatmap.png"),
        ("Bray-Curtis Heatmap", "braycurtis.png"),
        ("PCA", "pca.png"),
        ("PCoA", "pcoa.png")
    ]
    tabs = []
    for label, filename in image_files:
        img_path = os.path.join(output_dir, filename)
        if os.path.exists(img_path):
            with open(img_path, "rb") as f:
                encoded = base64.b64encode(f.read()).decode()
                img_data = f"data:image/png;base64,{encoded}"
                tabs.append(dcc.Tab(label=label, children=[
                    html.Img(src=img_data, style={"width": "100%", "marginBottom": "20px"})
                ]))
    return html.Div([
        html.H2("Visualization Results"),
        html.Br(),
        dbc.Button("Download All Results", id="download-button", color="success", className="me-2"),
        dbc.Button("Back to Home", href="/", color="secondary"),
        html.Br(),
        dcc.Tabs(tabs),

    ])

# ------------------------------
# Callback 1: Session Init on Page Load
# ------------------------------
@app.callback(
    Output("session-uid", "data"),
    Input("page-url", "href"),
    prevent_initial_call=False
)
def initiate_session(_):
    uid = str(uuid.uuid4())
    output_dir = os.path.join("output", uid)
    os.makedirs(output_dir, exist_ok=True)

    r.set(uid, "active", ex=SESSION_TTL)
    r.set(f"{uid}:methods", "", ex=SESSION_TTL)

    if os.path.exists("output"):
        for folder_name in os.listdir("output"):
            folder_path = os.path.join("output", folder_name)
            if os.path.isdir(folder_path) and not r.exists(folder_name):
                shutil.rmtree(folder_path)
                print(f"Deleted expired session folder: {folder_name}")

    return uid

# ------------------------------
# Callback 2: File Upload
# ------------------------------
@app.callback(
    Output("upload-status", "children"),
    Input("upload-matrix", "contents"),
    Input("upload-matrix", "filename"),
    Input("upload-metadata", "contents"),
    Input("upload-metadata", "filename"),
    State("session-uid", "data"),
    prevent_initial_call=True
)
def handle_upload(matrix_content, matrix_name, metadata_content, metadata_name, uid):
    if not uid:
        return "Session not initialized."

    upload_path = os.path.join("output", uid)
    saved = []

    def save_file(content, filename):
        content_type, content_string = content.split(",")
        decoded = base64.b64decode(content_string)
        save_path = os.path.join(upload_path, filename)
        with open(save_path, "wb") as f:
            f.write(decoded)
        return filename

    if matrix_content and matrix_name:
        save_file(matrix_content, "raw.csv")
        saved.append("matrix")

    if metadata_content and metadata_name:
        save_file(metadata_content, "metadata.csv")
        saved.append("metadata")

    return f"Uploaded: {', '.join(saved)}"

# ------------------------------
# Callback 3: Run Normalization
# ------------------------------
@app.callback(
    Output("output-area", "value"),
    Output("page-url", "pathname"),
    Input("run-button", "n_clicks"),
    State("method-checklist", "value"),
    State("session-uid", "data"),
    prevent_initial_call=True
)
def run_r_script(n_clicks, selected_methods, uid):
    if not selected_methods:
        return "No methods selected.", dash.no_update
    if not uid:
        return "Session ID missing.", dash.no_update

    output_dir = os.path.join("output", uid)
    method_arg = ",".join(selected_methods)

    r.set(f"{uid}:methods", method_arg, ex=SESSION_TTL)
    for key in r.keys(f"{uid}*"):
        r.expire(key, SESSION_TTL)

    print(f"[{uid}] Starting normalization...")
    try:
        print(f"[{uid}] Running methods.R with methods: {method_arg}")
        subprocess.run([
            "Rscript", "methods.R", method_arg, output_dir
        ], capture_output=True, text=True, encoding="utf-8", errors="replace", check=True)
        print(f"[{uid}] Normalization completed.")

        plot_scripts = [
            "alignment_score.R",
            "auc_heatmap.R",
            "bray_curtis_heatmap.R",
            "pca.R",
            "pcoa.R"
        ]

        for script in plot_scripts:
            script_path = os.path.join("plots", script)
            print(f"[{uid}] Running {script}...")
            subprocess.run(["Rscript", script_path, output_dir],
                           capture_output=True, text=True, encoding="utf-8", errors="replace", check=True)
            print(f"[{uid}] Completed {script}.")

        print(f"[{uid}] All visualization scripts finished.")
        return f"UID: {uid}\n\nNormalization complete.", "/visualize"

    except subprocess.CalledProcessError as e:
        print(f"[{uid}] Error occurred.")
        print(f"STDOUT:\n{e.stdout}")
        print(f"STDERR:\n{e.stderr}")
        return f"Error:\nSTDOUT:\n{e.stdout}\n\nSTDERR:\n{e.stderr}", dash.no_update


# ------------------------------
# Callback 4: Page Routing
# ------------------------------
@app.callback(
    Output("page-content", "children"),
    Input("page-url", "pathname"),
    State("session-uid", "data")
)
def display_page(pathname, uid):
    if pathname == "/visualize" and uid:
        return visualization_layout(uid)
    return main_layout()

# ------------------------------
# Callback 5: Download
# ------------------------------
@app.callback(
    Output("download", "data"),
    Input("download-button", "n_clicks"),
    State("session-uid", "data"),
    prevent_initial_call=True
)
def download_session_folder(n_clicks, uid):
    if not uid:
        return dash.no_update

    folder_path = os.path.join("output", uid)
    zip_path = os.path.join(folder_path, "results.zip")  # Save inside user's folder
    shutil.make_archive(zip_path.replace(".zip", ""), 'zip', folder_path)

    return dcc.send_file(zip_path)

# ------------------------------
# Run App
# ------------------------------
if __name__ == '__main__':
    app.layout = html.Div([
        dcc.Location(id='page-url', refresh=False),
        dcc.Store(id="session-uid", storage_type="session"),
        dcc.Download(id="download"),
        html.Div(id='page-content', children=main_layout())  # Render the main layout initially
    ])
    app.run_server(debug=True, host="127.0.0.1", port=8051)