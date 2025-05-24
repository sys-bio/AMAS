"""
This script provides a web interface for annotating SBML models using the AMAS (Automated Model Annotation System) framework.

Usage:
0. Navigate to the AMAS-v2/UI folder.

1. Run the application:
   python app.py

2. Open the web application in browser:
   Navigate to http://127.0.0.1:8050/ (default address) in web browser.

3. Workflow:
   a. Upload an SBML file for annotation. Only files in valid SBML format are supported.
   b. Choose whether to annotate all elements or only specific types (species, reactions, or genes).
   c. Adjust match score criteria, cutoff, minimum name length, and optimization settings as needed.
   d. Click "Recommend Annotations" to generate recommendations based on the selected parameters.
   e. Use the annotation table to make interactive updates to the recommendations.
   f. Download the annotated and updated SBML model to your local machine.

4. Stopping the Server:
   Press `Ctrl+C` in the terminal to stop the application.

"""

import dash
from dash import dcc, html, dash_table
from dash.dependencies import Input, Output, State
import pandas as pd
import base64
import subprocess
import os
import threading
import flask
server = flask.Flask(__name__)
app = dash.Dash(__name__, server=server)

# Paths
csv_output_path = "recommendations.csv"
amas_recommendation_path = "../AMAS/recommend_annotation.py"
amas_species_path = "../AMAS/recommend_species.py"
amas_reactions_path = "../AMAS/recommend_reactions.py"
amas_genes_path = "../AMAS/recommend_genes.py"
amas_update_path = "../AMAS/update_annotation.py"
current_model_path = "current_model.xml"
download_file_path = "new_model.xml"

# Global variable to hold AMAS process information
process_info = {"process": None, "output": ""}

# Function to read output in a separate thread
def read_process_output(process, params_used):
    process_info["output"] = params_used  # Initialize with parameters
    while True:
        line = process.stdout.readline()
        if not line and process.poll() is not None:
            break
        if line:
            process_info["output"] += line

# Global variable to hold update process information
update_process_info = {"process": None, "output": ""}

# Function to read output in a separate thread
def read_update_process_output(update_process, params_used):
    """Reads output from the update process while it's running."""
    update_process_info["output"] = params_used  # Initialize with parameters

    for line in iter(update_process.stdout.readline, ""):  # Read line-by-line
        if not line:
            break
        update_process_info["output"] += line
        print(line, end="")  # Debugging: Print output in real-time

    stdout, stderr = update_process.communicate()  # Ensure all remaining output is captured
    if stdout:
        update_process_info["output"] += stdout
    if stderr:
        update_process_info["output"] += f"\nError: {stderr}"

    print("Update process completed.")  # Debugging: Ensure process finishes
    print(update_process_info["output"])


# Define the layout of the app
app.layout = html.Div(
    style={"display": "flex", "flexDirection": "row", "width": "100%", "height": "95vh", "overflow": "hidden"},
    children=[
        # Left Section: Menu, Parameters, Messages
        html.Div(
            style={"width": "25%", "padding": "10px", "boxSizing": "border-box", "borderRight": "1px solid #ccc"},
            children=[
                html.H1("AMAS: Automated Model Annotation System", style={"fontSize": "20px", "textAlign": "center"}),

                # Menu Buttons
                html.Div([
                    dcc.Upload(id='upload-sbml', children=html.Button(html.Strong("Upload SBML Model"), style={"height": "50px", "fontSize": "12px"},
                    title = "Upload a model in SBML format."), multiple=False),
                    html.Button(html.Strong("Recommend Annotations"), id="annotation-button", n_clicks=0, style={"height": "50px", "fontSize": "12px"},
                    title = "Generate a table that contains the recommended annotations."),
                    html.Button(html.Strong("Update Annotations"), id="update-button", n_clicks=0, style={"height": "50px", "fontSize": "12px"},
                    title = "Update model annotation using instructions in table."),
                    html.Button(html.Strong("Download Annotated Model"), id="download-button", n_clicks=0, style={"height": "50px", "fontSize": "12px"},
                    title = "Download the annotated model to local."),
                    dcc.Download(id="download-component")
                ], style={"display": "flex", "gap": "10px", "marginBottom": "20px"}),

                # Annotation Type Selection
                html.Label(html.Strong("What to Annotate:"), 
                title = "To annotate all elements, or species/reactions/genes only."),
                dcc.RadioItems(
                    id="annotation-type",
                    options=[
                        {"label": "All", "value": "all"},
                        {"label": "Species", "value": "species"},
                        {"label": "Reactions", "value": "reactions"},
                        {"label": "Genes", "value": "genes"}
                    ],
                    value="all",
                    inline=True
                ),
                html.Br(),

                # Tax_id for gene annotation
                html.Label(html.Strong("Taxonomy:"), 
                title="Choose a taxonomy or enter its taxonomy ID for gene annotation."),
                dcc.Dropdown(
                    id="taxonomy-dropdown",
                    options=[
                        {"label": "Human", "value": "9606"},
                        {"label": "E coli", "value": "511145"},
                        {"label": "Mouse", "value": "10090"},
                    ],
                    value="9606",
                    placeholder="Select taxonomy for gene annotation",
                    style={"marginBottom": "10px", "width": "100%"},
                ),
                html.Label("Or enter its Taxonomy ID:"),
                dcc.Input(
                    id="taxonomy-id-input",
                    type="text",
                    placeholder="Enter taxonomy ID",
                    style={"marginBottom": "20px", "width": "30%"},
                ),
                html.Br(),                

                # Cross-reference fields
                html.Label(html.Strong("Cross-reference fields:"),
                title = "Select or enter cross-reference fields for gene annotation."),
                dcc.Dropdown(
                    id="crossref-dropdown",
                    options=[
                        {"label": "UniProt", "value": "uniprot"},
                        {"label": "HGNC", "value": "HGNC"},
                        {"label": "Ensembl", "value": "ensembl"},
                        {"label": "Ensembl.gene", "value": "ensembl.gene"},
                        {"label": "Ensembl.protein", "value": "ensembl.protein"},
                        {"label": "Refseq", "value": "refseq"},
                        {"label": "Refseq.rna", "value": "refseq.rna"},
                        {"label": "Refseq.protein", "value": "refseq.protein"},
                    ],
                    multi=True,
                    placeholder="Select or enter cross-reference fields",
                ),
                dcc.Input(
                    id="crossref-custom",
                    type="text",
                    placeholder="Enter fields ('ec,MIM')",
                ),
                html.Br(),
                html.Br(),

                # Parameters Section
                html.Div([
                    html.Label(html.Strong("Match Score Selection Criteria:"), 
                    title = "Decide how to recommend candidate annotations based on the match scores."),
                    dcc.RadioItems(
                        id='mssc-radio',
                        options=[
                            {'label': 'Top', 'value': 'top'},
                            {'label': 'Above', 'value': 'above'}
                        ],
                        value='top',  # Default value
                        inline=True
                    ),
                    html.Br(),
                    html.Label(html.Strong("Match Score Cutoff:"),
                    title = "The cutoff for match scores (applies to both top and above MSSC)"),
                    dcc.Slider(
                        id='cutoff-slider',
                        min=0.0,
                        max=1.0,
                        step=0.01,
                        value=0.0,  # Default value
                        marks={i / 10: f"{i / 10:.1f}" for i in range(11)}
                    ),
                    html.Br(),
                    html.Label(html.Strong("Minimum Length of Name:"),
                    title = "Only recommend annotations for those with longer displayed names."),
                    html.Br(),
                    dcc.Input(
                        id='min-len-input',
                        type='number',
                        value=0,  # Default value
                        min=0,
                        step=1,
                        style={"width": "10%"}
                    ),
                    html.Br(),
                    html.Br(),
                    html.Label("Optimize Predictions:",
                    title = "Whether to compare once-predicted annotations of species and reactions and iteratively updates them."),
                    dcc.RadioItems(
                        id='optimize-radio',
                        options=[
                            {'label': 'Yes', 'value': 'yes'},
                            {'label': 'No', 'value': 'no'}
                        ],
                        value='no',  # Default value
                        inline=True
                    ),
                    html.Br(),
                    html.Label("Convert annotation in <isDescribedBy> to <is>:",
                            title="Choose whether to convert annotations from <isDescribedBy> to <is>."),
                    dcc.RadioItems(
                        id="convert-annotation-radio",
                        options=[
                            {"label": "Yes", "value": "yes"},
                            {"label": "No", "value": "no"},
                            {"label": "All", "value": "all"}
                        ],
                        value="no",  # Default value
                        inline=True
                    ),
                ], style={"marginBottom": "20px"}),

                # Processing Messages
                html.Div(id="amas-output", style={
                    "whiteSpace": "pre-wrap", "border": "1px solid #ccc", "padding": "10px",
                    "maxHeight": "700px", "overflowY": "scroll"
                }, children="Please upload an SBML model to begin."),
                html.Br(),
                # Processing Messages for update
                html.Div(id="update-output", style={
                    "whiteSpace": "pre-wrap", "border": "1px solid #ccc", "padding": "10px",
                    "maxHeight": "700px", "overflowY": "scroll"
                }, children="")
            ]
        ),

        # Middle Section: SBML Model Display
        html.Div(
            id="sbml-container",
            style={"width": "30%", "padding": "10px", "boxSizing": "border-box", "borderRight": "1px solid #ccc"},
            children=[
                html.H2("SBML Model", style={"fontSize": "16px", "textAlign": "center"}),
                html.Div(id='sbml-content', style={
                    "whiteSpace": "pre-wrap", "border": "1px solid #ccc", "padding": "10px",
                    "maxHeight": "900px", "overflowY": "scroll", "fontSize": "14px"
                }, children="Please upload an SBML model to begin.")
            ]
        ),

        # Right Section: Table
        html.Div(
            id="table-container",
            style={"width": "45%", "padding": "10px", "boxSizing": "border-box"},
            children=[
                html.H2("Annotation Table", style={"fontSize": "18px", "textAlign": "center"}),
                html.H2("Please select instructions for updating the model in the UPDATE ANNOTATION column.", style={"fontSize": "14px", "textAlign": "left"}),
                html.H2("Shaded: Existing annotations; White: Recommended annotations.", style={"fontSize": "13px", "textAlign": "left"}),
                # Add Row Button
                html.Button("Add Row Below", id="add-row-button", n_clicks=0, style={"marginBottom": "0px", "fontSize": "11px", "width": "20%"}),

                dash_table.DataTable(
                    id="csv-table",
                    style_table={"maxHeight": "900px", "overflowY": "scroll", "marginTop": "20px"},
                    style_data={'whiteSpace': 'normal', 'height': 'auto', 'lineHeight': '15px'},
                    editable=True,
                    row_selectable="single",
                    filter_action="native", sort_action="native", sort_mode="multi",
                    columns=[
                        {"name": "Type", "id": "type"},
                        {"name": "ID", "id": "id"},
                        {"name": "Display Name", "id": "display name"},
                        {"name": "Annotation", "id": "annotation"},
                        {"name": "Annotation Label", "id": "annotation label"},
                        {"name": "Match Score", "id": "match score"},
                        {"name": "Existing", "id": "existing"}, 
                        {"name": "UPDATE ANNOTATION", "id": "UPDATE ANNOTATION", "presentation": "dropdown"}
                    ],
                    dropdown_conditional=[
                        {
                            "if": {"column_id": "UPDATE ANNOTATION", "filter_query": "{existing} = 1"},
                            "options": [
                                {"label": "keep", "value": "keep"},
                                {"label": "delete", "value": "delete"}
                            ],
                        },
                        {
                            "if": {"column_id": "UPDATE ANNOTATION", "filter_query": "{existing} = 0"},
                            "options": [
                                {"label": "ignore", "value": "ignore"},
                                {"label": "add", "value": "add"}
                            ],
                        },
                    ],
                    style_cell_conditional=[
                        {
                            'if': {
                                'column_id': 'existing'
                            },
                            'display': 'none'  # Hide the 'Existing' column by default
                        }
                    ],
                    style_data_conditional=[
                        # Highlight Annotation column with bold border
                        {
                            "if": {"column_id": "annotation"},
                            "border": "2px solid black",  # Bold border
                            # "backgroundColor": "rgba(230, 230, 255, 0.5)"  # Light blue background 
                        },
                        # Highlight Annotation Label column with bold border
                        {
                            "if": {"column_id": "annotation label"},
                            "border": "2px solid black",  # Bold border
                            # "backgroundColor": "rgba(230, 230, 255, 0.5)"  # Light blue background 
                        },
                        # Bold font for the UPDATE ANNOTATION column
                        {
                            "if": {"column_id": "UPDATE ANNOTATION"},
                            "fontWeight": "bold",  # Make text bold
                        },
                        # Highlight Existing column with light grey background
                        {
                            "if": {
                                "filter_query": "{existing} = 1"
                            },
                            "backgroundColor": "rgba(211, 211, 211, 0.5)",  # Light grey for entire row
                            "color": "black",
                        },
                    ]
                )
            ]
        ),

        # Interval for real-time AMAS output
        dcc.Interval(id="interval-component", interval=1000, n_intervals=0)
    ]
)

# Callback to manage workflow
@app.callback(
    [Output("amas-output", "children"),
     Output("sbml-content", "children"),
     Output("csv-table", "data"),
     Output("download-component", "data"),
     Output("interval-component", "disabled"),
     Output("update-output", "children")],
    [Input("upload-sbml", "contents"),
     Input("annotation-button", "n_clicks"),
     Input("update-button", "n_clicks"),
     Input("download-button", "n_clicks"),
     Input("interval-component", "n_intervals"),
     Input("add-row-button", "n_clicks")],
    [State("annotation-type", "value"),
     State("taxonomy-dropdown", "value"), 
     State("taxonomy-id-input", "value"),
     State("mssc-radio", "value"),
     State("cutoff-slider", "value"),
     State("optimize-radio", "value"),
     State("csv-table", "data"),
     State("csv-table", "selected_rows"),
     State("sbml-content", "children"),
     State("min-len-input","value"),
     State("convert-annotation-radio", "value"),
     State("crossref-dropdown", "value"),
     State("crossref-custom", "value")]
)
def manage_workflow(contents, annotation_clicks, update_clicks, download_clicks, n_intervals, add_row_clicks,
                    annotation_type, tax_dropdown, tax_input, mssc, cutoff, optimize, table_data, selected_rows, sbml_text,min_len,convert_annotation,crossref_dropdown,crossref_custom):
    ctx = dash.callback_context
    if not ctx.triggered:
        return "Please upload an SBML model to begin.", sbml_text, [], dash.no_update, True, ""

    button_id = ctx.triggered[0]["prop_id"].split(".")[0]

    # Handle SBML Upload
    if button_id == "upload-sbml" and contents:
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        with open(current_model_path, "wb") as f:
            f.write(decoded)
        model_content = decoded.decode("utf-8")
        return "Model uploaded successfully.\n Please adjust the parameters to generate recommendations.", model_content, [], dash.no_update, True, ""

    # Handle Annotation
    if button_id == "annotation-button" and annotation_clicks > 0:
        script_path = {
            "all": amas_recommendation_path,
            "species": amas_species_path,
            "reactions": amas_reactions_path,
            "genes": amas_genes_path
        }[annotation_type]

        if annotation_type in ["genes", "all"]:
            tax_id = tax_dropdown if tax_dropdown else tax_input
            if not tax_id:
                return "Please select a taxonomy or enter a taxonomy ID for gene annotation.", sbml_text, [], dash.no_update, True, ""

        else:
            tax_id = None
        
        command = [
            "python", script_path, current_model_path,
            "--mssc", mssc, "--cutoff", str(cutoff),
            "--outfile", csv_output_path
        ]

        params_used = f"Running AMAS {annotation_type} annotation with parameters:\n  --mssc: {mssc}\n  --cutoff: {cutoff}\n"
        if annotation_type in ["all", "genes"]:
            params_used += f"  --tax: {tax_id}\n"
            command += ["--tax", tax_id]

        if annotation_type == "all":
            params_used += f"--optimize: {optimize}\n"
            command += ["--save", "csv", "--optimize", optimize]
        else: 
            params_used += f" --min_len: {min_len}\n"
            command += ["--min_len", str(min_len)]

        print(command)
        process_info["process"] = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
        )
        threading.Thread(target=read_process_output, args=(process_info["process"], params_used)).start()
        return params_used, sbml_text, table_data, dash.no_update, False, "Please generate recommendations first."

    # Handle Add Row Below
    if button_id == "add-row-button" and add_row_clicks > 0 and selected_rows:
        selected_index = selected_rows[0]
        selected_row = table_data[selected_index]
        # Create a duplicate row with modified fields
        new_row = selected_row.copy()
        new_row["annotation"] = ""  
        new_row["annotation label"] = "" 
        new_row["existing"] = 0
        new_row["UPDATE ANNOTATION"] = "add" 
        table_data.insert(selected_index + 1, new_row)  # Insert new row below selected row
        return dash.no_update, sbml_text, table_data, dash.no_update, True, "Row added. \n Please update the model after confirming the instructions in the table."

    # Handle Updates
    if button_id == "update-button" and update_clicks > 0:
        pd.DataFrame(table_data).to_csv(csv_output_path, index=False)
        fields = crossref_dropdown if crossref_dropdown else []
        if crossref_custom:
            fields.append(crossref_custom)
        fields_str = ','.join(fields)
        update_process_info["output"] = f"Updating model with instructions in table...\n\n"
        command_update = [
            "python", amas_update_path, current_model_path, csv_output_path, download_file_path, "--convert", convert_annotation, "--fields", fields_str
        ]
        print(command_update)
        update_process_info["process"] = subprocess.Popen(
            command_update, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1
        )

        threading.Thread(target=read_update_process_output, args=(update_process_info["process"], update_process_info["output"]),daemon=True).start()
        with open(download_file_path, "r") as f:
            updated_model_content = f.read()
        return "", updated_model_content, dash.no_update, dash.no_update, False, update_process_info["output"]


    # Real-Time Output
    if button_id == "interval-component" and process_info["process"]:
        if process_info["process"].poll() is not None:
            df = pd.read_csv(csv_output_path)
            df["UPDATE ANNOTATION"] = df.apply(
            lambda row: "add" if row["existing"] == 0 and (
            str(row["annotation label"]).lower() == str(row["display name"]).lower() or 
            str(row["annotation label"]).lower() == str(row["id"]).lower()
            )
            else "keep" if row["existing"] == 1
            else "ignore",
            axis=1
            )
            process_info["process"] = None
            return process_info["output"], sbml_text, df.to_dict("records"), dash.no_update, True, "Please make instructions in the table, and then update the model."
        else:
            return dash.no_update, sbml_text, [], dash.no_update, False, "Please generate recommendations first."

    if button_id == "interval-component" and update_process_info["process"]:
        if update_process_info["process"].poll() is not None:
            with open(download_file_path, "r") as f:
                updated_model_content = f.read()
            update_process_info["process"] = None
            print(update_process_info["output"])
            return dash.no_update, updated_model_content, table_data, dash.no_update, False, update_process_info["output"]
        else:
            return dash.no_update, sbml_text, table_data, dash.no_update, False, dash.no_update

    # Handle Download
    if button_id == "download-button" and download_clicks > 0:
        return "Upload an new SBML model to annotate.", sbml_text, table_data, dcc.send_file(download_file_path), True, "Annotated model downloaded."

    return dash.no_update, sbml_text, table_data, dash.no_update, True, dash.no_update


if __name__ == "__main__":
    app.run_server(host='127.0.0.1', port=8050, debug=True)