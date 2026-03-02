from pathlib import Path
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

def extract_ca_coordinates_from_pdb(pdb_path: str):
    """
    Extract CA atom coordinates from a PDB file.

    Returns:
        list of (x, y, z)
    """
    coords = []

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                if atom_name == "CA":
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))

    return coords

def plot_with_protein_from_pdb(
    pdb_path,
    cofactor_coords,
    labels,
    names=None,
    title="Cofactors + Protein",
    cofactor_size=6,
    protein_size=2,
    protein_opacity=0.15,
    out_html=None,
):
    protein_coords = extract_ca_coordinates_from_pdb(pdb_path)

    # Cofactor dataframe
    df = pd.DataFrame(cofactor_coords, columns=["x", "y", "z"])
    df["label"] = labels
    if names is not None:
        df["name"] = names

    hover_cols = ["label"] + (["name"] if names is not None else [])

    fig = px.scatter_3d(
        df,
        x="x", y="y", z="z",
        color="label",
        hover_data=hover_cols,
        title=title,
    )
    fig.update_traces(marker=dict(size=cofactor_size))

    # Protein cloud
    if protein_coords:
        p = pd.DataFrame(protein_coords, columns=["x", "y", "z"])
        fig.add_trace(
            go.Scatter3d(
                x=p["x"],
                y=p["y"],
                z=p["z"],
                mode="markers",
                name="protein (CA)",
                marker=dict(
                    size=protein_size,
                    opacity=protein_opacity,
                ),
                hoverinfo="skip",
            )
        )

    fig.update_layout(
        margin=dict(l=0, r=0, t=40, b=0),
        scene=dict(xaxis_title="X", yaxis_title="Y", zaxis_title="Z"),
    )

    if out_html is not None:
        Path(out_html).parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(out_html, include_plotlyjs="cdn")

    return fig