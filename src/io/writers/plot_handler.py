import math
from pathlib import Path
from src.utils.constants import (
    FINGERPRINT_RADIUS
)
import numpy as np
from src.utils.geometric.calculate_geometric_centers import calculate_center
import plotly.graph_objects as go
from plotly.subplots import make_subplots


PALETTE = [
    "#4DA8E8", "#E8783C", "#3EC48A", "#C44EC4",
    "#E8C43C", "#E84E84", "#4EC4C4", "#C4844E",
]


def parse_pdb_ca(path):
    """
    Return Ca atoms from a PDB file as [{"x":…, "y":…, "z":…}, …].
    """
    atoms = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM  ", "HETATM")) and line[12:16].strip() == "CA":
                try:
                    atoms.append({
                        "x": float(line[30:38]),
                        "y": float(line[38:46]),
                        "z": float(line[46:54]),
                    })
                except ValueError:
                    pass
    return atoms

def _sphere_lines(cx, cy, cz, radius, n=40):
    """
    Wireframe sphere returned as Scatter3d-compatible x/y/z lists (None-separated).
    """
    xs, ys, zs = [], [], []
    # latitude rings
    for lat in np.linspace(-math.pi / 2, math.pi / 2, 9):
        t = np.linspace(0, 2 * math.pi, n)
        r = radius * math.cos(lat)
        xs.extend(cx + r * np.cos(t));      xs.append(None)
        ys.extend(cy + r * np.sin(t));      ys.append(None)
        zs.extend([cz + radius * math.sin(lat)] * n); zs.append(None)
    # longitude lines
    for lon in np.linspace(0, 2 * math.pi, 12, endpoint=False):
        t = np.linspace(-math.pi / 2, math.pi / 2, n)
        xs.extend(cx + radius * np.cos(t) * math.cos(lon)); xs.append(None)
        ys.extend(cy + radius * np.cos(t) * math.sin(lon)); ys.append(None)
        zs.extend(cz + radius * np.sin(t));                 zs.append(None)
    return xs, ys, zs


def _hex_to_rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))


def build_figure(protein_atoms, pockets, pred_per_pocket, sphere_radius=3.0):
    """
    """
    n = len(pockets)

    ready = []
    for i, cl in enumerate(pockets):

        class_names = list(pred_per_pocket[cl].keys())
        #print(class_names)

        pts   = pockets[cl]
        preds_per_cl = pred_per_pocket[cl]
        probs = [preds_per_cl[a] for a in preds_per_cl.keys()]
        #print(probs)
        ready.append({
            "label":  f"Cluster {cl}",
            "points": pts,
            "center": calculate_center(np.array(pockets[cl])),
            "probs":  probs,
            "color":  PALETTE[i % len(PALETTE)],
        })

    scene_y_bottom = 0.26
    bar_top        = 0.22
    col_w          = 1.0 / n

    fig = go.Figure()

    if protein_atoms:
        fig.add_trace(go.Scatter3d(
            x=[a["x"] for a in protein_atoms],
            y=[a["y"] for a in protein_atoms],
            z=[a["z"] for a in protein_atoms],
            mode="lines",
            line=dict(color="rgba(60,80,140,0.85)", width=5),
            name="Backbone",
            hoverinfo="skip",
        ))

    for i, cl in enumerate(ready):
        color   = cl["color"]
        label   = cl["label"]
        pts     = cl["points"]
        cen     = cl["center"]
        probs   = cl["probs"]
        top_i   = probs.index(max(probs))
        top_lbl = class_names[top_i]
        r, g, b = _hex_to_rgb(color)

        # Pocket scatter points
        fig.add_trace(go.Scatter3d(
            x=[p[0] for p in pts],
            y=[p[1] for p in pts],
            z=[p[2] for p in pts],
            mode="markers",
            marker=dict(
                size=4,
                color=color,
                opacity=0.90,
                line=dict(width=0),
            ),
            name=label,
            legendgroup=f"cl{i}",
            legendgrouptitle=dict(text="") if i == 0 else {},
            hovertemplate=(
                f"<b>{label}</b><br>"
                f"Predicted: <b>{top_lbl}</b> ({probs[top_i]*100:.0f}%)<br>"
                "x=%{x:.1f}  y=%{y:.1f}  z=%{z:.1f}<extra></extra>"
            ),
        ))

        sx, sy, sz = _sphere_lines(cen[0], cen[1], cen[2], sphere_radius)
        fig.add_trace(go.Scatter3d(
            x=sx, y=sy, z=sz,
            mode="lines",
            line=dict(color=f"rgba({r},{g},{b},0.25)", width=1),
            name=f"{label} sphere",
            legendgroup=f"cl{i}",
            showlegend=False,
            hoverinfo="skip",
        ))

        # Centre marker
        fig.add_trace(go.Scatter3d(
            x=[cen[0]], y=[cen[1]], z=[cen[2]],
            mode="markers",
            marker=dict(
                size=7,
                color="white",
                symbol="diamond",
                line=dict(color=color, width=3),
            ),
            name=f"{label} centre",
            legendgroup=f"cl{i}",
            showlegend=False,
            hovertemplate=(
                f"<b>{label} — centre</b><br>"
                f"x={cen[0]:.1f}  y={cen[1]:.1f}  z={cen[2]:.1f}"
                "<extra></extra>"
            ),
        ))

    for i, cl in enumerate(ready):
        probs  = cl["probs"]
        top_i  = probs.index(max(probs))
        color  = cl["color"]
        r, g, b = _hex_to_rgb(color)
        ax     = "" if i == 0 else str(i + 1)

        bar_colors = [
            color if j == top_i else f"rgba({r},{g},{b},0.18)"
            for j in range(len(class_names))
        ]

        fig.add_trace(go.Bar(
            x=class_names,
            y=[round(p * 100, 1) for p in probs],
            marker_color=bar_colors,
            marker_line_width=0,
            showlegend=False,
            name=cl["label"],
            xaxis=f"x{ax}",
            yaxis=f"y{ax}",
            hovertemplate="%{x}: <b>%{y:.1f}%</b><extra></extra>",
        ))

        x_domain = [i * col_w + 0.015, (i + 1) * col_w - 0.015]

        fig.update_layout(**{
            f"xaxis{ax}": dict(
                domain=x_domain,
                anchor=f"y{ax}",
                tickfont=dict(size=9, color="#6b7280"),
                showgrid=False,
                zeroline=False,
                tickangle=-30,
            ),
            f"yaxis{ax}": dict(
                domain=[0.0, bar_top],
                anchor=f"x{ax}",
                range=[0, 108],
                tickfont=dict(size=9, color="#6b7280"),
                gridcolor="rgba(0,0,0,0.06)",
                griddash="dot",
                zeroline=False,
                showline=False,
                title=dict(
                    text="%" if i == 0 else "",
                    font=dict(size=10, color="#9ca3af"),
                ),
            ),
        })

    annotations = [
        dict(
            text=f"<b>{ready[i]['label']}</b>",
            x=(i + 0.5) * col_w,
            y=bar_top + 0.002,
            xref="paper", yref="paper",
            showarrow=False,
            font=dict(size=11, color=ready[i]["color"]),
            xanchor="center", yanchor="bottom",
        )
        for i in range(n)
    ]

    fig.update_layout(
        scene=dict(
            domain=dict(x=[0, 1], y=[scene_y_bottom, 1.0]),
            bgcolor="white",
            xaxis=dict(
                showgrid=True, gridcolor="rgba(0,0,0,0.06)",
                zeroline=False, showticklabels=False,
                showspikes=False, title="",
                backgroundcolor="rgba(200,212,235,0.85)",
            ),
            yaxis=dict(
                showgrid=True, gridcolor="rgba(0,0,0,0.10)",
                zeroline=False, showticklabels=False,
                showspikes=False, title="",
                backgroundcolor="rgba(200,212,235,0.85)",
            ),
            zaxis=dict(
                showgrid=True, gridcolor="rgba(0,0,0,0.10)",
                zeroline=False, showticklabels=False,
                showspikes=False, title="",
                backgroundcolor="rgba(200,212,235,0.85)",
            ),
        ),
        paper_bgcolor="white",
        plot_bgcolor="white",
        font=dict(
            family="'Inter', 'Segoe UI', system-ui, sans-serif",
            color="#374151",
            size=11,
        ),
        height=860,
        margin=dict(l=10, r=10, t=60, b=10),
        legend=dict(
            x=1.01, y=0.98,
            xanchor="left",
            bgcolor="rgba(255,255,255,0.90)",
            bordercolor="rgba(0,0,0,0.08)",
            borderwidth=1,
            font=dict(size=12, color="#374151"),
            tracegroupgap=4,
        ),
        annotations=annotations,
    )

    return fig

def generate_html(
    protein_atoms,
    pockets,
    pred_per_pocket,
    sphere_radius=3.0,
    output_path="pocket_viewer.html",
    title="Protein Pocket Viewer",
):
    """
    Build and write a fully self-contained HTML file (no internet needed).
    """
    fig = build_figure(protein_atoms, pockets, pred_per_pocket, sphere_radius)
    fig.update_layout(
        title=dict(
            text=title,
            font=dict(size=18, color="#111827", family="'Inter','Segoe UI',sans-serif"),
            x=0.02,
            y=0.99,
        )
    )

    out = Path(output_path)
    fig.write_html(
        str(out),
        include_plotlyjs=True,   # self-contained – no CDN required
        full_html=True,
        config={
            "displaylogo": False,
            "scrollZoom": True,
            "modeBarButtonsToRemove": ["toImage"],
        },
    )
    return str(out.resolve())

def plot_protein(structure_path : Path, pockets : dict, pred_per_pocket : dict, output: Path, title: str = "Protein Plot"):
    """
    """
    protein_atoms = parse_pdb_ca(structure_path)
    output_path = output / "protein_plot.html"
    path = generate_html(
        protein_atoms = protein_atoms,
        pockets      = pockets,
        pred_per_pocket   = pred_per_pocket,
        sphere_radius = FINGERPRINT_RADIUS,
        output_path   = output_path,
        title         = title,
    )
