"""
app.py — SARS-CoV-2 Mutation Tracker
"""

import json
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

st.set_page_config(
    page_title="sars-cov-2 / mutations",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;500;600&family=Inter:wght@300;400;500;600&display=swap');

html, body, [class*="css"] { font-family: 'Inter', sans-serif; }
.stApp { background-color: #0d1117; color: #e6edf3; }

/* sidebar */
[data-testid="stSidebar"] { background-color: #090d13 !important; border-right: 1px solid #21262d; }
[data-testid="stSidebar"] * { color: #484f58 !important; font-family: 'JetBrains Mono', monospace !important; font-size: 11px !important; }
[data-testid="stSidebar"] .stSelectbox label,
[data-testid="stSidebar"] .stDateInput label { font-size: 9px !important; text-transform: uppercase; letter-spacing: 0.12em; color: #3d444d !important; }
[data-testid="stSidebar"] .stSelectbox > div > div { background: #0d1117 !important; border: 1px solid #21262d !important; border-radius: 4px !important; color: #484f58 !important; font-size: 10px !important; }
[data-testid="stSidebar"] .stRadio label { font-size: 10px !important; }

/* tabs → replaced by sidebar nav, but keep fallback styled */
.stTabs [data-baseweb="tab-list"] { background: #090d13; border-bottom: 1px solid #21262d; }
.stTabs [data-baseweb="tab"] { font-family: 'JetBrains Mono', monospace; font-size: 11px; color: #484f58; padding: 8px 16px; border-bottom: 2px solid transparent; background: transparent; }
.stTabs [aria-selected="true"] { color: #58a6ff !important; border-bottom-color: #1f6feb !important; background: transparent !important; }
.stTabs [data-baseweb="tab"]:hover { color: #8b949e !important; background: #161b22 !important; }

/* metrics */
[data-testid="stMetric"] { background: #161b22; border: 1px solid #21262d; border-radius: 6px; padding: 12px 16px; }
[data-testid="stMetricLabel"] { font-family: 'JetBrains Mono', monospace !important; font-size: 9px !important; text-transform: uppercase; letter-spacing: 0.12em; color: #3d444d !important; }
[data-testid="stMetricValue"] { font-family: 'JetBrains Mono', monospace !important; font-size: 22px !important; font-weight: 600 !important; color: #e6edf3 !important; }

/* dataframes */
[data-testid="stDataFrame"] { border: 1px solid #21262d !important; border-radius: 6px; }
[data-testid="stDataFrame"] * { font-family: 'JetBrains Mono', monospace !important; font-size: 11px !important; }

/* inputs */
.stTextInput > div > div > input { background: #161b22 !important; border: 1px solid #21262d !important; color: #8b949e !important; border-radius: 4px !important; font-family: 'JetBrains Mono', monospace !important; font-size: 11px !important; }
.stMultiSelect > div > div { background: #161b22 !important; border: 1px solid #21262d !important; }
.stRadio > label { font-family: 'JetBrains Mono', monospace !important; font-size: 11px !important; color: #484f58 !important; }

/* buttons */
.stButton > button, .stDownloadButton > button { background: #161b22 !important; color: #58a6ff !important; border: 1px solid #21262d !important; border-radius: 4px !important; font-family: 'JetBrains Mono', monospace !important; font-size: 11px !important; }
.stButton > button:hover, .stDownloadButton > button:hover { border-color: #1f6feb !important; background: #0d1e30 !important; }

/* alerts */
.stAlert { background: #161b22 !important; border: 1px solid #21262d !important; border-radius: 4px !important; font-family: 'JetBrains Mono', monospace !important; font-size: 11px !important; }

/* dividers */
hr { border-color: #21262d !important; }

/* headings */
h1, h2, h3 { font-family: 'JetBrains Mono', monospace !important; color: #e6edf3 !important; font-weight: 500 !important; }

/* spinner */
.stSpinner > div { border-top-color: #58a6ff !important; }

/* scrollbar */
::-webkit-scrollbar { width: 4px; height: 4px; }
::-webkit-scrollbar-track { background: #0d1117; }
::-webkit-scrollbar-thumb { background: #21262d; border-radius: 2px; }
</style>
""", unsafe_allow_html=True)

sys.path.insert(0, str(Path(__file__).parent))
from config import OUTPUT_DIR, KEY_MUTATIONS, VARIANT_DEFINITIONS
from visualization.plots import (
    timeline_plot, geographic_plotly, spike_protein_map,
    spike_rbd_detail, sweep_speed_chart, WHO_COLOURS,
)

PLOTLY_DARK = dict(
    paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
    font=dict(color="#484f58", family="JetBrains Mono, monospace", size=10),
    xaxis=dict(gridcolor="#161b22", linecolor="#21262d", tickcolor="#3d444d", tickfont=dict(size=9)),
    yaxis=dict(gridcolor="#161b22", linecolor="#21262d", tickcolor="#3d444d", tickfont=dict(size=9)),
    legend=dict(bgcolor="#161b22", bordercolor="#21262d", borderwidth=1, font=dict(size=9)),
    margin=dict(l=40, r=16, t=32, b=32),
)

def dark_fig(fig):
    fig.update_layout(**PLOTLY_DARK)
    return fig

def mono(text, color="#3d444d", size="9px", extra=""):
    return (f"<span style='font-family:JetBrains Mono,monospace;"
            f"font-size:{size};color:{color};{extra}'>{text}</span>")

def section_header(label, sub=""):
    s = f" <span style='color:#3d444d;font-size:10px'>— {sub}</span>" if sub else ""
    st.markdown(
        f"<div style='margin:16px 0 10px'>"
        f"<span style='font-family:JetBrains Mono,monospace;font-size:9px;"
        f"color:#3d444d;text-transform:uppercase;letter-spacing:0.14em'>{label}</span>{s}"
        f"</div>", unsafe_allow_html=True)

def event_card(date, name, color, desc):
    return (f"<div style='display:flex;gap:10px;padding:7px 0;"
            f"border-bottom:1px solid #21262d;align-items:flex-start'>"
            f"<span style='font-family:JetBrains Mono,monospace;font-size:9px;"
            f"color:#3d444d;min-width:46px;padding-top:2px'>{date}</span>"
            f"<div style='width:2px;border-radius:1px;background:{color};"
            f"align-self:stretch;flex-shrink:0;margin-top:3px'></div>"
            f"<div><div style='font-size:12px;font-weight:500;color:#c9d1d9'>{name}</div>"
            f"<div style='font-size:10px;color:#484f58;margin-top:1px'>{desc}</div>"
            f"</div></div>")

@st.cache_data(show_spinner=False)
def load_all_data():
    mut_table_path = OUTPUT_DIR / "mutation_table.parquet"
    if not mut_table_path.exists():
        with st.spinner("initialising pipeline..."):
            from pipeline.fetch_data import _generate_sample_metadata
            from pipeline.alignment import build_mutation_table
            from analysis.frequency_analysis import run_all_analyses
            meta_path = _generate_sample_metadata(n=5000)
            df = build_mutation_table(meta_path)
            run_all_analyses(df)
    data = {}
    def _load(name):
        p = OUTPUT_DIR / f"{name}.parquet"
        return pd.read_parquet(p) if p.exists() else pd.DataFrame()
    data["mutations"]     = _load("mutation_table")
    data["variant_freq"]  = _load("variant_frequency")
    data["mutation_freq"] = _load("mutation_frequency_time")
    data["convergent"]    = _load("convergent_evolution")
    data["hotspots"]      = _load("hotspots")
    data["sweeps"]        = _load("sweep_speed")
    data["geo"]           = _load("geographic_spread")
    for key in ["variant_freq", "mutation_freq"]:
        df = data[key]
        if not df.empty and "year_month" in df.columns:
            if not hasattr(df["year_month"].dtype, "freq"):
                df["year_month"] = pd.PeriodIndex(df["year_month"], freq="M")
            data[key] = df
    return data

def render_sidebar(data):
    st.sidebar.markdown(
        "<div style='padding:16px 0 12px'>"
        "<span style='font-family:JetBrains Mono,monospace;font-size:14px;"
        "font-weight:600;color:#e6edf3'>sars-cov-2</span>"
        "<span style='font-family:JetBrains Mono,monospace;font-size:14px;"
        "font-weight:300;color:#1f6feb'> / mutations</span>"
        "</div>", unsafe_allow_html=True)

    st.sidebar.markdown(
        "<div style='font-family:JetBrains Mono,monospace;font-size:9px;"
        "color:#3d444d;text-transform:uppercase;letter-spacing:0.14em;"
        "margin-bottom:4px'>// filters</div>", unsafe_allow_html=True)

    mutations_df = data["mutations"]
    date_range = None
    if not mutations_df.empty and "date" in mutations_df.columns:
        mutations_df["date"] = pd.to_datetime(mutations_df["date"])
        min_date = mutations_df["date"].min().date()
        max_date = mutations_df["date"].max().date()
        date_range = st.sidebar.date_input("date range",
            value=(min_date, max_date), min_value=min_date, max_value=max_date)

    region = "All"
    if not mutations_df.empty and "region" in mutations_df.columns:
        regions = ["All"] + sorted(mutations_df["region"].dropna().unique().tolist())
        region = st.sidebar.selectbox("region", regions)

    variant = "All"
    if not mutations_df.empty and "who_label" in mutations_df.columns:
        variants = ["All"] + sorted(mutations_df["who_label"].dropna().unique().tolist())
        variant = st.sidebar.selectbox("variant", variants)

    st.sidebar.divider()

    source_label = "sample data · 5k seqs"
    real_path = OUTPUT_DIR / "nextstrain_metadata.parquet"
    if real_path.exists():
        try:
            df_size = len(pd.read_parquet(real_path))
            if df_size > 5000:
                source_label = f"nextstrain · {df_size:,} seqs"
        except Exception:
            pass

    st.sidebar.markdown(
        f"<div style='font-family:JetBrains Mono,monospace;font-size:9px;"
        f"color:#3fb950;background:#0d2818;border:1px solid #1a4427;"
        f"border-radius:4px;padding:6px 10px;margin-top:4px'>● {source_label}</div>",
        unsafe_allow_html=True)

    return date_range, region, variant

def apply_filters(data, date_range, region, variant):
    df = data["mutations"].copy()
    if df.empty:
        return data
    df["date"] = pd.to_datetime(df["date"])
    if date_range and len(date_range) == 2:
        df = df[(df["date"] >= pd.Timestamp(date_range[0])) &
                (df["date"] <= pd.Timestamp(date_range[1]))]
    if region != "All" and "region" in df.columns:
        df = df[df["region"] == region]
    if variant != "All" and "who_label" in df.columns:
        df = df[df["who_label"] == variant]
    filtered = dict(data)
    filtered["mutations"] = df
    return filtered

def render_metrics(df):
    unique = set()
    for m in df["spike_mutations"].dropna():
        if isinstance(m, str):
            unique.update(json.loads(m))
    c1, c2, c3, c4, c5 = st.columns(5)
    with c1: st.metric("sequences", f"{len(df):,}")
    with c2: st.metric("mutations", f"{len(unique):,}")
    with c3: st.metric("variants", df["who_label"].nunique() if "who_label" in df.columns else 0)
    with c4: st.metric("countries", df["country"].nunique() if "country" in df.columns else 0)
    with c5:
        avg = df["rbd_mut_count"].mean() if "rbd_mut_count" in df.columns else 0
        st.metric("avg rbd muts", f"{avg:.1f}")

def tab_overview(data, filtered):
    render_metrics(filtered["mutations"])
    st.divider()

    col_left, col_right = st.columns([3, 2], gap="medium")

    with col_left:
        section_header("variant dominance", "monthly frequency · 2020–2026")
        vf = data["variant_freq"]
        if not vf.empty:
            fig = dark_fig(timeline_plot(vf, mode="variant"))
            fig.update_layout(
                title="", height=220,
                xaxis_title="", yaxis_title="frequency (%)",
            )
            st.plotly_chart(fig, width='stretch', key="chart_1")

    with col_right:
        section_header("evolutionary log")
        events = [
            ("2020-02", "D614G",   "#58a6ff", "universal backbone · sweeps globally"),
            ("2020-09", "Alpha",   "#1d9e75",  "N501Y + P681H · 50% advantage"),
            ("2020-10", "E484K",   "#f85149",  "convergent immune escape · 3 lineages"),
            ("2021-04", "Delta",   "#ef9f27",  "L452R + T478K · furin enhanced"),
            ("2021-11", "Omicron", "#d4537e",  "30+ spike muts · immune evasion"),
            ("2022-10", "XBB",     "#7f77dd",  "recombinant · F486P ACE2 escape"),
            ("2023-08", "JN.1",    "#9e6ac7",  "BA.2.86 descendant · L455S shift"),
        ]
        log_html = "<div style='border-top:1px solid #21262d'>"
        for date, name, color, desc in events:
            log_html += event_card(date, name, color, desc)
        log_html += "</div>"
        st.markdown(log_html, unsafe_allow_html=True)

def tab_timeline(data):
    section_header("mutation frequency")

    mode = st.radio("", ["variant dominance", "individual mutations"], horizontal=True)

    if mode == "individual mutations":
        all_muts = [m[0] for m in KEY_MUTATIONS]
        selected = st.multiselect("", all_muts,
            default=["D614G", "N501Y", "E484K", "L452R", "P681H"],
            label_visibility="collapsed")
        fig = dark_fig(timeline_plot(data["variant_freq"],
            mutation_freq=data["mutation_freq"], mode="mutation",
            selected_mutations=selected))
    else:
        fig = dark_fig(timeline_plot(data["variant_freq"], mode="variant"))

    fig.update_layout(title="", height=240, xaxis_title="", yaxis_title="frequency (%)")
    st.plotly_chart(fig, width='stretch', key="chart_2")

    st.divider()
    section_header("sweep speed", "days from 5% → 50% global frequency")

    sweeps = data["sweeps"]
    if not sweeps.empty:
        col1, col2 = st.columns([2, 1], gap="medium")
        with col1:
            fig2 = dark_fig(sweep_speed_chart(sweeps))
            fig2.update_layout(title="", height=240)
            st.plotly_chart(fig2, key="chart_3", width='stretch')
        with col2:
            st.dataframe(
                sweeps[["mutation","first_seen","sweep_days","peak_freq"]]
                    .rename(columns={"mutation":"mut","first_seen":"first",
                                     "sweep_days":"days","peak_freq":"peak %"}),
                width='stretch', hide_index=True)

def tab_geographic(data):
    section_header("geographic spread", "variant origins by country")
    geo = data["geo"]
    if geo.empty:
        st.info("no geographic data")
        return
    fig = dark_fig(geographic_plotly(geo))
    fig.update_layout(title="", height=300)
    st.plotly_chart(fig, width='stretch', key="chart_4")

    col1, col2 = st.columns(2, gap="medium")
    with col1:
        section_header("first detection")
        first = (geo.sort_values("first_detection")
                    .drop_duplicates(subset=["who_label"])
                    [["who_label","country","first_detection","sequences_count"]]
                    .rename(columns={"who_label":"variant","first_detection":"date","sequences_count":"seqs"}))
        st.dataframe(first, width='stretch', hide_index=True)
    with col2:
        section_header("top countries")
        top = (geo.groupby("country")["sequences_count"].sum().reset_index()
                  .sort_values("sequences_count", ascending=False).head(12)
                  .rename(columns={"sequences_count":"sequences"}))
        st.dataframe(top, width='stretch', hide_index=True)

def tab_spike(data):
    section_header("spike protein", "1,273 aa · mutational hotspot map")

    hs = data["hotspots"]
    fig = dark_fig(spike_protein_map(hs))
    fig.update_layout(title="", height=200)
    st.plotly_chart(fig, width='stretch', key="chart_5")

    col1, col2 = st.columns(2, gap="medium")
    with col1:
        section_header("rbd zoom", "positions 306–527")
        fig2 = dark_fig(spike_rbd_detail(hs))
        fig2.update_layout(title="", height=200)
        st.plotly_chart(fig2, width='stretch', key="chart_6")
    with col2:
        section_header("mutation categories")
        import plotly.graph_objects as go
        from config import FUNCTIONAL_COLOURS
        cats = {}
        for _, _, cat, _ in KEY_MUTATIONS:
            cats[cat] = cats.get(cat, 0) + 1
        fig3 = go.Figure(go.Pie(
            labels=list(cats.keys()), values=list(cats.values()),
            marker_colors=[FUNCTIONAL_COLOURS.get(c, "#3d444d") for c in cats],
            hole=0.55, textinfo="label+percent",
            textfont=dict(size=9, color="#8b949e", family="JetBrains Mono")))
        fig3.update_layout(showlegend=False, height=200,
            margin=dict(l=8,r=8,t=8,b=8),
            paper_bgcolor="#0d1117", plot_bgcolor="#0d1117")
        st.plotly_chart(fig3, width='stretch', key="chart_7")

    section_header("domain breakdown")
    if not hs.empty:
        ds = (hs.groupby("domain")
                .agg(mutations=("mutation_count","sum"),
                     positions=("position","count"),
                     top=("position", lambda x: ", ".join(x.nlargest(3).astype(str))))
                .sort_values("mutations", ascending=False).reset_index())
        st.dataframe(ds, width='stretch', hide_index=True)

def tab_convergent(data):
    section_header("convergent evolution",
                   "mutations that arose independently in multiple lineages")

    conv = data["convergent"]
    if conv.empty:
        st.info("no convergent evolution data")
        return

    import plotly.express as px
    fig = px.bar(conv.nlargest(12, "appearance_count"),
                 x="mutation", y="appearance_count",
                 color="appearance_count",
                 color_continuous_scale=[[0,"#0d2818"],[0.4,"#1a4427"],[1.0,"#3fb950"]],
                 labels={"mutation":"mutation","appearance_count":"lineage appearances"})
    fig.update_layout(coloraxis_showscale=False, showlegend=False,
                      title="", height=200,
                      xaxis_title="", yaxis_title="appearances",
                      **PLOTLY_DARK)
    fig.update_traces(marker_line_width=0)
    st.plotly_chart(fig, width='stretch', key="chart_8")

    st.markdown(
        "<div style='background:#161b22;border:1px solid #21262d;"
        "border-left:2px solid #3fb950;border-radius:0 4px 4px 0;"
        "padding:10px 14px;margin:8px 0'>"
        "<span style='font-family:JetBrains Mono,monospace;font-size:10px;color:#484f58'>"
        "when the same mutation appears independently across unrelated lineages it signals "
        "strong selective pressure — the virus keeps arriving at the same solution. "
        "<span style='color:#8b949e'>E484K</span> is the textbook case: "
        "Beta (South Africa) · Gamma (Brazil) · Iota (USA) · no shared ancestor."
        "</span></div>", unsafe_allow_html=True)

    col1, col2 = st.columns(2, gap="medium")
    with col1:
        section_header("top convergent mutations")
        st.dataframe(
            conv.nlargest(8,"appearance_count")[["mutation","appearance_count","lineages"]]
                .rename(columns={"mutation":"mutation","appearance_count":"appearances","lineages":"found in"}),
            width='stretch', hide_index=True)
    with col2:
        section_header("what this means")
        st.markdown(
            "<div style='font-family:JetBrains Mono,monospace;font-size:10px;"
            "color:#484f58;line-height:1.9;padding-top:4px'>"
            "high bar = high convergence<br>"
            "= strong selection pressure<br><br>"
            "almost always linked to:<br>"
            "<span style='color:#58a6ff'>→ ACE2 binding change</span><br>"
            "<span style='color:#f85149'>→ antibody escape</span><br>"
            "<span style='color:#d29922'>→ furin cleavage gain</span>"
            "</div>", unsafe_allow_html=True)

def tab_data(filtered):
    section_header("sequence database", "raw data · searchable · exportable")
    df = filtered["mutations"].copy()
    if df.empty:
        st.info("no sequences match current filters")
        return

    search = st.text_input("", placeholder="search sample id, country, or mutation  —  e.g. N501Y or India",
                            label_visibility="collapsed")
    if search:
        mask = (df["strain"].str.contains(search, case=False, na=False) |
                df["country"].str.contains(search, case=False, na=False) |
                df["spike_mutations"].str.contains(search, case=False, na=False))
        df = df[mask]

    show_cols = [c for c in ["strain","date","country","region","lineage",
                              "who_label","spike_mut_count","rbd_mut_count"] if c in df.columns]
    st.dataframe(
        df[show_cols].rename(columns={
            "strain":"sample_id","date":"date","country":"country",
            "region":"region","lineage":"pango","who_label":"who",
            "spike_mut_count":"spike_muts","rbd_mut_count":"rbd_muts"})
        .sort_values("date", ascending=False).head(500),
        width='stretch', hide_index=True)

    st.markdown(
        f"<div style='font-family:JetBrains Mono,monospace;font-size:9px;"
        f"color:#3d444d;margin-top:6px'>"
        f"showing {min(500,len(df)):,} of {len(df):,} sequences</div>",
        unsafe_allow_html=True)

    st.download_button("↓ export csv",
        df[show_cols].to_csv(index=False),
        "sars_cov2_mutations.csv", "text/csv")

def tab_about():
    section_header("about")
    st.markdown("""
<div style='font-family:JetBrains Mono,monospace;font-size:11px;color:#484f58;line-height:1.9;max-width:680px'>

<p style='color:#8b949e;margin-bottom:16px'>
a genomic surveillance pipeline tracking how SARS-CoV-2 spike protein mutations emerge,
spread geographically, and affect viral fitness. built on real sequence data from the
nextstrain open database.
</p>

<table style='border-collapse:collapse;width:100%'>
<tr style='border-bottom:1px solid #21262d'>
  <td style='padding:7px 20px 7px 0;color:#3d444d;white-space:nowrap'>fetch_data.py</td>
  <td style='padding:7px 0;color:#484f58'>downloads MN908947 reference · fetches nextstrain metadata</td>
</tr>
<tr style='border-bottom:1px solid #21262d'>
  <td style='padding:7px 20px 7px 0;color:#3d444d;white-space:nowrap'>alignment.py</td>
  <td style='padding:7px 0;color:#484f58'>parses amino-acid substitutions · builds mutation table</td>
</tr>
<tr style='border-bottom:1px solid #21262d'>
  <td style='padding:7px 20px 7px 0;color:#3d444d;white-space:nowrap'>frequency_analysis.py</td>
  <td style='padding:7px 0;color:#484f58'>variant frequency · sweep speed · convergent evolution · hotspots</td>
</tr>
<tr style='border-bottom:1px solid #21262d'>
  <td style='padding:7px 20px 7px 0;color:#3d444d;white-space:nowrap'>plots.py</td>
  <td style='padding:7px 0;color:#484f58'>all plotly figures · choropleth · spike protein map</td>
</tr>
<tr>
  <td style='padding:7px 20px 7px 0;color:#3d444d;white-space:nowrap'>app.py</td>
  <td style='padding:7px 0;color:#484f58'>this streamlit dashboard</td>
</tr>
</table>

<p style='color:#3d444d;margin-top:16px;font-size:9px'>
data · nextstrain.org/ncov/open · gisaid.org · ncbi.nlm.nih.gov/nuccore/MN908947<br>
sequences gratefully acknowledged from all contributing laboratories worldwide.
</p>

</div>
""", unsafe_allow_html=True)

def main():
    st.markdown(
        "<div style='padding:10px 0 6px'>"
        "<span style='font-family:JetBrains Mono,monospace;font-size:20px;"
        "font-weight:600;color:#e6edf3;letter-spacing:-0.02em'>sars-cov-2</span>"
        "<span style='font-family:JetBrains Mono,monospace;font-size:20px;"
        "font-weight:300;color:#1f6feb'> / mutations</span>"
        "<span style='font-family:JetBrains Mono,monospace;font-size:10px;"
        "color:#3d444d;margin-left:16px'>spike protein evolution · 2020–2026 · nextstrain open data</span>"
        "</div>", unsafe_allow_html=True)
    st.divider()

    with st.spinner("loading..."):
        data = load_all_data()

    date_range, region, variant = render_sidebar(data)
    filtered = apply_filters(data, date_range, region, variant)

    tabs = st.tabs([
        "overview", "timeline", "geographic",
        "spike protein", "convergence", "data", "about"
    ])
    with tabs[0]: tab_overview(data, filtered)
    with tabs[1]: tab_timeline(data)
    with tabs[2]: tab_geographic(data)
    with tabs[3]: tab_spike(data)
    with tabs[4]: tab_convergent(data)
    with tabs[5]: tab_data(filtered)
    with tabs[6]: tab_about()

if __name__ == "__main__":
    main()
