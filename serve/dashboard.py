"""
Haptophyte calcification - presence/absence dashboard (MVP).

Reads the dbt-built models from warehouse.duckdb and lets you explore which
gene families (orthogroups) are present across calcifying vs non-calcifying
species.

Run from the PROJECT ROOT:
    streamlit run serve/dashboard.py
"""

from pathlib import Path

import duckdb
import pandas as pd
import plotly.express as px
import streamlit as st

# warehouse.duckdb lives at the project root, one level up from serve/.
# Resolving relative to THIS file means it works no matter where you launch from
# (no more root-vs-genomics working-directory surprises).
DB_PATH = Path(__file__).resolve().parent.parent / "warehouse.duckdb"

CALC_COLOR = "#2a9d8f"
NONCALC_COLOR = "#adb5bd"
COLOR_MAP = {"Calcifying": CALC_COLOR, "Non-calcifying": NONCALC_COLOR}

st.set_page_config(page_title="Haptophyte Calcification", layout="wide")

if not DB_PATH.exists():
    st.error(
        f"Database not found at {DB_PATH}.\n\n"
        "Run the pipeline first: `python ingest\\parse_orthogroups.py`, "
        "then `cd genomics; dbt build`."
    )
    st.stop()


@st.cache_resource
def get_con():
    # read_only so the dashboard never locks the file or fights with dbt.
    return duckdb.connect(str(DB_PATH), read_only=True)


@st.cache_data
def run(sql: str) -> pd.DataFrame:
    return get_con().execute(sql).df()


# --------------------------------------------------------------------- header
st.title("Haptophyte calcification - gene family explorer")
st.caption(
    "Comparative genomics of 27 haptophyte species: which orthogroups "
    "are present in calcifying vs. non-calcifying species?"
)

# ------------------------------------------------------------------- overview
species = run(
    "select species_label, "
    "coalesce(nullif(species_name,''), species_label) as name, "
    "is_calcifying, clade "
    "from stg_species"
)
n_species = len(species)
n_calc = int(species.is_calcifying.sum())
n_noncalc = n_species - n_calc
n_og = int(run("select count(distinct orthogroup_id) as n from int_orthogroup_presence").n.iloc[0])

c1, c2, c3, c4 = st.columns(4)
c1.metric("Species", n_species)
c2.metric("Calcifying", n_calc)
c3.metric("Non-calcifying", n_noncalc)
c4.metric("Gene families", f"{n_og:,}")

st.divider()

# ----------------------------------------------------------- per-species coverage
st.subheader("Gene family coverage per species")
st.caption("How many orthogroups each species appears in, colored by calcification status.")

coverage = run(
    """
    select
        p.species_label,
        coalesce(nullif(s.species_name,''), p.species_label) as name,
        s.is_calcifying,
        count(*) filter (where p.is_present) as n_present
    from int_orthogroup_presence p
    join stg_species s using (species_label)
    group by 1, 2, 3
    order by is_calcifying desc, n_present desc
    """
)
coverage["status"] = coverage.is_calcifying.map({True: "Calcifying", False: "Non-calcifying"})

fig = px.bar(
    coverage,
    x="name",
    y="n_present",
    color="status",
    color_discrete_map=COLOR_MAP,
    labels={"name": "", "n_present": "Orthogroups present", "status": ""},
)
fig.update_layout(xaxis_tickangle=-45, height=450, legend_title="")
st.plotly_chart(fig, use_container_width=True)

st.divider()

# ------------------------------------------------------------------- explorer
st.subheader("Find and explore candidate gene families")

# Per-orthogroup presence counts (calc vs non-calc). Computed once, then
# filtered in pandas so the sliders respond instantly with no re-query.
candidates = run(
    """
    with totals as (
        select
            count(*) filter (where is_calcifying) as n_calc,
            count(*) filter (where not is_calcifying) as n_noncalc
        from (select distinct species_label, is_calcifying from int_orthogroup_presence)
    ),
    per_og as (
        select
            orthogroup_id,
            count(*) filter (where is_calcifying and is_present) as calc_present,
            count(*) filter (where not is_calcifying and is_present) as noncalc_present
        from int_orthogroup_presence
        group by 1
    )
    select
        orthogroup_id,
        calc_present,
        noncalc_present,
        round(calc_present::double / n_calc, 3)     as calc_frac,
        round(noncalc_present::double / n_noncalc, 3) as noncalc_frac,
        round(calc_present::double / n_calc
              - noncalc_present::double / n_noncalc, 3) as presence_diff
    from per_og cross join totals
    order by presence_diff desc, calc_present desc
    """
)
valid_ids = set(candidates.orthogroup_id)

# --- filter controls ---
st.markdown("**Filter by presence pattern**")
st.caption(
    f"Set the left slider to {n_calc} and the right to 0 to find gene families "
    "present in *every* calcifier and *no* non-calcifier - the strongest signal."
)
f1, f2 = st.columns(2)
with f1:
    min_calc = st.slider(
        "Minimum calcifying species present", 0, n_calc, n_calc
    )
with f2:
    max_noncalc = st.slider(
        "Maximum non-calcifying species present", 0, n_noncalc, 0
    )

filtered = candidates[
    (candidates.calc_present >= min_calc)
    & (candidates.noncalc_present <= max_noncalc)
]

st.caption(
    f"**{len(filtered):,}** of {len(candidates):,} orthogroups match "
    f"(present in \u2265{min_calc} calcifiers and \u2264{max_noncalc} non-calcifiers)."
)

col_l, col_r = st.columns([1, 2])

with col_l:
    st.markdown("**Matching orthogroups**")
    st.dataframe(
        filtered.head(500),
        height=430,
        use_container_width=True,
        hide_index=True,
    )

with col_r:
    options = filtered.orthogroup_id.head(500).tolist()
    if options:
        picked = st.selectbox("Pick a matching orthogroup", options, index=0)
    else:
        picked = None
        st.info("No orthogroups match the current filter. Loosen the sliders, or type an ID below.")

    typed = st.text_input("...or type any orthogroup ID (e.g. OG0000123)", "").strip()
    og = typed if typed else picked

    if not og:
        st.stop()
    if og not in valid_ids:
        st.warning(f"No orthogroup '{og}' found in the data.")
    else:
        row = candidates[candidates.orthogroup_id == og].iloc[0]
        m1, m2, m3 = st.columns(3)
        m1.metric("Calcifiers with gene", f"{int(row.calc_present)} / {n_calc}")
        m2.metric("Non-calcifiers with gene", f"{int(row.noncalc_present)} / {n_noncalc}")
        m3.metric("Presence difference", f"{row.presence_diff:+.0%}")

        detail = run(
            "select "
            "coalesce(nullif(s.species_name,''), p.species_label) as name, "
            "s.is_calcifying, p.is_present, p.gene_count "
            "from int_orthogroup_presence p "
            "join stg_species s using (species_label) "
            f"where p.orthogroup_id = '{og}'"
        )
        detail["status"] = detail.is_calcifying.map({True: "Calcifying", False: "Non-calcifying"})
        detail = detail.sort_values(["is_calcifying", "gene_count"], ascending=[False, False])

        fig2 = px.bar(
            detail,
            x="name",
            y="gene_count",
            color="status",
            color_discrete_map=COLOR_MAP,
            labels={"name": "", "gene_count": "Genes in family", "status": ""},
        )
        fig2.update_layout(xaxis_tickangle=-45, height=400, legend_title="")
        st.plotly_chart(fig2, use_container_width=True)

