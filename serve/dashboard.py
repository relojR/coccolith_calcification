"""
Haptophyte calcification - presence/absence dashboard.

Controls live in the sidebar; results render in the main canvas.

"""

from pathlib import Path

import duckdb
import pandas as pd
import plotly.express as px
import streamlit as st
import streamlit.components.v1 as components

DB_PATH = Path(__file__).resolve().parent.parent / "warehouse.duckdb"

CALC_COLOR = "#2a9d8f"
NONCALC_COLOR = "#adb5bd"
COLOR_MAP = {"Calcifying": CALC_COLOR, "Non-calcifying": NONCALC_COLOR}

st.set_page_config(page_title="Haptophyte Calcification", layout="wide")


st.markdown(
    """
    <style>
    [class*="st-key-section-"] {
        border: 1.5px solid #2a9d8f;
        border-radius: 10px;
        padding: 1rem 1.25rem;
        margin-bottom: 0.75rem;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

if not DB_PATH.exists():
    st.error(
        f"Database not found at {DB_PATH}.\n\n"
        "Run the pipeline first: `python ingest\\parse_orthogroups.py`, "
        "then `cd genomics; dbt build`."
    )
    st.stop()


@st.cache_resource
def get_con():
    return duckdb.connect(str(DB_PATH), read_only=True)


@st.cache_data
def run(sql: str) -> pd.DataFrame:
    return get_con().execute(sql).df()


def request_scroll():
    """Fired by the picker/search on_change so we only scroll on real user input."""
    st.session_state["scroll_to_detail"] = True


# header
st.title("Haptophyte calcification - gene family explorer")
st.caption(
    "Comparative genomics of 27 haptophyte species: which orthogroups "
    "are present in calcifying vs. non-calcifying species?"
)

# shared queries 
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

# sidebars
with st.sidebar:
    st.header("Controls")

    st.markdown("**Filter by presence pattern**")

    min_calc = st.slider("Minimum calcifying species present", 0, n_calc, n_calc)
    max_noncalc = st.slider("Maximum non-calcifying species present", 0, n_noncalc, 0)

    filtered = candidates[
        (candidates.calc_present >= min_calc)
        & (candidates.noncalc_present <= max_noncalc)
    ]
    st.caption(
        f"**{len(filtered):,}** of {len(candidates):,} orthogroups match."
    )

    st.divider()

    st.markdown("**Select an orthogroup**")
    options = filtered.orthogroup_id.head(500).tolist()
    if options:
        picked = st.selectbox(
            "Pick a matching orthogroup to see more details about it", options, index=0,
            on_change=request_scroll,
        )
    else:
        picked = None
        st.info("No matches. Loosen the sliders, or search for an existing ID below.")

    typed = st.text_input(
        "Or search any orthogroup ID (e.g. OG0000123)", "",
        on_change=request_scroll,
    ).strip()


og = typed if typed else picked

#  TABS
tab_finder, tab_overview = st.tabs(["Candidate Finder", "Dataset Overview"])

# Candidate Finder
with tab_finder:

    with st.container(key="section-matching"):
        st.subheader("Matching orthogroups")
        display_df = filtered.head(500).reset_index(drop=True)
        st.dataframe(
            display_df,
            height=300,
            use_container_width=True,
            hide_index=True,
        )

    with st.container(key="section-detail"):
        st.subheader("Selected orthogroup")
        if og is None:
            st.info("Adjust the sidebar controls to select an orthogroup.")
        elif og not in valid_ids:
            st.warning(f"No orthogroup '{og}' found in the data.")
        else:
            row = candidates[candidates.orthogroup_id == og].iloc[0]
            st.markdown(f"### {og}")
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
                detail, x="name", y="gene_count", color="status",
                color_discrete_map=COLOR_MAP,
                labels={"name": "", "gene_count": "Genes in family", "status": ""},
            )
            fig2.update_layout(xaxis_tickangle=-45, height=360, legend_title="")
            st.plotly_chart(fig2, use_container_width=True)

            members = run(
                "select "
                "coalesce(nullif(s.species_name,''), m.species_label) as species, "
                "s.is_calcifying, m.protein_id "
                "from stg_orthogroup_membership m "
                "join stg_species s using (species_label) "
                f"where m.orthogroup_id = '{og}' "
                "order by s.is_calcifying desc, species, protein_id"
            )
            members["status"] = members.is_calcifying.map({True: "Calcifying", False: "Non-calcifying"})

            st.markdown(f"**Member proteins** \u2014 {len(members):,} across "
                        f"{members.species.nunique()} species")
            st.dataframe(
                members[["species", "status", "protein_id"]],
                height=280,
                use_container_width=True,
                hide_index=True,
            )
            st.download_button(
                "Download protein IDs (CSV)",
                members[["species", "status", "protein_id"]].to_csv(index=False),
                file_name=f"{og}_members.csv",
                mime="text/csv",
            )

# Dataset Overview
with tab_overview:

    with st.container(key="section-metrics"):
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Total Species analyzed", n_species)
        c2.metric("Calcifying species", n_calc)
        c3.metric("Non-calcifying species", n_noncalc)
        c4.metric("Orthogroups", f"{n_og:,}")

    with st.container(key="section-coverage"):
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
            coverage, x="name", y="n_present", color="status",
            color_discrete_map=COLOR_MAP,
            labels={"name": "", "n_present": "Orthogroups present", "status": ""},
        )
        fig.update_layout(xaxis_tickangle=-45, height=450, legend_title="")
        st.plotly_chart(fig, use_container_width=True)


if st.session_state.get("scroll_to_detail"):

    st.session_state["scroll_nonce"] = st.session_state.get("scroll_nonce", 0) + 1
    nonce = st.session_state["scroll_nonce"]

    scroll_html = """
        <script>
        // scroll nonce __NONCE__
        setTimeout(function () {
            try {
                var el = window.parent.document.querySelector('.st-key-section-detail');
                if (el) { el.scrollIntoView({behavior: 'smooth', block: 'start'}); }
            } catch (e) {}
        }, 100);
        </script>
    """.replace("__NONCE__", str(nonce))

    components.html(scroll_html, height=0)
    st.session_state["scroll_to_detail"] = False