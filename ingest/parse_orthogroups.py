import duckdb
import pandas as pd


ORTHOGROUPS_TSV = "data/raw/Orthogroups_diamond.tsv"
DB_PATH = "warehouse.duckdb"


wide = pd.read_csv(ORTHOGROUPS_TSV, sep="\t", dtype=str)


long = wide.melt(
    id_vars="Orthogroup",
    var_name="species_label",
    value_name="protein_ids",
).dropna(subset=["protein_ids"])


long["protein_id"] = long["protein_ids"].str.split(r",\s*", regex=True)
long = long.explode("protein_id")

long = long.rename(columns={"Orthogroup": "orthogroup_id"})
long["protein_id"] = long["protein_id"].str.strip()
membership = long[["orthogroup_id", "species_label", "protein_id"]]
membership = membership[membership["protein_id"].notna() & (membership["protein_id"] != "")]


con = duckdb.connect(DB_PATH)
con.register("membership_df", membership)
con.execute("CREATE OR REPLACE TABLE raw_orthogroup_membership AS SELECT * FROM membership_df")

rows = con.execute("SELECT count(*) FROM raw_orthogroup_membership").fetchone()[0]
ogs = con.execute("SELECT count(DISTINCT orthogroup_id) FROM raw_orthogroup_membership").fetchone()[0]
spp = con.execute("SELECT count(DISTINCT species_label) FROM raw_orthogroup_membership").fetchone()[0]
con.close()

print(f"Loaded {rows:,} protein rows across {ogs:,} orthogroups and {spp} species.")