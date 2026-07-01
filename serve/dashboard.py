import streamlit as st, duckdb
df = duckdb.connect("warehouse.duckdb").execute("select * from stg_orthogroups limit 20").df()
st.title("Haptophyte calcification — skeleton")
st.dataframe(df)