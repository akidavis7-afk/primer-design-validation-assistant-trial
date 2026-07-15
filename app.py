from __future__ import annotations

import io
from time import perf_counter

import pandas as pd
import streamlit as st
import yaml

from src.core import analyze_pair, read_fasta_text

MAX_PRIMER_PAIRS = 5000
MAX_TARGET_BASES = 5_000_000


@st.cache_data(show_spinner=False, max_entries=20)
def analyze_primers_cached(target_text: str, primer_csv: bytes) -> pd.DataFrame:
    sequence = read_fasta_text(target_text)
    primer_df = pd.read_csv(io.BytesIO(primer_csv))
    required = {"pair_id", "forward", "reverse"}
    missing = sorted(required.difference(primer_df.columns))
    if missing:
        raise ValueError(f"Primer CSV is missing required columns: {', '.join(missing)}")
    if len(primer_df) > MAX_PRIMER_PAIRS:
        raise ValueError(f"Primer CSV contains more than {MAX_PRIMER_PAIRS} pairs.")

    rows = []
    for row in primer_df.itertuples(index=False):
        rows.append(
            {
                "pair_id": getattr(row, "pair_id"),
                **analyze_pair(sequence, str(getattr(row, "forward")), str(getattr(row, "reverse"))),
            }
        )
    return pd.DataFrame(rows)


with open("configs/default.yaml", encoding="utf-8") as handle:
    cfg = yaml.safe_load(handle) or {}

st.set_page_config(page_title=cfg.get("app_title", "Primer Design & Validation Assistant"), layout="wide")
st.title(cfg.get("app_title", "Primer Design & Validation Assistant"))
st.caption("Batch primer QC, binding-site checks, and downloadable reports.")

with st.form("primer_validation_form"):
    target_text = st.text_area("Target sequence or FASTA", height=180)
    primer_file = st.file_uploader("Primer CSV with pair_id, forward, reverse", type=["csv"])
    submitted = st.form_submit_button("Validate primers", type="primary")

if submitted:
    sequence = read_fasta_text(target_text)
    if not sequence or primer_file is None:
        st.error("Provide a target sequence and a primer CSV file.")
    elif len(sequence) > MAX_TARGET_BASES:
        st.error(f"Target sequence exceeds {MAX_TARGET_BASES:,} bases.")
    else:
        started = perf_counter()
        try:
            with st.spinner("Checking primer pairs..."):
                result = analyze_primers_cached(target_text, primer_file.getvalue())
            st.session_state["primer_result"] = result
            st.session_state["primer_elapsed"] = perf_counter() - started
        except Exception as exc:
            st.exception(exc)

result = st.session_state.get("primer_result")
if result is not None:
    st.success(f"Analysis completed in {st.session_state.get('primer_elapsed', 0.0):.2f} seconds.")
    st.dataframe(result, use_container_width=True, hide_index=True)
    st.download_button(
        "Download report",
        result.to_csv(index=False).encode("utf-8"),
        "primer_report.csv",
        "text/csv",
    )
