# Primer Design & Validation Assistant

Validates PCR, sequencing and mutagenesis primer pairs using a reusable CLI + Streamlit interface.

[Click here to view the Live Interactive Web App Demo](https://primer-design-validation-assistant-trial-eh2wj8grt6dkfdfuhtx69.streamlit.app/)

## Features

- GC percentage, length and Wallace-rule melting temperature
- Homopolymer and GC-clamp warnings
- Self-complementarity and pair-dimer heuristics
- Forward/reverse binding search
- Amplicon size validation
- Mutagenesis mismatch summary
- CSV and HTML-ready outputs

## Local setup (Windows CMD)

```cmd
py -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
pytest -q
python cli.py --target examples/target.fasta --primers examples/primers.csv --output outputs
streamlit run app.py
```

Open `http://localhost:8501`.

## Docker

```cmd
docker build -t primer-assistant .
docker run --rm -p 8501:8501 primer-assistant
```

Personalize with `configs/lab_profile.yaml` and target-specific synthetic or public files under `examples/<lab_name>/`.

## Streamlit performance

The app uses a submit form, `st.cache_data`, immutable upload bytes, and session-state result persistence. This prevents the analysis from running again when a widget changes. The sequence-comparison apps also use a direct linear comparison for closely related equal-length sequences and reserve global alignment for likely indels or larger differences.
