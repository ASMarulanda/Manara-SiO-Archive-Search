# Manara SiO(v=0) Archive Search Pipeline

This folder contains a complete, reproducible Python pipeline for identifying **all ALMA archival observations** associated with the **Manara (PPVII) YSO survey** that include **any SiO(v=0) rotational transition** from  
**J = 1→0** up to **J = 20→19**.

The script automates:

1. **Manara sample preprocessing**  
2. **ALMA archive querying (ALminer/TAP)**  
3. **Spectral-window filtering for SiO lines**  
4. **Generation of per-SPW and per-MOUS tables**  
5. **Export to CSV and LaTeX for thesis/paper use**

The output defines the **SiO-sensitive subset** of the Manara survey and feeds directly into the next pipeline stages (download, calibration, SPW extraction).

---

## Contents

```
manara_sio_archive_search.py   # Main pipeline script
manara_refined_sample.csv      # Clean Manara sample (auto-generated)
manara_alma_obscore.csv        # Combined ALMA ObsCore table (auto-generated)
manara_sio_spw_matches.csv     # SPW-level SiO matches (auto-generated)
manara_sio_mous_summary.csv    # MOUS-level SiO matches (auto-generated)
manara_sio_spw_matches.tex     # SPW-level table in LaTeX format
manara_sio_mous_summary.tex    # MOUS-level table in LaTeX format
```

---

## Pipeline Overview

### 1. **Load and clean the Manara catalogue**

The script downloads:

```
PP7-Surveys_2022-10-19_PPVII_website.tsv
```

Then it:

- removes irrelevant columns, fixes numeric types  
- applies quality cuts:  
  - must have **L\***, **M\***, and **distance**  
- classifies sources as:
  - **LM** (low-mass, M⋆ < 1 M⊙)  
  - **HM** (high-mass)  
- keeps only entries with **ALMA** in the “Observatory” column  

This produces the refined sample:

```
manara_refined_sample.csv
```

---

## Installation

We recommend running inside a virtual environment:

```bash
python3 -m venv venv
source venv/bin/activate
```

Install required packages:

```bash
pip install pandas numpy seaborn matplotlib astroquery alminer
```

---

## Running the Pipeline

From inside this folder:

```bash
python3 manara_sio_archive_search.py
```

The script will automatically:

1. Download & preprocess the Manara table  
2. Filter the sample  
3. Query ALMA via ALminer  
4. Build the ObsCore table  
5. Find SiO SPWs  
6. Build per-SPW and per-MOUS tables  
7. Export CSV + LaTeX  

---

## Output Summary

After running, you should see:

| File | Description |
|------|-------------|
| **manara_refined_sample.csv** | Clean sample of Manara YSOs with ALMA coverage |
| **manara_alma_obscore.csv** | All ALMA archival rows found for the refined sample |
| **manara_sio_spw_matches.csv / .tex** | One row per SPW/SiO match |
| **manara_sio_mous_summary.csv / .tex** | One row per MOUS ID with aggregated SiO info |
- This does not affect the SiO filtering—unresolved names simply produce no ObsCore rows.

---

## License

MIT License.
