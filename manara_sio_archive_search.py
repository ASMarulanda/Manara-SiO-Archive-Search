#!/usr/bin/env python3
"""
Manara survey: SiO(v=0) archive search and per-MOUS summary.

Pipeline steps:

1) Load the Manara (PPVII) survey table from the public URL.
2) Clean and filter the sample:
     - keep only rows with valid L*, M*, distance
     - classify sources as LM / HM using M* < 1 M_sun
     - keep only entries with ALMA coverage in the Observatory field
3) For each unique ALMA source, query the ALMA archive via ALminer
   (ObsCore / TAP) and combine all results into a single ObsCore table.
   This table is cached to disk as 'manara_alma_obscore.csv'.
   Name resolution is robust:
     - try many name variants (spaces, underscores, A/B suffixes, etc.)
     - for each variant, try alminer.target() and alminer.keysearch()
4) For every ObsCore row, check whether its [min_freq_GHz, max_freq_GHz]
   range contains any SiO(v=0) transition from J=1–0 to J=20–19.
5) Build:
     - manara_sio_spw_matches.csv / .tex  (per SPW, one row per SiO line)
     - manara_sio_mous_summary.csv / .tex (per MOUS, listing all SiO lines)

Only MOUS IDs with at least one SiO(v=0) line appear in the final tables.
This matches the “SiO-positive subset” branch in the Manara flowchart.
"""

from __future__ import annotations
from typing import Dict, List

import os
import warnings

import numpy as np
import pandas as pd
import alminer

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------
# 1. User configuration
# --------------------------------------------------------------------

# Source of the Manara / PPVII survey table
MANARA_URL = "http://ppvii.org/chapter/15/PP7-Surveys_2022-10-19_PPVII_website.tsv"

# Output / cache filenames
REFINED_SAMPLE_CSV = "manara_refined_sample.csv"
ALMA_OBSCORE_CSV   = "manara_alma_obscore.csv"

SPW_CSV_OUT  = "manara_sio_spw_matches.csv"
SPW_TEX_OUT  = "manara_sio_spw_matches.tex"
MOUS_CSV_OUT = "manara_sio_mous_summary.csv"
MOUS_TEX_OUT = "manara_sio_mous_summary.tex"

# ALminer / ALMA query settings
TAP_SERVICE = "ESO"           # ESO mirror
SEARCH_RADIUS_ARCMIN = 2.0    # cone-search radius (arcmin)

# SiO(v=0) transitions in GHz
SIO_V0_TRANSITIONS_GHZ: Dict[str, float] = {
    "J=1-0":   43.423864,
    "J=2-1":   86.846960,
    "J=3-2":  130.268610,
    "J=4-3":  173.688310,
    "J=5-4":  217.104980,
    "J=6-5":  260.518200,
    "J=7-6":  303.927030,
    "J=8-7":  347.331000,
    "J=9-8":  390.728730,
    "J=10-9": 434.120450,
    "J=11-10":477.506120,
    "J=12-11":520.885480,
    "J=13-12":564.258560,
    "J=14-13":607.625260,
    "J=15-14":650.985560,
    "J=16-15":694.339440,
    "J=17-16":737.686780,
    "J=18-17":781.027470,
    "J=19-18":824.361490,
    "J=20-19":867.688720,
}

# Expected ObsCore columns from ALminer
COL_SOURCE   = "Source"        # Manara source name
COL_PROJECT  = "project_code"
COL_BAND     = "band_list"
COL_MIN_FREQ = "min_freq_GHz"
COL_MAX_FREQ = "max_freq_GHz"
COL_ANG_RES  = "ang_res_arcsec"
COL_MOUS     = "MOUS_id"       # ALminer usually calls this 'MOUS_id'


# --------------------------------------------------------------------
# 2. Helper: generate name variants for ALMA queries
# --------------------------------------------------------------------

def generate_name_variants(src: str) -> List[str]:
    """
    Generate a wide range of possible ALMA naming variants for a source.

    ALMA names are highly inconsistent across projects, so we try:
      - exact name
      - different separators (space, underscore, hyphen, none)
      - lower/upper/title case
      - without A/B/C component suffixes
      - truncated to first two tokens
      - remove digits, hyphens, underscores

    Returns a list of unique candidate names (order is not guaranteed).
    """
    variants = set()

    # Original
    variants.add(src)

    # Normalize spacing
    tokens = src.split()
    if tokens:
        variants.add("".join(tokens))          # "DPTauA"
        variants.add("_".join(tokens))         # "DP_Tau_A"
        variants.add("-".join(tokens))         # "DP-Tau-A"

    # Case variations
    variants.add(src.lower())
    variants.add(src.upper())
    variants.add(src.title())

    # Remove trailing component labels (A, B, C)
    if tokens and tokens[-1].upper() in ["A", "B", "C"]:
        base = " ".join(tokens[:-1])          # "DP Tau"
        variants.add(base)
        variants.add(base.replace(" ", ""))
        variants.add(base.replace(" ", "_"))
        variants.add(base.replace(" ", "-"))

    # Truncated (first two tokens)
    if len(tokens) >= 2:
        first_two = " ".join(tokens[:2])      # "DP Tau" from "DP Tau A"
        variants.add(first_two)
        variants.add(first_two.replace(" ", ""))
        variants.add(first_two.replace(" ", "_"))
        variants.add(first_two.replace(" ", "-"))

    # Remove digits (some ALMA entries drop catalog numbering)
    stripped = "".join([ch for ch in src if not ch.isdigit()])
    variants.add(stripped.strip())

    # Remove hyphens entirely
    variants.add(src.replace("-", ""))

    # Remove underscores entirely
    variants.add(src.replace("_", ""))

    # Drop empty strings if any
    variants = {v for v in variants if v.strip()}

    return list(variants)


# --------------------------------------------------------------------
# 3. Manara sample: load and refine
# --------------------------------------------------------------------

def load_and_refine_manara() -> pd.DataFrame:
    """Load PPVII/Manara table, clean numeric columns, filter sample."""
    print(f"Loading Manara survey from:\n  {MANARA_URL}")
    df = pd.read_csv(MANARA_URL, delimiter="\t")

    # Drop irrelevant unnamed columns
    cols_to_drop = [c for c in df.columns if c.startswith("Unnamed:")]
    if cols_to_drop:
        df = df.drop(columns=cols_to_drop)

    # Select relevant columns
    selected_cols = [
        "Region", "Source", "RA", "Dec", "Disk",
        "Observatory", "Freq_GHz", "Wavelength",
        "Beam_arcsec", "RMS_mJy/beam", "Survey_Ref",
        "Flux_mJy", "Flux_mJy_Standardized", "Mdust_Mearth_Standardized",
        "Notes", "Object", "SpT_xs", "Teff_xs", "Av_xs",
        "Lstar_xs", "Mstar_PPVII", "logMacc_PPVII",
        "notes_Macc_PPVII", "dist_PPVII",
    ]
    df = df[selected_cols].copy()

    # Numeric fields
    num_vars = [
        "Teff_xs", "Av_xs", "Lstar_xs",
        "Mstar_PPVII", "logMacc_PPVII", "dist_PPVII",
    ]
    for col in num_vars:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df[num_vars] = df[num_vars].replace(-99, np.nan)

    # Keep sources with complete L*, M*, distance
    df_ref = df.dropna(subset=["Lstar_xs", "Mstar_PPVII", "dist_PPVII"]).copy()

    # Tag LM / HM
    df_ref["Mstar_PPVII_tag"] = np.where(df_ref["Mstar_PPVII"] < 1.0, "LM", "HM")

    print("\nRefined sample size after L, M, d cuts:")
    print(f"  Total rows: {len(df_ref)}")
    print("  LM/HM counts:")
    print(df_ref["Mstar_PPVII_tag"].value_counts())

    # Keep only ALMA observations
    alma_filter = [
        "[ALMA, ALMA]", "[ALMA, ]", "[, ALMA]",
        "[ALMA, pre-ALMA]", "[pre-ALMA, ALMA]",
    ]
    df_alma = df_ref[df_ref["Observatory"].isin(alma_filter)].copy()

    print("\nRows with ALMA coverage:")
    print(f"  {len(df_alma)} rows, {df_alma['Source'].nunique()} unique sources")

    df_alma.to_csv(REFINED_SAMPLE_CSV, index=False)
    print(f"Saved refined ALMA sample to: {REFINED_SAMPLE_CSV}")

    return df_alma


# --------------------------------------------------------------------
# 4. Query ALMA via ALminer for all Manara ALMA sources
# --------------------------------------------------------------------

def query_alma_for_manara(df_alma: pd.DataFrame) -> pd.DataFrame:
    """
    For each unique source in df_alma, query ALMA via ALminer.

    For each Manara source name:
      - generate a list of possible ALMA naming variants
      - for each variant:
          * try alminer.target()
          * if empty, try alminer.keysearch({'target_name': [variant]})
      - if any variant returns rows, adopt that result and move on
      - if all variants fail, record failure and continue

    Combine all ObsCore rows into a single table.
    """
    unique_sources = sorted(df_alma["Source"].unique())
    print(f"\nQuerying ALMA for {len(unique_sources)} unique Manara sources...")

    all_obs = []
    stats = []

    for idx, src in enumerate(unique_sources, start=1):
        print(f"\n[{idx}/{len(unique_sources)}] Source: {src}")
        candidates = generate_name_variants(src)
        result = None
        total_rows_for_src = 0

        for cand in candidates:
            print(f"    Trying name variant: '{cand}'")

            # --- Attempt 1: alminer.target() ---
            try:
                r = alminer.target(
                    [cand],
                    search_radius=SEARCH_RADIUS_ARCMIN,
                    tap_service=TAP_SERVICE,
                    point=False,
                    public=True,
                    published=None,
                    print_query=False,
                    print_targets=False,
                )
                if r is not None and len(r) > 0:
                    print(f"      SUCCESS with target(): {cand} ({len(r)} rows)")
                    r = r.copy()
                    r[COL_SOURCE] = src  # preserve original Manara name
                    result = r
                    total_rows_for_src = len(r)
                    break
            except Exception as exc:
                print(f"      target() FAILED for '{cand}': {exc}")

            # --- Attempt 2: alminer.keysearch() ---
            try:
                r = alminer.keysearch({'target_name': [cand]})
                if r is not None and len(r) > 0:
                    print(f"      SUCCESS with keysearch(): {cand} ({len(r)} rows)")
                    r = r.copy()
                    r[COL_SOURCE] = src
                    result = r
                    total_rows_for_src = len(r)
                    break
            except Exception as exc2:
                print(f"      keysearch() FAILED for '{cand}': {exc2}")

        # After trying all variants
        if result is None or total_rows_for_src == 0:
            print(f"    NO ALMA RESULTS FOUND after trying all name variants for '{src}'")
            stats.append({
                "Source": src,
                "success": False,
                "n_results": 0,
                "error_msg": "all variants failed",
            })
            continue

        # Successful retrieval
        all_obs.append(result)
        stats.append({
            "Source": src,
            "success": True,
            "n_results": total_rows_for_src,
            "error_msg": "",
        })

    if not all_obs:
        raise RuntimeError("No ALMA ObsCore rows returned for any Manara source.")

    obs_df = pd.concat(all_obs, ignore_index=True)
    print(f"\nTotal ObsCore rows retrieved (all sources combined): {len(obs_df)}")

    # Save stats and ObsCore combined table
    stats_df = pd.DataFrame(stats)
    stats_df.to_csv("manara_alma_query_stats.csv", index=False)
    print("Saved query stats to: manara_alma_query_stats.csv")

    obs_df.to_csv(ALMA_OBSCORE_CSV, index=False)
    print(f"Saved combined ObsCore table to: {ALMA_OBSCORE_CSV}")

    return obs_df


def load_or_query_obscore(df_alma: pd.DataFrame) -> pd.DataFrame:
    """
    If a cached ObsCore CSV exists, load it. Otherwise, run the ALMA queries.
    """
    if os.path.exists(ALMA_OBSCORE_CSV):
        print(f"\nFound existing ObsCore cache: {ALMA_OBSCORE_CSV}")
        obs = pd.read_csv(ALMA_OBSCORE_CSV)
        # Ensure the Source column is present
        if COL_SOURCE not in obs.columns:
            obs[COL_SOURCE] = "UNKNOWN"
        print(f"Loaded {len(obs)} rows from cache.")
        return obs
    else:
        return query_alma_for_manara(df_alma)


# --------------------------------------------------------------------
# 5. SiO filtering: per-SPW and per-MOUS tables
# --------------------------------------------------------------------

def find_sio_matches_per_spw(obs_df: pd.DataFrame) -> pd.DataFrame:
    """
    For each ObsCore row (SPW/product), check which SiO transitions lie within
    [min_freq_GHz, max_freq_GHz]. Return a dataframe with one row per
    (SPW, transition) match.
    """
    # Ensure numeric min/max
    for col in (COL_MIN_FREQ, COL_MAX_FREQ):
        if col not in obs_df.columns:
            raise KeyError(f"ObsCore table missing column '{col}'.")
        obs_df[col] = pd.to_numeric(obs_df[col], errors="coerce")

    rows: List[dict] = []

    for _, row in obs_df.iterrows():
        nu_min = row[COL_MIN_FREQ]
        nu_max = row[COL_MAX_FREQ]

        if not np.isfinite(nu_min) or not np.isfinite(nu_max):
            continue
        if nu_min > nu_max:
            nu_min, nu_max = nu_max, nu_min

        for trans, nu_sio in SIO_V0_TRANSITIONS_GHZ.items():
            if nu_min <= nu_sio <= nu_max:
                rows.append({
                    "Source":        row.get(COL_SOURCE, "UNKNOWN"),
                    "Project":       row.get(COL_PROJECT, ""),
                    "ALMA_Band":     row.get(COL_BAND, ""),
                    "min_freq_GHz":  nu_min,
                    "max_freq_GHz":  nu_max,
                    "SiO_transition": trans,
                    "SiO_freq_GHz":   nu_sio,
                    "ang_res_arcsec": row.get(COL_ANG_RES, np.nan),
                    "MOUS_ID":       row.get(COL_MOUS, ""),
                })

    if not rows:
        print("No SiO(v=0) transitions found in any SPW.")
        return pd.DataFrame()

    matches = pd.DataFrame(rows)
    print(f"\nTotal SiO-covering SPW rows (Manara): {len(matches)}")

    # Quick per-source summary
    print("\nSiO-covering SPWs per source:")
    print(matches["Source"].value_counts())

    return matches


def build_per_mous_table(spw_df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate SPW-level matches to MOUS level:
    one row per MOUS_ID with lists of transitions and frequencies.
    """
    def join_unique(series: pd.Series) -> str:
        return ", ".join(sorted(series.astype(str).unique()))

    grouped = (
        spw_df.groupby("MOUS_ID")
        .agg({
            "Source":         join_unique,
            "Project":        join_unique,
            "ALMA_Band":      join_unique,
            "SiO_transition": join_unique,
            "SiO_freq_GHz":   join_unique,
        })
        .reset_index()
    )

    grouped = grouped.rename(columns={
        "SiO_transition": "SiO_transitions",
        "SiO_freq_GHz":   "SiO_freqs_GHz",
    })

    print(f"\nTotal unique MOUS IDs with SiO coverage (Manara): {len(grouped)}")

    return grouped.sort_values(["Source", "Project", "MOUS_ID"])


def save_tables(spw_df: pd.DataFrame, mous_df: pd.DataFrame) -> None:
    """Save CSV + LaTeX tables for SPW-level and MOUS-level summaries."""
    # CSV
    spw_df.to_csv(SPW_CSV_OUT, index=False)
    mous_df.to_csv(MOUS_CSV_OUT, index=False)
    print(f"\nSaved CSV tables:\n  {SPW_CSV_OUT}\n  {MOUS_CSV_OUT}")

    # LaTeX (for your thesis / paper)
    with open(SPW_TEX_OUT, "w") as f:
        f.write(
            spw_df.to_latex(
                index=False,
                float_format="%.6f",
                columns=[
                    "Source", "Project", "ALMA_Band",
                    "min_freq_GHz", "max_freq_GHz",
                    "SiO_transition", "SiO_freq_GHz",
                    "ang_res_arcsec", "MOUS_ID",
                ],
            )
        )

    with open(MOUS_TEX_OUT, "w") as f:
        f.write(
            mous_df.to_latex(
                index=False,
                float_format="%.6f",
                columns=[
                    "MOUS_ID", "Source", "Project", "ALMA_Band",
                    "SiO_transitions", "SiO_freqs_GHz",
                ],
            )
        )

    print(f"Saved LaTeX tables:\n  {SPW_TEX_OUT}\n  {MOUS_TEX_OUT}")


# --------------------------------------------------------------------
# 6. Main
# --------------------------------------------------------------------

def main() -> None:
    print("\n==== Manara SiO(v=0) Archive Search ====\n")

    # Step 1: refine Manara sample
    df_alma = load_and_refine_manara()

    # Step 2: ALMA ObsCore table (load cache or query)
    obs = load_or_query_obscore(df_alma)

    # Step 3: SiO frequency filtering
    sio_spw = find_sio_matches_per_spw(obs)
    if sio_spw.empty:
        print("\nNo SiO-positive SPWs found; nothing more to do.")
        return

    # Step 4: per-MOUS table
    mous_table = build_per_mous_table(sio_spw)

    # Step 5: save outputs
    save_tables(sio_spw, mous_table)


if __name__ == "__main__":
    main()

