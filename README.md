# Revio_get_readstat — Per‑well HiFi Read Statistics from PacBio Revio BAMs ✨🧬

A lightweight Bash workflow that aggregates per‑well HiFi read statistics from PacBio Revio runs. It scans each well folder for `hifi_reads/*.hifi_reads.bam`, extracts per‑read fields via `samtools`, derives helper files (lengths and per‑read Phred Q), and writes compact, memory‑efficient per‑well summaries using `gawk`. 🚀

---

## ✨ Highlights

- 🔎 Auto‑discovers wells by pattern `*/hifi_reads/*.hifi_reads.bam` from the run folder
- 🧮 Computes key HiFi metrics per well: total bases/reads, average length, N50, median per‑read Q, mean passes
- ⚙️ Fast & lean: streaming with `samtools`, histogram summaries in `gawk` (low RAM)
- 🧱 Robust filtering: primary alignments only (`-F 0x900`); requires `rq`; `np` optional
- 🗂️ Produces clean artifacts in `out_by_well/`


---

## 🗺️ Overview

The pipeline runs in three stages:

1) Per‑well extraction (append across multiple BAMs of the same well)  
   - Filter: `-F 0x900` (exclude secondary `0x100` and supplementary `0x800`)  
   - Require `rq` (per‑read accuracy); optional `np` (passes; default 0 if missing)  
   - Output per read: `<length>\t<rq>\t<np>` → `out_by_well/<WELL>.readstats`

2) Helper files  
   - `<WELL>_HiFi_Length.txt`: one read length per line (from column 1 of `.readstats`)  
   - `<WELL>.hifi_reads_Phred.txt`: per‑read Phred Q computed from `rq`

3) Per‑well summary via histogram‑based `gawk`  
   - `<WELL>_stat.txt` with header:  
     `Polymerase_Read_Bases, Polymerase_Reads, Polymerase_N50, Ave_Polymerase_Length, Subreads_Bases, Subreads_Reads, Ave_Subreads_Length, Subreads_N50, HiFi_Bases, HiFi_Reads, Ave_HiFi_Length, HiFi_N50, HiFi_MeanPass, HiFi_Quality`  
   - Polymerase/Subreads are placeholders here; only HiFi metrics are computed.

---

## 🧩 File discovery and naming

- Run the script from the Revio “run folder” that contains per‑well subfolders.  
  Pattern searched: `*/hifi_reads/*.hifi_reads.bam`  
- The well name is the first path component before `/hifi_reads/` (i.e., the directory name).  
- All BAMs for a well are appended into a single `<WELL>.readstats` file.

---

## 🧰 Requirements

- Linux with Bash 4+
- `samtools` (e.g., 1.19.x; configurable via `SAMTOOLS`)
- `gawk`
- Coreutils (`cut`, `mkdir`, etc.)

Optional/HPC:
- `NSLOTS` (SGE) to auto‑set thread count

---

## ⚙️ Setup

- Ensure `samtools` and `gawk` are installed.  
- Either set an absolute path, `SAMTOOLS=/abs/path/to/samtools`, or rely on `samtools` in `PATH`.

---

## 🚀 Quick start

From the run folder that contains well directories with `hifi_reads` subfolders:

```bash
# Example structure
# r84135_20251003_081955/
# ├── WELL_A/hifi_reads/xxx.hifi_reads.bam
# ├── WELL_B/hifi_reads/yyy.hifi_reads.bam
# └── ...

# Threads: NSLOTS > T > 16
bash Revio_get_readstat.sh
```

Environment variables:
- Threads: `T="${NSLOTS:-${T:-16}}"` → used by `samtools view -@ T`  
- `SAMTOOLS`: optional absolute path (fallback: `samtools` in PATH)

Artifacts are written under `out_by_well/`:

```
out_by_well/
├── WELL_A.readstats
├── WELL_A_HiFi_Length.txt
├── WELL_A.hifi_reads_Phred.txt
├── WELL_A_stat.txt
└── ...
```

---

## 📦 Outputs

- `<WELL>.readstats` (TSV)  
  - Col1: read length (bp)  
  - Col2: `rq` (per‑read accuracy, 0–1)  
  - Col3: `np` (passes; `0` if missing)

- `<WELL>_HiFi_Length.txt`  
  - One integer length per line (from `.readstats` col1)

- `<WELL>.hifi_reads_Phred.txt`  
  - One floating‑point Phred value per line (see formula below)

- `<WELL>_stat.txt` (two rows: header + values)  
  - Header:  
    `Polymerase_Read_Bases, Polymerase_Reads, Polymerase_N50, Ave_Polymerase_Length, Subreads_Bases, Subreads_Reads, Ave_Subreads_Length, Subreads_N50, HiFi_Bases, HiFi_Reads, Ave_HiFi_Length, HiFi_N50, HiFi_MeanPass, HiFi_Quality`  
  - Values: Polymerase/Subreads placeholders; HiFi metrics filled

---

## 🧮 Calculations and formulas

### 1) Per‑read Phred Q from `rq`

Let `acc` be the per‑read accuracy from `rq`. Values that reach 1.0 are clamped to avoid $\log(0)$.

- Clamp:

$$
\mathrm{acc}' =
\begin{cases}
1-\varepsilon, & \mathrm{acc} \ge 1 \\
\mathrm{acc}, & \text{otherwise}
\end{cases}
\quad (\varepsilon = 10^{-12})
$$

- Per‑read accuracy Phred:

$$
Q = -10\,\log_{10}\bigl(1-\mathrm{acc}'\bigr)
$$

- Rounding for the per‑well median histogram:

$$
Q_i = \max\bigl(0,\, \mathrm{round}(Q)\bigr)
$$

### 2) HiFi length metrics

Totals:

$$
B = \sum_{r=1}^{R} L_r,\quad R = \text{reads with }\mathtt{rq}\text{ present}
$$

Average HiFi length:

$$
\mathrm{AveHiFiLength} = \left\lfloor \dfrac{\sum_{r=1}^{R} L_r}{R} \right\rfloor
$$

### 4) HiFi mean passes

Only reads with `np` present contribute to the mean (rounded to nearest integer):

$$
\mathrm{HiFi\_MeanPass} = \mathrm{round}\!\left( \dfrac{1}{n_P} \sum_{i=1}^{n_P} P_i \right)
$$

---

## 🧱 Filters and tag handling

- `samtools view -F 0x900` → primary alignments only (drop secondary `0x100`, supplementary `0x800`)  
- Records must carry `rq` to enter `.readstats`  
- `np` is optional; written as `0` per read if missing, but the mean uses only reads where `np` exists

---

## ⚡ Performance notes

- Streaming decode with `samtools -@ T`  
- Histogram‑based `gawk` keeps memory small  
- `shopt -s nullglob` avoids literal globs on empty matches

---

## ⚠️ Limitations and caveats

- Polymerase/Subreads columns in `_stat.txt` are placeholders in this workflow  
- Only reads with `rq` are counted; missing‑`rq` reads are ignored  
- `Ave_HiFi_Length` is integer‑truncated (not rounded)  
- If a well has no reads, `_stat.txt` contains zeros/hyphens and `Q0` by design  
- The search pattern is fixed to `*/hifi_reads/*.hifi_reads.bam` from the current directory

---

## 🛟 Troubleshooting

- “No inputs processed” → check you are in the correct run folder and the pattern exists  
- “Command not found” → ensure `samtools` and `gawk` are installed and on `PATH` (or set `SAMTOOLS`)  
- Empty `_stat.txt` → confirm any records had `rq`  
- Extremely high Q → `rq=1.0` is clamped with $\varepsilon=10^{-12}$ to avoid $\log(0)$

---

## 🔁 Reproducibility tips

- Record `samtools --version` and the exact script commit  
- Pin your environment (modules/containers) across runs  
- Keep `_stat.txt` alongside raw BAMs for traceability

---

## 🙏 Acknowledgments

- PacBio Revio data model; `samtools`, `gawk`, and GNU coreutils

