# Norovirus sequence and metadata retrieval

A reproducible protocol for assembling a norovirus VP1/genome sequence set from NCBI GenBank, and for selecting a single genotype lineage (here GII.4) from it by local BLAST against a reference.

The protocol is written so that a reader with no prior knowledge of the project can regenerate the dataset from scratch. Every step lists its inputs, the exact command, and the expected output. Sequence counts are given for the retrieval date recorded in `logs/retrieval_date.txt`; GenBank grows continuously, so your counts will differ.

---

## Contents
Part A:
1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Directory layout](#3-directory-layout)
4. [Step 1 — Define the GenBank queries](#4-step-1--define-the-genbank-queries)
5. [Step 2 — Download sequences and metadata](#5-step-2--download-sequences-and-metadata)
6. [Step 3 — Merge and deduplicate](#6-step-3--merge-and-deduplicate)
7. [Step 4 — Build a local BLAST database](#7-step-4--build-a-local-blast-database)
8. [Step 5 — Select the GII.4 lineage](#8-step-5--select-the-gii4-lineage)
9. [Step 6 — Extract the selected sequences](#9-step-6--extract-the-selected-sequences)
10. [Step 7 — Parse metadata and flag chronic infection](#10-step-7--parse-metadata-and-flag-chronic-infection)
11. [Step 8 — Quality control](#11-step-8--quality-control)
12. [Outputs](#12-outputs)
13. [Adapting to another genotype](#13-adapting-to-another-genotype)
Part B-E:
    
F. [Known limitations](#14-known-limitations)
G. [Citation and licence](#15-citation-and-licence)

---
### Part A

## 1. Overview

```
           NCBI GenBank (Entrez)
                    │
        ┌───────────┴───────────┐
        │                       │
  Query A: VP1 / ORF2 /    Query B: complete
  capsid, 1.6–10 kb        genomes, 7–8 kb
        │                       │
        └───────────┬───────────┘
                    │  merge, deduplicate by accession
                    ▼
          all_norovirus.fasta  +  all_norovirus_metadata.tsv
                    │
                    │  makeblastdb
                    ▼
            local nucleotide BLAST database
                    │
                    │  blastn, query = GII.4 reference (MK073894)
                    │  keep hits with identity >= 80% to the reference
                    ▼
            gii4_accessions.txt
                    │
                    │  seqkit grep
                    ▼
   gii4.fasta  +  gii4_metadata.tsv  (+ chronic-infection flag)
```

Two Entrez queries are used rather than one. The first targets annotated VP1/ORF2/capsid records across a wide length range; the second targets complete genomes, which are sometimes annotated in a way the first query misses. The union is taken so that no complete genome is lost to annotation inconsistency.

Genotype selection is done by similarity to a reference rather than by the GenBank genotype field, because that field is free text, frequently absent, and not reliably assigned. The reference used is **MK073894**, the same GII.4 reference used by Nextstrain, so the resulting set is comparable to Nextstrain norovirus builds.

---

## 2. Requirements

| Tool | Version used | Purpose |
| --- | --- | --- |
| NCBI Entrez Direct (EDirect) | 22.x | query and download from GenBank |
| NCBI BLAST+ | 2.15.0 | `makeblastdb`, `blastn` |
| seqkit | 2.8.x | FASTA subsetting and statistics |
| Python | 3.12 | metadata parsing |
| pandas | 2.x | metadata parsing |

Install with conda:

```bash
conda create -n noro -c conda-forge -c bioconda \
    entrez-direct blast seqkit python=3.12 pandas
conda activate noro
```

Register an NCBI API key (free, raises the rate limit from 3 to 10 requests/second) and export it before running any download step:

```bash
export NCBI_API_KEY=your_key_here
```

No compilation or project-specific installation is required.

---

## 3. Directory layout

```
.
├── data/
│   ├── raw/          # untouched Entrez downloads
│   ├── interim/      # merged, deduplicated
│   └── final/        # analysis-ready sequence sets
├── db/               # local BLAST databases (gitignored)
├── ref/              # reference sequences
├── scripts/
├── logs/
└── README.md
```

Create it with `mkdir -p data/{raw,interim,final} db ref scripts logs`.

`db/` and `data/raw/` should be in `.gitignore`; they are large and fully regenerable from the commands below.

---

## 4. Step 1 — Define the GenBank queries

### Query A — VP1 / ORF2 / capsid records

Human-readable form:

```
norovirus[Organism] AND (VP1 OR ORF2 OR capsid protein) AND (cds OR genome)
```

Fully qualified form, as submitted (length filter 1,600–10,000 nt):

```
"Norovirus"[Organism] AND (VP1[All Fields] OR ORF2[All Fields] OR capsid protein[All Fields])
AND (cds[All Fields] OR genome[All Fields])
AND (viruses[filter] AND biomol_genomic[PROP] AND ("1600"[SLEN] : "10000"[SLEN]))
```

Returned **10,112 sequences**.

The lower bound of 1,600 nt is set just below the length of a complete VP1 ORF (~1,620 nt in GII.4), so that records containing a full capsid gene are retained and short partial-capsid fragments are excluded. The upper bound of 10,000 nt comfortably exceeds the ~7.5 kb norovirus genome and removes chimeric or vector-containing records.

> The prose description of this search mentions a 1,600–8,000 nt filter; the query actually submitted used 1,600–10,000. The qualified string above is authoritative.

### Query B — complete genomes

Human-readable form:

```
"Norovirus"[Organism] AND (complete[Title] AND genome[Title])
```

Fully qualified form, as submitted:

```
"Norovirus"[Organism] AND (complete[Title] AND genome[Title]) AND ("7000"[SLEN] : "8000"[SLEN])
```

Returned **1,778 sequences**.

This query is a safety net. Some complete genomes lack the `VP1`, `ORF2` or `capsid` strings in the fields Query A searches, and would otherwise be dropped.

---

## 5. Step 2 — Download sequences and metadata

Run both queries through EDirect. Record the retrieval date — GenBank is not versioned, and the date is the only thing that makes the download reproducible.

```bash
date -u +%Y-%m-%d > logs/retrieval_date.txt
```

**Query A:**

```bash
esearch -db nucleotide -query '"Norovirus"[Organism] AND (VP1[All Fields] OR ORF2[All Fields] OR capsid protein[All Fields]) AND (cds[All Fields] OR genome[All Fields]) AND (viruses[filter] AND biomol_genomic[PROP] AND ("1600"[SLEN] : "10000"[SLEN]))' \
  | efetch -format fasta > data/raw/queryA_vp1.fasta
```

**Query B:**

```bash
esearch -db nucleotide -query '"Norovirus"[Organism] AND (complete[Title] AND genome[Title]) AND ("7000"[SLEN] : "8000"[SLEN])' \
  | efetch -format fasta > data/raw/queryB_genomes.fasta
```

**Metadata for both**, as a tab-separated table. Fields are pulled from the GenBank source feature, so missing values are expected and are emitted as empty cells:

```bash
for q in A B; do
  case $q in
    A) QUERY='"Norovirus"[Organism] AND (VP1[All Fields] OR ORF2[All Fields] OR capsid protein[All Fields]) AND (cds[All Fields] OR genome[All Fields]) AND (viruses[filter] AND biomol_genomic[PROP] AND ("1600"[SLEN] : "10000"[SLEN]))' ;;
    B) QUERY='"Norovirus"[Organism] AND (complete[Title] AND genome[Title]) AND ("7000"[SLEN] : "8000"[SLEN])' ;;
  esac
  esearch -db nucleotide -query "$QUERY" \
    | efetch -format gpc \
    | xtract -insd INSDSeq_accession-version INSDSeq_length INSDSeq_create-date \
             INSDSeq_definition source collection_date country host strain isolate note \
    > data/raw/query${q}_metadata.tsv
done
```

Add a header row once, then concatenate:

```bash
printf 'accession\tlength\tcreate_date\tdefinition\tsource\tcollection_date\tcountry\thost\tstrain\tisolate\tnote\n' \
  > data/interim/all_norovirus_metadata.tsv
cat data/raw/queryA_metadata.tsv data/raw/queryB_metadata.tsv \
  >> data/interim/all_norovirus_metadata.tsv
```

Downloads of this size take roughly 15–40 minutes depending on the connection. If `efetch` is interrupted, re-run it — EDirect does not resume partial downloads, and a truncated FASTA will fail silently downstream.

---

## 6. Step 3 — Merge and deduplicate

The two queries overlap substantially; complete genomes satisfying Query A appear in both. Deduplicate by accession, keeping the first occurrence:

```bash
cat data/raw/queryA_vp1.fasta data/raw/queryB_genomes.fasta \
  | seqkit rmdup -n \
  > data/interim/all_norovirus.fasta

seqkit stats data/interim/all_norovirus.fasta | tee logs/all_norovirus_stats.txt
```

`-n` deduplicates on the full sequence name. Use `-s` instead to deduplicate on sequence content, which additionally collapses identical sequences submitted under different accessions — only do this if downstream analysis assumes unique sequences, since it discards genuine independent isolates that happen to be identical.

Deduplicate the metadata table on the same key:

```bash
awk -F'\t' 'NR==1 || !seen[$1]++' data/interim/all_norovirus_metadata.tsv \
  > data/interim/all_norovirus_metadata.dedup.tsv
mv data/interim/all_norovirus_metadata.dedup.tsv data/interim/all_norovirus_metadata.tsv
```

---

## 7. Step 4 — Build a local BLAST database

Fetch the reference. **MK073894** is the Nextstrain norovirus GII.4 reference:

```bash
efetch -db nucleotide -id MK073894 -format fasta > ref/noro_gii4_refseq.fasta
```

Build the database from the merged set:

```bash
makeblastdb \
  -in data/interim/all_norovirus.fasta \
  -title "norovirus GenBank all" \
  -dbtype nucl \
  -parse_seqids \
  -out db/noro_db_all
```

`-parse_seqids` is worth including; it lets you pull sequences straight out of the database later with `blastdbcmd`.

---

## 8. Step 5 — Select the GII.4 lineage

BLAST the reference against the database. The reference is the **query**; every GenBank record is a **subject**.

Save the full hit table once, then filter it. Columns are, in order: `qaccver`, `saccver`, `pident`, `evalue`, `qcovs`.

```bash
blastn \
  -db db/noro_db_all \
  -query ref/noro_gii4_refseq.fasta \
  -outfmt "6 qaccver saccver pident evalue qcovs" \
  -max_target_seqs 70000 \
  -num_threads 8 \
  > data/interim/gii4_blast_hits.tsv
```

```bash
awk -F'\t' '$3 >= 80 {print $2}' data/interim/gii4_blast_hits.tsv \
  | sort -u \
  > data/final/gii4_accessions.txt

wc -l data/final/gii4_accessions.txt
```

The selection criterion is **percent identity ≥ 80** across the aligned region (`pident`, column 3). Writing the hit table to disk before filtering is the only change from the one-liner: the same table then serves the threshold sweep and the divergent-tail check below, without re-running the search each time.

`-max_target_seqs 70000` is well above the database size, so the search is not truncated and every record that aligns at all is reported. That is the right setting here — see the note below on why a low value breaks the complement check.

Two points on the threshold itself:

- **80% nucleotide identity is permissive for a GII.4 set.** GII genogroup members routinely align to a GII.4 reference above 80%, so the set will contain non-GII.4 GII genotypes. If the intent is a GII.4 lineage specifically, check how the count behaves as the threshold moves — a genotype boundary should show up as a shoulder in the distribution rather than a smooth slope:

  ```bash
  for t in 80 85 88 90 92 95; do
    n=$(awk -F'\t' -v t=$t '$3 >= t {print $2}' data/interim/gii4_blast_hits.tsv | sort -u | wc -l)
    echo -e "${t}\t${n}"
  done
  ```

- **No coverage filter is applied.** A short fragment matching 300 nt of the reference at 95% identity passes `pident >= 80` as readily as a full-length genome does. Add `&& $5 >= 90` if the downstream analysis assumes near-complete VP1 coverage; keep it off if partial records are wanted deliberately.

The original output file was named `*.fasta` but contains a list of accessions, one per line. It is renamed `gii4_accessions.txt` here to avoid confusion.

### Inspecting the divergent tail

The complement of the selection — hits that align to the reference but fall below the threshold. These are the records the GII.4 set deliberately excludes: other genogroups caught by the broad organism query, recombinants, and mis-annotated entries.

```bash
awk -F'\t' '$3 < 80' data/interim/gii4_blast_hits.tsv \
  > data/interim/gii4_low_identity.tsv     # identity < 80: the excluded ones

wc -l data/interim/gii4_low_identity.tsv
```

Run against the saved hit table, this is a genuine check: the two filters partition the table, so the row counts must sum to its total.

```bash
total=$(wc -l < data/interim/gii4_blast_hits.tsv)
kept=$(awk -F'\t' '$3 >= 80' data/interim/gii4_blast_hits.tsv | wc -l)
low=$(wc -l < data/interim/gii4_low_identity.tsv)
echo "$total rows = $kept kept + $low excluded"
test "$total" -eq "$(( kept + low ))" && echo OK || echo "MISMATCH"
```

**Do not re-run `blastn` with a low `-max_target_seqs` for this check.** BLAST applies that cutoff during the search and keeps the best-scoring subjects, so the rows a low cap discards are exactly the low-identity rows the check is looking for. A second search capped at 1,400 against a database of ~11,000 records returns a short or empty `< 80` list no matter what the selection filter did — the result is consistent with a working filter and with a truncated search equally, so it cannot distinguish them. Filtering the one saved table twice avoids the problem entirely and costs nothing.

**Records with no hit at all are invisible here.** A sequence too divergent to produce any alignment never appears in BLAST output, so no identity filter can catch it. Get those by difference:

```bash
cut -f2 data/interim/gii4_blast_hits.tsv | sort -u > /tmp/hit_accessions.txt
grep '^>' data/interim/all_norovirus.fasta | cut -d' ' -f1 | tr -d '>' | sort -u \
  | comm -23 - /tmp/hit_accessions.txt \
  > data/interim/no_blast_hit.txt
```

Between `gii4_low_identity.tsv` and `no_blast_hit.txt` you have the full set of records the GII.4 reference does not describe well. Spot-check a handful of each against their GenBank `definition` field before deciding what the cause is — in practice the low-identity tail is mostly GI and non-GII.4 GII genotypes picked up by the broad organism query, while the no-hit set is more often short or heavily ambiguous records.

---

## 9. Step 6 — Extract the selected sequences

```bash
seqkit grep -f data/final/gii4_accessions.txt data/interim/all_norovirus.fasta \
  > data/final/gii4.fasta

seqkit stats data/final/gii4.fasta | tee logs/gii4_stats.txt
```

Subset the metadata to match:

```bash
head -1 data/interim/all_norovirus_metadata.tsv > data/final/gii4_metadata.tsv
grep -F -f data/final/gii4_accessions.txt data/interim/all_norovirus_metadata.tsv \
  >> data/final/gii4_metadata.tsv
```

Confirm the two agree:

```bash
test "$(grep -c '^>' data/final/gii4.fasta)" -eq "$(( $(wc -l < data/final/gii4_metadata.tsv) - 1 ))" \
  && echo "FASTA and metadata counts agree" \
  || echo "MISMATCH — investigate before proceeding"
```

---

## 10. Step 7 — Parse metadata and flag chronic infection

Records relating to chronic, persistent or within-host infection are identified by keyword search across the free-text metadata fields (`definition`, `note`, `strain`, `isolate`, `source`).

Keywords used:

```
chronic, immuno, immune, immunity, long-term, persistent,
within-host, intra-host, shedder
```

```bash
python scripts/flag_chronic.py \
  --metadata data/final/gii4_metadata.tsv \
  --out data/final/gii4_metadata_flagged.tsv
```

`scripts/flag_chronic.py`:

```python
#!/usr/bin/env python3
"""Flag GenBank records whose free-text fields suggest chronic or
within-host norovirus infection. Keyword matching only — the output is a
candidate list for manual review, not a classification."""

import argparse
import pandas as pd

KEYWORDS = [
    "chronic", "immuno", "immune", "immunity", "long-term",
    "persistent", "within-host", "intra-host", "shedder",
]
TEXT_FIELDS = ["definition", "note", "strain", "isolate", "source"]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--metadata", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    df = pd.read_csv(args.metadata, sep="\t", dtype=str).fillna("")
    present = [c for c in TEXT_FIELDS if c in df.columns]
    blob = df[present].agg(" ".join, axis=1).str.lower()

    for kw in KEYWORDS:
        df[f"kw_{kw.replace('-', '_')}"] = blob.str.contains(kw, regex=False)

    df["chronic_candidate"] = blob.apply(
        lambda s: any(kw in s for kw in KEYWORDS)
    )
    df.to_csv(args.out, sep="\t", index=False)

    n = int(df["chronic_candidate"].sum())
    print(f"{n} of {len(df)} records flagged as chronic candidates")


if __name__ == "__main__":
    main()
```

Keyword matching is deliberately permissive. `immuno` matches `immunocompromised` but also `immunoassay`; `persistent` appears in unrelated phrasing. Treat the flagged subset as a shortlist to read through by hand, and keep a record of which accessions survive that review in `data/final/chronic_confirmed.txt`. Many chronic-infection series are described only in the associated publication, not in the GenBank record, so this step has poor recall as well as poor precision.

---

## 11. Step 8 — Quality control

Before using `gii4.fasta` downstream:

```bash
# Length distribution — a long tail below ~1.5 kb suggests the length
# filter is not doing what you expect
seqkit fx2tab -nl data/final/gii4.fasta | awk '{print $NF}' | sort -n | uniq -c

# Ambiguous bases — sequences heavy in N will distort alignment and
# phylogenetic inference
seqkit fx2tab -n -B N data/final/gii4.fasta | awk -F'\t' '$NF > 5'

# Duplicate accessions
grep '^>' data/final/gii4.fasta | cut -d' ' -f1 | sort | uniq -d
```

A reasonable default is to drop sequences with more than 5% ambiguous bases:

```bash
seqkit fx2tab -n -B N data/final/gii4.fasta \
  | awk -F'\t' '$NF <= 5 {print $1}' \
  | cut -d' ' -f1 \
  > data/final/gii4_pass_qc.txt
seqkit grep -f data/final/gii4_pass_qc.txt data/final/gii4.fasta \
  > data/final/gii4_qc.fasta
```

Check the results with norovirus typing tools: 1) https://www.genomedetective.com/app/typingtool/nov/ 2) https://mpf.rivm.nl/mpf/typingtool/norovirus/

---

## 12. Outputs

| File | Contents |
| --- | --- |
| `data/interim/all_norovirus.fasta` | merged, deduplicated norovirus sequences from both queries |
| `data/interim/all_norovirus_metadata.tsv` | GenBank metadata for the above |
| `data/interim/gii4_blast_hits.tsv` | raw BLAST output, reference vs. all |
| `data/interim/gii4_low_identity.tsv` | hits below the identity threshold, for inspection |
| `data/interim/no_blast_hit.txt` | records producing no alignment to the reference |
| `data/final/gii4_accessions.txt` | accessions passing the identity threshold |
| `data/final/gii4.fasta` | GII.4 sequence set |
| `data/final/gii4_metadata.tsv` | metadata for the GII.4 set |
| `data/final/gii4_metadata_flagged.tsv` | the above, with chronic-infection keyword flags |
| `data/final/gii4_qc.fasta` | GII.4 set after ambiguity filtering |
| `logs/retrieval_date.txt` | GenBank retrieval date — cite this |

---

## 13. Adapting to another genotype

Only Step 5 is genotype-specific. To build a different lineage set, replace the reference and re-run Steps 4–6:

| Genotype | Suggested reference |
| --- | --- |
| GII.4 | MK073894 (Nextstrain reference) |
| GII.2 | choose a RefSeq or Nextstrain reference for the target genotype |
| GI.1 | M87661 (Norwalk virus prototype) |

Identity thresholds are not transferable between genotypes. Norovirus genotype boundaries are conventionally defined on VP1 amino-acid distance, not whole-genome nucleotide identity, so an 80% nucleotide cutoff is a pragmatic proxy rather than a genotype definition. Where genotype assignment matters, confirm the set with a dedicated typing tool — the RIVM Human Calicivirus Typing Tool or a Nextclade norovirus dataset — and report the agreement between the two approaches.

There's a snakemake pipeline to carry out this step too.


### Part B

Alignment and Phylogeny

Mafft based alignment and then iqtree based ML phylogeny. We also used beast to generate discrete phylodynamics analysis.


### Part C

Mutation profile


### Part D
Deep node finding and antigenicity prediction of the ancestral sequences

### Part E
Antigenic cartography


---

## Known limitations

- **GenBank is a moving target.** Record counts are not reproducible without the retrieval date. Cite the date in `logs/retrieval_date.txt` in any write-up.
- **BLAST similarity is not genotyping.** A single-reference identity threshold will include recombinants whose capsid is not GII.4 and exclude divergent but genuine GII.4 sequences. Recombination between ORF1 and ORF2 is common in norovirus, so whole-genome identity to a GII.4 reference can misrepresent the capsid genotype. Where the capsid genotype is what matters, BLAST against the VP1 region alone.
- **Metadata is sparse and inconsistent.** `collection_date` is missing or year-only for a large fraction of records, and `country` is free text. Both need normalisation before use in a time-scaled analysis.
- **The chronic-infection flag is a keyword heuristic**, with both false positives and false negatives, and is not a substitute for reading the source publications.
- **No ORF-level annotation is performed here.** This protocol delivers whole records; extracting the VP1 ORF itself requires a separate alignment or annotation step.

---

## Citation and licence

If you use this protocol, please cite the repository and state the GenBank retrieval date.

Underlying tools to cite: NCBI Entrez Direct (Kans, *Entrez Direct: E-utilities on the Unix Command Line*), BLAST+ (Camacho *et al.* 2009, *BMC Bioinformatics*), seqkit (Shen *et al.* 2016, *PLoS ONE*).

Sequence data are from NCBI GenBank and are subject to NCBI's terms of use. Code in this repository is released under the MIT licence unless stated otherwise.
