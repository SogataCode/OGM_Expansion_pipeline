OGM REPEAT EXPANSION DISTANCE PIPELINE
======================================

OVERVIEW
--------
One-command workflow to measure Bionano OGM single-molecule inter-label
lengths for a reference locus, deduplicate molecule IDs across multiple
assembly contigs, fit the original DBSCAN + Gaussian mixture algorithm, and
save plots of distance and reference-relative expansion length.

ENTRYPOINT
----------
bash expansion_distribution.sh \
  --input sample/folder \
  --output output/dir \
  --prefix saphyr186 \
  --locus chr4:3934999-3953000

--prefix is optional:
bash expansion_distribution.sh --input sample/folder --output output/dir \
  --locus chr4:3934999-3953000

INPUT DISCOVERY
---------------
--input accepts (a) the Bionano 'contigs' directory itself or (b) any parent
containing exactly one suitable 'contigs' directory at arbitrary depth.
The expected internal layout is:
  .../output/contigs/
    annotation/exp_refineFinal1_merged.xmap
    annotation/exp_refineFinal1_merged_q.cmap
    annotation/exp_refineFinal1_merged_r.cmap
    exp_refineFinal1/alignmol/merge/exp_refineFinal1_contigNNN.xmap
    exp_refineFinal1/alignmol/merge/exp_refineFinal1_contigNNN_q.cmap
    exp_refineFinal1/alignmol/merge/exp_refineFinal1_contigNNN_r.cmap

The tool stops if zero or multiple suitable contigs directories are found.
The locus chromosome number must match RefContigID in the Bionano reference
CMAP, and locus coordinates must use the same reference genome build.

OUTPUT DIRECTORY
----------------
--output DIR --prefix SAMPLE --locus chr4:START-END creates:
  DIR/SAMPLE_chr4_START_END/
Without --prefix it creates DIR/chr4_START_END/

Files:
  reference_interval.tsv  Reference label IDs, positions, reference distance
  manifest.tsv            Selected contigs and matched assembly label IDs
  contigNNN.distances.tsv Individual molecule distance from each contig
  all_distances.tsv       All contig measurements, including duplicates
  duplicated_molecules.tsv  Duplicate molecule records for audit
  merged_distances.tsv    One row per unique molecule ID
  gmm/allele_summary.tsv  Mean/SD/weight and relative expansion bp
  gmm/model_summary.tsv   Model BIC/AIC and counts
  gmm/molecule_assignment.tsv Molecule component, posterior, outlier status
  gmm/distance_distribution_density.png
  gmm/distance_distribution_count.png
  gmm/expansion_distribution_density.png
  gmm/expansion_distribution_count.png

METRICS
-------
Inter-label distance: observed molecule distance between mapped flanking
reference label IDs (in bp, truncated to integer for compatibility with
original aln_modify.py output).
Reference-relative expansion (bp) = molecule distance - reference inter-label
distance. For a GMM component, expansion_bp = Gaussian mean - reference
distance. Negative values are valid and are not clipped to zero.
This measurement is a difference relative to the reference and is NOT the
absolute physical repeat tract length, nor evidence of motif identity.

HOW IT WORKS
------------
1. Read annotation XMAP/CMAP to identify every assembly contig whose alignment
   spans the requested locus. XMAP Alignment pairs are (Reference SiteID,
   Query SiteID). Require the same reference flanking labels in all contigs;
   do not pool measurements of different reference intervals.
2. Cross-reference each assembly label to alignmol reference CMAP by SiteID.
   Coordinate values are taken from the alignmol CMAP, since annotation and
   alignmol CMAP values may differ due to map processing.
3. For each contig, read alignmol XMAP, require exactly one molecule label
   corresponding to each flanking reference label, then fetch those exact
   molecule SiteIDs from query CMAP and compute absolute position difference.
4. Merge across contigs by QryContigID. If duplicate IDs have identical
   integer distance, keep one record; if they conflict, write
   conflicting_molecules.tsv and stop (no arbitrary tie breaking).
5. Call extract_gaussian_alleles() supplied by user unchanged. That original
   function runs DBSCAN and selects a GMM among >=2 components; a fitted
   component must NOT be assumed to correspond to a validated biological
   allele. DBSCAN outliers are excluded from fitting, retained in assignments,
   and shown as red crosses on graphs.
6. Draw both normalized-density and molecule-count histograms with weighted
   Gaussian components and combined Gaussian mixture. Count curves are scaled
   by retained-molecule count times histogram bin width.

REQUIREMENTS / INSTALL
----------------------
Python 3.10+ recommended, POSIX bash, packages in requirements.txt:
  python -m pip install -r requirements.txt
Default Python can be overridden with: PYTHON=/path/to/python bash ...

VALIDATION AND LIMITATIONS
--------------------------
- Tested against synthetic miniature XMAP/CMAP fixtures. Actual Bionano Solve
  output naming, schemas and labeling may vary; verify the manifest for your
  first real sample, and inspect selected molecules for alignment anomalies.
- Expects Saphyr De novo pipeline exp_refineFinal1 naming/layout.
- Source molecule IDs assumed unique within the one Bionano analysis/run.
  Do NOT combine unrelated runs without incorporating the run ID.
- A single contig may contain duplicate or ambiguous mappings; these are
  excluded rather than guessed. There is no automatic allele phasing.
- A shared reference flanking pair is essential for pooled absolute lengths.
- The original GMM implementation may fail when KneeLocator returns None,
  DBSCAN removes too many molecules, or clusters are too small; original
  statistical behavior is preserved intentionally.
- The original GMM starts from at least 2 Gaussian components, so 2 components
  is NOT proof of two resolvable alleles. Review raw histograms and model fit.
- Intermediate results are overwritten if the command is rerun with the same
  output directory and locus; use a unique --prefix or output folder.
- The pipeline does not upload data; all processing happens locally.

TROUBLESHOOTING
---------------
No matching contigs folder: point --input to a directory containing the
pipeline output/contigs tree, and check annotation and alignmol/merge exist.
Multiple contigs folders: narrow the --input directory.
No shared flanking label: inspect annotation XMAP and repeat reference build;
not all assembly contigs necessarily have equivalent label mappings.
Conflicting molecules: inspect conflicting_molecules.tsv and source XMAPs.
GMM error: inspect merged_distances.tsv and original DBSCAN/GMM assumptions.

SCRIPTS
-------
expansion_distribution.sh    Main CLI entrypoint
ogm_pipeline.py             Contig discovery, locus resolution, distances, dedup
finish_gmm.py               Connects shell entrypoint to GMM visualization
run_gmm.py                  Plots and GMM tables using existing function
extract_gaussian_alleles.py User-supplied original statistical algorithm
