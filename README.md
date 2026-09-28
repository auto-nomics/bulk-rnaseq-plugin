# autonomics bulk-rnaseq plugin

Manifest plugin family for the bulk RNA-seq wrapper nodes: `limma_voom`
(limma + voom + eBayes differential expression) and `wgcna` (the official
WGCNA pipeline). One plugin — the two nodes share the pinned image and bind
no panels.

Layered on top of the published `autonomics/deseq2:1.50.2` image
(`sha256:8b2e2a78d87293e6cae6dbed2e283dff1dd8461a7f993cb347ca9810698b2b3e`).
The base image already ships R 4.5.3 and the pinned CRAN snapshot used by
the autonomics build pipeline; this stage adds:

- `limma` 3.66.0 (Bioconductor 3.22) — voom + eBayes + contrasts
- `edgeR` 4.8.2 (Bioconductor 3.22) — TMM normalization + GLM
- `dynamicTreeCut`, `preprocessCore`, `impute` (Bioconductor 3.22) +
  `fastcluster`, `matrixStats`, `doParallel`, `foreach`, `Rcpp`, `Hmisc`
  from CRAN
- `WGCNA` 1.74 (CRAN, sha256-verified)

The two runner scripts are baked into `/opt/autonomics/`:

- `limma_voom_runner.R` — consumes a `gene_id` + sample count matrix and a
  `sample_id` + covariate metadata table, fits
  `limma::voom + lmFit + eBayes` (with optional `makeContrasts`), applies
  TMM (default) or quantile normalization, and publishes
  `results.tsv`, `voom_weights.tsv`, `normalized_expression.tsv`,
  `contrast_summary.tsv` and a JSON run report.
- `wgcna_runner.R` — consumes a `gene_id` + sample expression matrix plus an
  optional sample-metadata table, runs `WGCNA::blockwiseModules` with signed
  adjacency, TOM, dynamic tree cut, module merge and signed kME, and
  publishes `soft_threshold.tsv`, `adjacency_stats.tsv`, `tom_stats.tsv`,
  `modules.tsv`, `module_eigengenes.tsv` and a JSON run report. TOM is
  kept inside the container — only summary statistics land on VFS, so the
  node never materializes a large TOM dataframe in the agent runtime.

## Layout

```text
bulk-rnaseq/
  manifest.toml          # node contract: two [[nodes]] (limma_voom, wgcna)
  Dockerfile             # image provenance (build + push via GHCR)
  limma_voom_runner.R    # baked into /opt/autonomics/
  wgcna_runner.R         # baked into /opt/autonomics/
  test_smoke.sh          # build + version smoke + push + digest
  README.md
```

## Build, smoke, publish

```bash
./test_smoke.sh
```

The script builds the image, verifies the exact R and package versions in an
isolated Podman container, pushes it to GHCR, and prints the immutable
manifest digest. The default build uses `--no-cache`; set `BUILD_FLAGS`
explicitly when reusing local layers during development.

`manifest.toml` pins the published digest:

Current immutable manifest digest (1.0.1):
`sha256:23071886af3864a0338242987980727753923a0a67f362aef7366d1277c39734`.

## Migration parity

Migrated from the hardcoded wrappers
`crates/node-bundles/nodes-io/src/limma_voom_container.rs` and
`wgcna_container.rs` (kind suffix `_container` dropped, artifact prefixes
renamed to match). The compiled container spec is byte-identical to the
legacy `container_spec` for image, command, outputs, env, resources and the
hardened defaults — asserted by
`crates/container-plugin/tests/bulk_rnaseq_migration.rs`.

Two legacy `validate()` rules the manifest DSL cannot express moved into
`limma_voom_runner.R` with the legacy error messages (this is the 1.0.1
image change; 1.0.0 lacked them):

- `contrast_levels` must contain exactly two entries for categorical outcome;
- covariates must be duplicate-free, non-reserved (outcome, `sample_id`),
  simple R identifiers.

The list-valued params (`contrast_levels`, `covariates`) travel as
comma-joined strings — byte-identical to the legacy `join(",")` env channel
for every legacy-legal input.

The legacy wrapper's optional second wgcna input (sample metadata) is not
declarable in the v0 DSL (no optional-port flag); the manifest declares only
the wired expression input and the runner treats a missing
`AUTONOMICS_INPUT1` as no metadata, exactly as before.

## Baselines

The real-Podman DAG baselines moved with the wrapper into the plugin
checkout; the recorded pasilla baselines remain:

- categorical `condition_treated_vs_untreated`, TMM, covariate `type`:
  14,599 input genes, 8,066 genes after expression filtering, 7 samples;
- continuous `duration`, quantile normalization, no contrast levels:
  14,599 input genes, 7,919 genes after expression filtering, 7 samples;
- both use limma 3.66.0, edgeR 4.8.2, R 4.5.3 and explicit BH correction.

The WGCNA baseline uses a deterministic 60-gene by 12-sample matrix with
optional sample metadata. It returns three modules after signed adjacency,
TOM, dynamic tree cut and merge; module eigengenes retain the metadata
columns. The runner supports Pearson and bicor correlations. The full TOM is
not published to VFS.

## Network-restricted environments

The build requires network access to download `limma`, `edgeR`, the WGCNA
closure, and Bioconductor's `preprocessCore` / `impute`. In offline
environments the same image can be built with the CRAN snapshot pinned by
the `DESEQ2_CRAN_SNAPSHOT` env var and the local Posit snapshot mirror
configured for the deployment.
