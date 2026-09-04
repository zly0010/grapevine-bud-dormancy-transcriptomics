# GSE337039 raw FASTQ branch

Frozen tool versions: fastp 0.24.3, STAR 2.7.11b, featureCounts 2.0.8. fastp and STAR are recorded in the existing software_versions.tsv. Its featureCounts field is blank; the existing per-run forward count header explicitly records featureCounts v2.0.8. Read length is 86 bp, single-end, with 60 run accessions in metadata/validation/GSE337039_sample_manifest.tsv.

The manifest records public SRA run URLs. The actual STAR rescue script downloads SRA objects from `https://sra-pub-run-odp.s3.amazonaws.com/sra/<RUN>/<RUN>` using aria2c, validates them with vdb-validate and falls back to prefetch. fasterq-dump creates `<RUN>.fastq`; fastp creates `<RUN>.clean.fastq.gz`. No public raw reads are redistributed here.

Reference assembly: NCBI Vitis vinifera 12X, GCF_000003745.3. The script downloads `GCF_000003745.3_12X_genomic.fna.gz` and `GCF_000003745.3_12X_genomic.gtf.gz` from `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/745/GCF_000003745.3_12X/`. featureCounts uses the decompressed same-assembly GTF (`genes.gtf`), feature `exon`, attribute `gene_id`. This is distinct from the release-62 gene-to-GO mapping resource.

The user must provide storage and construct the STAR index locally; `04b_rescue_gse337039_star.sh` does this when no complete index exists. It uses sjdbOverhang 85, genomeSAindexNbases 13, genomeSAsparseD 2 and a 9,000,000,000-byte index-generation RAM setting. No index is included.

Within the chosen work directory: `reference/ncbi_12X_refseq/`, `star_sra/`, `star_fastq/`, `star_clean_fastq/`, `qc/fastp_star/`, `star_align/`, `star_counts/`, `logs_star/`, `metadata_star/`. QC includes per-run fastp JSON/HTML, STAR reports, featureCounts summaries and MultiQC. The future runner removes each successful sample's intermediate raw/alignment files by its original logic.

Both unstranded (`-s 0`) and forward (`-s 1`) count outputs are generated for assessment. Final analysis uses `<RUN>.forward.txt`. The DESeq2 consumer's `auto` branch checks strandedness before selecting; the frozen run used forward counts. The count/VST tables in the archive are the existing frozen outputs.

Sequential raw counting syntax (documentation only):

```bash
mkdir -p data/gse337039_raw_work
bash scripts/extended_pipeline/04b_rescue_gse337039_star.sh data/gse337039_raw_work metadata/validation/GSE337039_sample_manifest.tsv
```

Then the existing 05b R script takes `<work_dir> <manifest> <mapping> <results_dir> [auto|forward|unstranded]`; use `forward` to reproduce the frozen choice. The optional `04c_run_gse337039_full_parallel.sh` launches two workers and invokes 05b with `auto`. Its workers can both attempt an index build if no index exists; build the index first using the sequential branch before using this historical parallel runner. No raw command was run in the packaging task.

Other executable dependencies actually referenced: SRA Toolkit (prefetch, fasterq-dump, vdb-validate), samtools, pigz, curl, aria2c, MultiQC, GNU core utilities and R. Version gaps remain explicitly unfilled where no reliable record exists.
