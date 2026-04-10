# GPBLUP User Guide

This guide summarizes how to run each program. For full details, use the
program-specific documentation listed below.

## End-to-End Workflow

1. ReadFR: convert Illumina FinalReport to numeric GENO with initial QC
2. popQC: population-level SNP and animal QC
3. phimpute: impute missing SNPs in the valid SNP set
4. relgeno: compute G and G inverse from SNP data
5. relped: compute A and A inverse from pedigree data

## ReadFR

Purpose: parse Illumina FinalReport files and generate GENO output.

```bash
ReadFR <parameter_file>
```

Docs: [ReadFR/README.md](ReadFR/README.md), [ReadFR/READFR_USER_MANUAL.md](ReadFR/READFR_USER_MANUAL.md)

## popQC

Purpose: SNP and animal QC pipeline.

```bash
popQC <parameter_file>
```

Docs: [popQC/README.md](popQC/README.md), [popQC/PIPELINE_GUIDE.md](popQC/PIPELINE_GUIDE.md)

## relped

Purpose: pedigree renumbering and A/Ainv matrix outputs.

```bash
relped <parameter_file>
```

Docs: [relped/README.md](relped/README.md)

## relgeno

Purpose: genomic relationship matrix processing.

```bash
relgeno <parameter_file>
```

Docs: [relgeno/RELGENO_PARAMS.md](relgeno/RELGENO_PARAMS.md)

## phimpute

Purpose: phasing and imputation.

```bash
phimpute <parameter_file>
```

Docs: [PHimpute/README.md](PHimpute/README.md)

## Additional References

- Algorithm background: [docs/ALGORITHM_BACKGROUND.md](docs/ALGORITHM_BACKGROUND.md)
- Model proof (sketches): [docs/MODEL_PROOF.md](docs/MODEL_PROOF.md)