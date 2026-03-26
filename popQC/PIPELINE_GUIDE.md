# popQC Pipeline Guide

## Overview

popQC implements a comprehensive quality control pipeline for genomic SNP data. 
This guide describes the workflow and how to use popQC in your analysis pipeline.

## Workflow

```
ReadFR Output (GENO file)
        |
        v
    popQC
    / | | \
   /  |  |  \
  /   |  |   \
Calc Call Calc Allele Calc HWE Detect Outliers
Rates Frequencies P-values
  |      |         |        |
  v      v         v        v
   Generate QC Reports
        |
        v
  Filtered Genotypes
        |
        v
Quality Control Complete
```

## Step-by-Step Usage

### 1. Generate ReadFR Output

First, run ReadFR to generate genotype data:

```bash
cd /home/dhlee/GPBLUP/ReadFR
./ReadFR parameter.txt
# Output: GENO_*.geno file
```

### 2. Prepare popQC Parameter File

Create a parameter file for popQC:

```bash
cd /home/dhlee/GPBLUP/popQC
cp test/parameter_example.txt parameter.txt
```

Edit `parameter.txt` to specify:
- GENO_FILE: Path to ReadFR output
- PED_FILE: Pedigree file
- MAP_FILE: Marker map file
- QC Thresholds: Call rates, MAF ranges, HWE p-values, etc.

### 3. Compile popQC

```bash
cd /home/dhlee/GPBLUP/popQC
make clean
make
```

### 4. Run popQC

```bash
cd /home/dhlee/GPBLUP/popQC
./popQC parameter.txt
```

### 5. Review QC Reports

popQC generates several output files:
- `popQC_Results_animal_qc.txt` - Animal-level QC summary
- `popQC_Results_snp_qc.txt` - SNP-level QC summary
- `popQC_Results_statistics.txt` - Population statistics
- `popQC_Results_filtered.geno` - Filtered genotype file

## Configuration Options

### QC Thresholds

| Parameter | Default | Recommended Range | Description |
|-----------|---------|-------------------|-------------|
| MIN_CALL_RATE_ANIMAL | 0.95 | 0.90 - 0.99 | Minimum call rate per animal |
| MIN_CALL_RATE_SNP | 0.95 | 0.90 - 0.99 | Minimum call rate per SNP |
| MIN_MAF | 0.01 | 0.01 - 0.05 | Minimum minor allele frequency |
| MAX_MAF | 0.50 | 0.45 - 0.50 | Maximum minor allele frequency |
| MIN_HWE_PVALUE | 1e-6 | 1e-4 - 1e-8 | HWE test threshold |

### Output Options

```
DETAILED_REPORTS = yes    # Generate detailed per-sample/SNP reports
SAVE_FILTERED_GENO = yes  # Save filtered genotype data
OUTPUT_PREFIX = popQC_Results  # Output file prefix
```

## Troubleshooting

### Issue: "Parameter file not found"
- **Solution**: Check file path, use absolute or relative path from popQC directory

### Issue: "Cannot open GENO file"
- **Solution**: Ensure ReadFR was executed first, check file path in parameter file

### Issue: "PED file format error"
- **Solution**: Verify PED file format matches ReadFR specification

### Issue: Low call rate for all animals
- **Solution**: Check QC threshold values, may be too stringent

## QC Best Practices

1. **Start with lenient thresholds** - Gradually tighten to understand data quality
2. **Examine outliers** - Investigate animals/SNPs failing QC
3. **Check for batch effects** - Ensure samples from different dates pass QC
4. **Validate filtered data** - Confirm expected number of animals/SNPs retained
5. **Document decisions** - Keep record of thresholds used and samples filtered

## Integration with Downstream Analysis

After QC filtering:

```
popQC Output (filtered GENO)
        |
        v
  GPBLUP Analysis
        |
        v
Genomic Predictions / Variance Components
```

## Performance Considerations

- **Large datasets**: popQC is optimized for datasets with:
  - 1,000s to 100,000s of animals
  - 10,000s to 1,000,000s of SNPs
- **Memory**: Requires sufficient RAM for genotype array
- **Runtime**: Scales with number of animals × number of SNPs

## Next Steps

See [README.md](README.md) for additional information.
For metric definitions, see [QC_METRICS_GUIDE.md](QC_METRICS_GUIDE.md).
