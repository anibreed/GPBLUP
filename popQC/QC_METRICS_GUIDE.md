# QC Metrics Guide

This document defines the quality control metrics calculated by popQC.

## Animal-Level Metrics

### Call Rate (CR)
- **Definition**: Proportion of SNPs with valid genotype calls for an animal
- **Formula**: (Number of non-missing SNPs) / (Total number of SNPs)
- **Range**: 0.0 - 1.0
- **Threshold**: Typically >= 0.95
- **Use**: Filter animals with low genotyping coverage

### Heterozygosity
- **Definition**: Proportion of heterozygous SNPs in an animal
- **Formula**: (Number of heterozygous SNPs) / (Total SNPs with calls)
- **Expected**: Depends on population structure
- **Use**: Detect potential sample mix-ups or contamination

## SNP-Level Metrics

### Missing Rate
- **Definition**: Proportion of animals with missing genotype at a SNP
- **Formula**: (Number of animals with missing genotype) / (Total animals)
- **Range**: 0.0 - 1.0
- **Threshold**: Typically <= 0.05 (5%)
- **Use**: Filter low-quality SNPs

### Allele Frequency
- **Definition**: Frequency of each allele in the population
- **Calculation**: Count alleles across all animals
- **Minor Allele Frequency (MAF)**: Frequency of less common allele
- **Range MAF**: 0.0 - 0.5
- **Threshold**: Typically MAF >= 0.01 (1%)

### Hardy-Weinberg Equilibrium (HWE)
- **Definition**: Test deviation from expected genotype frequencies
- **Test**: Chi-square test comparing observed vs expected frequencies
- **Expected**: p-value > threshold (e.g., 1e-6)
- **Use**: Detect genotyping errors or population stratification
- **Formula**: 
  - Expected frequency: p² + 2pq + q² = 1
  - Where p, q are allele frequencies

## Population-Level Metrics

### Linkage Disequilibrium (LD)
- **Definition**: Non-random association of alleles at nearby loci
- **Measure**: r² or D' values
- **Use**: SNP selection, population structure analysis

### Population Structure
- **Definition**: Genetic differentiation between subpopulations
- **Measures**: F-statistics, Principal Component Analysis (PCA)
- **Use**: Account for population stratification in analysis

## Missing Data Patterns

### Type I: Random Missing
- Missing data randomly distributed across animals and SNPs
- Suggests technical issues (e.g., poor sample quality)

### Type II: SNP-Specific Missing
- High missing rate for specific SNPs
- Suggests probe/assay design issues

### Type III: Animal-Specific Missing
- High missing rate for specific animals
- Suggests poor DNA quality or low genotyping efficiency

## Output Report Descriptions

### Animal QC Summary
- List of all animals with calculated metrics
- Call rate, heterozygosity, status (PASS/FAIL)
- Reasons for failure if applicable

### SNP QC Summary
- List of all SNPs with calculated metrics
- Missing rate, allele frequencies, HWE test results
- Status (PASS/FAIL) for each SNP

### QC Statistics
- Overall call rate (population level)
- Distribution of call rates across animals
- Distribution of missing rates across SNPs
- Summary of filtered samples and SNPs

## References

- Hardy-Weinberg Equilibrium: https://en.wikipedia.org/wiki/Hardy%E2%80%93Weinberg_principle
- Population Genetics: Hartl & Clark "Principles of Population Genetics"
