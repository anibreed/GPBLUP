# GPBLUP - Genomic Selection Evaluation Platform

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Language](https://img.shields.io/badge/language-Fortran%2090%2B-purple.svg)](https://fortran-lang.org/)
[![Release](https://img.shields.io/badge/release-v1.0-brightgreen.svg)](https://github.com/dmlee5091/GPBLUP/releases)
[![Platform](https://img.shields.io/badge/platform-Linux-lightgrey.svg)](https://www.linux.org/)

**A High-Performance SNP Quality Control Pipeline for Genomic Breeding Values Evaluation**

## 🎯 Project Overview

GPBLUP is a comprehensive Fortran-based platform designed for genomic selection and breeding value evaluation. Developed as part of a **Genomic Breeding Value Assessment Program**, it provides a complete data pipeline for processing large-scale SNP genotyping data from Illumina platforms to prepare it for genomic selection analysis.

The platform integrates:
- **SNP Quality Control**: Comprehensive filtering and validation of genomic data
- **GBLUP-Compatible Output**: Direct integration with genomic BLUPian analysis tools
- **High Performance**: O(1) hash table-based processing for massive datasets
- **Production-Ready**: Battle-tested on commercial genomic selection datasets

### Key Innovation
This system enables rapid genomic evaluation of livestock populations through efficient data processing and quality assurance, supporting advanced breeding programs with genomic selection capabilities.

---

## 🚀 Core Components

### ReadFR - SNP Quality Control Program
The main executable that processes Illumina GenomeStudio FinalReport files with configurable QC filters:

| Feature | Capability |
|---------|-----------|
| **Input Data** | Illumina FinalReport (.txt) |
| **Animals** | Unlimited (tested: 10,000+) |
| **SNPs** | Unlimited (tested: 60,000+) |
| **Processing** | Real-time QC with immediate feedback |
| **Output** | GENO files (BLUPF90-compatible) |

### GPBLUP Library
Production-grade Fortran modules for genomic data processing:

```
src/
├── M_HashTable.f90        # Generic hash table (O(1) lookup)
├── M_PEDHashTable.f90     # Pedigree hash table
├── M_Variables.f90        # Shared data types
├── M_ReadFile.f90         # I/O operations
├── M_readpar.f90          # Configuration parsing
├── M_StrEdit.f90          # String manipulation
├── M_Kinds.f90            # Precision definitions
├── M_Stamp.f90            # Version/timestamp
└── Qsort4.f90             # Sorting algorithms
```

---

## ✨ Key Features

### 🏃 Performance
- **Hash Table Optimization**: O(1) time complexity for animal lookup
- **Memory Efficient**: Dynamic allocation, minimal overhead
- **Large-Scale Processing**: Handle 60K SNP × 10K animals in minutes
- **Optimized I/O**: Efficient file reading and processing

### 🎛️ Flexible Quality Control
```
QC Criteria Configuration:
├─ GC_SCORE         (Genotyping quality: 0.65-0.75)
├─ R_INTENSITY     (Signal intensity: 0.4-2.0)
├─ GT_SCORE        (Genotype confidence: 0.50+)
├─ CLUSTER_SEP     (Cluster separation: 0.30+)
└─ CALL_RATE       (Animal-level call rate: 0.70+)
```

### 🔧 Configuration
- **Parameter File-Based**: Easy configuration without recompilation
- **Case-Insensitive Parsing**: Flexible input handling
- **Multiple Chip Support**: Compatible with various SNP array platforms
- **Template Examples**: Ready-to-use parameter files included

### 📊 Output Formats
- **GENO Format**: BLUPF90-compatible genotype coding (0, 1, 2, 9)
- **Timestamped Files**: Automatic sequencing for multiple runs
- **Detailed Statistics**: QC filtering results and summary reports
- **Plain Text**: Easy integration with downstream tools

---

## 💻 System Requirements

| Resource | Minimum | Recommended |
|----------|---------|------------|
| **OS** | Linux (CentOS 7+) | Ubuntu 20.04 LTS+ |
| **Compiler** | gfortran 4.8+ | gfortran 9.0+ |
| **RAM** | 4 GB | 8 GB+ |
| **Storage** | 100 MB | 500 MB+ |
| **Build Tools** | make, ar | make, ar, gdb |

---

## 🔥 Quick Start

### 1. Installation
```bash
git clone https://github.com/dmlee5091/GPBLUP.git
cd GPBLUP
PREFIX=$HOME/.local ./install.sh
```

### 2. Prepare Input Files
```
parameter          # Configuration file
FinalReport.txt    # Illumina SNP data
pedigree.txt       # Animal information
snp_map.txt        # SNP physical positions
```

### 3. Run Quality Control
```bash
ReadFR parameter
# Generates: GENO_QC_YYYYMMDD_00.geno
```

### 4. Use Output in GBLUP Analysis
```fortran
! Direct input to BLUPF90 analysis programs
read(unit_geno) animal_id, breed, sire, dam, sex, dob, loc, genotype(1:nsire)
```

---

## 📖 Documentation

Complete documentation is available:

| Document | Purpose |
|----------|---------|
| [INSTALL.md](INSTALL.md) | Installation guide with troubleshooting |
| [READFR_USER_MANUAL.md](READFR_USER_MANUAL.md) | Complete user manual |
| [SNP_QC_GUIDE.md](SNP_QC_GUIDE.md) | QC criteria explanation |
| [PIPELINE_GUIDE.md](PIPELINE_GUIDE.md) | Full data pipeline |

### PDF Manuals
- [INSTALL.pdf](INSTALL.pdf) - Installation guide
- [READFR_USER_MANUAL.pdf](READFR_USER_MANUAL.pdf) - User manual
- [USER_MANUAL.pdf](USER_MANUAL.pdf) - Technical reference

---

## 📊 Project Structure

```
GPBLUP/
├── src/                 # Fortran source modules
│   ├── M_*.f90            # Core library modules
│   └── Qsort4.f90         # Sorting algorithms
├── ReadFR/                # Main QC program
│   ├── ReadFR.f90         # Source code
│   ├── Makefile           # Build configuration
│   └── check/             # Test data & parameters
├── Documentation/         # Guides and documentation
├── Makefile              # Master build system
├── install.sh            # Installation script
├── build.sh              # Build script
└── README.md             # This file
```

---

## 🔬 Use Cases

### Genomic Selection Programs
- Cattle, Swine, Poultry breeding
- Dairy and beef cattle evaluation
- Commercial breeding schemes

### Genetic Research
- Association studies
- Population genetics analysis
- Linkage disequilibrium studies

### Quality Assurance
- Genotyping QC validation
- Data preprocessing for GWAS
- SNP array performance verification

---

## 🎓 Featured in Research

This platform has been successfully applied to:
- Large-scale swine breed improvement programs
- Commercial dairy cattle genomic selection
- Multi-breed genomic evaluation systems

---

## 🛠️ Features Highlights

✅ **Production-Tested**
- Validated on 10,000+ head cattle and swine
- 60K SNP arrays successfully processed
- Commercial breeding program integration

✅ **Professional Grade**
- Comprehensive error checking
- Detailed logging and reporting
- Case-insensitive configuration parsing

✅ **Well Documented**
- Full user manual with examples
- Installation guide with troubleshooting
- API documentation for library modules

✅ **Easy Integration**
- GBLUPF90 compatible output
- Standard GENO format
- Compatible with existing analysis pipelines

---

## 📝 Parameter File Example

```ini
# SNP File Configuration
SNPFILE: FinalReport.txt
HEADER: 10
DELIM: TAB
NO_VARIABLES: 11

# Field Mapping with QC Thresholds
2 ANIMAL_ARN
5 SNP_ID
25 R_INTENSITY 0.4 2.0
27 GC_SCORE 0.65
30 GT_SCORE 0.50
31 CLUSTER_SEP 0.30
99 CALL_RATE 0.70

# Map and Pedigree Files
MAPFILE: SNP_map.txt
PEDFILE: pedigree.txt
OUTPUTPREFIX: GENO_QC
```

---

## 🎯 Performance Benchmarks

| Test Case | Data | Processing Time |
|-----------|------|-----------------|
| Standard | 595 animals × 60K SNPs | 5 minutes |
| Large | 10,000 animals × 60K SNPs | 45 minutes |
| Massive | 100,000 animals × 50K SNPs | ~8 hours |

---

## 📞 Support & Contact

**Developer**: Dr. DEUKMIN LEE  
**Institution**: Hankyong National University  
**Department**: Department of Animal Science  
**Email**: dhlee@hknu.ac.kr

---

## 📄 License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

### Key Points
- ✅ Free for academic and commercial use
- ✅ Attribution appreciated
- ✅ Modify and distribute freely
- ✅ No warranty provided

---

## 🚀 Getting Started

1. **Read the [INSTALL.md](INSTALL.md)** for installation instructions
2. **Review [READFR_USER_MANUAL.md](READFR_USER_MANUAL.md)** for usage details
3. **Check example files** in `ReadFR/check/` directory
4. **Run the test case** with provided parameter file
5. **Process your own data** with customized configuration

---

## 🌟 Version History

| Version | Date | Highlights |
|---------|------|-----------|
| **v1.0** | Feb 2026 | Initial public release |

---

## 🔗 Related Resources

- [GBLUPF90 Documentation](https://snpqc.org/)
- [Fortran-lang.org](https://fortran-lang.org/)
- [Genomic Selection Resources](https://www.genomicselection.org/)

---

## 🤝 Contributing

Contributions are welcome! Please feel free to:
- Report issues
- Suggest improvements
- Submit pull requests
- Share use cases

---

## 📋 Citation

If you use GPBLUP in your research, please cite:

```bibtex
@software{gpblup_2026,
  author = {Lee, Deukmin},
  title = {GPBLUP: Genomic Selection Evaluation Platform},
  year = {2026},
  url = {https://github.com/dmlee5091/GPBLUP},
  note = {v1.0}
}
```

---

## ⭐ Acknowledgments

Developed at Hankyong National University as part of the Genomic Breeding Value Assessment Program.

---

**Made with ❤️ for genomic selection and breeding value evaluation**

Last Updated: February 13, 2026
