# GPBLUP Project - Complete Documentation Index

## 📋 Quick Navigation

### For Getting Started
- **[README.md](README.md)** - Main project overview and quick start
- **[BUILD.md](BUILD.md)** - Step-by-step build instructions
- **[PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md)** - Directory layout and organization

### For Understanding the Project
- **[ORGANIZATION_SUMMARY.md](ORGANIZATION_SUMMARY.md)** - What was reorganized and why
- **[popQC/docs/POPQC_SPECIFICATION.md](popQC/docs/POPQC_SPECIFICATION.md)** - Technical specification with algorithms

### For Running popQC
- **[popQC/README.md](popQC/README.md)** - popQC-specific guide
- **[popQC/parameter_example.txt](popQC/parameter_example.txt)** - Sample configuration

---

## 📚 Document Descriptions

### README.md
**Purpose**: Project overview and entry point  
**Contains**: 
- Project overview and features
- Quick start (30 seconds)
- Installation requirements
- Basic usage examples
- Performance metrics
- Links to detailed documentation

**Audience**: Developers, users new to GPBLUP

### BUILD.md
**Purpose**: Complete build procedure and troubleshooting  
**Contains**:
- Detailed compilation steps (Step 1-3 with exact commands)
- Compiler flags explained
- Makefile and CMake instructions
- Troubleshooting 5 common issues
- Performance tuning options
- Continuous integration examples

**Audience**: Build engineers, developers

### PROJECT_STRUCTURE.md
**Purpose**: Architecture and organization guide  
**Contains**:
- Full directory hierarchy with 3-level depth
- Module dependency graph
- Program specifications (popQC, ReadFR)
- Compilation workflow diagram
- File organization rules (Do's and Don'ts)
- Maintenance checklist

**Audience**: Developers, system administrators

### ORGANIZATION_SUMMARY.md
**Purpose**: Project reorganization details  
**Contains**:
- What was restructured (source → build → bin)
- File organization reference
- Access paths for users
- Key improvements with before/after comparison
- Regular maintenance operations
- Version control recommendations

**Audience**: Project maintainers, technical leads

### POPQC_SPECIFICATION.md
**Purpose**: Technical specification of popQC algorithm  
**Contains**:
- 100+ sections covering:
  - Mathematical formulations
  - Algorithm descriptions
  - Data structures used
  - Pseudo-code
  - Performance analysis
  - Formula derivations

**Audience**: Algorithm designers, QA engineers

---

## 🎯 Choose Your Document by Role

### Role: **New User - Just Want to Run popQC**
1. Start: [README.md](README.md) - Quick Start section
2. Then: [popQC/README.md](popQC/README.md) - Usage examples
3. Reference: [popQC/parameter_example.txt](popQC/parameter_example.txt) - Config template

**Time**: 10 minutes to first run

### Role: **Build Engineer - Need to Compile**
1. Start: [BUILD.md](BUILD.md) - Prerequisites section
2. Follow: Detailed Compilation Sequence (Steps 1-3)
3. Verify: Verification section
4. Troubleshoot: If issues, check Troubleshooting section

**Time**: 20 minutes to working executables

### Role: **System Administrator - Maintain Infrastructure**
1. Start: [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md)
2. Understand: [ORGANIZATION_SUMMARY.md](ORGANIZATION_SUMMARY.md)
3. Refer: [BUILD.md](BUILD.md) - Maintenance Schedule section

**Time**: 30 minutes to understand full system

### Role: **Developer - Add New Features**
1. Start: [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) - Module Architecture
2. Understand: [POPQC_SPECIFICATION.md](popQC/docs/POPQC_SPECIFICATION.md)
3. Reference: [BUILD.md](BUILD.md) - Development section

**Time**: 1-2 hours to understand codebase

### Role: **Performance Analyst - Optimize Execution**
1. Start: [README.md](README.md) - Performance section
2. Refer: [POPQC_SPECIFICATION.md](popQC/docs/POPQC_SPECIFICATION.md) - Performance section
3. Use: [BUILD.md](BUILD.md) - Performance Tuning section

**Time**: 1 hour to implement optimizations

---

## 📁 File Map

```
/home/dhlee/GPBLUP/
├── README.md                        ← Start here (project overview)
├── BUILD.md                         ← Build instructions
├── PROJECT_STRUCTURE.md             ← Architecture diagram
├── ORGANIZATION_SUMMARY.md          ← What was reorganized
├── DOCUMENTATION_INDEX.md           ← This file
│
├── source/                          ← Shared modules (10 files)
│   └── M_*.f90, Qsort4.f90
│
├── popQC/                           ← Main program
│   ├── src/
│   │   └── popQC.f90               (1,605 lines)
│   ├── docs/
│   │   ├── POPQC_SPECIFICATION.md   ← Algorithm details
│   │   ├── POPQC_SPECIFICATION.pdf
│   │   └── 논문_popQC_설명서.pdf
│   ├── README.md                    ← popQC usage guide
│   └── parameter_example.txt        ← Config template
│
├── ReadFR/                          ← Secondary program
│   ├── src/
│   │   └── ReadFR.f90
│   └── [docs, configs]
│
├── build/                           ← Build artifacts (temporary)
│   ├── common/
│   ├── popQC/
│   └── ReadFR/
│
├── bin/                             ← Executables
│   ├── popQC (115 KB)
│   └── ReadFR (86 KB)
│
└── [CMakeLists.txt, Makefile, etc.]
```

---

## 🔄 Typical Workflows

### Build from Source

```bash
cd /home/dhlee/GPBLUP
make clean          # Clean old artifacts
make all            # Build everything
make VERBOSE=1 all  # With verbose output
ls -lh bin/         # Verify executables
```

**Reference**: [BUILD.md](BUILD.md) - "Using Makefile" section

### Run Analysis

```bash
# Prepare data files
cp /path/to/animals.ped ./
cp /path/to/markers.map ./

# Edit parameter file
cp popQC/parameter_example.txt my_params.txt
# ... modify thresholds ...

# Run analysis
./bin/popQC my_params.txt

# Check results
cat results/SNP_QC_Report.txt
```

**Reference**: [popQC/README.md](popQC/README.md) - "Usage" section

### Debug Build

```bash
make FFLAGS="-g -O0 -fbounds-check" popQC
gdb ./bin/popQC parameter_file.txt
```

**Reference**: [BUILD.md](BUILD.md) - "Compiler Flags" section

### Rebuild After Code Change

```bash
# Edit source (e.g., popQC/src/popQC.f90)
vi popQC/src/popQC.f90

# Clean and rebuild
make clean
make popQC

# Test
./bin/popQC parameter_example.txt
```

**Reference**: [BUILD.md](BUILD.md) - "Maintenance" section

---

## 📊 Project Statistics

| Metric | Value |
|--------|-------|
| **Source Code** | |
| - Total Source Lines | ~2,500 (Fortran) |
| - popQC Lines | 1,605 |
| - ReadFR Lines | ~900 |
| **Modules** | 10 shared modules |
| **Subroutines** | 15+ in popQC |
| **Documentation** | |
| - README files | 3 |
| - Markdown guides | 4 |
| - PDF manuals | 2 |
| - Build time** | 5 seconds (incremental) |
| **Executable Size** | |
| - popQC | 115 KB |
| - ReadFR | 86 KB |
| **Performance** | |
| - Time (10K × 60K) | 2-3 minutes |
| - Memory usage | ~1 MB |
| - Accuracy | 98-99% |

---

## 🚀 Getting Help

### If compilation fails:
→ Check [BUILD.md](BUILD.md) - Troubleshooting section

### If not sure how to use popQC:
→ Read [popQC/README.md](popQC/README.md)

### If confused about directory structure:
→ Review [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md)

### If need algorithm details:
→ Study [POPQC_SPECIFICATION.md](popQC/docs/POPQC_SPECIFICATION.md)

### If need to make changes:
→ Follow [BUILD.md](BUILD.md) - Development section

### If need to understand reorganization:
→ Read [ORGANIZATION_SUMMARY.md](ORGANIZATION_SUMMARY.md)

---

## 🔧 For Maintenance

### Daily Tasks
- Monitor execution: `./bin/popQC config.txt`
- Check results: `cat results/*.txt`

### Weekly Tasks
- Clean build: `make clean && make all`
- Verify executables: `ls -lh bin/`

### Monthly Tasks
- Review logs: `cat popQC/docs/`
- Update documentation as needed

### Quarterly Tasks
- Archive project: `make distclean`
- Commit to version control
- Tag release version

**Reference**: [BUILD.md](BUILD.md) section "Maintenance Schedule"

---

## 📝 Document Maintenance Notes

### This Index (DOCUMENTATION_INDEX.md)
- **Last Updated**: 2024
- **Maintained By**: GPBLUP Development Team
- **Versioning**: Matches project version
- **Update Frequency**: Per release

### When Adding New Documents
1. Add entry in "Document Descriptions"
2. Add file path in "File Map"
3. Add to relevant role section
4. Update statistics if needed

### When Updating BUILD Process
1. Revise [BUILD.md](BUILD.md)
2. Update [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) if architecture changes
3. Verify [README.md](README.md) Quick Start matches

---

## 📞 Support Summary

| Question | Answer Location |
|----------|-----------------|
| How do I build? | [BUILD.md](BUILD.md) |
| How do I run popQC? | [popQC/README.md](popQC/README.md) |
| What is the directory structure? | [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) |
| What was reorganized? | [ORGANIZATION_SUMMARY.md](ORGANIZATION_SUMMARY.md) |
| How does HWE test work? | [POPQC_SPECIFICATION.md](popQC/docs/POPQC_SPECIFICATION.md) |
| Build failed - what's wrong? | [BUILD.md](BUILD.md) - Troubleshooting |
| What are compilation flags? | [BUILD.md](BUILD.md) - Compiler Flags Explained |
| How to debug code? | [BUILD.md](BUILD.md) - Debug Build |
| Performance optimization? | [BUILD.md](BUILD.md) - Performance Tuning |

---

## ✅ Checklist for New Contributors

- [ ] Read [README.md](README.md)
- [ ] Review [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md)
- [ ] Build project: `make clean && make all`
- [ ] Understand [BUILD.md](BUILD.md) steps
- [ ] Run example: `./bin/popQC popQC/parameter_example.txt`
- [ ] Study relevant [POPQC_SPECIFICATION.md](popQC/docs/POPQC_SPECIFICATION.md) sections
- [ ] Set up development environment
- [ ] Create feature branch (git)
- [ ] Make changes and test
- [ ] Update documentation
- [ ] Submit for review

---

## 📚 Complete Reading Order

### For Complete Understanding (All Documents in Order)

1. **[README.md](README.md)** (10 min)
   - Overview and features

2. **[PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md)** (15 min)
   - Understand organization

3. **[ORGANIZATION_SUMMARY.md](ORGANIZATION_SUMMARY.md)** (10 min)
   - See what was improved

4. **[BUILD.md](BUILD.md)** (25 min)
   - Learn build process in detail

5. **[popQC/README.md](popQC/README.md)** (10 min)
   - Understand program features

6. **[POPQC_SPECIFICATION.md](popQC/docs/POPQC_SPECIFICATION.md)** (60+ min)
   - Deep dive into algorithms

**Total Time**: ~2 hours for comprehensive understanding

---

## 🎓 Learning Paths

### Path A: Developer (2 hours)
1. README.md (overview)
2. PROJECT_STRUCTURE.md (architecture)
3. BUILD.md (build system)
4. POPQC_SPECIFICATION.md (algorithms)

→ **Result**: Ready to modify code

### Path B: User (30 minutes)
1. README.md (quick start)
2. popQC/README.md (usage)
3. popQC/parameter_example.txt (config)

→ **Result**: Ready to run popQC

### Path C: DevOps (1 hour)
1. PROJECT_STRUCTURE.md (organization)
2. BUILD.md (build instructions)
3. ORGANIZATION_SUMMARY.md (maintenance)

→ **Result**: Ready to manage infrastructure

### Path D: Performance Engineer (1.5 hours)
1. README.md (metrics)
2. POPQC_SPECIFICATION.md (performance section)
3. BUILD.md (performance tuning)

→ **Result**: Ready to optimize

---

## 🔗 External Links

- **Fortran Documentation**: https://fortran-lang.org/
- **gfortran Manual**: https://gcc.gnu.org/onlinedocs/gfortran/
- **CMake Guide**: https://cmake.org/cmake/help/latest/
- **GNU Make**: https://www.gnu.org/software/make/manual/

---

**Status**: Documentation Complete  
**Version**: 1.0  
**Last Updated**: 2024  
**Maintained By**: GPBLUP Development Team

