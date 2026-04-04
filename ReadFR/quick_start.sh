#!/bin/bash
#
set -e

# Quick Start Template (placeholders only)
# Replace all paths before use.

PIPELINE_ROOT="/path/to/your/ReadFR/Pipeline"
READFR_BIN="/path/to/ReadFR"
PED_FILE="/path/to/pedigree.txt"
MAP_V1="/path/to/MAP_V1.txt"
MAP_V2="/path/to/MAP_V2.txt"

echo "========================================="
echo "  Genome QC Pipeline - Quick Start"
echo "========================================="
echo ""

echo "Step 1: Creating directory structure..."
mkdir -p ${PIPELINE_ROOT}/{Raw_Data/{ChipV1,ChipV2},QC_Results/{ChipV1,ChipV2},Merged_Data,Population_QC,Imputation,Final_Results}
echo "✓ Directory structure created"
echo ""

echo "Step 2: Creating parameter file templates..."

cat > ${PIPELINE_ROOT}/parameter_v1 << EOF
COMMENT PED file
PEDFile
${PED_FILE}

COMMENT SNP FinalReport
SNPFile
[FINALREPORT_FILE]
ANIMAL-ARN 2
SNP_Name 1
Chr 10
Position 11
Allele1-AB 13
Allele2-AB 14
GC_Score 27
R-Intensity 25
GT_Score 29
Cluster_Sep 30

COMMENT MAP V1
MAPFile
${MAP_V1}
EOF

cat > ${PIPELINE_ROOT}/parameter_v2 << EOF
COMMENT PED file
PEDFile
${PED_FILE}

COMMENT SNP FinalReport
SNPFile
[FINALREPORT_FILE]
ANIMAL-ARN 2
SNP_Name 1
Chr 10
Position 11
Allele1-AB 13
Allele2-AB 14
GC_Score 27
R-Intensity 25
GT_Score 29
Cluster_Sep 30

COMMENT MAP V2
MAPFile
${MAP_V2}
EOF

echo "✓ Parameter templates created"
echo ""

echo "Step 3: Creating batch processing scripts..."

cat > ${PIPELINE_ROOT}/batch_qc_v1.sh << 'EOFBASH'
#!/bin/bash
# Chip V1 QC Processing (template)

PARAM_TEMPLATE="../parameter_v1"
INPUT_DIR="Raw_Data/ChipV1"
OUTPUT_DIR="QC_Results/ChipV1"
LOG_DIR="QC_Results/ChipV1/logs"
READFR="/path/to/ReadFR"

mkdir -p ${OUTPUT_DIR} ${LOG_DIR}

total_files=$(ls ${INPUT_DIR}/*.txt 2>/dev/null | wc -l)
current=0

echo "Processing Chip V1 files (Total: ${total_files})"

for finalreport in ${INPUT_DIR}/*.txt; do
    current=$((current + 1))
    basename=$(basename ${finalreport} .txt)
    
    echo "[${current}/${total_files}] ${basename}"
    
    sed "s|\[FINALREPORT_FILE\]|${finalreport}|g" ${PARAM_TEMPLATE} > temp_param.txt
    
    ${READFR} temp_param.txt > ${LOG_DIR}/${basename}.log 2>&1
    
    if [ -f QC_PASSED_GENO.txt ]; then
        mv QC_PASSED_GENO.txt ${OUTPUT_DIR}/${basename}_GENO.txt
        echo "  ✓ GENO created"
    else
        echo "  ✗ No GENO output"
    fi
    
    rm -f temp_param.txt
done

echo "Chip V1 QC completed!"
EOFBASH

sed 's/V1/V2/g; s/v1/v2/g' ${PIPELINE_ROOT}/batch_qc_v1.sh > ${PIPELINE_ROOT}/batch_qc_v2.sh
chmod +x ${PIPELINE_ROOT}/batch_qc_*.sh

echo "✓ Batch processing scripts created"
echo ""

cat > ${PIPELINE_ROOT}/merge_geno_files.sh << 'EOFMERGE'
#!/bin/bash
# Merge all GENO files (template)

cd Merged_Data

OUTPUT="ALL_INDIVIDUALS_GENO.txt"

first_file=$(find ../QC_Results -name "*_GENO.txt" | head -1)
head -1 "$first_file" > ${OUTPUT}
find ../QC_Results -name "*_GENO.txt" -exec tail -n +2 {} \; >> ${OUTPUT}

total=$(tail -n +2 ${OUTPUT} | wc -l)
echo "Merged ${total} individuals into ${OUTPUT}"

awk 'NR>1 {print $1}' ${OUTPUT} | sort | uniq -d > duplicates.txt
if [ -s duplicates.txt ]; then
    echo "Warning: Found $(wc -l < duplicates.txt) duplicate IDs"
else
    echo "No duplicates found"
    rm duplicates.txt
fi
EOFMERGE

chmod +x ${PIPELINE_ROOT}/merge_geno_files.sh

echo "✓ Merge script created"
echo ""

echo "========================================="
echo "  Setup completed successfully!"
echo "========================================="
