# popQC (Population Quality Control) Algorithms Detailed Specification

본 문서는 `popQC` 프로그램에서 수행하는 각 QC 단계별 세부 알고리즘과 입력 모수의 활용 및 계산 로직을 설명합니다.

---

## 1. SNP-level Quality Control (SNP 기반 품질 관리)

### 1.1 Call Rate (SNP 결측률 계산)
*   **목적**: 대다수 개체에서 전형(Genotyping)이 되지 않은 SNP를 제거하여 데이터 신뢰도를 확보함.
*   **입력 모수 (Parameter File)**: `SNPQC % MIN_CALLRATE` (예: 0.90)
*   **알고리즘**:
    1.  전체 개체($N$)를 대상으로 특정 SNP $j$의 결측값(Missing value, 보통 5 또는 -9) 개수를 산출($M_j$).
    2.  $CallRate_j = \frac{N - M_j}{N}$
    3.  $CallRate_j < MIN\_CALLRATE$ 인 경우 해당 SNP를 분석에서 제외(`snp_qc_status = 1`).

### 1.2 MAF (Minor Allele Frequency, 희귀 대립유전자 빈도)
*   **목적**: 집단 내 출현 빈도가 너무 낮은 대립유전자를 가진 SNP는 통계적 유의성이 낮으므로 필터링함.
*   **입력 모수**: `SNPQC % MIN_MAF` (예: 0.01)
*   **알고리즘**:
    1.  대립유전자 개수 계산: SNP $j$에서 Allele 1의 개수($n_{A1}$)와 Allele 2의 개수($n_{A2}$)를 계산.
    2.  빈도 산출: $Freq(A1) = \frac{n_{A1}}{n_{A1} + n_{A2}}$, $Freq(A2) = 1 - Freq(A1)$
    3.  $MAF_j = \min(Freq(A1), Freq(A2))$
    4.  $MAF_j < MIN\_MAF$ 인 경우 필터링(`snp_qc_status = 1`).

### 1.3 HWE (Hardy-Weinberg Equilibrium, 유전적 평형 검정)
*   **목적**: 인위적 선택이 없는 집단에서 기대되는 유전자형 빈도와 실제 관측값을 비교하여 오차(오류) 검출.
*   **입력 모수**: `SNPQC % HWE_PVALUE` (예: 1e-6), `SNPQC % HWE_METHOD` (BONFERRONI 또는 FDR)
*   **알고리즘**:
    1.  **Chi-square Test**: 자유도 1의 카이제곱 검정을 수행.
        *   $\chi^2 = \sum \frac{(Observed - Expected)^2}{Expected}$
    2.  **Breed-stratified Analysis (품종별 층화)**: `animal_breeds` 정보를 활용하여 품종별로 HWE p-value를 계산 후 통합.
    3.  **Multiple Testing Correction**:
        *   **Bonferroni**: $Threshold = \frac{HWE\_PVALUE}{N_{SNPs}}$
        *   **FDR (Benjamini-Hochberg)**: p-value를 정렬하여 위양성(False Discovery) 비율을 제어.
    4.  $P_{HWE} < Threshold$ 인 경우 기술적 오류(Genotyping Error)로 간주하여 필터링.

---

## 2. Animal-level Quality Control (개체 기반 품질 관리)

### 2.1 Sex Concordance (성별 일치 여부 확인)
*   **목적**: 가계부(PED) 상의 성별과 유전체(X-염색체) 데이터로부터 추정된 성별이 일치하는지 확인.
*   **입력 데이터**: `MAP` 파일의 Chromosome 23(또는 X), `PED` 파일의 `Sex` 필드.
*   **알고리즘**:
    1.  **Heterozygosity Calculation**: X 염색체 SNP들에 대해 각 개체의 이형접합성($F_{het}$)을 계산.
    2.  **Sex Inference**:
        *   수컷(Male): X 염색체가 하나뿐이므로 $F_{het}$이 매우 낮아야 함 (0에 수렴).
        *   암컷(Female): X 염색체가 두 개이므로 상염색체와 유사한 $F_{het}$을 가짐.
    3.  **Comparison**: 추정된 성별과 PED 파일의 성별을 비교하여 불일치(Mismatch) 시 보고(`sex_qc_status = 1`).

### 2.2 Mendelian Inconsistency (멘델 유전 모순 검사)
*   **목적**: 부모-자손 간의 유전형 전이 법칙에 어긋나는 오류(Parentage Error)를 검출.
*   **입력 데이터**: PED($Sire, Dam, Animal$), GENO matrix.
*   **알고리즘**:
    1.  **Opposite Homozygotes Check**:
        *   부모가 `AA`이고 자손이 `BB`인 경우(또는 그 반대)는 절대 발생할 수 없음.
    2.  **Inconsistency Rate**: 전체 SNP 중 모순이 발생하는 SNP의 비율을 계산.
    3.  **Handling**: 모순율이 임계값을 초과하는 개체는 부모 정보가 잘못되었거나(Parentage mismatch) 샘플 오염으로 간주.

---

## 3. 알고리즘 적용 시 고려사항
*   **Monomorphic SNPs**: 전체 집단에서 모든 개체가 동일한 유전자형(`AA`만 존재 등)을 가진 경우, 통계 분석의 수치적 불안정성을 방지하기 위해 가장 먼저 필터링(`snp_monomorphic = 1`).
*   **Efficiency**: 대용량 유전체(수십만 개의 SNP) 처리를 위해 `M_PEDHashTable` 및 효율적인 메모리 관리 기법을 사용함.
