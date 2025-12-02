
# Sponge metagenomes → MAGs pipeline (Dendrilla, Myxilla, Mycale)

This document summarizes the **full binning pipeline** from **quality filtering** of reads to **refined MAGs with DAS Tool** for the Antarctic sponge metagenomes (Dendrilla, Myxilla, Mycale).  
All scripts are written assuming a Linux HPC environment and a project rooted under `/datos2/mmoreno/data/01_metaChile/`.

For GitHub readability, long scripts are wrapped in collapsible sections (`<details>`).

---

## 0. Directory layout (example)

```text
01_metaChile/
├── 00_raw_reads/           # raw FASTQ (illumina)
├── 01_spongebin/           # cleaned metagenome reads (post-QC, post-host)
├── 02_spongeBCs/           # extra BC reads (e.g. BC-E117, BC-E123)
├── 03_spades_coassembly/   # SPAdes co-assemblies by sponge species
├── 04_bowtie_depth/        # mapping, BAMs, depth, binning
└── ...                     # other analysis
```

Adjust paths as needed.

---

## 1. Quality filtering + host removal (per metagenome)

This step is a **template**. Replace file names and references with your real ones. It illustrates:

- adapter/quality trimming with `fastp`
- host removal with `bowtie2`-based filtering (sponge host)
- keeping clean reads for binning

<details>
<summary><code>01_qc_hostfilter.sh</code></summary>

```bash
#!/usr/bin/env bash
set -euo pipefail

# Example for one sample GM2034-1
RAW_DIR="/datos2/mmoreno/data/01_metaChile/00_raw_reads"
OUT_DIR="/datos2/mmoreno/data/01_metaChile/01_spongebin"
HOST_IDX="/datos2/mmoreno/data/01_metaChile/host_refs/sponge_host_idx"  # bowtie2 index prefix

SAMPLE="GM2034-1_JAO"

mkdir -p "${OUT_DIR}"

# 1) Quality + adapter trimming
fastp \
  -i "${RAW_DIR}/${SAMPLE}.1.fastq.gz" \
  -I "${RAW_DIR}/${SAMPLE}.2.fastq.gz" \
  -o "${OUT_DIR}/${SAMPLE}.clean.1.fastq.gz" \
  -O "${OUT_DIR}/${SAMPLE}.clean.2.fastq.gz" \
  --detect_adapter_for_pe \
  --qualified_quality_phred 15 \
  --length_required 50 \
  --thread 16 \
  --html "${OUT_DIR}/${SAMPLE}.fastp.html" \
  --json "${OUT_DIR}/${SAMPLE}.fastp.json"

# 2) Host removal (keep only unmapped reads)
bowtie2 \
  -x "${HOST_IDX}" \
  -1 "${OUT_DIR}/${SAMPLE}.clean.1.fastq.gz" \
  -2 "${OUT_DIR}/${SAMPLE}.clean.2.fastq.gz" \
  --very-sensitive \
  --threads 16 \
  --un-conc-gz "${OUT_DIR}/${SAMPLE}_noHost.fastq.gz" \
  -S /dev/null

# Resulting clean, host-filtered reads (paired):
#  ${OUT_DIR}/${SAMPLE}_noHost.1.fastq.gz
#  ${OUT_DIR}/${SAMPLE}_noHost.2.fastq.gz
```
</details>

You can generalize this with a loop over all samples (GM2034-1, GM2034-2, …) and the BC samples (BC-E117, BC-E123, etc.).

---

## 2. Co-assembly by sponge species (SPAdes)

Co-assemble **3 individuals por especie**:

- **Dendrilla**: GM2034-1, GM2034-2, GM2034-3  
- **Myxilla**:  GM2034-4, GM2034-6, BC-E123  
- **Mycale**:   GM2034-5, GM2034-15, BC-E117  

<details>
<summary><code>02_spades_coassembly.sh</code></summary>

```bash
#!/usr/bin/env bash
set -euo pipefail

SPADES="/datos2/mmoreno/software/SPAdes-3.15.3-Linux/bin/spades.py"

BASE="/datos2/mmoreno/data/01_metaChile/03_spades_coassembly"
BIN1="/datos2/mmoreno/data/01_metaChile/01_spongebin"
BC="/datos2/mmoreno/data/01_metaChile/02_spongeBCs"

THREADS=30
MEM=500  # GB, ajusta según el servidor

mkdir -p "${BASE}"

########################################
# 1) DENDRILLA: E1, E4, E19
#    GM2034-1, GM2034-2, GM2034-3
########################################
echo "============================================================"
echo "[INFO] Co-ensamble Dendrilla (E1, E4, E19)"
echo "============================================================"

"${SPADES}" -k 97,119,121 --only-assembler \
  --pe1-1 "${BIN1}/GM2034-1_JAO.1.fastq" --pe1-2 "${BIN1}/GM2034-1_JAO.2.fastq" \
  --pe2-1 "${BIN1}/GM2034-2_JAO.1.fastq" --pe2-2 "${BIN1}/GM2034-2_JAO.2.fastq" \
  --pe3-1 "${BIN1}/GM2034-3_JAO.1.fastq" --pe3-2 "${BIN1}/GM2034-3_JAO.2.fastq" \
  -t "${THREADS}" -m "${MEM}" \
  -o "${BASE}/spades_Dendrilla"

########################################
# 2) MYXILLA: E22, E25, E123
#    GM2034-4, GM2034-6, BC-E123
########################################
echo "============================================================"
echo "[INFO] Co-ensamble Myxilla (E22, E25, E123)"
echo "============================================================"

"${SPADES}" -k 97,119,121 --only-assembler \
  --pe1-1 "${BIN1}/GM2034-4_JAO.1.fastq" --pe1-2 "${BIN1}/GM2034-4_JAO.2.fastq" \
  --pe2-1 "${BIN1}/GM2034-6_JAO.1.fastq" --pe2-2 "${BIN1}/GM2034-6_JAO.2.fastq" \
  --pe3-1 "${BC}/BC-E123_noHost.1.fastq"   --pe3-2 "${BC}/BC-E123_noHost.2.fastq" \
  -t "${THREADS}" -m "${MEM}" \
  -o "${BASE}/spades_Myxilla"

########################################
# 3) MYCALE: E24, E93, E117
#    GM2034-5, GM2034-15, BC-E117
########################################
echo "============================================================"
echo "[INFO] Co-ensamble Mycale (E24, E93, E117)"
echo "============================================================"

"${SPADES}" -k 97,119,121 --only-assembler \
  --pe1-1 "${BIN1}/GM2034-5_JAO.1.fastq"  --pe1-2 "${BIN1}/GM2034-5_JAO.2.fastq" \
  --pe2-1 "${BIN1}/GM2034-15_JAO.1.fastq" --pe2-2 "${BIN1}/GM2034-15_JAO.2.fastq" \
  --pe3-1 "${BC}/BC-E117_noHost.1.fastq"  --pe3-2 "${BC}/BC-E117_noHost.2.fastq" \
  -t "${THREADS}" -m "${MEM}" \
  -o "${BASE}/spades_Mycale"

echo "[INFO] Co-ensambles SPAdes completados."
```
</details>

Copiar los contigs finales a la carpeta de binning:

```bash
cd /datos2/mmoreno/data/01_metaChile/04_bowtie_depth

cp ../03_spades_coassembly/spades_Dendrilla/contigs.fasta  contigs.dendr.fasta
cp ../03_spades_coassembly/spades_Myxilla/contigs.fasta    contigs.myxilla.fasta
cp ../03_spades_coassembly/spades_Mycale/contigs.fasta     contigs.mycale.fasta
```

---

## 3. Mapping: Bowtie2 (metagenomes → contigs)

Mapear cada metagenoma contra el contigs de la especie correspondiente.  
**Importante:** máximo **30 hilos totales** (1 bowtie2 a la vez con 30 threads).

<details>
<summary><code>03_bowtie_mapping_sponge.sh</code></summary>

```bash
#!/bin/bash
set -euo pipefail

# Ruta de trabajo
#   /datos2/mmoreno/data/01_metaChile/04_bowtie_depth
THREADS=30

echo "========================================"
echo "[INFO] Iniciando mapeos Bowtie2 (máx ${THREADS} hilos)"
echo "========================================"

########################################
# 1) DENDRILLA: GM2034-1, GM2034-2, GM2034-3
########################################
echo "============================================================"
echo "[INFO] Mapeos Dendrilla (E1, E4, E19)"
echo "============================================================"

bowtie2 --threads ${THREADS} --reorder --very-sensitive-local \
  -x dendr_idx \
  -1 ../01_spongebin/GM2034-1_JAO.1.fastq \
  -2 ../01_spongebin/GM2034-1_JAO.2.fastq \
  --no-unal -S GM2034-1_Dendrilla.sam &> GM2034-1_Dendrilla.log

bowtie2 --threads ${THREADS} --reorder --very-sensitive-local \
  -x dendr_idx \
  -1 ../01_spongebin/GM2034-2_JAO.1.fastq \
  -2 ../01_spongebin/GM2034-2_JAO.2.fastq \
  --no-unal -S GM2034-2_Dendrilla.sam &> GM2034-2_Dendrilla.log

bowtie2 --threads ${THREADS} --reorder --very-sensitive-local \
  -x dendr_idx \
  -1 ../01_spongebin/GM2034-3_JAO.1.fastq \
  -2 ../01_spongebin/GM2034-3_JAO.2.fastq \
  --no-unal -S GM2034-3_Dendrilla.sam &> GM2034-3_Dendrilla.log


########################################
# 2) MYXILLA: GM2034-4, GM2034-6, BC-E123
########################################
echo "============================================================"
echo "[INFO] Mapeos Myxilla (E22, E25, E123)"
echo "============================================================"

bowtie2 --threads ${THREADS} --reorder --very-sensitive-local \
  -x myxilla_idx \
  -1 ../01_spongebin/GM2034-4_JAO.1.fastq \
  -2 ../01_spongebin/GM2034-4_JAO.2.fastq \
  --no-unal -S GM2034-4_Myxilla.sam &> GM2034-4_Myxilla.log

bowtie2 --threads ${THREADS} --reorder --very-sensitive-local \
  -x myxilla_idx \
  -1 ../01_spongebin/GM2034-6_JAO.1.fastq \
  -2 ../01_spongebin/GM2034-6_JAO.2.fastq \
  --no-unal -S GM2034-6_Myxilla.sam &> GM2034-6_Myxilla.log

bowtie2 --threads ${THREADS} --reorder --very-sensitive-local \
  -x myxilla_idx \
  -1 ../02_spongeBCs/BC-E123_noHost.1.fastq \
  -2 ../02_spongeBCs/BC-E123_noHost.2.fastq \
  --no-unal -S BC-E123_Myxilla.sam &> BC-E123_Myxilla.log


########################################
# 3) MYCALE: GM2034-5, GM2034-15, BC-E117
########################################
echo "============================================================"
echo "[INFO] Mapeos Mycale (E24, E93, E117)"
echo "============================================================"

bowtie2 --threads ${THREADS} --reorder --very-sensitive-local \
  -x mycale_idx \
  -1 ../01_spongebin/GM2034-5_JAO.1.fastq \
  -2 ../01_spongebin/GM2034-5_JAO.2.fastq \
  --no-unal -S GM2034-5_Mycale.sam &> GM2034-5_Mycale.log

bowtie2 --threads ${THREADS} --reorder --very-sensitive-local \
  -x mycale_idx \
  -1 ../01_spongebin/GM2034-15_JAO.1.fastq \
  -2 ../01_spongebin/GM2034-15_JAO.2.fastq \
  --no-unal -S GM2034-15_Mycale.sam &> GM2034-15_Mycale.log

bowtie2 --threads ${THREADS} --reorder --very-sensitive-local \
  -x mycale_idx \
  -1 ../02_spongeBCs/BC-E117_noHost.1.fastq \
  -2 ../02_spongeBCs/BC-E117_noHost.2.fastq \
  --no-unal -S BC-E117_Mycale.sam &> BC-E117_Mycale.log

echo "============================================================"
echo "[INFO] Todos los mapeos finalizaron."
echo "============================================================"
```
</details>

> Nota: los índices `dendr_idx`, `myxilla_idx`, `mycale_idx` ya deben existir (construidos previamente con `bowtie2-build contigs.*.fasta idxname`).

---

## 4. SAM → sorted BAM (+ index)

Convertir todos los `.sam` a `.sorted.bam` e indexar, usando `samtools`.

<details>
<summary><code>04_sam2bam_sorted.sh</code></summary>

```bash
#!/bin/bash
set -euo pipefail

# Número de hilos para samtools (≤ 30)
THREADS=10

echo "======================================="
echo "[INFO] Convirtiendo SAM → BAM sorted"
echo "======================================="

for SAM in *.sam; do
    BASE="${SAM%.sam}"
    RAW_BAM="${BASE}.bam"
    SORTED_BAM="${BASE}.sorted.bam"

    echo "[INFO] Procesando ${SAM} ..."

    # SAM → BAM
    samtools view -@ ${THREADS} -b "${SAM}" -o "${RAW_BAM}"

    # Ordenar BAM
    samtools sort -@ ${THREADS} -o "${SORTED_BAM}" "${RAW_BAM}"

    # Indexar BAM ordenado
    samtools index -@ ${THREADS} "${SORTED_BAM}"

    # Eliminar BAM intermedio
    rm "${RAW_BAM}"

    echo "[OK] ${SORTED_BAM} + ${SORTED_BAM}.bai listo."
    echo ""
done

echo "======================================="
echo "[INFO] Conversión completada."
echo "======================================="
```
</details>

---

## 5. MetaBAT2: depth + binning (multi -m, por especie)

Generar tabla de coberturas con `jgi_summarize_bam_contig_depths` y correr MetaBAT2 con varios `-m` por especie.

<details>
<summary><code>05_metabat_multi_m.sh</code></summary>

```bash
#!/usr/bin/env bash
set -euo pipefail

# ============================================
# MetaBAT2 en co-ensambles de esponjas
# ============================================

THREADS=30
MINLEN_LIST=("1500" "2000" "2500")

run_metabat_species () {
    local SPECIES="$1"        # Dendrilla / Myxilla / Mycale
    local CONTIGS="$2"        # contigs.<especie>.fasta
    local BAM_PATTERN="$3"    # *Dendrilla.sorted.bam, etc.

    echo "--------------------------------------------------------"
    echo "[INFO] Procesando especie: ${SPECIES}"
    echo "  Contigs : ${CONTIGS}"
    echo "  BAM patt: ${BAM_PATTERN}"
    echo "--------------------------------------------------------"

    local OUT_ROOT="metabat2_${SPECIES}"
    mkdir -p "${OUT_ROOT}/logs" "${OUT_ROOT}/bins"

    mapfile -t BAMS < <(ls -1 ${BAM_PATTERN} 2>/dev/null | sort || true)

    if (( ${#BAMS[@]} == 0 )); then
        echo "[WARN] No se encontraron BAMs para ${SPECIES} con patrón ${BAM_PATTERN}. Saltando."
        return 0
    fi

    echo "[INFO] BAMs usados para ${SPECIES}:"
    printf '  %s\n' "${BAMS[@]}"
    echo ""

    local JGI_DEPTH="${OUT_ROOT}/${SPECIES}_coverage_jgi.tsv"

    echo "[INFO] Generando tabla de cobertura con jgi_summarize_bam_contig_depths..."
    jgi_summarize_bam_contig_depths \
        --outputDepth "${JGI_DEPTH}" \
        --minContigLength 1000 \
        "${BAMS[@]}"

    local DEPTH_FOR_METABAT="${JGI_DEPTH}"

    for M in "${MINLEN_LIST[@]}"; do
        local OUT_PREFIX="${OUT_ROOT}/bins/metabat2_${SPECIES}_m${M}"
        echo "[INFO] Ejecutando MetaBAT2 para ${SPECIES} con -m ${M}..."
        metabat2 \
            -i "${CONTIGS}" \
            -a "${DEPTH_FOR_METABAT}" \
            -o "${OUT_PREFIX}" \
            -m "${M}" \
            -t "${THREADS}" \
            --seed 1 \
            &> "${OUT_ROOT}/logs/metabat2_${SPECIES}_m${M}.log"
    done

    echo "[OK] MetaBAT2 completado para ${SPECIES}."
    echo ""
}

run_metabat_species "Dendrilla" "contigs.dendr.fasta"   "*Dendrilla.sorted.bam"
run_metabat_species "Myxilla"   "contigs.myxilla.fasta" "*Myxilla.sorted.bam"
run_metabat_species "Mycale"    "contigs.mycale.fasta"  "*Mycale.sorted.bam"

echo "========================================================"
echo "[FIN] MetaBAT2 completado para Dendrilla, Myxilla, Mycale."
echo "========================================================"
```
</details>

---

## 6. MaxBin2 (por especie)

Usando la versión que ya adaptamos, basándonos en `run_MaxBin.pl` y listas de BAMs.

<details>
<summary><code>06_maxbin_sponge.sh</code></summary>

```bash
#!/usr/bin/env bash
set -euo pipefail

MAXBIN="/media/disk01/software/MaxBin-2.2.7/run_MaxBin.pl"

THREADS=30
MIN_CONTIG=1500
PROB=0.9

echo "============================================================"
echo "[INFO] Ejecutando MaxBin2 para Dendrilla, Myxilla y Mycale"
echo "============================================================"

# 1) DENDRILLA
echo "------------------------------------------------------------"
echo "[INFO] Procesando DENDRILLA"
echo "------------------------------------------------------------"

DENDR_BAMS=$(ls *Dendrilla.sorted.bam)

echo "${DENDR_BAMS}" | tr ' ' '\n' > list_depth_dendr.txt
mkdir -p maxbin_dendr

${MAXBIN} \
   -thread ${THREADS} \
   -min_contig_length ${MIN_CONTIG} \
   -prob_threshold ${PROB} \
   -verbose \
   -contig contigs.dendr.fasta \
   -abund_list list_depth_dendr.txt \
   -out maxbin_dendr/bin

# 2) MYXILLA
echo "------------------------------------------------------------"
echo "[INFO] Procesando MYXILLA"
echo "------------------------------------------------------------"

MYX_BAMS=$(ls *Myxilla.sorted.bam)
echo "${MYX_BAMS}" | tr ' ' '\n' > list_depth_myxilla.txt
mkdir -p maxbin_myxilla

${MAXBIN} \
   -thread ${THREADS} \
   -min_contig_length ${MIN_CONTIG} \
   -prob_threshold ${PROB} \
   -verbose \
   -contig contigs.myxilla.fasta \
   -abund_list list_depth_myxilla.txt \
   -out maxbin_myxilla/bin

# 3) MYCALE
echo "------------------------------------------------------------"
echo "[INFO] Procesando MYCALE"
echo "------------------------------------------------------------"

MYC_BAMS=$(ls *Mycale.sorted.bam)
echo "${MYC_BAMS}" | tr ' ' '\n' > list_depth_mycale.txt
mkdir -p maxbin_mycale

${MAXBIN} \
   -thread ${THREADS} \
   -min_contig_length ${MIN_CONTIG} \
   -prob_threshold ${PROB} \
   -verbose \
   -contig contigs.mycale.fasta \
   -abund_list list_depth_mycale.txt \
   -out maxbin_mycale/bin

echo "============================================================"
echo "[FIN] MaxBin2 completado para las 3 especies."
echo "============================================================"
```
</details>

---

## 7. CONCOCT (por especie)

CONCOCT requiere:

1. Cortar contigs en fragmentos (`cut_up_fasta.py`)  
2. Generar tabla de coberturas (`concoct_coverage_table.py`)  
3. Correr `concoct`  
4. Hacer `merge_cutup_clustering.py`  
5. Extract bins finales (`extract_fasta_bins.py`)

<details>
<summary><code>07_concoct_sponge.sh</code></summary>

```bash
#!/bin/bash
set -euo pipefail

# Ejecutar en: /datos2/mmoreno/data/01_metaChile/04_bowtie_depth
# Activar entorno:
#   mamba activate concoct-1p1

THREADS=30

echo "============================================================"
echo "[INFO] Corriendo CONCOCT para Dendrilla, Myxilla y Mycale"
echo "============================================================"

#########################
# 1) DENDRILLA
#########################
CONTIGS_DENDR="contigs.dendr.fasta"
BED_DENDR="dendr_contigs_10K.bed"
CUTUP_DENDR="dendr_contigs_10K.fa"
COV_DENDR="dendr_coverage_table.tsv"
OUT_DENDR="concoct_output_dendr"
MERGED_DENDR="${OUT_DENDR}/clustering_merged_dendr.csv"
BINS_DENDR="bins_dendr_concoct"

DENDR_BAMS=$(ls *Dendrilla.sorted.bam)

cut_up_fasta.py "${CONTIGS_DENDR}" \
  -c 10000 -o 0 --merge_last \
  -b "${BED_DENDR}" > "${CUTUP_DENDR}"

concoct_coverage_table.py "${BED_DENDR}" ${DENDR_BAMS} > "${COV_DENDR}"

concoct \
  --threads ${THREADS} \
  --composition_file "${CUTUP_DENDR}" \
  --coverage_file "${COV_DENDR}" \
  -b "${OUT_DENDR}/"

merge_cutup_clustering.py \
  "${OUT_DENDR}/clustering_gt1000.csv" > "${MERGED_DENDR}"

mkdir -p "${BINS_DENDR}"

extract_fasta_bins.py "${CONTIGS_DENDR}" \
  "${MERGED_DENDR}" \
  --output_path "${BINS_DENDR}/"

#########################
# 2) MYXILLA
#########################
CONTIGS_MYX="contigs.myxilla.fasta"
BED_MYX="myxilla_contigs_10K.bed"
CUTUP_MYX="myxilla_contigs_10K.fa"
COV_MYX="myxilla_coverage_table.tsv"
OUT_MYX="concoct_output_myxilla"
MERGED_MYX="${OUT_MYX}/clustering_merged_myxilla.csv"
BINS_MYX="bins_myxilla_concoct"

MYX_BAMS=$(ls *Myxilla.sorted.bam)

cut_up_fasta.py "${CONTIGS_MYX}" \
  -c 10000 -o 0 --merge_last \
  -b "${BED_MYX}" > "${CUTUP_MYX}"

concoct_coverage_table.py "${BED_MYX}" ${MYX_BAMS} > "${COV_MYX}"

concoct \
  --threads ${THREADS} \
  --composition_file "${CUTUP_MYX}" \
  --coverage_file "${COV_MYX}" \
  -b "${OUT_MYX}/"

merge_cutup_clustering.py \
  "${OUT_MYX}/clustering_gt1000.csv" > "${MERGED_MYX}"

mkdir -p "${BINS_MYX}"

extract_fasta_bins.py "${CONTIGS_MYX}" \
  "${MERGED_MYX}" \
  --output_path "${BINS_MYX}/"

#########################
# 3) MYCALE
#########################
CONTIGS_MYC="contigs.mycale.fasta"
BED_MYC="mycale_contigs_10K.bed"
CUTUP_MYC="mycale_contigs_10K.fa"
COV_MYC="mycale_coverage_table.tsv"
OUT_MYC="concoct_output_mycale"
MERGED_MYC="${OUT_MYC}/clustering_merged_mycale.csv"
BINS_MYC="bins_mycale_concoct"

MYC_BAMS=$(ls *Mycale.sorted.bam)

cut_up_fasta.py "${CONTIGS_MYC}" \
  -c 10000 -o 0 --merge_last \
  -b "${BED_MYC}" > "${CUTUP_MYC}"

concoct_coverage_table.py "${BED_MYC}" ${MYC_BAMS} > "${COV_MYC}"

concoct \
  --threads ${THREADS} \
  --composition_file "${CUTUP_MYC}" \
  --coverage_file "${COV_MYC}" \
  -b "${OUT_MYC}/"

merge_cutup_clustering.py \
  "${OUT_MYC}/clustering_gt1000.csv" > "${MERGED_MYC}"

mkdir -p "${BINS_MYC}"

extract_fasta_bins.py "${CONTIGS_MYC}" \
  "${MERGED_MYC}" \
  --output_path "${BINS_MYC}/"

echo "============================================================"
echo "[FIN] CONCOCT completado para Dendrilla, Myxilla y Mycale."
echo "============================================================"
```
</details>

---

## 8. DAS Tool (integrando MetaBAT2 + MaxBin2 + CONCOCT)

DAS Tool combina las predicciones de varios binners para obtener un set final de MAGs por especie.

<details>
<summary><code>08_dastool_sponge.sh</code></summary>

```bash
#!/usr/bin/env bash
set -euo pipefail

# Ejecutar desde: /datos2/mmoreno/data/01_metaChile/04_bowtie_depth

THREADS=30

# Asegúrate de que estos estén en el PATH o usa rutas absolutas
DASTOOL="DAS_Tool"
FASTA2BIN="Fasta_to_Scaffolds2Bin.sh"

echo "============================================================"
echo "[INFO] Corriendo DAS Tool para Dendrilla, Myxilla y Mycale"
echo "============================================================"

run_dastool_species () {
    local SPECIES="$1"      # Dendrilla / Myxilla / Mycale
    local CONTIGS="$2"      # contigs.<sp>.fasta
    local META_DIR="$3"     # metabat2_<sp>/bins
    local MAXB_DIR="$4"     # bins_<sp>
    local CONC_DIR="$5"     # bins_<sp>_concoct

    echo "------------------------------------------------------------"
    echo "[INFO] Especie: ${SPECIES}"
    echo "  Contigs      : ${CONTIGS}"
    echo "  MetaBAT2 dir : ${META_DIR}"
    echo "  MaxBin2 dir  : ${MAXB_DIR}"
    echo "  CONCOCT dir  : ${CONC_DIR}"
    echo "------------------------------------------------------------"

    local META_TSV="${SPECIES}_metabat2_scaffolds2bin.tsv"
    local MAXB_TSV="${SPECIES}_maxbin2_scaffolds2bin.tsv"
    local CONC_TSV="${SPECIES}_concoct_scaffolds2bin.tsv"

    # 1) MetaBAT2 (.fa)
    ${FASTA2BIN} -e fa -i "${META_DIR}" > "${META_TSV}"

    # 2) MaxBin2 (.fa)
    ${FASTA2BIN} -e fa -i "${MAXB_DIR}" > "${MAXB_TSV}"

    # 3) CONCOCT (.fa)
    ${FASTA2BIN} -e fa -i "${CONC_DIR}" > "${CONC_TSV}"

    local OUT_PREFIX="DASTool_${SPECIES}"

    ${DASTOOL} \
      -i "${META_TSV},${MAXB_TSV},${CONC_TSV}" \
      -l metabat2,maxbin2,concoct \
      -c "${CONTIGS}" \
      -o "${OUT_PREFIX}" \
      --threads ${THREADS} \
      --search_engine diamond \
      --write_bins 1

    echo "[OK] DAS Tool completado para ${SPECIES}."
    echo ""
}

run_dastool_species \
  "Dendrilla" \
  "contigs.dendr.fasta" \
  "metabat2_Dendrilla/bins" \
  "bins_dendr" \
  "bins_dendr_concoct"

run_dastool_species \
  "Myxilla" \
  "contigs.myxilla.fasta" \
  "metabat2_Myxilla/bins" \
  "bins_myxilla" \
  "bins_myxilla_concoct"

run_dastool_species \
  "Mycale" \
  "contigs.mycale.fasta" \
  "metabat2_Mycale/bins" \
  "bins_mycale" \
  "bins_mycale_concoct"

echo "============================================================"
echo "[FIN] DAS Tool completado para las 3 esponjas."
echo "============================================================"
```
</details>

---

## 9. Notas finales

- Después de DAS Tool, los bins refinados quedan, por ejemplo, en:
  - `DASTool_Dendrilla_DASTool_bins/`
  - `DASTool_Myxilla_DASTool_bins/`
  - `DASTool_Mycale_DASTool_bins/`

- Siguiente paso natural: **CheckM2** para evaluar completitud/contaminación de MAGs y filtrar HQ/MQ para anotación funcional y mapeo de metatranscriptomas.
