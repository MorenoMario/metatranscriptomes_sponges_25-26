# Read recruitment of sponge metagenomes to pooled MAGs (Linn)

This document describes the pipeline used to:

1. Renombrar MAGs por especie de esponja (Dendrilla / Myxilla / Mycale).
2. Construir un **pool de MAGs** (`MAGs_all`).
3. Mapear los 9 metagenomas (3 × 3 esponjas) contra ese pool.
4. Calcular **coverage y breadth por MAG × muestra** con `coverm genome`.
5. Analizar presencia por host (`host_specific`, `multi_host`, `shifted_host`) en R.

> 🔧 **Servidor:** Linn  
> 🧬 **Contexto:** metagenomas de esponjas antárticas Dendrilla, Myxilla y Mycale.  
> 📁 **Base:** `/datos2/mmoreno/data/01_metaChile/`

---

## 0. Requerimientos

- `bowtie2`
- `samtools`
- `coverm >= 0.7.0` (en `coverm_env` u otro env conda)
- `Fasta_to_Contig2Bin.sh` (de DAS Tool / binning helper scripts)
- `R` con `tidyverse` (instalado vía conda: `r-base`, `r-tidyverse`)

---

## 1. Renombrar MAGs con tags por especie

**Objetivo:**  
Agregar tags al nombre de cada bin para saber de qué esponja proviene:

- `dend_` → MAGs de **Dendrilla**
- `myx_`  → MAGs de **Myxilla**
- `myc_`  → MAGs de **Mycale**

Los bins de DAS Tool provienen de:

- `DASTool_Dendrilla_thr0p1_DASTool_bins/`
- `DASTool_Myxilla_thr0p1_DASTool_bins/`
- `DASTool_Mycale_thr0p1_DASTool_bins/`

y se copian a un directorio común `MAGs_all/` con los nuevos nombres.

> 🔁 Este paso se ejecutó originalmente sobre discos de trabajo intermedios
> (`/media/mmorenos/disk3/005_MAGs_3Sponges/thr0p1`) y luego se copió
> `MAGs_all` a Linn en:  
> `/datos2/mmoreno/data/01_metaChile/MAGs_all`.

<details>
<summary><code>01_rename_bins_sponges.sh</code></summary>

```bash
#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Script: 01_rename_bins_sponges.sh
# Descripción:
#   - Renombra bins de DAS Tool añadiendo prefijos por esponja:
#       Dendrilla -> dend_
#       Myxilla   -> myx_
#       Mycale    -> myc_
#   - Copia los bins renombrados a MAGs_all/
# ============================================================

# Directorio base (ajústalo si cambias de carpeta)
BASE_DIR="/media/mmorenos/disk3/005_MAGs_3Sponges/thr0p1"
cd "$BASE_DIR"

OUT_DIR="MAGs_all"
mkdir -p "$OUT_DIR"

rename_bins () {
    local INDIR="$1"
    local TAG="$2"

    echo "============================================================"
    echo "[INFO] Procesando bins en: $INDIR  -> tag: $TAG"
    echo "============================================================"

    for f in "$INDIR"/*.fa; do
        [ -e "$f" ] || continue  # por si no hay .fa

        base=$(basename "$f")
        new="${TAG}_${base}"

        # Copiar (si prefieres mover, cambia cp por mv)
        cp "$f" "$OUT_DIR/$new"

        echo "  $f  -->  $OUT_DIR/$new"
    done
}

# Dendrilla -> dend
rename_bins "DASTool_Dendrilla_thr0p1_DASTool_bins" "dend"

# Myxilla -> myx
rename_bins "DASTool_Myxilla_thr0p1_DASTool_bins" "myx"

# Mycale -> myc
rename_bins "DASTool_Mycale_thr0p1_DASTool_bins" "myc"

echo "============================================================"
echo "[OK] Todos los bins renombrados y copiados a $OUT_DIR"
echo "============================================================"
```
</details>

---

## 2. Tabla contig → bin (AllMAGs_contig2bin.tsv)

**Objetivo:**  
Crear una tabla que indique a qué MAG pertenece cada contig del pool.

- Input: `MAGs_all/*.fa` (renombrados con `dend_`, `myx_`, `myc_`)
- Output: `MAGs_all/AllMAGs_contig2bin.tsv` (2 columnas: `contig_id`, `bin_id`)

> 🖥️ Ejecutado en **Linn**.

<details>
<summary><code>02_run_Fasta_to_Contig2Bin.sh</code></summary>

```bash
#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Script: 02_run_Fasta_to_Contig2Bin.sh
# Descripción:
#   - Genera tabla contig -> MAG (contig2bin) para el pool de MAGs_all
# ============================================================

MAG_DIR="/datos2/mmoreno/data/01_metaChile/MAGs_all"
OUT_TSV="${MAG_DIR}/AllMAGs_contig2bin.tsv"

Fasta_to_Contig2Bin.sh \
  -e fa \
  -i "${MAG_DIR}" \
  > "${OUT_TSV}"

echo "[OK] Generado contig2bin:"
echo "     ${OUT_TSV}"
```
</details>

---

## 3. Mapping de metagenomas al pool de MAGs (Bowtie2)

**Objetivo:**  
Mapear los 9 metagenomas (3 por esponja) contra el **pool completo de MAGs** (`MAGs_all`) y generar BAM ordenados e indexados.

- MAGs pool: `/datos2/mmoreno/data/01_metaChile/MAGs_all`
- Índice Bowtie2: `/datos2/mmoreno/data/01_metaChile/ref_allMAGs/AllMAGs.*`
- BAMs: `/datos2/mmoreno/data/01_metaChile/map_allMAGs/*_allMAGs.sorted.bam`

Metagenomas usados (coinciden con los co-ensambles por especie):

- **Dendrilla**: GM2034-1, GM2034-2, GM2034-3
- **Myxilla**: GM2034-4, GM2034-6, BC-E123
- **Mycale**:  GM2034-5, GM2034-15, BC-E117

> 🖥️ Ejecutado en **Linn** (máx 32 hilos).

<details>
<summary><code>03_map_all_metagenomes_to_allMAGs.sh</code></summary>

```bash
#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Script: 03_map_all_metagenomes_to_allMAGs.sh
# Descripción:
#   - Construye índice Bowtie2 del pool de MAGs (MAGs_all)
#   - Mapea los 9 metagenomas (3x Dendrilla, 3x Myxilla, 3x Mycale)
#     contra ese pool
#   - Genera BAMs ordenados + indexados por muestra
#   - Ejecutar en servidor Linn.
# ============================================================

THREADS=32

# ----- MAGs (pool) -----
BASE_MAG="/datos2/mmoreno/data/01_metaChile"
MAG_DIR="${BASE_MAG}/MAGs_all"
REF_DIR="${BASE_MAG}/ref_allMAGs"
BAM_DIR="${BASE_MAG}/map_allMAGs"

# ----- Reads usados en los co-ensambles -----
BIN1="${BASE_MAG}/01_spongebin"
BC="${BASE_MAG}/02_spongeBCs"

mkdir -p "${REF_DIR}" "${BAM_DIR}"

echo "============================================================"
echo "[INFO] Construyendo índice de Bowtie2 del pool de MAGs"
echo "       MAG_DIR : ${MAG_DIR}"
echo "       REF_DIR : ${REF_DIR}"
echo "============================================================"

# 1) Construir índice de Bowtie2 del pool de MAGs
cat "${MAG_DIR}"/*.fa > "${REF_DIR}/AllMAGs.fna"
bowtie2-build "${REF_DIR}/AllMAGs.fna" "${REF_DIR}/AllMAGs"

echo "============================================================"
echo "[INFO] Iniciando mapeos de metagenomas contra AllMAGs"
echo "       BAM_DIR : ${BAM_DIR}"
echo "============================================================"

map_sample () {
    local SAMPLE_ID="$1"
    local R1="$2"
    local R2="$3"

    echo "------------------------------------------------------------"
    echo "[INFO] Mapeando muestra: ${SAMPLE_ID}"
    echo "       R1: ${R1}"
    echo "       R2: ${R2}"
    echo "------------------------------------------------------------"

    bowtie2 -x "${REF_DIR}/AllMAGs" \
        -1 "${R1}" -2 "${R2}" \
        -p "${THREADS}" \
        --very-sensitive \
        | samtools view -bS - \
        | samtools sort -@ "${THREADS}" -o "${BAM_DIR}/${SAMPLE_ID}_allMAGs.sorted.bam"

    samtools index "${BAM_DIR}/${SAMPLE_ID}_allMAGs.sorted.bam"

    echo "[OK] ${SAMPLE_ID} -> ${BAM_DIR}/${SAMPLE_ID}_allMAGs.sorted.bam"
}

############################
# DENDRILLA: GM2034-1/2/3
############################
map_sample "GM2034-1_Dendrilla" \
    "${BIN1}/GM2034-1_JAO.1.fastq" \
    "${BIN1}/GM2034-1_JAO.2.fastq"

map_sample "GM2034-2_Dendrilla" \
    "${BIN1}/GM2034-2_JAO.1.fastq" \
    "${BIN1}/GM2034-2_JAO.2.fastq"

map_sample "GM2034-3_Dendrilla" \
    "${BIN1}/GM2034-3_JAO.1.fastq" \
    "${BIN1}/GM2034-3_JAO.2.fastq"

############################
# MYXILLA: GM2034-4, GM2034-6, BC-E123
############################
map_sample "GM2034-4_Myxilla" \
    "${BIN1}/GM2034-4_JAO.1.fastq" \
    "${BIN1}/GM2034-4_JAO.2.fastq"

map_sample "GM2034-6_Myxilla" \
    "${BIN1}/GM2034-6_JAO.1.fastq" \
    "${BIN1}/GM2034-6_JAO.2.fastq"

map_sample "BC-E123_Myxilla" \
    "${BC}/BC-E123_noHost.1.fastq" \
    "${BC}/BC-E123_noHost.2.fastq"

############################
# MYCALE: GM2034-5, GM2034-15, BC-E117
############################
map_sample "GM2034-5_Mycale" \
    "${BIN1}/GM2034-5_JAO.1.fastq" \
    "${BIN1}/GM2034-5_JAO.2.fastq"

map_sample "GM2034-15_Mycale" \
    "${BIN1}/GM2034-15_JAO.1.fastq" \
    "${BIN1}/GM2034-15_JAO.2.fastq"

map_sample "BC-E117_Mycale" \
    "${BC}/BC-E117_noHost.1.fastq" \
    "${BC}/BC-E117_noHost.2.fastq"

echo "============================================================"
echo "[OK] Todos los mapeos contra AllMAGs finalizaron."
echo "============================================================"
```
</details>

---

## 4. Coverage y breadth por MAG × muestra (coverM)

**Objetivo:**  
Calcular:

- `Covered Fraction` (breadth, fracción del MAG cubierta)
- `Mean` (profundidad media)

para cada combinación **MAG × muestra** a partir de los BAM del paso 3.

- Input:
  - MAGs: `/datos2/mmoreno/data/01_metaChile/MAGs_all/*.fa`
  - BAMs: `/datos2/mmoreno/data/01_metaChile/map_allMAGs/*_allMAGs.sorted.bam`
- Output:
  - `/datos2/mmoreno/data/01_metaChile/AllMAGs_coverage_coverm.tsv`

> 🧪 Umbrales de alineamiento usados:
> - `--min-read-percent-identity 95`
> - `--min-read-aligned-percent 75`

<details>
<summary><code>04_coverm_allMAGs_coverage.sh</code></summary>

```bash
#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Script: 04_coverm_allMAGs_coverage.sh
# Descripción:
#   - Calcula coverage y breadth por MAG (genome) para todos
#     los metagenomas mapeados al pool AllMAGs.
#   - Requiere coverM instalado.
#   - Ejecutar en servidor Linn.
# ============================================================

THREADS=32

BASE="/datos2/mmoreno/data/01_metaChile"

MAG_DIR="${BASE}/MAGs_all"          # dend_*.fa, myx_*.fa, myc_*.fa
BAM_DIR="${BASE}/map_allMAGs"       # *_allMAGs.sorted.bam
OUT_TSV="${BASE}/AllMAGs_coverage_coverm.tsv"

echo "============================================================"
echo "[INFO] Ejecutando coverM genome"
echo "       MAG_DIR : ${MAG_DIR}"
echo "       BAM_DIR : ${BAM_DIR}"
echo "       OUT_TSV : ${OUT_TSV}"
echo "============================================================"

coverm genome \
  --bam-files "${BAM_DIR}"/*_allMAGs.sorted.bam \
  --genome-fasta-directory "${MAG_DIR}" \
  --genome-fasta-extension fa \
  --methods "covered_fraction" "mean" \
  --min-read-percent-identity 95 \
  --min-read-aligned-percent 75 \
  --threads "${THREADS}" \
  --output-file "${OUT_TSV}"

echo "============================================================"
echo "[OK] Matriz de coverage generada:"
echo "     ${OUT_TSV}"
echo "============================================================"
```
</details>

---

## 5. Clasificación de MAGs por host (R + tidyverse)

**Objetivo:**  
A partir de `AllMAGs_coverage_coverm.tsv`:

1. Crear una tabla larga **MAG × muestra** con:
   - `covered_fraction` (breadth),
   - `mean_depth`,
   - `host_sample` (Dendrilla / Myxilla / Mycale),
   - `mag_host` (host de origen del MAG, según prefijo `dend_`, `myx_`, `myc_`),
   - `present_sample` (TRUE/FALSE según umbrales).
2. Resumir presencia por **MAG × host** (Dendrilla / Myxilla / Mycale).
3. Clasificar cada MAG en:
   - `host_specific`   → sólo se detecta en su host de origen.
   - `multi_host`      → se detecta en ≥2 hosts.
   - `shifted_host`    → ensamblado en un host pero se ve principalmente en otro.
   - `no_host`         → sin reclutamiento consistente.

**Umbrales usados:**

- Presencia a nivel de **muestra**:
  - `covered_fraction ≥ 0.5` (50% del MAG cubierto)
  - `mean_depth ≥ 2×`
- Presencia a nivel de **host**:
  - `present_in_host = TRUE` si **≥ 2 muestras** de esa esponja pasan el filtro.

**Outputs:**

- `AllMAGs_coverage_long.tsv`  
- `AllMAGs_host_summary.tsv`  
- `AllMAGs_host_presence_pattern.tsv`

<details>
<summary><code>05_analyze_AllMAGs_coverage.R</code></summary>

```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# ============================================================
# Script: 05_analyze_AllMAGs_coverage.R
# Descripción:
#   - Lee AllMAGs_coverage_coverm.tsv (output de coverM genome)
#   - Pasa a formato largo: MAG x muestra
#   - Anota host de la muestra y host de origen del MAG
#   - Define presencia por MAG x host usando thresholds
#   - Resume patrones de presencia (host_specific, multi_host, etc.)
# ============================================================

BASE <- "/datos2/mmoreno/data/01_metaChile"
COV_FILE <- file.path(BASE, "AllMAGs_coverage_coverm.tsv")

OUT_LONG    <- file.path(BASE, "AllMAGs_coverage_long.tsv")
OUT_HOST    <- file.path(BASE, "AllMAGs_host_summary.tsv")
OUT_PATTERN <- file.path(BASE, "AllMAGs_host_presence_pattern.tsv")

message("============================================================")
message("[INFO] Leyendo: ", COV_FILE)
message("============================================================")

cov <- readr::read_tsv(COV_FILE, show_col_types = FALSE)

# Asegurarnos de que la primera columna se llame 'genome'
if (ncol(cov) < 2) {
  stop("La tabla de coverM parece tener menos de 2 columnas. Revisa el archivo: ", COV_FILE)
}

id_col <- colnames(cov)[1]
message("[INFO] La primera columna en el archivo es: ", id_col, " -> la renombro a 'genome'")
colnames(cov)[1] <- "genome"

# 1) Separar Covered Fraction y Mean en tablas largas
# CoverM usa nombres del tipo:
# 'GM2034-1_Dendrilla_allMAGs.sorted Mean'
# 'GM2034-1_Dendrilla_allMAGs.sorted Covered Fraction'

# Covered Fraction
cov_cf <- cov %>%
  dplyr::select(genome, dplyr::ends_with(" Covered Fraction")) %>%
  tidyr::pivot_longer(
    cols = -genome,
    names_to = "sample",
    values_to = "covered_fraction"
  ) %>%
  dplyr::mutate(
    # eliminar el sufijo " Covered Fraction"
    sample = stringr::str_remove(sample, " Covered Fraction$")
  )

# Mean depth
cov_mean <- cov %>%
  dplyr::select(genome, dplyr::ends_with(" Mean")) %>%
  tidyr::pivot_longer(
    cols = -genome,
    names_to = "sample",
    values_to = "mean_depth"
  ) %>%
  dplyr::mutate(
    # eliminar el sufijo " Mean"
    sample = stringr::str_remove(sample, " Mean$")
  )

# 2) Unir cf + mean
cov_long <- cov_cf %>%
  dplyr::left_join(cov_mean, by = c("genome", "sample"))

# 3) Anotar host de la muestra y host de origen del MAG
cov_long <- cov_long %>%
  dplyr::mutate(
    host_sample = dplyr::case_when(
      stringr::str_detect(sample, "Dendrilla") ~ "Dendrilla",
      stringr::str_detect(sample, "Myxilla")   ~ "Myxilla",
      stringr::str_detect(sample, "Mycale")    ~ "Mycale",
      TRUE ~ "unknown"
    ),
    mag_host = dplyr::case_when(
      stringr::str_starts(genome, "dend_") ~ "Dendrilla",
      stringr::str_starts(genome, "myx_")  ~ "Myxilla",
      stringr::str_starts(genome, "myc_")  ~ "Mycale",
      TRUE ~ "unknown"
    )
  )

# 4) Definir presencia por MAG x muestra
CF_THRESH    <- 0.5  # breadth mínima (50% del MAG cubierto)
DEPTH_THRESH <- 2    # cobertura media mínima (2x)

cov_long <- cov_long %>%
  dplyr::mutate(
    present_sample = dplyr::if_else(
      !is.na(covered_fraction) & !is.na(mean_depth) &
        covered_fraction >= CF_THRESH & mean_depth >= DEPTH_THRESH,
      TRUE, FALSE
    )
  )

# Guardar tabla larga (MAG x muestra)
readr::write_tsv(cov_long, OUT_LONG)
message("[OK] Tabla larga guardada en: ", OUT_LONG)

# 5) Resumen por MAG x host_sample
host_summary <- cov_long %>%
  dplyr::filter(host_sample != "unknown") %>%
  dplyr::group_by(genome, mag_host, host_sample) %>%
  dplyr::summarise(
    mean_covered_fraction = mean(covered_fraction, na.rm = TRUE),
    mean_depth            = mean(mean_depth, na.rm = TRUE),
    n_samples             = dplyr::n(),
    n_present             = sum(present_sample, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    # presente en un host si ≥2 muestras de esa esponja pasan el filtro
    present_in_host = n_present >= 2
  )

readr::write_tsv(host_summary, OUT_HOST)
message("[OK] Resumen por MAG x host guardado en: ", OUT_HOST)

# 6) Patrón de presencia por MAG (qué hosts lo tienen)
pattern <- host_summary %>%
  dplyr::group_by(genome, mag_host) %>%
  dplyr::summarise(
    hosts_present = paste(
      sort(host_sample[present_in_host]),
      collapse = ";"
    ),
    n_hosts_present = sum(present_in_host),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    host_category = dplyr::case_when(
      n_hosts_present == 0                              ~ "no_host",
      n_hosts_present == 1 & hosts_present == mag_host ~ "host_specific",
      n_hosts_present == 1 & hosts_present != mag_host ~ "shifted_host",
      n_hosts_present > 1                              ~ "multi_host",
      TRUE                                             ~ "unknown"
    )
  )

readr::write_tsv(pattern, OUT_PATTERN)
message("[OK] Patrón de presencia por MAG guardado en: ", OUT_PATTERN)

message("============================================================")
message("[FIN] Análisis de coverage AllMAGs completado.")
message("============================================================")
```
</details>

---

## 6. Resumen conceptual

Con estos pasos tienes:

- Un **pool de MAGs** etiquetados por host (`dend_`, `myx_`, `myc_`).
- Reclutamiento de lecturas de los **9 metagenomas** contra ese pool.
- Una matriz de coverage/breadth por **MAG × muestra** (`AllMAGs_coverage_coverm.tsv`).
- Tablas derivadas para interpretar **patrones de host**:
  - `AllMAGs_coverage_long.tsv` (MAG × muestra).
  - `AllMAGs_host_summary.tsv` (MAG × host con medias de breadth y depth).
  - `AllMAGs_host_presence_pattern.tsv` (categoría: `host_specific`, `multi_host`, `shifted_host`, `no_host`).

Esto te permite, por ejemplo:

- Detectar MAGs ensamblados en Myxilla pero **presentes y activos** en Dendrilla (caso *Nitrosopumilus*).
- Definir un subconjunto de MAGs **multi-host** para usarlos como referencia en el mapeo de metatranscriptomas.
- Contrastar simbiontes **especialistas de host** vs **generalistas entre esponjas**.
