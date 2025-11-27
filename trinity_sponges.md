# Pipeline Trinity para metatranscriptomas de esponjas

## 1. Contexto general

- FASTQ (ya filtrados de rRNA):  
  `/media/mmorenos/disk3/000_rnaMeta_raw/03_sortMe/00_sponge_reads`

- Ensamblajes de Trinity + temporales pesados:  
  `/media/mmorenos/clau1/trinity_assemblies`      # outputs finales  
  `/media/mmorenos/clau1/trinity_tmp`             # temporales (sort, bowtie2, etc.)

Especies y réplicas:

- **Dendrilla** → E42, E54, E64  
- **Mycale sp.** → E44, E68, E95  
- **Myxilla sp.** → E22, E39, E53  

Ejemplos de archivos FASTQ:

```
1_E42_R1_non-rRNA.fastq   1_E42_R2_non-rRNA.fastq
2_E54_R1_non-rRNA.fastq   2_E54_R2_non-rRNA.fastq
3_E64_R1_non-rRNA.fastq   3_E64_R2_non-rRNA.fastq
4_E44_R1_non-rRNA.fastq   4_E44_R2_non-rRNA.fastq
5_E68_R1_non-rRNA.fastq   5_E68_R2_non-rRNA.fastq
6_E95_R1_non-rRNA.fastq   6_E95_R2_non-rRNA.fastq
7_E39_R1_non-rRNA.fastq   7_E39_R2_non-rRNA.fastq
8_E53_R1_non-rRNA.fastq   8_E53_R2_non-rRNA.fastq
11_E22_R1_non-rRNA.fastq  11_E22_R2_non-rRNA.fastq
```

---

## 2. Entorno conda de Trinity

```
conda create -n trinity_env -c conda-forge -c bioconda trinity bowtie2 samtools coreutils
conda activate trinity_env
```

---

## 3. `samples_file` para Trinity

### Estructura:

```
condición   réplica   left.fq   right.fq
```

### `Dendrilla.samples.txt`

```
Dendrilla	E42	1_E42_R1_non-rRNA.fastq	1_E42_R2_non-rRNA.fastq
Dendrilla	E54	2_E54_R1_non-rRNA.fastq	2_E54_R2_non-rRNA.fastq
Dendrilla	E64	3_E64_R1_non-rRNA.fastq	3_E64_R2_non-rRNA.fastq
```

### `Mycale.samples.txt`

```
Mycale	E44	4_E44_R1_non-rRNA.fastq	4_E44_R2_non-rRNA.fastq
Mycale	E68	5_E68_R1_non-rRNA.fastq	5_E68_R2_non-rRNA.fastq
Mycale	E95	6_E95_R1_non-rRNA.fastq	6_E95_R2_non-rRNA.fastq
```

### `Myxilla.samples.txt`

```
Myxilla	E22	11_E22_R1_non-rRNA.fastq	11_E22_R2_non-rRNA.fastq
Myxilla	E39	7_E39_R1_non-rRNA.fastq	7_E39_R2_non-rRNA.fastq
Myxilla	E53	8_E53_R1_non-rRNA.fastq	8_E53_R2_non-rRNA.fastq
```

---

## 4. Script principal Trinity usando disco externo

Guardar como: `run_trinity_by_species_clau1.sh`

```
#!/usr/bin/env bash
set -euo pipefail

THREADS=20
MEM="60G"

RAW_BASE="/media/mmorenos/disk3/000_rnaMeta_raw/03_sortMe/00_sponge_reads"
OUT_BASE="/media/mmorenos/clau1/trinity_assemblies"
TMP_BASE="/media/mmorenos/clau1/trinity_tmp"
TRINITY_BIN="Trinity"

mkdir -p "${OUT_BASE}"
mkdir -p "${TMP_BASE}"

export TMPDIR="${TMP_BASE}"
export TRINITY_MAX_SORT_BUFFER="2G"

cd "${RAW_BASE}"
SAMPLES_DIR="${RAW_BASE}/samples_files"

${TRINITY_BIN} --seqType fq --max_memory "${MEM}" --CPU "${THREADS}"   --samples_file "${SAMPLES_DIR}/Dendrilla.samples.txt"   --SS_lib_type RF --output "${OUT_BASE}/Trinity_Dendrilla_out"

${TRINITY_BIN} --seqType fq --max_memory "${MEM}" --CPU "${THREADS}"   --samples_file "${SAMPLES_DIR}/Mycale.samples.txt"   --SS_lib_type RF --output "${OUT_BASE}/Trinity_Mycale_out"

${TRINITY_BIN} --seqType fq --max_memory "${MEM}" --CPU "${THREADS}"   --samples_file "${SAMPLES_DIR}/Myxilla.samples.txt"   --SS_lib_type RF --output "${OUT_BASE}/Trinity_Myxilla_out"
```

---

## 5. Errores comunes y soluciones

### ❌ Error: `cannot open file: samples_files/Dendrilla.samples.txt`

**Causa:** ruta relativa incorrecta.

**Soluciones:**

Usar ruta absoluta:

```
--samples_file /media/mmorenos/disk3/000_rnaMeta_raw/03_sortMe/00_sponge_reads/samples_files/Dendrilla.samples.txt
```

O moverse al directorio donde están los FASTQ antes de correr:

```
cd /media/mmorenos/disk3/000_rnaMeta_raw/03_sortMe/00_sponge_reads
```

---

### ❌ Error: `No space left on device`

Ocurre en pasos como:

```
sort ... readsToComponents.out > readsToComponents.out.sort
```

**Soluciones:**

- Enviar temporales al disco externo:

```
export TMPDIR=/media/mmorenos/clau1/trinity_tmp
```

- Reducir buffer del sort:

```
export TRINITY_MAX_SORT_BUFFER=2G
```

---

### ❌ Error en `ReadsToTranscripts` (Chrysalis)

Trinity muestra:

```
Trinity run failed... ReadsToTranscripts...
```

**Solución aplicada: ejecutar ReadsToTranscripts manualmente**

```
ReadsToTranscripts -i both.fa  -f bundled_iworm_contigs.fasta  -o readsToComponents.out  -t 10 -max_mem_reads 20000000 -strand -p 5
```

Luego marcar como completado:

```
touch readsToComponents.out.ok
```

Y si se necesita:

```
touch readsToComponents.out.sort.ok
```

Relanzar Trinity:

```
Trinity ... --output Trinity_Dendrilla_out
```

---

## 6. Archivos finales esperados

```
/media/mmorenos/clau1/trinity_assemblies/Trinity_Dendrilla_out/Trinity.fasta
/media/mmorenos/clau1/trinity_assemblies/Trinity_Mycale_out/Trinity.fasta
/media/mmorenos/clau1/trinity_assemblies/Trinity_Myxilla_out/Trinity.fasta
```

---

## 7. Próximos pasos

- TransDecoder para predecir ORFs  
- Cuantificación con Salmon  
- Repartición de TPM por ORF  
- Agregado por KO  

