# Anvi'o workflow – MAGs + Metatranscriptomas

Pipeline en pasos, listo para pegar en GitHub. Ajusta rutas, nombres de archivos y número de hilos según tu servidor.

---

## 1. Generar bases de datos de contigs para cada MAG

Asumiendo que `list_genomes` contiene una lista de archivos FASTA (uno por MAG):

```bash
for sample in `cat list_genomes`; do
  anvi-gen-contigs-database -f $sample -o ${sample}.db -T 40 -n "sponges_MAGs_metatrans_db"
done &> out_anvi_gen.txt
```

Salida: un archivo `*.db` por MAG (base de datos de contigs de anvi'o) y el log en `out_anvi_gen.txt`.

---

## 2. Generar `GENOMES.db` (genomes storage)

Asumiendo que `genome_path` es el archivo `external-genomes` con la ruta a cada `.db`:

```bash
anvi-gen-genomes-storage   --external-genomes genome_path   --gene-caller prodigal   -o sponges_MAGs_metatrans-GENOMES.db   &> out_genomes-storage.txt
```

Salida: `sponges_MAGs_metatrans-GENOMES.db` y log `out_genomes-storage.txt`.

---

## 3. De-replicar genomas con pyANI

Usando pyANI a modo de clustering / de-replicación:

```bash
anvi-dereplicate-genomes   -e genome_path   --skip-fasta-report   --program pyANI   -o pyANI_0.95   --similarity-threshold 0.85   --representative-method length   --force-overwrite   -T 10
```

Salida: directorio `pyANI_0.95` con los resultados de de-replicación y lista de genomas representativos.

---

## 4. Anotación funcional con NCBI COGs

Asumiendo que `genome_path` contiene los nombres de las bases de datos de contigs (`*.db`) a anotar:

```bash
for sample in `cat genome_path`; do
  anvi-run-ncbi-cogs -c $sample -T 10     --cog-data-dir /home/mmorenos/miniconda3/envs/anvio-8/lib/python3.10/site-packages/anvio/data/misc/COG/
done
```

Salida: tablas de COGs asociadas a cada base de datos de contigs.

---

## 5. Anotación funcional con KEGG KOfams

Asumiendo que `list_genomes` contiene los nombres de las bases de datos (`*.db`) o genomas a anotar:

```bash
for sample in `cat list_genomes`; do
  anvi-run-kegg-kofams -c $sample -T 10
done
```

Salida: anotaciones KEGG KOfams asociadas a cada MAG.

---

## 6. Matriz de funciones (COG20) a través de genomas

Una vez anotados los genomas, generar una matriz de funciones (COG20) para todos los MAGs:

```bash
anvi-script-gen-function-matrix-across-genomes   -e genome_path2   --annotation-source COG20_FUNCTION   --output-file-prefix MY-GENOMES
```

Notas:
- `genome_path2` es un archivo `external-genomes` (o similar) que apunta a los genomas/DBs que quieres incluir.
- Salida: archivos `MY-GENOMES*` con la matriz de funciones por genoma.

---

## 7. Estimación de metabolismo (KEGG módulos)

Generar una matriz de presencia/copia de módulos metabólicos:

```bash
anvi-estimate-metabolism   -e genome_path2   --matrix-format   --add-copy-number   --include-zeros   --include-metadata
```

Notas:
- Usa el mismo `genome_path2` que en el paso anterior, para mantener consistencia en los genomas incluidos.
- Salida: matrices de metabolismo (módulos KEGG) listas para análisis en R / Python.

---
