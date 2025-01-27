#DNA_ambiental
##Juan Alfredo Hernández-García


#1.- Bucle para Quality Control con **FASTQC**
    
    #corremos el análisis de cada una de las muestras:

    for file in raw/*.fastq.gz; do
    fastqc -t 16 "$file" -o 01_qc/
    done
    #También podemos correr el análisis concatenando todas las secuencias
    cat raw/*.fastq.gz > all.fq.gz
    fastqc --nogroup -f fastq all.fq.gz

#2.- Realizamos el trimeado de las secuencias obtenidas

    #Ejecutar Trimmomatic en todas las muestras
    #Crear un directorio para almacenar las lecturas filtradas
    mkdir -p ../02_trim


    # Ejecutar Trimmomatic en todas las muestras
    # Crear un directorio para almacenar las lecturas filtradas
    mkdir -p ../02_trim

#se renombraron los archivos para poder correr el bucle de trimmomatic
rename 's/_Custom_1/_R1/; s/_Custom_2/_R2/' *.fastq.gz

# Bucle para procesar cada par de lecturas (R1 y R2)

## Crear un directorio para almacenar las lecturas filtradas
mkdir -p ../02_trim
out=../02_trim
for r1 in *_R1.fastq.gz; do
    # Obtener el nombre del archivo R2 (reemplazando _R1 por _R2)
    r2=${r1/_R1/_R2}

    # Verificar si el archivo R2 existe
    if [[ -f "$r2" ]]; then
        # Extraer el nombre base de la muestra (eliminando _R1.fastq.gz)
        base=$(basename "$r1" _R1.fastq.gz)

        # Ejecutar Trimmomatic
        echo "Procesando muestra: $base"
        trimmomatic PE -phred33 \
            "$r1" "$r2" \
            "$out/${base}_R1_paired.fastq.gz" "$out/${base}_R1_unpaired.fastq.gz" \
            "$out/${base}_R2_paired.fastq.gz" "$out/${base}_R2_unpaired.fastq.gz" \
            ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
    else
        echo "¡Advertencia! No se encontró el archivo R2 para $r1"
    fi
done
#verificamos la calidad de los archivos trimeados

    cat 02_trim/*_paired.fastq.gz > 02_trim/all.fq.gz
    fastqc --nogroup -f fastq 02_trim/all.fq.gz


# Crear un directorio para guardar los resultados ensamblados
    mkdir -p 03_assembled

# Lista de muestras (nombres base sin _R1 o _R2)
    samples=("B10" "B11" "B12" "B13" "B14" "B6" "B7" "B8" "B9")

# Bucle para procesar cada muestra
´for sample in "${samples[@]}"; do
    # Definir los archivos de entrada (R1 y R2)
    R1="02_trim/${sample}_1_R1_paired.fastq.gz"
    R2="02_trim/${sample}_1_R2_paired.fastq.gz"

    # Definir el archivo de salida (prefijo)
    output_prefix="03_assembled/${sample}_assembled"

    # Ejecutar PEAR para ensamblar las lecturas pareadas
    pear -f "$R1" -r "$R2" -o "$output_prefix"

    # Mensaje de confirmación
    echo "Procesada muestra: $sample"
    done
    echo "¡Proceso de ensamblaje completado!"
#Filtrado por longitud y calidad con vsearch
#  vsearch para filtrar secuencias
    mkdir -p 05_filtered
    for sample in B10 B11 B12 B13 B14 B6 B7 B8 B9; do
    vsearch --fastq_filter 03_assembled/${sample}_assembled.assembled.fastq \
            --fastq_maxee 1.0 \
            --fastq_minlen 200 \
            --fastq_maxlen 500 \
            --fastaout 05_filtered/${sample}_filtered.fasta
    done
## vsearch para dereplicar secuencias, para esto es necesario modifcar los headers
# Directorio de salida
    mkdir -p 06_dereplicated

# Bucle para procesar cada muestra
    for sample in B10 B11 B12 B13 B14 B6 B7 B8 B9; do
    input_file="05_filtered/${sample}_filtered.fasta"
    simplified_file="06_dereplicated/${sample}_simplified.fasta"
    dereplicated_file="06_dereplicated/${sample}_dereplicated.fasta"
    mapping_file="06_dereplicated/${sample}_mapping.txt"

    if [[ -f "$input_file" ]]; then
        echo "Procesando muestra: $sample"

        # Crear archivo de mapeo y simplificar encabezados
        awk '/^>/ {print ">seq" ++i "\t" $0; next} {print}' "$input_file" > "$mapping_file"
        awk '/^>/ {print ">seq" ++i; next} {print}' "$input_file" > "$simplified_file"

        # Dereplicar
        vsearch --derep_fulllength "$simplified_file" \
                --output "$dereplicated_file" \
                --sizeout

        echo "Muestra $sample procesada correctamente."
    else
        echo "Error: El archivo $input_file no existe. Saltando muestra $sample."
    fi
done

echo "¡Proceso de dereplicación completado!"

echo "¡Proceso de dereplicación completado!"
done
#5.- Clustering de OTUS
#Clustering: Usar vsearch o usearch para agrupar secuencias en OTUs (Unidades Taxonómicas Operativas).
# Crear el directorio de salida si no existe
mkdir -p 07_otus

# Lista de muestras
samples=("B10" "B11" "B12" "B13" "B14" "B6" "B7" "B8" "B9")

# Bucle para procesar cada muestra
for sample in "${samples[@]}"; do
    # Definir el archivo de entrada (secuencias dereplicadas)
    input_file="06_dereplicated/${sample}_dereplicated.fasta"

    # Definir los archivos de salida
    centroids_file="07_otus/${sample}_otus.fasta"
    uc_file="07_otus/${sample}_clusters.uc"

    # Verificar si el archivo de entrada existe
    if [[ -f "$input_file" ]]; then
        echo "Procesando muestra: $sample"

        # Ejecutar vsearch para clustering
        vsearch --cluster_size "$input_file" \
                --id 0.97 \
                --centroids "$centroids_file" \
                --uc "$uc_file"

        echo "Muestra $sample procesada correctamente."
    else
        echo "Error: El archivo $input_file no existe. Saltando muestra $sample."
    fi
    done

    echo "¡Proceso de clustering completado!"

#08.- Asignación Taxonómica
# Formatear la base de datos (solo necesitas hacer esto una vez)

    makeblastdb -in coi_metazoa.fasta -dbtype nucl -out coi_metazoa

# Crear el directorio de salida si no existe
    mkdir -p 08_taxonomy

# Define la lista de archivos FASTA
files=(
  "07_otus/B10_otus.fasta"
  "07_otus/B12_otus.fasta"
  "07_otus/B14_otus.fasta"
  "07_otus/B7_otus.fasta"
  "07_otus/B9_otus.fasta"
  "07_otus/B11_otus.fasta"
  "07_otus/B13_otus.fasta"
  "07_otus/B6_otus.fasta"
  "07_otus/B8_otus.fasta"
)

# Itera sobre cada archivo en la lista
    for file in "${files[@]}"; do
  # Extrae el nombre base del archivo sin la extensión
      base_name=$(basename "$file" .fasta)

  # Define el archivo de salida para los resultados de BLAST
      output_file="08_taxonomy/blast_results_${base_name}.txt"

  # Ejecuta el comando BLAST con la base de datos formateada
      blastn -query "$file" -db coi_metazoa -out "$output_file" -outfmt 6

  # Opcional: Imprime un mensaje para indicar que el archivo ha sido procesado
      echo "Procesado: $file -> $output_file"
    done
