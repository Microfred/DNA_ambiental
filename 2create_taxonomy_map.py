import os
import pandas as pd

# Directorio donde se encuentran los archivos de resultados de BLAST
results_dir = "08_taxonomy"

# Archivo de referencia que contiene la taxonomía
reference_file = "coi_metazoa.fasta"

# Diccionario para mapear subject_id a taxonomía
taxonomy_dict = {}

# Leer el archivo de referencia para extraer la taxonomía
with open(reference_file, "r") as fasta_file:
    for line in fasta_file:
        if line.startswith(">"):  # Encabezado de la secuencia
            parts = line.strip().split()  # Dividir el encabezado en partes
            subject_id = parts[0][1:]  # Extraer el ID de referencia (eliminar ">")
            taxonomy = " ".join(parts[1:])  # Unir el resto como taxonomía
            taxonomy_dict[subject_id] = taxonomy  # Mapear subject_id a taxonomía

# Lista para almacenar todos los datos
data = []

# Iterar sobre todos los archivos en el directorio de resultados
for filename in os.listdir(results_dir):
    if filename.startswith("blast_results_") and filename.endswith(".txt"):
        filepath = os.path.join(results_dir, filename)

        # Leer el archivo de resultados de BLAST
        with open(filepath, "r") as file:
            for line in file:
                # Dividir la línea en columnas (formato outfmt 6)
                columns = line.strip().split("\t")

                # Extraer la información relevante
                otu_id = columns[0]  # ID de la OTU
                subject_id = columns[1]  # ID de la secuencia de referencia
                identity = float(columns[2])  # Porcentaje de identidad
                e_value = float(columns[10])  # Valor de e-value

                # Obtener la taxonomía del diccionario
                taxonomy = taxonomy_dict.get(subject_id, "Unknown")  # Si no se encuentra, usar "Unknown"

                # Agregar los datos a la lista
                data.append([otu_id, taxonomy, identity, e_value])

# Crear un DataFrame con los datos
columns = ["OTU_ID", "Taxonomy", "Identity", "E-value"]
df = pd.DataFrame(data, columns=columns)

# Guardar la tabla de OTUs en un archivo CSV
output_file = "08_taxonomy/otu_table_full_taxonomy.csv"
df.to_csv(output_file, index=False)

print(f"Tabla de OTUs guardada en: {output_file}")
