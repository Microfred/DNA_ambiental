import os
import pandas as pd

# Directorio donde se encuentran los archivos de resultados de BLAST
results_dir = "08_taxonomy"

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

                # Aquí puedes agregar más procesamiento para obtener la taxonomía
                # Por ejemplo, si el subject_id contiene información taxonómica
                taxonomy = subject_id  # Esto es un ejemplo, ajusta según tu caso

                # Agregar los datos a la lista
                data.append([otu_id, taxonomy, identity, e_value])

# Crear un DataFrame con los datos
columns = ["OTU_ID", "Taxonomy", "Identity", "E-value"]
df = pd.DataFrame(data, columns=columns)

# Guardar la tabla de OTUs en un archivo CSV
output_file = "08_taxonomy/otu_table.csv"
df.to_csv(output_file, index=False)

print(f"Tabla de OTUs guardada en: {output_file}")
