# Importar bibliotecas necesarias
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.ordination import pcoa

# Cargar los datos desde el archivo Excel
df = pd.read_excel("abundancias.xlsx", index_col="index")

# Mostrar las primeras filas del DataFrame
print("Datos cargados:")
print(df.head())

# Paso 1: Calcular frecuencias en porcentajes
df_percent = df.div(df.sum(axis=1), axis=0) * 100
print("\nFrecuencias en porcentajes:")
print(df_percent.head())

# Paso 2: Generar gráficos de frecuencias

# 1. Heatmap de porcentajes
plt.figure(figsize=(12, 8))
sns.heatmap(df_percent, cmap="YlOrRd", annot=True, fmt=".2f", linewidths=0.5)
plt.title("Heatmap de abundancia relativa (%)")
plt.xlabel("Taxones")
plt.ylabel("Muestras")
plt.tight_layout()
plt.show()

# 2. Gráfico de barras (porcentaje por muestra)
plt.figure(figsize=(12, 6))
df_percent.sum(axis=1).plot(kind="bar", color="skyblue")
plt.title("Abundancia relativa total por muestra (%)")
plt.xlabel("Muestras")
plt.ylabel("Abundancia relativa (%)")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# 3. Gráfico de barras (porcentaje por taxón)
plt.figure(figsize=(12, 6))
df_percent.sum().plot(kind="bar", color="lightgreen")
plt.title("Abundancia relativa total por taxón (%)")
plt.xlabel("Taxones")
plt.ylabel("Abundancia relativa (%)")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()

# Paso 3: Calcular y visualizar diversidad alfa y beta

# 1. Diversidad alfa (índice de Shannon)
shannon = alpha_diversity("shannon", df.values, ids=df.index)
print("\nDiversidad alfa (Shannon):")
print(shannon)

# Gráfico de barras para la diversidad alfa
plt.figure(figsize=(8, 5))
shannon.plot(kind="bar", color="orange")
plt.title("Diversidad alfa (Índice de Shannon) por muestra")
plt.xlabel("Muestras")
plt.ylabel("Índice de Shannon")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# 2. Diversidad beta (matriz de distancia de Bray-Curtis)
bray_curtis = beta_diversity("braycurtis", df.values, ids=df.index)

# Convertir a DataFrame para visualización
bray_curtis_df = pd.DataFrame(bray_curtis.data, index=df.index, columns=df.index)

# Heatmap de la matriz de distancia
plt.figure(figsize=(8, 6))
sns.heatmap(bray_curtis_df, cmap="YlOrRd", annot=True, fmt=".2f", linewidths=0.5)
plt.title("Matriz de distancia de Bray-Curtis")
plt.tight_layout()
plt.show()

# Análisis de coordenadas principales (PCoA)
pcoa_result = pcoa(bray_curtis)
sns.scatterplot(x=pcoa_result.samples["PC1"], y=pcoa_result.samples["PC2"], hue=df.index, s=100)
plt.title("PCoA basado en Bray-Curtis")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.legend(title="Muestras")
plt.tight_layout()
plt.show()
