# Importar bibliotecas necesarias
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
from mpl_toolkits.mplot3d import Axes3D  # Para gráficos 3D

# Cargar los datos desde el archivo Excel
df = pd.read_excel("abundancias.xlsx", index_col="index")

# Calcular la matriz de distancia de Bray-Curtis
bray_curtis = beta_diversity("braycurtis", df.values, ids=df.index)

# Realizar el PCoA
pcoa_result = pcoa(bray_curtis)

# Extraer las coordenadas de las dos primeras componentes principales (PC1, PC2, PC3)
co
