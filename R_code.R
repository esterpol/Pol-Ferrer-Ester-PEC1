# Cargamos las librerías necesarias
library(readxl, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(S4Vectors, quietly = TRUE)
library(SummarizedExperiment, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(tidyr, quietly = TRUE)

# Cargamos el dataset
dataset_path <- "TIO2+PTYR-human-MSS+MSIvsPD.XLSX" 
data <- read_excel(dataset_path)

# Examinamos el dataset
head(data)
colnames(data)

# Seleccionamos columnas relevantes y preparar el dataset
data_abundances <- data %>% select(SequenceModifications, starts_with("M"), starts_with("T"))
head(data_abundances)

# Creamos metadatos de grupo
sample_names <- c("M1_1_MSS", "M1_2_MSS", "M5_1_MSS", "M5_2_MSS", "T49_1_MSS", "T49_2_MSS", 
                  "M42_1_PD", "M42_2_PD", "M43_1_PD", "M43_2_PD", "M64_1_PD", "M64_2_PD")
groups <- c("MSS", "MSS", "MSS", "MSS", "MSS", "MSS", "PD", "PD", "PD", "PD", "PD", "PD")
metadata <- DataFrame(Sample = sample_names, Group = groups)

# Creamos el objeto SummarizedExperiment
abundances_matrix <- as.matrix(data_abundances %>% select(-SequenceModifications))
SE <- SummarizedExperiment(assays = list(counts = abundances_matrix), colData = metadata)

# Guardamos el objeto para análisis posterior
save(SE, file = "dataset_metabolomics.Rda")

# Extraemos la matriz de abundancias del objeto SummarizedExperiment
abundances_matrix_boxplot <- assay(SE, "counts")

# Normalizamos logarítmicamente para reducir el rango de valores
abundances_matrix_boxplot <- log2(abundances_matrix_boxplot + 1)

# Convertimos la matriz en un dataframe y añadir los grupos
plot_data <- as.data.frame(t(abundances_matrix_boxplot))  # Transponemos para tener una columna por fosfopéptido
plot_data$Group <- colData(SE)$Group  # Añadir los grupos como columna

# Convertimos abundances_matrix en formato largo
plot_data_long <- plot_data %>%
  pivot_longer(cols = starts_with("V"), names_to = "Fosfopeptido", values_to = "Abundancia")

# Creamos el boxplot con múltiples fosfopéptidos
ggplot(plot_data_long, aes(x = Group, y = Abundancia, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribución de Abundancias Normalizadas de Fosfopéptidos por Grupo",
       x = "Grupo",
       y = "Abundancia de Fosfopéptidos (log2)") +
  theme(legend.position = "none")

# Filtramos filas con baja variación o constantes
abundances_matrix_pca <- abundances_matrix[apply(abundances_matrix, 1, function(x) {
  sd(x) > 1e-5  # Filtra filas con desviación estándar significativa (ajusta si es necesario)
}), ]

# Normalizamos logarítmicamente para reducir el rango de valores
abundances_matrix_pca <- log2(abundances_matrix_pca + 1)

# Realizamos el PCA en la matriz filtrada
pca <- prcomp(t(abundances_matrix_pca), scale = TRUE)
pca_data <- data.frame(pca$x, Group = colData(SE)$Group)

# Graficamos el PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA de las Abundancias de Fosfopéptidos",
       x = "PC1",
       y = "PC2")

# Extraemos y normalizamos la matriz de abundancias (si no está ya normalizada)
abundances_matrix <- assay(SE, "counts")
abundances_matrix <- log2(abundances_matrix + 1)

# Creamos dataframe para realizar comparaciones
long_data <- as.data.frame(t(abundances_matrix))  # Transponer para tener una columna por fosfopéptido
long_data$Group <- colData(SE)$Group  # Añadir la información de grupo
long_data <- pivot_longer(long_data, cols = starts_with("V"), names_to = "Fosfopeptido", values_to = "Abundancia")

# Realizamos prueba t para cada fosfopéptido y obtener p-values
results <- long_data %>%
  group_by(Fosfopeptido) %>%
  summarise(p_value = t.test(Abundancia ~ Group)$p.value)

# Aplicamos corrección de p-values para múltiples pruebas (FDR)
results <- results %>%
  mutate(adj_p_value = p.adjust(p_value, method = "fdr"))

# Filtramos fosfopéptidos significativos
significant_results <- results %>% filter(adj_p_value < 0.05)

# Ordenamos los resultados por adj_p_value de menor a mayor
ordered_results <- significant_results %>% arrange(adj_p_value)

# Mostramos fosfopéptidos con diferencias significativas entre MSS y PD ordenados
print(ordered_results)

