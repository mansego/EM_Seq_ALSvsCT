import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ks_2samp

print("🔄 Cargando archivo TSV grande con Polars...")
df = pl.read_csv("fragment_analysis/stats/combined_fragment_lengths.tsv", separator="\t")
print(f"✅ Datos cargados: {df.shape[0]:,} fragmentos")

# === Estadísticas resumen ===
print("📊 Calculando estadísticas por grupo (ALS / Control)...")
summary = (
    df.groupby("group")
    .agg([
        pl.col("length").mean().alias("mean_length"),
        pl.col("length").median().alias("median_length"),
        pl.col("length").std().alias("sd_length"),
        (pl.col("length") < 150).mean().alias("prop_short"),
        (pl.col("length") > 250).mean().alias("prop_long")
    ])
)
summary.write_csv("fragment_analysis/stats/summary_stats.tsv", separator="\t")
print("✅ Estadísticas guardadas en 'summary_stats.tsv'")

# === KS test ===
print("🧪 Ejecutando test de Kolmogorov-Smirnov entre ALS y Control...")
als_lengths = df.filter(pl.col("group") == "ALS").get_column("length").to_numpy()
ctrl_lengths = df.filter(pl.col("group") == "Control").get_column("length").to_numpy()
ks_result = ks_2samp(als_lengths, ctrl_lengths)
with open("fragment_analysis/stats/ks_test.txt", "w") as f:
    f.write(str(ks_result))
print("✅ Resultado guardado en 'ks_test.txt'")
print("KS Test:", ks_result)

# === Muestreo para gráficos ===
print("📉 Realizando muestreo de 1 millón de fragmentos para visualización...")
sample = df.sample(n=1_000_000, seed=42)

# === Clasificación por categoría ===
print("📦 Clasificando fragmentos como cortos, largos o medios...")
sample = sample.with_columns([
    pl.when(pl.col("length") < 150).then("Short (<150 bp)")
      .when(pl.col("length") > 250).then("Long (>250 bp)")
      .otherwise("Medium")
      .alias("category")
])

# Convertimos a pandas para gráficos con seaborn
sample_pd = sample.to_pandas()
summary_pd = summary.to_pandas()

# === Visualización ===
sns.set(style="whitegrid")

# --- Densidad ---
print("📈 Generando gráfico de densidad de longitud de fragmentos...")
plt.figure(figsize=(8, 6))
sns.kdeplot(data=sample_pd, x="length", hue="group", fill=True, common_norm=False, alpha=0.5)
plt.axvline(166, linestyle="--", color="gray")
plt.axvline(320, linestyle="--", color="gray")
plt.title("cfDNA Fragment Length Distribution")
plt.xlabel("Fragment Length (bp)")
plt.ylabel("Density")
plt.tight_layout()
plt.savefig("fragment_analysis/plots/fragment_length_density.png")
plt.close()
print("✅ Gráfico guardado: fragment_length_density.png")

# --- Boxplot ---
print("📦 Generando boxplot para fragmentos cortos y largos...")
short_long_pd = sample_pd[sample_pd["category"] != "Medium"]
plt.figure(figsize=(6, 5))
sns.boxplot(data=short_long_pd, x="group", y="length", hue="category")
plt.title("Short vs. Long cfDNA Fragments by Group")
plt.tight_layout()
plt.savefig("fragment_analysis/plots/short_vs_long_boxplot.png")
plt.close()
print("✅ Gráfico guardado: short_vs_long_boxplot.png")

# --- Barras proporción ---
print("📊 Generando gráfico de barras con proporciones de fragmentos...")
melted = summary_pd.melt(id_vars="group", value_vars=["prop_short", "prop_long"],
                         var_name="category", value_name="proportion")
plt.figure(figsize=(6, 5))
sns.barplot(data=melted, x="group", y="proportion", hue="category")
plt.title("Proportion of Short/Long cfDNA Fragments")
plt.tight_layout()
plt.savefig("fragment_analysis/plots/short_long_ratio.png")
plt.close()
print("✅ Gráfico guardado: short_long_ratio.png")

print("🎉 Análisis completo.")
