# 1. Prepare DATA
import pandas as pd
from pathlib import Path

def run_prepare_data(data_dir="data/crudos_por_cohorte", output_tsv="data/espectro_all_chromosomes.tsv"):
    DATA_DIR = Path(data_dir)
    columnas = ["chr", "k", "phase", "N", "power", "sample"]
    archivos = sorted(DATA_DIR.glob("espectro_*.tsv"))
    
    if not archivos:
        print(f"No se encontraron archivos en {DATA_DIR}")
        return

    # Usar dtype={'chr': str} para evitar DtypeWarning en la columna 0
    dfs = [pd.read_csv(a, sep="\t", usecols=columnas, dtype={'chr': str}) for a in archivos]
    espectro_all = pd.concat(dfs, ignore_index=True)
    
    Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
    espectro_all.to_csv(output_tsv, sep="\t", index=False)
    print(f"✅ Prepare Data completado. Archivo guardado en: {output_tsv}")

# 2. FASE 0
import polars as pl
from pathlib import Path

def write_matrix(sub: pl.DataFrame, out_path: Path):
    mat = sub.pivot(
        values="power",
        index="sample",
        on="k",
        aggregate_function="mean",
    ).sort("sample")
    mat = mat.rename({"sample": "Muestra"})
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mat.write_csv(out_path, separator="\t")
    print(f"✅ {out_path}  ({mat.height} muestras × {mat.width - 1} ks)")

def write_matrix_all(sub: pl.DataFrame, out_path: Path):
    sub = sub.with_columns(
        (pl.col("chr") + pl.lit("_") + pl.col("k").cast(pl.Utf8)).alias("k_chr")
    )
    mat = sub.pivot(
        values="power",
        index="sample",
        on="k_chr",
        aggregate_function="mean",
    ).sort("sample")
    mat = mat.rename({"sample": "Muestra"})
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mat.write_csv(out_path, separator="\t")
    print(f"✅ {out_path}  ({mat.height} muestras × {mat.width - 1} ks)")

def run_fase_0(input_tsv="data/espectro_all_chromosomes.tsv", out_dir="results/FASE0"):
    OUTDIR = Path(out_dir)
    
    if not Path(input_tsv).exists():
        print(f"Archivo {input_tsv} no encontrado.")
        return

    df = pl.read_csv(input_tsv, separator="\t", schema_overrides={"chr": pl.Utf8})
    df = df.select(["chr", "k", "power", "sample"])
    
    df = df.with_columns(
        pl.col("k").cast(pl.Int64),
        pl.col("power").cast(pl.Float64),
        pl.col("chr").str.strip_chars(),
    )
    
    write_matrix_all(df, OUTDIR / "ALL" / "data.tsv")
    
    for chrom in sorted(df["chr"].unique().to_list()):
        sub = df.filter(pl.col("chr") == chrom)
        write_matrix(sub, OUTDIR / f"CHR_{chrom}" / "data.tsv")
        
    print("\n🎉 Fase 0 completada.")

# 3. FASE 1
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

def find_elbow_point(x_vals, y_vals):
    x = np.array(x_vals, dtype=float)
    y = np.array(y_vals, dtype=float)
    
    x_norm = (x - x.min()) / (x.max() - x.min() + 1e-12)
    y_norm = (y - y.min()) / (y.max() - y.min() + 1e-12)
    
    p1 = np.array([x_norm[0], y_norm[0]])
    p2 = np.array([x_norm[-1], y_norm[-1]])
    line_vec = p2 - p1
    line_len = np.linalg.norm(line_vec)
    line_unitvec = line_vec / (line_len + 1e-12)
    
    points = np.column_stack((x_norm, y_norm))
    vec_p1_to_p = points - p1
    
    proj_len = np.dot(vec_p1_to_p, line_unitvec)
    proj_vec = np.outer(proj_len, line_unitvec)
    perp_vec = vec_p1_to_p - proj_vec
    distances = np.linalg.norm(perp_vec, axis=1)
    
    best_idx = np.argmax(distances)
    return int(x[best_idx])

class DenseAutoencoder(nn.Module):
    def __init__(self, input_dim, latent_dim):
        super(DenseAutoencoder, self).__init__()
        hidden_dim = max(latent_dim * 2, min(input_dim // 2, 256))
        
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, latent_dim)
        )
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, input_dim)
        )

    def forward(self, x):
        z = self.encoder(x)
        out = self.decoder(z)
        return out

def run_fase_1(base_in_dir="results/FASE0", base_out_dir="results/FASE1", max_latent=100, device="cpu"):
    torch.manual_seed(42)
    np.random.seed(42)
    
    chr_folders = sorted(
        glob.glob(os.path.join(base_in_dir, "CHR_*")),
        key=lambda p: (len(os.path.basename(p)), os.path.basename(p))
    )
    
    print(f"Cromosomas encontrados a procesar: {len(chr_folders)}")
    
    for chr_path in chr_folders:
        chr_name = os.path.basename(chr_path)
        file_path = os.path.join(chr_path, "data.tsv")
        
        if not os.path.exists(file_path):
            continue
            
        out_dir = os.path.join(base_out_dir, chr_name)
        os.makedirs(out_dir, exist_ok=True)
        
        print(f"\n==================== Procesando {chr_name} ====================")
        
        df = pd.read_csv(file_path, sep="\t", index_col=0)
        X_raw = df.values
        n_samples, n_features = X_raw.shape
        
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_raw)
        
        max_eval_k = min(max_latent, n_features, n_samples - 1)
        k_range = list(range(1, max_eval_k + 1))
        
        # 1. PCA
        pca_full = PCA(n_components=max_eval_k)
        pca_full.fit(X_scaled)
        pca_var_cum = np.cumsum(pca_full.explained_variance_ratio_) * 100.0
        elbow_pca = find_elbow_point(k_range, pca_var_cum)
        print(f"[{chr_name}] Codo PCA detectado en k={elbow_pca} (Varianza: {pca_var_cum[elbow_pca-1]:.2f}%)")
        
        # 2. Autoencoder
        tensor_x = torch.tensor(X_scaled, dtype=torch.float32)
        dataset = TensorDataset(tensor_x, tensor_x)
        loader = DataLoader(dataset, batch_size=64, shuffle=True)
        
        ae_mse_losses = []
        criterion = nn.MSELoss()
        
        for l_dim in k_range:
            model = DenseAutoencoder(input_dim=n_features, latent_dim=l_dim).to(device)
            optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-5)
            
            model.train()
            for epoch in range(35):
                for batch_x, _ in loader:
                    batch_x = batch_x.to(device)
                    optimizer.zero_grad()
                    preds = model(batch_x)
                    loss = criterion(preds, batch_x)
                    loss.backward()
                    optimizer.step()
            
            model.eval()
            with torch.no_grad():
                full_x = tensor_x.to(device)
                recon = model(full_x)
                final_mse = criterion(recon, full_x).item()
                ae_mse_losses.append(final_mse)
                
        elbow_ae = find_elbow_point(k_range, ae_mse_losses)
        print(f"[{chr_name}] Codo Autoencoder detectado en z={elbow_ae} (MSE: {ae_mse_losses[elbow_ae-1]:.4f})")
        
        # Guardar métricas
        metrics_df = pd.DataFrame({
            'dim': k_range,
            'pca_cum_variance_pct': pca_var_cum,
            'ae_reconstruction_mse': ae_mse_losses
        })
        metrics_df.to_csv(os.path.join(out_dir, "elbow_metrics.tsv"), sep="\t", index=False)
        
        # 3. Gráficas
        fig1, ax1 = plt.subplots(figsize=(8, 5))
        ax1.plot(k_range, pca_var_cum, color='#2B6CB0', lw=2.2, label='Varianza acumulada (%)')
        ax1.axvline(elbow_pca, color='#E53E3E', linestyle='--', lw=1.5, label=f'Codo PCA (k={elbow_pca})')
        ax1.scatter([elbow_pca], [pca_var_cum[elbow_pca-1]], color='#E53E3E', s=80, zorder=5)
        ax1.set_title(f'PCA: Varianza Explicada - {chr_name}', fontweight='bold')
        ax1.set_xlabel('Número de Componentes Principales')
        ax1.set_ylabel('Varianza Explicada Acumulada (%)')
        ax1.grid(True, linestyle=':', alpha=0.6)
        ax1.legend(loc='lower right')
        plt.tight_layout()
        fig1.savefig(os.path.join(out_dir, "grafica_pca_codo.png"), dpi=200)
        fig1.savefig(os.path.join(out_dir, "grafica_pca_codo.pdf"), dpi=200)
        plt.close(fig1)

        fig2, ax2 = plt.subplots(figsize=(8, 5))
        ax2.plot(k_range, ae_mse_losses, color='#805AD5', lw=2.2, label='Error Reconstrucción (MSE)')
        ax2.axvline(elbow_ae, color='#DD6B20', linestyle='--', lw=1.5, label=f'Codo AE (z={elbow_ae})')
        ax2.scatter([elbow_ae], [ae_mse_losses[elbow_ae-1]], color='#DD6B20', s=80, zorder=5)
        ax2.set_title(f'Autoencoder: Pérdida vs Dimensión Latente - {chr_name}', fontweight='bold')
        ax2.set_xlabel('Dimensión del Espacio Latente (z)')
        ax2.set_ylabel('MSE de Reconstrucción')
        ax2.grid(True, linestyle=':', alpha=0.6)
        ax2.legend(loc='upper right')
        plt.tight_layout()
        fig2.savefig(os.path.join(out_dir, "grafica_autoencoder_codo.png"), dpi=200)
        fig2.savefig(os.path.join(out_dir, "grafica_autoencoder_codo.pdf"), dpi=200)
        plt.close(fig2)

        fig3, ax_p = plt.subplots(figsize=(9, 5.5))
        ax_ae = ax_p.twinx()

        l1 = ax_p.plot(k_range, pca_var_cum, color='#2B6CB0', lw=2.2, label='PCA: Varianza acumulada (%)')
        pt1 = ax_p.scatter([elbow_pca], [pca_var_cum[elbow_pca-1]], color='#E53E3E', s=70, zorder=5)
        ax_p.axvline(elbow_pca, color='#E53E3E', linestyle=':', lw=1.2, alpha=0.8)

        l2 = ax_ae.plot(k_range, ae_mse_losses, color='#805AD5', lw=2.2, label='AE: Reconstrucción (MSE)')
        pt2 = ax_ae.scatter([elbow_ae], [ae_mse_losses[elbow_ae-1]], color='#DD6B20', s=70, zorder=5)
        ax_ae.axvline(elbow_ae, color='#DD6B20', linestyle=':', lw=1.2, alpha=0.8)

        ax_p.set_title(f'Comparativa y Codos: PCA vs Autoencoder - {chr_name}', fontweight='bold')
        ax_p.set_xlabel('Dimensión Reducida (k componentes / z latente)')
        ax_p.set_ylabel('PCA: Varianza Acumulada (%)', color='#2B6CB0')
        ax_ae.set_ylabel('Autoencoder: MSE Reconstrucción', color='#805AD5')

        ax_p.tick_params(axis='y', labelcolor='#2B6CB0')
        ax_ae.tick_params(axis='y', labelcolor='#805AD5')
        ax_p.grid(True, linestyle=':', alpha=0.6)

        lines = l1 + l2
        labels = [l.get_label() for l in lines] + [f'Codo PCA ({elbow_pca})', f'Codo AE ({elbow_ae})']
        ax_p.legend(lines + [pt1, pt2], labels, loc='center right')

        plt.tight_layout()
        fig3.savefig(os.path.join(out_dir, "grafica_traslapada_codos.png"), dpi=200)
        fig3.savefig(os.path.join(out_dir, "grafica_traslapada_codos.pdf"), dpi=200)
        plt.close(fig3)

    print("\nProcesamiento y generación de gráficas completada para todos los cromosomas.")

# 4. ORQUESTADOR
import torch

# Diccionario de configuración central
CONFIG = {
    "RUN_PREPARE_DATA": True,
    "RUN_FASE_0": True,
    "RUN_FASE_1": True,
    "MAX_LATENT_DIM": 100,
    "DEVICE": "cuda" if torch.cuda.is_available() else "cpu"
}

if __name__ == "__main__":
    if CONFIG["RUN_PREPARE_DATA"]:
        print("--- Ejecutando Prepare Data ---")
        run_prepare_data(
            data_dir="data/crudos_por_cohorte",
            output_tsv="data/espectro_all_chromosomes.tsv"
        )
        
    if CONFIG["RUN_FASE_0"]:
        print("\n--- Ejecutando Fase 0 ---")
        run_fase_0(
            input_tsv="data/espectro_all_chromosomes.tsv",
            out_dir="results/FASE0"
        )
        
    if CONFIG["RUN_FASE_1"]:
        print("\n--- Ejecutando Fase 1 ---")
        run_fase_1(
            base_in_dir="results/FASE0",
            base_out_dir="results/FASE1",
            max_latent=CONFIG["MAX_LATENT_DIM"],
            device=CONFIG["DEVICE"]
        )
    
    print("\nOrquestador finalizado.")
