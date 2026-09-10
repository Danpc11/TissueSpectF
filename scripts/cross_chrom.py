#!/usr/bin/env python3
"""cross_chrom.py -- ¿un mismo periodo se modula a la vez en varios cromosomas?

    python3 scripts/cross_chrom.py \
        --results-dir results_bp250 --dataset GSE162694 \
        --bin-mb 0.25 [--self-test]

LA PREGUNTA
-----------
Un espectro global por concatenación no sirve: pegar los cromosomas uno tras
otro pone el último gen de chr1 junto al primero de chr2, y esas 23
discontinuidades inyectan artefactos espectrales propios.

Lo que sí es una pregunta global: para una escala dada, ¿las muestras con
potencia alta en chr20 son las mismas que la tienen alta en chr16? Si la
potencia a 3 Mb sube y baja de forma coordinada entre cromosomas, hay un
proceso que modula esa escala en todo el genoma. Si cada cromosoma varía por su
cuenta, lo que hay es 23 fenómenos locales.

Se mide con la potencia POR MUESTRA, no con la del consenso: el consenso ya
promedió las muestras y con él no hay covarianza que medir.

LO QUE NO ES IDENTIFICABLE
--------------------------
Un factor que afecta a TODOS los cromosomas por igual es indistinguible de la
profundidad de secuenciación: en los datos son la misma estructura. Así que se
descuenta la primera componente principal --que lo absorbe, sea técnico o
biológico-- y se mide el acoplamiento RESIDUAL. La pregunta contestable es más
estrecha: ¿hay un SUBCONJUNTO de cromosomas acoplado entre sí por encima de lo
que mueve a todos?

Un proceso genómico que modulara los 23 cromosomas por igual sería invisible.
No es una limitación de esta implementación: ese caso no se puede separar de lo
técnico sin un control externo.

EL NULO
-------
El problema: la potencia espectral de una muestra depende de su profundidad de
secuenciación, de su calidad de RNA y de su lote. Eso correlaciona TODOS los
cromosomas de una muestra entre sí, sin que exista ninguna organización a la
escala en cuestión. Un nulo que permutara muestras libremente destruiría esa
covarianza técnica y declararía significativo cualquier par de cromosomas.

Por eso el estadístico se calcula sobre la potencia ESTANDARIZADA dentro de
cada (muestra, cromosoma) --lo que ya hace power_normalised-- y el nulo permuta
las muestras DENTRO de cada cromosoma por separado, conservando la
distribución marginal de cada uno y destruyendo sólo el emparejamiento entre
cromosomas. Eso aísla el acoplamiento, que es la afirmación.

`--self-test` comprueba que el nulo separa los casos antes de leer datos: un
acoplamiento plantado debe salir, y una covarianza puramente técnica no.
"""

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd


def coupling(M, n_null=2000, seed=42):
    """Acoplamiento medio entre columnas de M (muestras x cromosomas)."""
    ok = np.isfinite(M).all(axis=1)
    M = M[ok]
    if M.shape[0] < 8 or M.shape[1] < 2:
        return {"n": M.shape[0], "obs": np.nan, "null": np.nan, "p": np.nan}

    def stat(X):
        # QUÉ SE PUEDE Y QUÉ NO SE PUEDE IDENTIFICAR
        #
        # Un factor que afecta a TODOS los cromosomas por igual es
        # indistinguible de la profundidad de secuenciación: en los datos son
        # la misma estructura. Tres intentos previos fallaron por ignorarlo --
        # el caso de acoplamiento plantado y el de escala técnica salían
        # idénticos, porque lo son.
        #
        # Así que se quita la primera componente principal, que absorbe
        # cualquier efecto global sea técnico o biológico, y se mide el
        # acoplamiento RESIDUAL. Eso responde una pregunta más estrecha y
        # contestable: ¿hay un subconjunto de cromosomas acoplado entre sí más
        # de lo que el azar explica, por encima de lo que mueve a todos?
        #
        # Lo que se pierde: un proceso genuinamente genómico que modulara los
        # 23 cromosomas por igual sería invisible. No es una limitación de esta
        # implementación, es que ese caso no es identificable sin un control
        # externo de la parte técnica.
        # RANGOS por columna, no valores. Multiplicar por la profundidad de
        # secuenciación hace que una muestra tenga valores más dispersos en
        # TODOS sus cromosomas: la correlación media sigue en cero pero la cola
        # de |r| se infla, y el nulo --que permuta y destruye ese
        # emparejamiento-- no la reproduce. El rango elimina la
        # heterocedasticidad por construcción, así que esto es Spearman.
        Z = np.apply_along_axis(
            lambda v: np.argsort(np.argsort(v)).astype(float), 0, X)
        Z = Z - Z.mean(axis=0, keepdims=True)
        sd = Z.std(axis=0, keepdims=True)
        sd[sd == 0] = 1
        Z = Z / sd
        if Z.shape[1] >= 3:
            u, sv, vt = np.linalg.svd(Z, full_matrices=False)
            Z = Z - np.outer(u[:, 0] * sv[0], vt[0])
            sd = Z.std(axis=0, keepdims=True)
            sd[sd == 0] = 1
            Z = Z / sd
        C = np.corrcoef(Z, rowvar=False)
        iu = np.triu_indices(C.shape[1], 1)
        v = C[iu]
        v = v[np.isfinite(v)]
        # El máximo de los pares y no la media: el acoplamiento de un
        # subconjunto sube unos pocos pares y deja el resto en cero, así que
        # promediar lo diluye entre los 253 pares de 23 cromosomas.
        return float(np.quantile(np.abs(v), 0.95)) if len(v) else np.nan

    obs = stat(M)
    rng = np.random.default_rng(seed)
    null = []
    for _ in range(n_null):
        # Permutar las muestras DENTRO de cada cromosoma, por separado. Conserva
        # la distribución marginal de cada columna y rompe sólo el
        # emparejamiento, que es lo que se está probando.
        P = np.column_stack([rng.permutation(M[:, j]) for j in range(M.shape[1])])
        null.append(stat(P))
    null = np.array([v for v in null if np.isfinite(v)])
    if not len(null):
        return {"n": M.shape[0], "obs": obs, "null": np.nan, "p": np.nan}
    return {"n": M.shape[0], "obs": obs, "null": float(np.median(null)),
            "null_q95": float(np.quantile(null, 0.95)),
            "p": (1 + int((null >= obs).sum())) / (len(null) + 1)}


def self_test():
    """El nulo tiene que separar acoplamiento real de covarianza técnica."""
    rng = np.random.default_rng(0)
    n_s, n_c = 120, 20
    ok = True

    # CASO 1: acoplamiento en un SUBCONJUNTO. Cinco cromosomas comparten un
    # factor latente y los otros quince no. Es la afirmación identificable.
    lat = rng.normal(size=n_s)
    M = rng.normal(size=(n_s, n_c))
    M[:, :5] += 1.1 * lat[:, None]
    r = coupling(M, n_null=400, seed=1)
    print(f"  subconjunto acoplado    obs={r['obs']:.3f} nulo={r['null']:.3f} "
          f"p={r['p']:.4f}", end="")
    if r["p"] <= 0.05:
        print("   PASS")
    else:
        print("   FAIL: no detecta acoplamiento real")
        ok = False

    # CASO 2: sin acoplamiento. Columnas independientes -- no debe salir.
    M = rng.normal(size=(n_s, n_c))
    r = coupling(M, n_null=400, seed=2)
    print(f"  independiente           obs={r['obs']:.3f} nulo={r['null']:.3f} "
          f"p={r['p']:.4f}", end="")
    if r["p"] > 0.05:
        print("   PASS")
    else:
        print("   FAIL: falso positivo")
        ok = False

    # CASO 3: covarianza puramente TÉCNICA. Cada muestra tiene un factor de
    # escala global --profundidad, calidad de RNA-- que multiplica todos sus
    # cromosomas. Tras estandarizar dentro de cromosoma eso no debe producir
    # acoplamiento, y es el caso que un nulo mal hecho declararía positivo.
    depth = rng.lognormal(0, 0.6, size=n_s)
    M = rng.normal(size=(n_s, n_c)) * depth[:, None]
    r = coupling(M, n_null=400, seed=3)
    print(f"  sólo escala por muestra obs={r['obs']:.3f} nulo={r['null']:.3f} "
          f"p={r['p']:.4f}", end="")
    if r["p"] > 0.05:
        print("   PASS")
    else:
        print("   FAIL: confunde escala técnica con acoplamiento")
        ok = False

    # CASO 4: acoplamiento GLOBAL. Documentado como NO identificable: es la
    # misma estructura que el caso 3 y debe salir negativo. Si algún día sale
    # positivo, el estadístico dejó de descontar el efecto global y el caso 3
    # empezará a dar falsos positivos.
    lat = rng.normal(size=n_s)
    M = rng.normal(size=(n_s, n_c)) + 1.1 * lat[:, None]
    r = coupling(M, n_null=400, seed=4)
    print(f"  global (no identifiable) obs={r['obs']:.3f} nulo={r['null']:.3f} "
          f"p={r['p']:.4f}", end="")
    if r["p"] > 0.05:
        print("   PASS (esperado: se descuenta con la PC1)")
    else:
        print("   FAIL: reporta como acoplamiento algo indistinguible de lo técnico")
        ok = False

    print("\nself-test", "PASSED" if ok else "FAILED")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir")
    ap.add_argument("--dataset")
    ap.add_argument("--bin-mb", type=float, default=0.25,
                    help="ancho de bin del eje; 0 para el eje de gen")
    ap.add_argument("--n-null", type=int, default=2000)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--out")
    a = ap.parse_args()

    if a.self_test:
        sys.exit(0 if self_test() else 1)
    if not (a.results_dir and a.dataset):
        ap.error("--results-dir y --dataset son obligatorios sin --self-test")

    print("Validando el nulo antes de tocar los datos:")
    if not self_test():
        sys.exit("el nulo no separa los casos; nada de lo de abajo significa nada")
    print()

    pat = os.path.join(a.results_dir, a.dataset, "spectra", "spectra_samples_*.tsv")
    files = sorted(glob.glob(pat))
    if not files:
        sys.exit(f"sin espectros por muestra en {pat}")
    frames = []
    for f in files:
        cond = os.path.basename(f)[len("spectra_samples_"):-len(".tsv")]
        d = pd.read_csv(f, sep="\t", low_memory=False,
                        usecols=lambda c: c in {"chr", "k", "N", "sample",
                                                "power_normalised", "period"},
                        dtype={"chr": "string", "sample": "string"})
        d["condition"] = cond
        frames.append(d)
    sp = pd.concat(frames, ignore_index=True)
    sp["chr"] = sp["chr"].astype(str)
    # chrY fuera: su cobertura es del 4% y su potencia no es comparable
    sp = sp[sp.chr != "Y"]
    sp["period_mb"] = sp.period * a.bin_mb if a.bin_mb > 0 else np.nan
    if a.bin_mb <= 0:
        sys.exit("el eje de gen no da un período en Mb comparable entre "
                 "cromosomas; usá un árbol de bins")

    print(f"{a.dataset}: {sp['sample'].nunique()} muestras, "
          f"{sp.chr.nunique()} cromosomas, "
          f"{sp.period_mb.min():.2f}-{sp.period_mb.max():.1f} Mb")

    br = np.exp(np.linspace(np.log(1), np.log(100), 21))
    mid = np.sqrt(br[1:] * br[:-1])
    q = sp[(sp.period_mb >= br[0]) & (sp.period_mb <= br[-1])].copy()
    q["b"] = np.clip(np.digitize(q.period_mb, br) - 1, 0, len(mid) - 1)

    rows = []
    print(f"\n{'período (Mb)':>13} {'chr':>4} {'muestras':>9} "
          f"{'|r| obs':>8} {'|r| nulo':>9} {'p':>8}")
    for b, g in q.groupby("b"):
        # una potencia por (muestra, cromosoma) en esta banda
        m = (g.groupby(["sample", "chr"])["power_normalised"]
               .mean().unstack("chr"))
        # estandarizar dentro de cromosoma: la escala por muestra
        # --profundidad, calidad-- correlaciona todo sin organización alguna
        m = (m - m.mean()) / m.std()
        m = m.dropna(axis=1, how="all").dropna()
        if m.shape[1] < 3 or m.shape[0] < 8:
            continue
        r = coupling(m.to_numpy(), n_null=a.n_null, seed=42 + int(b))
        r.update(period_mb=round(float(mid[b]), 2), n_chr=m.shape[1])
        rows.append(r)
        star = " *" if (np.isfinite(r["p"]) and r["p"] <= 0.05) else ""
        print(f"{mid[b]:13.2f} {m.shape[1]:4d} {r['n']:9d} "
              f"{r['obs']:8.3f} {r['null']:9.3f} {r['p']:8.4f}{star}")

    if rows and a.out:
        pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
        print(f"\nescrito {a.out}")

    print("\nQué significa un p pequeño: las muestras se ordenan igual en varios "
          "cromosomas a esa escala. No dice que el período sea la causa, ni "
          "distingue un proceso genómico de un efecto de lote que sobreviva la "
          "estandarización por cromosoma.")


if __name__ == "__main__":
    main()
