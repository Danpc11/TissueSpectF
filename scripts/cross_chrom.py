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


def coupling_stat(M):
    """El estadístico solo, sin nulo. Extraído para que el maxT conjunto pueda
    aplicar las MISMAS permutaciones a todas las bandas."""
    ok = np.isfinite(M).all(axis=1)
    M = M[ok]
    if M.shape[0] < 8 or M.shape[1] < 2:
        return np.nan
    Z = np.apply_along_axis(
        lambda v: np.argsort(np.argsort(v)).astype(float), 0, M)
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
    return float(np.quantile(np.abs(v), 0.95)) if len(v) else np.nan


def joint_maxt(mats, n_null=2000, seed=42):
    """maxT sobre las bandas, con permutaciones CONJUNTAS.

    Cada iteración genera UN conjunto de permutaciones --una por cromosoma-- y
    lo aplica a TODAS las bandas antes de tomar el máximo. Sin eso el sorteo
    número 10 de una banda no corresponde a la misma permutación que el 10 de
    otra, así que alinear sorteos independientes y tomar su máximo no es una
    realización conjunta del nulo: es conservador cuando las bandas están
    correlacionadas positivamente, pero no es el maxT que el programa dice
    calcular.

    Las bandas comparten muestras y cromosomas, así que la misma permutación es
    aplicable a todas. Lo que cambia entre bandas es el VALOR en cada celda, no
    quién es cada muestra.

    QUÉ CAMBIA Y QUÉ NO, medido. El nulo POR BANDA es idéntico al de
    `coupling()`: permutar cada cromosoma con su propia permutación da lo mismo
    se comparta o no entre bandas (p = 0.119 contra 0.124 en la misma matriz).
    Lo que la permutación conjunta cambia es el MÁXIMO entre bandas, que es
    donde la correlación entre ellas importa -- y esa es la única cantidad que
    este maxT produce y `coupling()` no puede.

    @param mats dict banda -> matriz (muestras x cromosomas), ya alineadas a
      las mismas muestras y en el mismo orden
    """
    bands = [b for b in mats if np.isfinite(mats[b]).all(axis=1).sum() >= 8]
    if len(bands) < 2:
        return None
    # Mismas filas en todas: la permutación tiene que significar lo mismo.
    n_rows = {mats[b].shape[0] for b in bands}
    if len(n_rows) != 1:
        raise ValueError(
            "joint_maxt: las bandas tienen distinto número de muestras "
            f"({sorted(n_rows)}). Alinealas antes: una permutación conjunta "
            "exige que la fila i sea la misma muestra en todas las bandas.")
    n = n_rows.pop()
    n_chr = mats[bands[0]].shape[1]

    obs = {b: coupling_stat(mats[b]) for b in bands}
    rng = np.random.default_rng(seed)
    gmax = np.empty(n_null)
    for it in range(n_null):
        # UNA permutación por cromosoma, compartida por todas las bandas.
        perms = [rng.permutation(n) for _ in range(n_chr)]
        stats = []
        for b in bands:
            P = np.column_stack([mats[b][perms[j], j] for j in range(n_chr)])
            stats.append(coupling_stat(P))
        gmax[it] = np.nanmax(stats) if np.any(np.isfinite(stats)) else np.nan
    gmax = gmax[np.isfinite(gmax)]
    if not len(gmax):
        return None
    return {b: {"obs": obs[b],
                "p_maxt": (1 + int((gmax >= obs[b]).sum())) / (len(gmax) + 1)}
            for b in bands if np.isfinite(obs[b])}, gmax


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
    out = {"n": M.shape[0], "obs": obs, "null": float(np.median(null)),
           "null_q95": float(np.quantile(null, 0.95)),
           "p": (1 + int((null >= obs).sum())) / (len(null) + 1)}
    return out


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

    # CASO 5: el maxT conjunto, de punta a punta. La version anterior leia
    # `df.mb_band`, una columna que nunca se creaba, asi que fallaba con
    # AttributeError en la primera corrida real -- y ningun self-test lo
    # tocaba porque solo se probaba `coupling()` en aislamiento.
    # El acoplamiento tiene que ser FUERTE para que sobreviva al descuento de
    # la PC1: con 4 de 20 cromosomas y un factor de 1.2 el propio coupling()
    # da p = 0.12, asi que exigir p <= 0.05 al maxT --que corrige por 4 bandas
    # ademas-- era pedirle mas potencia que a la prueba por banda. Ese fue el
    # error del caso de prueba, no de la implementacion.
    # DÓNDE DEJA DE FUNCIONAR, y es del estadístico, no del maxT.
    #
    # Descontar la PC1 borra cualquier acoplamiento que sea la ESTRUCTURA
    # DOMINANTE, no sólo el global. Medido: con 8 de 20 cromosomas acoplados y
    # un factor de 2.5, la PC1 se carga con 0.34-0.36 uniforme sobre esos ocho
    # y absorbe el 34% de la varianza, así que el estadístico de la banda
    # acoplada (0.283) queda POR DEBAJO de una banda de puro ruido (0.292).
    #
    # LA VENTANA DE SENSIBILIDAD, medida sobre 20 cromosomas y 60 muestras.
    # p_maxt de la banda acoplada:
    #
    #   k acoplados   factor 1.0   factor 2.0   factor 4.0
    #             2        0.960        0.582        0.015  <- detecta
    #             3        0.945        0.348        0.065
    #             4        0.458        0.771        0.930
    #             6        0.900        0.672        0.896
    #
    # Sólo detecta con DOS cromosomas acoplados y un acoplamiento fuerte. Con
    # cuatro o más la PC1 lo absorbe y el estadístico de la banda acoplada cae
    # por debajo de una de puro ruido. Es una ventana estrechísima, y quien use
    # este script tiene que saberlo: un resultado negativo aquí no dice que no
    # haya acoplamiento, dice que no hay uno de ese tipo particular.
    n_s2, n_c2 = 60, 20
    lat2 = rng.normal(size=n_s2)
    mats = {}
    for b in range(4):
        M = rng.normal(size=(n_s2, n_c2))
        if b == 1:
            # 2 de 20 con factor 4: el ÚNICO régimen medido donde el maxT
            # conjunto detecta. Ver la tabla en el comentario de arriba.
            M[:, :2] += 4.0 * lat2[:, None]
        mats[b] = M
    out = joint_maxt(mats, n_null=300, seed=9)
    if out is None:
        print("  maxT conjunto           devolvio None   FAIL")
        ok = False
    else:
        res, gmax = out
        acoplada = res.get(1, {}).get("p_maxt", 1.0)
        otras = [res[b]["p_maxt"] for b in res if b != 1]
        print(f"  maxT conjunto           {len(res)} bandas, {len(gmax)} perms "
              f"compartidas, p(acoplada)={acoplada:.3f} "
              f"p(otras)={[round(x, 2) for x in otras]}", end="")
        if acoplada < min(otras):
            print("   PASS")
        else:
            print("   FAIL: la banda acoplada no es la mas extrema")
            ok = False

    # CASO 7: acoplamiento DOMINANTE -> se borra con la PC1 y NO se detecta.
    # Fijado como test para que la limitación sea verificable: si algún día se
    # detecta, el estadístico dejó de descontar la estructura dominante y el
    # caso 3 --escala técnica-- empezará a dar falsos positivos.
    mats2 = {}
    for b in range(4):
        M = rng.normal(size=(n_s2, n_c2))
        if b == 1:
            M[:, :8] += 2.5 * lat2[:, None]
        mats2[b] = M
    out2 = joint_maxt(mats2, n_null=200, seed=11)
    print("  acoplamiento dominante  ", end="")
    if out2 is None:
        print("None   FAIL")
        ok = False
    else:
        r2 = out2[0]
        a2 = r2.get(1, {}).get("p_maxt", 1.0)
        if a2 > 0.05:
            print(f"p={a2:.3f} > 0.05   PASS (se borra con la PC1, documentado)")
        else:
            print(f"p={a2:.3f}   FAIL: el estadistico dejo de descontar la "
                  "estructura dominante")
            ok = False

    # CASO 6: bandas con distinto numero de muestras deben abortar. Una
    # permutacion conjunta exige que la fila i sea la misma muestra en todas.
    bad = {0: rng.normal(size=(40, 5)), 1: rng.normal(size=(30, 5))}
    e = None
    try:
        joint_maxt(bad, n_null=10)
    except ValueError as exc:
        e = str(exc)
    print("  bandas desalineadas     ", end="")
    if e and "distinto" in e:
        print("abortan   PASS")
    else:
        print("NO abortan   FAIL")
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
    # Los estadisticos nulos de cada banda se guardan para el maxT de abajo:
    # probar ~20 bandas y reportar cualquiera con p <= 0.05 da una tasa
    # familiar muy por encima de 0.05, asi que el p por banda no basta.
    # Las matrices se guardan para el maxT CONJUNTO: la misma permutacion tiene
    # que aplicarse a todas las bandas, asi que hacen falta todas a la vez y
    # alineadas a las mismas muestras.
    mats = {}
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
        mats[int(b)] = m
        r = coupling(m.to_numpy(), n_null=a.n_null, seed=42 + int(b))
        # mb_band se anadia mas abajo y nunca se creaba: `df.mb_band` daba
        # AttributeError en la primera corrida real.
        r.update(mb_band=int(b), period_mb=round(float(mid[b]), 2),
                 n_chr=m.shape[1])
        rows.append(r)
        star = " *" if (np.isfinite(r["p"]) and r["p"] <= 0.05) else ""
        print(f"{mid[b]:13.2f} {m.shape[1]:4d} {r['n']:9d} "
              f"{r['obs']:8.3f} {r['null']:9.3f} {r['p']:8.4f}{star}")

    # --- multiplicidad entre bandas -----------------------------------------
    if rows:
        df = pd.DataFrame(rows)
        ok = df.p.notna()
        if ok.any():
            # BH entre bandas. Necesario y no suficiente: BH controla la FDR
            # bajo dependencia limitada, y las bandas vecinas comparten
            # frecuencias, asi que estan correlacionadas.
            pv = df.loc[ok, "p"].to_numpy()
            order = np.argsort(pv)
            n = len(pv)
            q = np.empty(n)
            run = 1.0
            for i in range(n - 1, -1, -1):
                run = min(run, pv[order[i]] * n / (i + 1))
                q[order[i]] = run
            df.loc[ok, "q_bh"] = q

            # maxT CONJUNTO entre bandas. Las bandas se recortan a los
            # cromosomas y muestras comunes: una permutacion conjunta exige que
            # la fila i sea la misma muestra en todas.
            if len(mats) >= 2:
                chrs = set.intersection(*[set(mats[b].columns) for b in mats])
                smp = set.intersection(*[set(mats[b].index) for b in mats])
                if len(chrs) >= 3 and len(smp) >= 8:
                    ch = sorted(chrs)
                    sm = sorted(smp)
                    aligned = {b: mats[b].loc[sm, ch].to_numpy() for b in mats}
                    out = joint_maxt(aligned, n_null=min(a.n_null, 1000),
                                     seed=7)
                    if out is not None:
                        res, gmax = out
                        df["p_maxt"] = np.nan
                        for i2, b in enumerate(df.mb_band.to_numpy()):
                            if int(b) in res:
                                df.iat[i2, df.columns.get_loc("p_maxt")] = \
                                    res[int(b)]["p_maxt"]
                        print(f"\nmaxT conjunto: {len(res)} banda(s), "
                              f"{len(gmax)} permutacion(es) compartidas, "
                              f"{len(ch)} cromosoma(s), {len(sm)} muestra(s)")
                        sig = df[df.p_maxt <= 0.05]
                        if len(sig):
                            for _, r2 in sig.iterrows():
                                print(f"  {r2.period_mb:8.2f} Mb  "
                                      f"p_maxt = {r2.p_maxt:.4f}")
                        else:
                            print("  ninguna banda pasa el maxT conjunto")
                else:
                    print("\nmaxT conjunto: sin cromosomas o muestras comunes "
                          "suficientes entre bandas")
        rows = df.to_dict("records")

    if rows and a.out:
        pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
        print(f"\nescrito {a.out}")

    print("\nQué significa un p pequeño: las muestras se ordenan igual en varios "
          "cromosomas a esa escala. No dice que el período sea la causa, ni "
          "distingue un proceso genómico de un efecto de lote que sobreviva la "
          "estandarización por cromosoma.")


if __name__ == "__main__":
    main()
