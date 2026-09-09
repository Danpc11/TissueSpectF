# PIPELINE.md — la corrida completa, eje por eje

Tres ejes espectrales en paralelo, cada uno en su propio árbol, comparables sólo
por resultado y nunca frecuencia a frecuencia.

| eje | `--grid-axis` | árboles | qué permite |
|---|---|---|---|
| rango de gen | `gene` (default) | `interim_ncbi` / `results_ncbi` | reproduce lo ya corrido |
| 100 kb | `bp --bin-size 100000` | `interim_bp100` / `results_bp100` | resolución de TAD |
| 250 kb | `bp --bin-size 250000` | `interim_bp250` / `results_bp250` | cobertura casi completa |

**Por qué tres y no uno.** En el eje de rango de gen un periodo es un conteo de
genes, y los genes no están equiespaciados: la densidad génica varía más de diez
veces a lo largo de un cromosoma **y varía con el estado de cromatina que un
resultado espectral querría explicar** — las regiones densas son cromatina
abierta, replicación temprana, GC alto; las dispersas son LADs, replicación
tardía, desiertos génicos. El factor rango→Mb no es constante, así que
`periodo = 30 genes` no nombra una distancia y ningún mecanismo se le puede
asignar: TADs 0.1–1 Mb, dominios de replicación 0.4–0.8 Mb, LADs 0.1–10 Mb.

**Una banda que aparece a una anchura de bin y no a la otra es propiedad del
binning, no del genoma.** De ahí que se corran las dos.

---

## 0. Antes de nada

```bash
cd /scratch/home/dperez/GPIB/gene_notes/TissueSpectF
mkdir -p logs

make test
./tsf selfcheck
```

Si alguno falla, para. No tiene sentido gastar horas sobre un árbol que falla en
datos sintéticos con un pico inyectado conocido.

```bash
./tsf check --config config/project.R --geo-dir data
```

---

## Eje de 100 kb

### 1. Ingest — aquí se arma la malla

```bash
./tsf ingest GSE130970 GSE135251 GSE142530 GSE162694 GSE276114 \
  --config config/project.R \
  --geo-dir data --interim-dir interim_bp100 --results-dir results_bp100 \
  --grid-axis bp --bin-size 100000 --bin-aggregate mean \
  --cores 24 --force \
  --log logs/bp100.log
```

**PUNTO DE PARADA.** Tres líneas deciden si vale seguir, y una de ellas cambia
un parámetro del paso 4:

```
Grid axis: base pairs, 100,000 bp bins. Occupied bins per chromosome: median X%
Binned N gene(s) into M occupied bin(s) by mean (median K gene(s) per bin)
Grid coverage: median X%, range A-B%
```

- **Ocupación de bins.** Se espera ~80%. Si sale bajo el 50%, el bin es
  demasiado fino para esta densidad génica y conviene saltar directo a 250 kb.
- **Genes por bin.** Con mediana 1 la agregación casi no actúa y el eje es casi
  biyectivo con los genes. Con mediana 4 o más, cada posición es un promedio y
  la estructura corta se suavizó.
- **`Identifier map`** debe decir 100% con la anotación de NCBI.
- **`Controles`** debe reportar 6, 8, 11 y 31 por cohorte. Si GSE135251 dice 10,
  el `exclude_samples` de las dos muestras con fibrosis incidental no se aplicó.
- **La malla debe ser idéntica** en las cinco cohortes, o nada aguas abajo es
  comparable.

`--bin-aggregate mean` y no `sum`: son preguntas distintas. La media es la
actividad promedio de la región y es robusta a cuántos genes contenga; la suma
es la salida transcripcional total y **crece con la densidad génica**, que en
este eje es exactamente el confusor que el eje existe para quitar.

### 2. Spectra — ~1 min

```bash
./tsf spectra GSE130970 GSE135251 GSE142530 GSE162694 GSE276114 \
  --config config/project.R \
  --geo-dir data --interim-dir interim_bp100 --results-dir results_bp100 \
  --cores 24 --log logs/bp100.log
```

Con ~2,480 bins en chr1 en vez de 2,066 genes habrá **más** frecuencias que las
9,578 del eje de gen. Eso empeora la multiplicidad; el paso 4 lo compensa.

Confirma la columna de la que depende el piso técnico:

```bash
head -1 results_bp100/GSE135251/spectra/spectra_samples_F3.tsv \
  | tr '\t' '\n' | grep -x coverage
```

### 3. maxT y condition — ~80 y ~19 min

```bash
./tsf maxt GSE130970 GSE135251 GSE142530 GSE162694 GSE276114 \
  --config config/project.R \
  --geo-dir data --interim-dir interim_bp100 --results-dir results_bp100 \
  --cores 24 --maxt-b 1000 --seed 42 --primary-scheme all \
  --log logs/bp100.log
```

`--primary-scheme all` porque el comentario del propio config lo pide para
cualquier afirmación sobre periodicidad: con `full` el pico sólo tiene que
superar una permutación que destruye toda la autocorrelación local; con `all`
tiene que sobrevivir también los esquemas de bloque, es decir ser algo más que
correlación local.

24 núcleos y no 48: el eje de paralelismo son los cromosomas.

```bash
./tsf condition GSE130970 GSE135251 GSE142530 GSE162694 GSE276114 \
  --config config/project.R \
  --geo-dir data --interim-dir interim_bp100 --results-dir results_bp100 \
  --cores 24 --condition-b 1000 --seed 42 \
  --log logs/bp100.log
```

Corre **después** de maxT: antes, `p_stouffer` se omite y el log lo dice en cada
condición. No cambia `n_fwer` ni `n_fdr` —esas salen de `p_condition` sola— pero
pierdes una línea de evidencia independiente y `stability --criterion
consistency` deja de ser posible.

### 4. Consensus — la etapa larga

**Decide el piso biológico antes de correr.** En el eje de gen eran 10 genes.
Aquí `--min-period-biological 10` son **10 bins = 1 Mb**, un piso mucho más
agresivo que descarta todo lo sub-megabase, incluidos los TADs. Si el interés
está en escala de TAD:

```bash
./tsf consensus GSE130970 GSE135251 GSE142530 GSE162694 GSE276114 \
  --config config/project.R \
  --geo-dir data --interim-dir interim_bp100 --results-dir results_bp100 \
  --cores 48 \
  --n-null 1000 --n-boot 200 --n-contrast 1000 --seed 42 \
  --min-period auto --period-margin 2 --min-period-biological 3 \
  --log logs/bp100.log
```

`--min-period-biological 3` = 300 kb. Es una **preespecificación**: se fija antes
de ver ningún resultado y se deja escrita. Ajustarla después de mirar qué
componentes sobreviven es un procedimiento distinto con garantías distintas,
por muchas banderas que se pongan.

Aquí sí escalan los 48 núcleos: el eje son los sorteos del nulo.

Dos líneas que buscar en el log:

```
null carries a maxT-based score as well: peaks are tested with the statistic
they were selected by
```
Si no aparece, maxT no llegó y estás en la ruta de rango.

```
family N -> M frequencies: rank-1 BH diagnostic (conservative) ...
```
El diagnóstico de rango 1 es **conservador**: supone cero empates. La línea
siguiente cuenta cuántos p-valores caen en el piso de permutación y da el q
alcanzable real. Léelo antes de decidir subir `--n-contrast`: si los scores
observados no son extremos respecto al nulo, más sorteos no abren la ruta.

Salidas: `consensus_spectrum_<cond>.tsv`, `signature_<cond>.tsv`,
`condition_contrast.tsv`.

### 5. Differential — la firma característica por condición

```bash
./tsf differential GSE130970 GSE135251 GSE162694 GSE276114 \
  --config config/project.R \
  --geo-dir data --interim-dir interim_bp100 --results-dir results_bp100 \
  --cores 24 --period-bins \
  --stage-order Controles,F0,F1,F2,F3,F4 \
  --log logs/bp100.log
```

GSE142530 fuera: aporta una sola condición y no hay nada que contrastar.

**Esta etapa contesta una pregunta distinta de todo lo anterior.** `consensus`,
`condition` y `stability` preguntan **detección** —¿esta frecuencia es más fuerte
que un nulo donde se barajaron las posiciones?— y eso exige que el componente
destaque *dentro* de una condición. `differential` pregunta **comparación**: ¿la
potencia difiere *entre* condiciones? Una diferencia no necesita destacar en
ningún lado para existir, y como el espectro es una transformada lineal del
vector ordenado de expresión, la expresión diferencial implica diferencia
espectral por construcción.

Sale en `<dataset>/differential/differential_spectrum.tsv`. Filtra
`test == "one_vs_rest"` para la firma por condición, y `test == "trend"` para la
progresión monótona por estadio.

**Lee `effect` junto a `q`.** Con decenas de muestras un corrimiento trivial
alcanza significancia. El efecto mínimo detectable a 80% de potencia es ~0.85 sd
con m≈2000 y ~0.75 sd con 40 bandas: la multiplicidad es la restricción
vinculante, y `--period-bins` es lo que hace detectable un efecto de ρ≈0.3.

Nota sobre el rango: `--period-bins` usa `FINGERPRINT_PERIOD_BREAKS`, de 10 a
500 **en unidades del eje**. Con bins de 100 kb eso es de 1 Mb a 50 Mb.

### 6. Clean, stability, peaks, compare, window

```bash
./tsf run --from clean --to compare GSE130970 GSE135251 GSE142530 GSE162694 GSE276114 \
  --config config/project.R \
  --geo-dir data --interim-dir interim_bp100 --results-dir results_bp100 \
  --cores 24 --seed 42 --log logs/bp100.log

./tsf window GSE130970 GSE135251 GSE142530 GSE162694 GSE276114 \
  --config config/project.R \
  --geo-dir data --interim-dir interim_bp100 --results-dir results_bp100 \
  --cores 24 --log logs/bp100.log
```

`window` **después** de `stability`, o no tiene picos estables que colocar y
reporta `no stable peak to place in the window`, que no informa de nada. Los 24
archivos `window_chr*.tsv` se escriben igual, pero la tabla que cruza picos con
la ventana necesita `stability`.

### 7. Reference — la base de espectros y el control

```bash
for f in amplitude period_bins_genomic band_ratios expression_baseline; do
  ./tsf reference GSE130970 GSE135251 GSE142530 GSE162694 GSE276114 \
    --config config/project.R \
    --geo-dir data --interim-dir interim_bp100 --results-dir results_bp100_$f \
    --cores 24 --target class_id --seed 42 \
    --features $f --n-features 200 \
    --log logs/bp100_feat_$f.log
done

grep "Out-of-cohort accuracy" logs/bp100_feat_*.log
```

`--target class_id`, no `condition`. `condition` guarda la etiqueta **cruda**,
antes de aplicar el mapeo del vocabulario —ese mapeo sólo construye `class_id`—
así que con `condition` una fusión de clases no tiene efecto y `Controles`
vuelve a partirse en tres clases de una cohorte cada una, imposibles de aprender
por leave-one-cohort-out.

`--n-features 200` **idéntico en las cuatro**: cambiarlo en una invalida la
comparación.

**`expression_baseline` es el que interpreta las otras tres.** Es la expresión
génica en la misma malla, sin transformada. Si clasifica mejor que el espectro,
la transformada descarta información y hay que decirlo; si clasifica peor, el
espectro es una compresión real.

### 8. Clasificadores

```bash
pip install pycatch22 sktime

python3 scripts/classify_spectra.py \
  --results-dir results_bp100 \
  --datasets GSE130970,GSE135251,GSE162694,GSE276114 \
  --n-bins 200 \
  --out results_bp100/classifier_comparison.tsv
```

200 bandas y no 40: la clasificación **no corrige por multiplicidad** —alimenta
un modelo, no prueba m hipótesis— así que la resolución es información, y 40
puntos son demasiado cortos para que los kernels dilatados tengan rango de
escalas.

`minirocket` es el **techo alcanzable**, no el modelo final: da exactitud alta y
ninguna interpretabilidad. Si con un método que las revisiones consideran estado
del arte tampoco supera el baseline de expresión, la información no está en el
espectro — y cerrar la pregunta así es mucho más difícil de objetar que cerrarla
con un centroide.

---

## Eje de 250 kb

Idéntico, cambiando tres cosas en todos los comandos:

```
--bin-size 250000
--interim-dir interim_bp250
--results-dir results_bp250
--log logs/bp250.log
```

Y `--min-period-biological`: 3 bins a 250 kb son 750 kb, no 300 kb. Para
mantener el mismo piso físico usa `--min-period-biological 2` (500 kb) o
acéptalo distinto y decláralo.

---

## Eje de rango de gen

Ya corrido en `results_ncbi`. Para reproducirlo, los mismos comandos **sin**
`--grid-axis` ni `--bin-size` — `gene` es el default — con
`--interim-dir interim_ncbi --results-dir results_ncbi` y
`--min-period-biological 10`.

---

## Empaquetar los resultados

```bash
OUT=tsf_resultados_$(date +%Y%m%d); mkdir -p "$OUT"
cp logs/*.log "$OUT/" 2>/dev/null
for tree in results_ncbi results_bp100 results_bp250; do
  [ -d "$tree" ] || continue
  for d in "$tree"/*/; do
    ds=$(basename "$d")
    [ -f "$d/differential/differential_spectrum.tsv" ] &&
      cp "$d/differential/differential_spectrum.tsv" "$OUT/${tree}_${ds}_differential.tsv"
    [ -f "$d/consensus/condition_contrast.tsv" ] &&
      cp "$d/consensus/condition_contrast.tsv" "$OUT/${tree}_${ds}_contrast.tsv"
  done
  {
    echo -e "dataset\tcondition\tchr\tN\tk\tperiod\tmedian_power_normalised\tprevalence\tplv\tp_null_fwer\tp_null_fwer_maxt\tq_null"
    for f in "$tree"/*/consensus/consensus_spectrum_*.tsv; do
      [ -f "$f" ] || continue
      ds=$(echo "$f" | cut -d/ -f2)
      cn=$(basename "$f" .tsv | sed 's/consensus_spectrum_//')
      awk -F'\t' -v d="$ds" -v c="$cn" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
        {print d"\t"c"\t"$h["chr"]"\t"$h["N"]"\t"$h["k"]"\t"$h["period"]"\t"\
         $h["median_power_normalised"]"\t"$h["prevalence"]"\t"$h["plv"]"\t"\
         $h["p_null_fwer"]"\t"$h["p_null_fwer_maxt"]"\t"$h["q_null"]}' "$f"
    done
  } | gzip > "$OUT/${tree}_espectros.tsv.gz"
done
for r in results_*_*/reference results_ncbi/reference; do
  [ -d "$r" ] || continue
  tag=$(echo "$r" | cut -d/ -f1)
  cp "$r/confusion_matrix.tsv" "$OUT/${tag}_confusion.tsv" 2>/dev/null
  cp "$r/out_of_cohort_predictions.tsv" "$OUT/${tag}_predictions.tsv" 2>/dev/null
done
tar czf "$OUT.tar.gz" "$OUT" && du -sh "$OUT.tar.gz"
```

---

## Limpieza

```bash
make clean-dry                                  # lista, no borra
Rscript scripts/clean_results.R --results-dir results_bp100
```

Pide teclear el nombre del árbol. `results_dir` apunta a trabajo real —un
`consensus` son decenas de minutos y no es reproducible desde lo que queda— así
que reporta el conteo y el tamaño antes de preguntar, y un stdin cerrado cuenta
como negativa para que un job no borre nada por accidente.

**No uses `make clean` con los árboles nuevos:** lee `config/project.R`, que no
apunta a `results_bp100`.

---

## Lo que no cambia con ninguna corrida

**Sin negativos.** La calibración acota con qué frecuencia se rechaza a un
miembro verdadero. No acota con qué frecuencia se **acepta** a alguien de fuera
del dominio, porque no hubo negativos en la validación. Para una referencia de
hígado eso significa otros tejidos: `config/datasets/GTEx.R` y
`scripts/gtex_subset.R` están listos y sin correr.

**Condiciones estructuralmente excluidas.** La alineación de fase necesita ~11
muestras: GSE130970 F4 (n=2), Normal_histology (6), GSE162694 F3 (8),
GSE130970 F2 (9), Control_disease_cohort (8). Cinco de 21, por diseño y no por
ausencia de señal.

**GSE142530 no puede confirmar nada.** Once muestras y la condición son las
once, así que no hay subconjunto que sortear y no existe nulo. Se salta como
hold-out en toda validación cruzada.

**`N` termina en el último gen anotado, no al final del cromosoma.** Decisión
deliberada, porque la anotación cargada no trae longitudes de cromosoma. En
acrocéntricos —13, 14, 15, 21, 22— el brazo corto es casi todo heterocromatina
sin genes, así que el eje empieza tarde y `N` subestima el cromosoma. Sin efecto
para periodos de pocas Mb; con efecto para los de decenas.

**Un gen que cruza un límite de bin cuenta sólo donde empieza.** Asignación por
posición de inicio, o sea por promotor, que es defendible —es donde la cromatina
decide— pero es una decisión y va declarada.

**`differential` no es una prueba de detección.** Establece que la potencia
difiere entre estadios, no que esas bandas contengan periodicidad detectable
dentro de ninguna condición. Son afirmaciones independientes.
