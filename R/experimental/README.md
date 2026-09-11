# R - módulos experimentales

Módulos **no cargados por el core** y que **ninguna etapa del pipeline llama**.
Se cargan a mano:

```r
source("R/experimental/background.R")
source("R/experimental/pdm.R")
```

Están aquí en vez de en `R/` porque cargarlos con el core los haría parecer
parte activa del método cuando no alteran ningún resultado. Un lector del
repositorio no puede distinguir un módulo que corre de uno que sólo existe si
los dos se cargan igual.

## Pendiente que bloquea el eje bp

Las bandas de cobertura **no se calibran con `--grid-axis bp`**, y por eso una
consulta sale `UNCALIBRATED_COVERAGE`. Dos razones medidas:

- `mask_grid_genes()` opera sobre `chrom_idx`, cuyas filas en bp son bins, así
  que simula perder bins completos. Una consulta real pierde genes **dentro**
  de los bins, y eso es lo que `bin_min_coverage` decide.
- No se puede re-normalizar: `sinh(mean(asinh(v)))` es la media geométrica, no
  la suma. Para un bin de 1000/100/10/1 TPM el total es 1111 y la inversión da
  33.2 — razón 0.03.

El arreglo exige que ingest conserve los **conteos por gen** y que, por cada
máscara, se rehaga la cadena entera: conteos → quitar genes → TPM → asinh →
máscara común → cobertura y agregación por bin → espectro → huella. Es un
cambio en lo que ingest guarda, no un parche en la calibración.

Hasta entonces, el eje bp sirve para `differential` y `consensus` --que no usan
umbrales de rechazo-- y no para clasificar consultas nuevas.

### Qué le pasa a una consulta parcial en bp

| cobertura | decisión | qué significa |
|---|---|---|
| 90–100% | se clasifica con el umbral completo | la única que da una identificación |
| 50–90% | `UNCALIBRATED_COVERAGE` | **no hay umbral para esa banda**, no es un rechazo |
| <50% | `LOW_COVERAGE` | por debajo del mínimo, no se clasifica |

`UNCALIBRATED_COVERAGE` no significa "no coincide". Significa que no existe un
umbral calibrado para esa cobertura, así que la similitud calculada no se puede
convertir en una decisión. Reportarlo como si fuera un rechazo, o al revés
tomar la clase más cercana ignorando la decisión, son las dos lecturas
equivocadas.

**Una consulta parcial no es una identificación definitiva.** Con eje bp, sólo
las de 90% o más lo son, y esa restricción tiene que aparecer junto a cualquier
resultado de clasificación.

## `background.R` — fondo 1/f

`peaks_over_background()`, `mark_over_background()`, `pooled_background_null()`,
`fisher_g_test()`.

Medido y documentado: la g de Fisher tiene **100% de falsos positivos** sobre
ruido 1/f al 26% de cobertura, y el fondo por muestra los controla en 1.3% con
99% de potencia a amplitud 8. Pero **el pipeline no corrige el fondo 1/f**:
ninguna etapa llama a estas funciones. Afirmar lo contrario sería falso.

Para integrarlo haría falta decidir dónde entra --¿antes de maxT, como
prefiltro? ¿después, como columna adicional?-- y que el nulo de permutación
use el mismo fondo en el observado y en los sorteos.

## `pdm.R` — minimización de dispersión de fase y GLS ponderado

`pdm_theta()`, `pdm_spectrum()`, `gls_weighted()`, `counts_to_weights()`.

**El PDM pierde contra el GLS en los cuatro regímenes probados**, incluida la
onda cuadrada donde la teoría predecía lo contrario: 97% contra 52% con sd 2.5.
La implementación es correcta --sin ruido theta es 0 exacto en el periodo real--
así que el resultado es del método: el GLS es el estimador de máxima
verosimilitud bajo ruido gaussiano, y agrupar en bins de fase tira la
información de dónde cae cada punto. Se conserva por el resultado negativo, que
descarta el supuesto de forma con un número.

**El GLS ponderado sí gana**, 76% contra 62% en el límite de detección. No está
integrado porque no usa FFT --la FFT exige pesos uniformes-- y a O(n·m) es del
orden de un segundo por espectro: sirve para un análisis dirigido, no para las
mil permutaciones de maxT.

**El GLS ponderado sí gana**, 76% contra 62% en el límite de detección. No está
integrado porque no usa FFT --la FFT exige pesos uniformes-- y a O(n·m) es del
orden de un segundo por espectro: sirve para un análisis dirigido, no para las
mil permutaciones de maxT.
