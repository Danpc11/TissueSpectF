Hallazgo estadístico principal
1. El consenso observado y su nulo no utilizan exactamente la misma estadística
El consenso observado se calcula con:
consensus_spectrum(sp, maxt = inp$maxt[[cond]])
Por tanto, su prevalencia significa “fracción de muestras donde la frecuencia fue significativa mediante maxT”.
El nulo se calcula en null_consensus_scores() sin resultados maxT:
consensus_spectrum(sub, n_boot = 0L)
En ese caso, la prevalencia significa “frecuencia por encima del percentil 95 dentro de cada espectro”.
Así, el score observado y el score nulo combinan definiciones diferentes de prevalencia. Compararlos directamente no constituye una prueba de permutación válida, aunque probablemente sea conservador.
Hay dos soluciones posibles:
- usar prevalencia por ranking tanto en observado como en nulo para decidir signature_class;
- o construir una tabla maxT conjunta y pasar al nulo los p-valores correspondientes a las muestras seleccionadas.
Prefiero la primera: mantener dos scores explícitos.
consensus_score_rank
consensus_score_maxt
Usar consensus_score_rank para la prueba permutada y consensus_score_maxt como evidencia confirmatoria secundaria.
2. El nulo devuelve un único máximo global
Cada permutación conserva sólo:
max(cs$consensus_score)
Esto controla aproximadamente el máximo sobre todas las frecuencias y cromosomas, lo cual es muy estricto. En selfcheck, los umbrales nulos quedaron alrededor de 0.93–0.94, y ninguna condición fue confirmada.
No es necesariamente incorrecto, pero conviene distinguir:
- nulo global maxT: controla errores sobre toda la firma;
- nulo por frecuencia o cromosoma: permite localizar componentes;
- q ajustado por BH sobre p-valores empíricos: opción intermedia.
Guardaría la distribución nula por (chr, k), calcularía un p-valor empírico por componente y después aplicaría BH. El máximo global podría conservarse como criterio “muy estricto”.
3. missing_blocks no garantiza bloques disjuntos
En mask_grid_genes(), los intervalos eliminados se seleccionan independientemente y pueden solaparse. Si se solapan, se eliminan menos genes que drop_total.
La cobertura real se calcula posteriormente, por lo que las bandas no quedan falseadas, pero el nombre y comentario “disjoint intervals” no siempre son ciertos.
Debe seleccionar nuevos intervalos solamente entre posiciones aún conservadas o repetir hasta alcanzar el número objetivo de genes eliminados.
4. La validación por cobertura continúa siendo global, no por clase
Cada banda obtiene un solo umbral para todas las condiciones. Sin embargo, las clases pueden tener dispersiones muy distintas. La calibración de cobertura completa sí contempla umbrales por clase, pero las bandas parciales no.
Cuando haya suficientes aciertos, generaría:
band × predicted_class × threshold
y usaría el umbral global de la banda sólo como fallback.
