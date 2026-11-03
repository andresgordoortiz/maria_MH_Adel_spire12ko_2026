# SPIRE y FMN2: revisión desde los FASTQ originales

La inferencia del manuscrito se actualizó el 3 de noviembre de 2026; véase
[la auditoría de GO, fondos y corrección](PUBLICATION_INFERENCE_REVIEW.md).
Los apartados de reconstrucción siguientes documentan las salidas anteriores.
La referencia de la revisión de muestras es el
[informe de revisión SPIRE](../reports/spire_review/REPORT.md).
Las conclusiones anteriores se conservan en
[el archivo histórico](../legacy/pre_fastq_2026-11-02/docs/CONCLUSIONS.md).

## Las etiquetas antiguas no coinciden con las de los FASTQ

La correspondencia entre los seis perfiles es unívoca usando 4.316 genes variables,
un evento por gen y excluyendo Spire1/2. Coinciden los mejores emparejamientos por
RMSE y Spearman:

| Etiqueta antigua | Perfil de los FASTQ originales |
|---|---|
| Cont_a | Cont_1 |
| Cont_b | Cont_2 |
| Cont_c | DKO_4 |
| DKO_a | DKO_5 |
| DKO_b | DKO_6 |
| DKO_c | Cont_7 |

Esto identifica una inconsistencia entre las etiquetas históricas y las de los
FASTQ. No establece en que paso ocurrio ni sustituye un genotipado de ADN.
La nueva asignación es consistente con la evidencia de RNA en los loci diana:
9 de 13 eventos exónicos observados de Spire1 y los 6 de Spire2 tienen soporte
VAST>=10 en los tres controles y <10 en los tres DKO. Otros cuatro exones de
Spire1 conservan soporte en DKO; por eso la expresión total y la cobertura de
los exones afectados deben revisarse por separado.

## Cambia mucho la lista de eventos, aunque el total sea parecido

Se mantienen los criterios: N>=10 y PSI completa en las seis muestras,
FPR empirico betAS<=0,05 y |delta PSI ajustado|>=10 puntos porcentuales.
El campo FDR de betAS es un FPR por evento; no es FDR corregida por pruebas múltiples.

| Medida | SPIRE anterior | SPIRE desde FASTQ |
|---|---:|---:|
| Eventos evaluados | 289.583 | 287.393 |
| Eventos seleccionados | 1.436 | 1.547 |
| Genes afectados con nombre | 945 | 1.000 |
| EX / IR / Alt5 / Alt3 seleccionados | 537 / 241 / 293 / 365 | 551 / 271 / 330 / 395 |
| Uso aumentado / disminuido | 711 / 725 | 784 / 763 |
| Retención de intrones aumentada / disminuida | 127 / 114 | 139 / 132 |

En los 286.899 eventos comunes, con anotación identica, 358 mantienen selección
y dirección, 102 invierten dirección, 967 dejan de seleccionarse y 1.080 pasan
a seleccionarse. Otros 9 eventos seleccionados antes ya no son evaluables;
7 seleccionados ahora no se evaluaban antes. Por tanto, la lista antigua de
candidatos no debe seguir usandose como lista actual. Los cambios incluyen
reprocesamiento y reasignacion de condiciones; no se atribuyen exclusivamente
al etiquetado. Las listas completas y las transiciones se exportan en
`reports/spire_review/tables/`.

No se observa un predominio general de aumento de retención de intrones.
Spire1 MmuINT0151614 sigue seleccionado y aumenta a +91,7 pp de delta ajustado
(Pdiff=1; PSI cruda 4,49-6,93% en controles frente a 76-100% en DKO).
Es un evento del locus manipulado y requiere esa cautela al interpretarlo;
no demuestra por si solo autoregulacion. Spire2 queda sin eventos evaluables
tras el filtro que exige N10 también en los DKO.

## Expresión y heterogeneidad

En las salidas DESeq2 suministradas, contraste DKO/Cont, hay 111 genes con
padj<=0,05: 40 aumentan y 71 disminuyen. De ellos, 50 cumplen ademas
|log2FC|>=1. Se conservan los log2FC exportados con shrinkage y sus p-valores;
no se ha vuelto a ejecutar DESeq2.

- Spire2: log2FC=-7,749; padj=0,000382. Los controles tienen aproximadamente
  1.959-2.460 cuentas normalizadas y los DKO 0-26.
- Spire1: log2FC=-0,177; padj=0,207. El soporte residual de algunos exones
  coexiste con la pérdida casi completa de otros.
- Fmn2: log2FC=+0,111; padj=0,419.

El investigador confirmó el 2 de noviembre de 2026 que las seis muestras SPIRE
proceden de animales independientes (tres por condición). La asignación coherente
de condiciones no elimina la heterogeneidad. DKO_4 y
Cont_7 siguen separados de las demas muestras en las PCA de expresión y splicing,
que excluyen los loci Spire1/2. No aparecen dos grupos compactos. STAR tiene
86,26-88,98% de alineamiento unico; ese QC no explica por si solo los ejes de
variación. Cont_7 lleva la etiqueta de fecha 14mars y las otras cinco muestras,
12mars; no se ha supuesto un lote biológico a partir de esos nombres.
No se ha descartado ninguna muestra ni aplicado una correccion de lote.

## Funcion y comparación con FMN2

SPIRE conserva la sensibilidad de GO al fondo: 106 términos BP con el fondo
OrgDb completo y ninguno con el universo de genes evaluados N10. También quedan
sin términos CC/MF significativos con ese universo. Los temas obtenidos solo
con el fondo amplio no deben presentarse como vias confirmadas.

La comparación actual usa 48.323 eventos en 4.229 genes evaluados en ambos
experimentos: 11 eventos seleccionados comunes (7 concordantes y 4 opuestos),
frente a 18 antes, y 22 genes afectados comunes, frente a 25 antes. Hay más
solapamiento que bajo un nulo uniforme de genes (22 frente a 3,44 esperados;
Fisher p=6,22e-13), que no corrige longitud de gen ni número de eventos por gen.
La concordancia de dirección no supera de forma convincente el 50%
(prueba de signo exploratoria p=0,549); la correlación de deltas medios crudos
en todos los eventos comunes es 0,0065.

Ktn1 se mantiene como evento concordante. Huwe1, Plag1, Pla2g4c y Zmynd11
aparecen entre los concordantes actuales. Rev3l ya no cumple la selección
compartida. Las diferencias de controles, cobertura y replicacion impiden
tratar los knockouts como experimentos intercambiables. FMN2 conserva dos pools
por grupo y controles heterocigotos. Ver las listas, PSI individuales, GO y
características de secuencia en [el informe comparativo](../reports/comparison/REPORT.md).
Ninguno de los siete eventos concordantes compartidos presenta separacion
completa de las PSI individuales en ambos experimentos; siguen siendo candidatos.

No hay términos GO enriquecidos en los conjuntos compartidos ni en cada estudio
al usar ese mismo universo común, en BP/CC/MF.

Se mantiene la asociacion descriptiva entre mayor inclusión exónica y flancos
intronicos más largos en ambos estudios. En SPIRE, las medianas de los eventos
EX-up frente a ndiff son 2.118 frente a 1.444 nt aguas arriba y 1.838 frente a
1.376,5 nt aguas abajo; ambas comparaciones pasan BH<0,01. Son características
genómicas de esos eventos, no intrones que se alarguen tras el knockout.
En SPIRE el mismo patrón de intrones flanqueantes más largos aparece también
en EX-down; por tanto, no distingue de forma específica los aumentos de inclusión.

## RBP: mapas recalculados con las condiciones corregidas

rMAPS ha terminado EX, IR, Alt5 y Alt3: 114 motivos bajo 91 etiquetas RBP,
con todos los eventos de fondo preparados. La auditoría valida 741.000 p-valores
de ventanas y 456 pares PDF/PNG por estudio; los recuentos independientes de
las ventanas más extremas coinciden con los resultados nativos.

El criterio de p crudo<=0,05 en la misma etiqueta, región y dirección en ambos
estudios selecciona 82/91 etiquetas para IR-up y 84/91 para IR-down; para EX-up
son 86/91 y EX-down 82/91. Es una interseccion demasiado amplia para identificar
un regulador común específico. Un mismo nombre puede ganar con variantes de
motivo y ventanas diferentes.

Los mínimos puntuales de SPIRE corresponden a motivos RBM4 en EX-down e IR-up,
y RBM42 en Alt5-up y Alt3-up. Su soporte es escaso: cuatro eventos con hits entre
271 en el ejemplo EX, tres entre 138 en IR, tres entre 83 en Alt5 y tres entre
114 en Alt3. En la ventana binaria Alt5, el p nativo es 5,96e-12 y el diagnóstico
condicional exacto 0,000719; en Alt3, 2,51e-10 frente a 0,00113. Ninguno de esos
diagnósticos corrige todas las búsquedas entre motivos y ventanas. Por tanto,
RBM4/RBM42 son etiquetas de motivos candidatas, no reguladores causales establecidos.

Las figuras 2H y S6A–E y las extensiones están en la
[galería RBP actual](../reports/rbp/index.html). Los recuentos, ventanas y limites
se detallan en [su informe](../reports/rbp/REPORT.md). La expresión o abundancia
de un RBP no se deduce de un enriquecimiento de su motivo.

## Alcance y archivos

Los FASTQ no se han vuelto a alinear aquí: se utilizan los resultados nuevos
suministrados de RNA-seq, DGE, VAST y betAS. Se han validado PSI, anotaciones,
orden de eventos, contraste, muestras y estadísticas contra los checkpoints,
y regenerado figuras, GO y características Matt/GC con las entradas nuevas.
Los originales y los resultados numericos de FMN2 se conservan.

Matt/GC y rMAPS están completos con las listas corregidas; la Figura 2 integra
los seis paneles de splicing y Matt. Se han retirado solo las secuencias generadas
del intento interrumpido, conservando su procedencia y la ejecución completa.
Los originales y los informes anteriores siguen archivados.

La conclusión computacional es una alteración selectiva y bidireccional del
splicing, con algunos genes y características de secuencia compartidos con FMN2.
No se establece una respuesta direccional global común, una vía GO robusta al
fondo evaluado ni un RBP causal común. La heterogeneidad entre muestras, el FPR
por evento de betAS y la calibración de los motivos limitan la inferencia.

## Figuras para el manuscrito

La versión actual está en `results/publication/index.html`: dos figuras SPIRE
y siete suplementos; una figura FMN2 y cinco suplementos; una comparación
compacta. GO usa genes evaluables, mapeo directo/alias único, exclusión de
los loci manipulados y BH sobre BP/CC/MF. La unión de división contiene 127
eventos en 79 genes, sin sobrerrepresentación: 47 genes de división frente
a 47,4 esperados. La presencia regional exacta de motivos pasa BH para
PCBP2/EX y PPRC1/Alt5 en SPIRE; ninguna asociación pasa en FMN2. Los mapas
nativos de mínimos asintóticos quedan como referencia. Las familias de
pruebas, métodos, resultados y límites están documentados en el paquete.
