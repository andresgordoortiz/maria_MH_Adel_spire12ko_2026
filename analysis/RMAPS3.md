# rMAPS3 y Matt: SPIRE con los grupos de los FASTQ originales

SPIRE usa Cont_1/2/7 frente a DKO_4/5/6. Matt está completo y validado;
rMAPS ha terminado los cuatro tipos y ha pasado la auditoría. La galería conjunta
está en `reports/rbp/index.html`. FMN2 conserva sus mapas originales.

## Reproducir rMAPS

Ejecutar desde la raíz, con las entradas de la ejecución completa de SPIRE
y destinos nuevos. El driver comprueba hashes, manifiestos, geometria de
coordenadas y unicidad de eventos antes de iniciar los calculos.

```bash
python3 -B tests/rmaps3_checks.py
PYTHON_CPU_COUNT=6 pixi run --manifest-path ../rMAPS3/pixi.toml --locked python analysis/run_rmaps3.py \
  reports/spire_final spire_fastq reports/spire_final/rmaps3_outputs ../rMAPS3
pixi run --manifest-path ../rMAPS3/pixi.toml --locked python -B tests/rmaps3_output_checks.py spire_fastq
bash scripts/run_local_R.sh analysis/render_rmaps3.R
```

Se reutiliza el entorno Pixi aislado existente; no se instalan paquetes en las
bibliotecas R/betAS ni del sistema. La ejecución usa cinco trabajadores de motivos
(`PYTHON_CPU_COUNT=6`) y un hilo BLAS/OpenMP por proceso. La cantidad de trabajadores
no cambia la definicion del contraste ni los eventos de fondo.

El contador corregido se mantiene en el repositorio vecino `../rMAPS3`, rama
`fix/motif-window-counting`, commit `1076240`. Los hashes de las fuentes, referencias,
entradas y coordenadas convertidas, los comandos y versiones quedan registrados
en `provenance.json`. El renderer contrasta el dataset y los hashes de entrada
con cada análisis para impedir combinar las listas nuevas con mapas antiguos.

## Métodos y comprobaciones

Se ejecutan EX, IR, Alt5 y Alt3 sobre mm10, 114 motivos bajo 91 etiquetas RBP,
extensiones exon/intron de 50/250 nt, ventanas de 50 nt y pasos de 1 nt,
con Mann–Whitney unilateral. Se conservan todos los eventos ndiff preparados
como fondo. Las coordenadas VASTDB cerradas de base 1 se convierten a intervalos
semiabiertos de base 0; las pruebas cubren ambas hebras y los sitios alternativos.
Up/down indica mayor/menor uso del evento en KO, no expresión de un RBP.

La auditoría comprueba 741.000 p-valores de ventanas por dataset, mínimos regionales,
orden y contenido de las secuencias, y 456 pares PDF/PNG. Recuenta de forma
independiente la ventana más extrema de cada tipo; para recuentos binarios
exporta también la probabilidad condicional exacta como diagnóstico.

Las figuras conjuntas incluyen 2H y S6A–E, junto con extensiones de sitios
alternativos y recuentos. Mantienen las selecciones, regiones, orden y colores
del análisis anterior. El suavizado LOESS usa span 0,75 porque el 0,3 depositado
era singular con solo cinco u ocho posiciones regionales. Las curvas describen
puntuaciones regionales, no posiciones de union experimentales.

Los mínimos entre motivos y ventanas son p-valores crudos exploratorios.
Los recuentos escasos, eventos dependientes, variantes de motivos y búsquedas
múltiples limitan su especificidad. Un motivo común no demuestra union,
abundancia ni un regulador causal común. Las tablas exportadas conservan los
p-valores completos; la saturacion de las figuras solo afecta a la visualizacion.

## Matt y archivos

Matt ya se ejecuto con las listas nuevas y el contenedor del trabajo original,
GENCODE vM10, GRCm38/mm10 y un fondo ndiff de hasta 10.000 eventos, semilla 42.
Los inputs y referencias se montan en modo lectura. Resultados y hashes están
en `data/processed/spire_fastq/matt/`; las 73 figuras de características/GC están
en `reports/spire_final/`.

```bash
python3 -B analysis/validate_matt.py reports/spire_final data/processed/spire_fastq/matt
```

Los resultados anteriores de SPIRE, las figuras conjuntas anteriores y sus
métodos se conservan en `legacy/pre_fastq_2026-11-02/`. Las auditorias de la
finalizacion actual están en `audit/spire_reanalysis_2026-11-02/rmaps_completion/`.
Los ejemplos antiguos de soporte de motivos siguen ligados a sus entradas
históricas; no se reutilizan como evidencia de las nuevas listas.
