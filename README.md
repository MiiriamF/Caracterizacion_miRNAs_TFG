# Caracterización miRNAs
Programas para el análisis de los motivos, número de copias y secuencia de los miRNAs predichos con MirMachine a partir del archivo de salida mirna.json obtenido al ejecutar la herramienta anot de srnatoolbox. Se ha ejecutado en un servidor remoto de Linux y con un entorno virtual conda.

## Análisis motivos y parálogos: 
- Extraccion_mirna_json.py: extracción del estado de cada motivo, así como de la secuencia de todos los miRNAs predichos en el mirna.json para las especies deseadas. Requiere 4 parámatros.
  - Lista de especies: lista de los números de acceso del NCBI para cada especie que se desea analizar.
  - Archivo de salida: nombre del archivo de salida de tipo json en el que se  recogerá la información deseada para los miRNAs predichos en las especies listadas.
  - Correo: correo electrónico que permita acceder a la base de datos del NCBI y determinar el análisis del linaje para cada especie.

- Analizar_motivos_paralogos.py: agrupación de las especies que presentan cada motivo en cada familia de miRNAs, elaboración de un árbol taxonómico con el módulo Phylo de BioPython y determinación del nodo del ancestro común para la aparición de cada motivo. Asimismo, se contabilizan las copias de cada miRNA en especies parasitiformes y el resto. Requiere 5 parámetros.
  - Archivo de salida de "Extraccion_mirna_json.py".
  - Archivo de salida 1 (nodos): archivo de salida de tipo json que mostrará cada familia de miRNAs y los motivos encontrados, así como su estado, el nodo en el que aparece y el número de copias con el mismo estado para ese motivo.
  - Archivo de salida 2 (conteo): archivo de salida de tipo json que va a recoger el conteo total de especies pertenecientes a cada nodo del árbol taxonómico, así como conteo en las diferentes clasificaciones taxonómicas para cada familia de miRNAs.
  - Archivo de salida 3 (conteo especies): archivo de salida de tipo json en el que se muestra el número de acceso del NCBI de las especies que presentan cada familia de miRNAs.
  - Archivo de salida 4 (parálogos): archivo de salida de tipo json que contendrá el número total de copias, la media aritmética y el número máximo de copias en una especie para cada familia de miRNAs en especies parasitiformes y el resto.

- Interpretar_motivos.py: se busca expresar la proporción de copias que presentan cada motivo en tanto por uno. Requiere 3 parámetros.
  - Archivo de salida 1 de "Analizar_motivos_paralogos.py".
  - Archivo de salida 2 de "Analizar_motivos_paralogos.py".
  - Archivo de salida: archivo de tipo json que va a recoger cada familia de miRNAs, motivos y si se encuentra presente o no, así como el nodo del ancestro común y el tanto por uno de copias que presentan ese estado del número total de copias de esa familia de miRNAs.

## Análisis alineamiento múltiple: 
- Multifasta.py: extracción de las secuencias y números de acceso del NCBI a partir del archivo de salida de "Extraccion_mirna_json.py" para generar un archivo MULTIFASTA (.fasta) para cada famlia de miRNAs. Requiere 1 parámetro.
  - Archivo de salida de "Extraccion_mirna_json.py".

- Muscle.py: ejecución de MUSCLE v3.8.1551 para obtener archivos _muscle.afa y posterior alineamiento con el módulo AlignIO de BioPython (_alineamiento.txt). Asimismo, se recoge en un único documento, consensos.txt, las secuencias consenso de todas las familias de miRNAs. Requiere 1 parámetro.
  - Archivo de salida de "Extraccion_mirna_json.py".

- Alineamiento.py: comparación de cada posición de cada secuencia recogida en los archivos del alineamiento múlitple con la secuencia consenso de esa familia de miRNAs. Requiere 2 parámetros.
  - Archivo de salida de "Extraccion_mirna_json.py".
  - Archivo de salida: archivo de tipo json que mostrará para cada familia de  miRNAs las posiciones en las que varía, el número de acceso de NCBI de la especie a la que pertenece esa copia y el cambio que se produce.

- Resultados_alineamiento.py: generación de un árbol taxonómico y obtención del nodo del ancestro común para cada tipo de cambio en cada posición de cada familia de miRNAs. Requiere 3 parámetros.
  - Archivo de salida de "Extraccion_mirna_json.py".
  - Archivo de salida de "Alineamiento.py".
  - Archivo de salida: nombre del archivo de salida de tipo json, el cual mostrará para cada familia de miRNAs, las posiciones en las que se dan cambios, el tipo de modificación, el nodo del ancestro común para esa mutación y el número de copias en las que se producce.

- Interpretar_alineamientos.py: se busca expresar los resultados en tanto por uno, número de copias en las que se da cada tipo de cambio respecto a las copias totales en esa familia de miRNAs. Requiere 3 parámetros.
  - Archivo de salida de "Resultados_Alineamiento.py".
  - Archivo de salida 2 de "Analizar_motivos_paralogos.py".
  - Archivo de salida: archivo de tipo json en el que se recogerá para cada familia de miRNAs las posiciones en las que hay modificación, el tipo de mutación, el nodo en el que aparece y el tanto por uno de copias que lo presentan.
