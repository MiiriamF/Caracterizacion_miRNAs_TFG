import subprocess
import sys
import os
import copy
import json
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade

def interpretacion(nombre_json, resultados_json, nodos_json):
    with open(nombre_json, "r") as f:
        rutas = json.load(f)

        #Buscar nodo en el que hay modificaciones en cada miRNA
        tree = Phylo.read("arbol_taxonomico.xml", "phyloxml")
        with open(resultados_json, "r") as f:
            resultado = json.load(f)
            nodos = {}
            for mir,inf in resultado.items():
                nodos[mir] = {}
                for pos,esp in inf.items():
                    cambios = {}
                    nodos[mir][pos] = {}
                    for gcf, tipo in esp.items():
                        if tipo not in cambios.keys():
                            cambios[tipo] = []
                        cambios[tipo].append(gcf)
                    for k in cambios.keys():
                        nodo_comun = tree.common_ancestor(cambios[k])
                        nodos[mir][pos][k] = [nodo_comun.name,len(cambios[k])]

    with open(nombre_json, "w") as f:
        json.dump(rutas, f, indent=4)

    with open(resultados_json, "w") as f:
        json.dump(resultado, f, indent=4)

    with open(nodos_json, "w") as f:
        json.dump(nodos, f, indent=4)
    
interpretacion(sys.argv[1], sys.argv[2], sys.argv[3])
