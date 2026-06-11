import subprocess
import sys
import os
import copy
import json
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade

def interpretacion(nombre_json, nodos_json, conteo_json, conteo_especies_json, paralogos_json):
    with open(nombre_json, "r") as f:
        rutas = json.load(f)
        dic_arbol = {}
        conteo = {}
        conteo_esp = {}
        conteo["organismos"] = {}
        conteo["mirnas"] = {}
        conteo["mirnas_por_especies"] = {}
        conteo["mirnas_alin"] = {}
        conteo["mirnas_alin_especies"] = {}
        paralogos = {}

        #Creación del árbol taxonómico
        for e,inf in rutas.items():
            dic_arbol[e] = rutas[e]["linaje"]

        root = Clade(name="Root")
        tree = Tree(root=root)

        for gcf, linaje_str in dic_arbol.items():
            niveles = [n.strip() for n in linaje_str.split(";") if n.strip()]
            
            actual = root
            for nivel in niveles:
                if nivel not in conteo["organismos"].keys():
                    conteo["organismos"][nivel] = 0
                conteo["organismos"][nivel] += 1
                encontrado = next((c for c in actual.clades if c.name == nivel), None)
                
                if encontrado:
                    actual = encontrado
                else:
                    nueva_rama = Clade(name=nivel)
                    actual.clades.append(nueva_rama)
                    actual = nueva_rama
            
            hoja = Clade(name=gcf)
            actual.clades.append(hoja)

        Phylo.write(tree, "arbol_taxonomico.xml", "phyloxml")
        #Abrir con Phylo.io

        #Buscar nodo en el que hay presencia de motivos en cada miRNA
        tree = Phylo.read("arbol_taxonomico.xml", "phyloxml")
        nodos = {}
        cambios = {}
        for e,inf in rutas.items():
            #Estudio de los motivos
            for mirna, inf_2 in rutas[e]["motifs"].items():
                if mirna.count("_") >= 1:
                    mirna_familia = mirna.split("_",1)[0]
                elif mirna.endswith("a") or mirna.endswith("b") or mirna.endswith("c"):
                    mirna_familia = mirna[:-1]
                else: 
                    mirna_familia = mirna
                if  mirna_familia.count("-") > 1:
                    mirna_familia = mirna.split("-")[0] + "-" + mirna.split("-")[1]
                if mirna_familia.count("Bantam") == 1:
                    mirna_familia = "Bantam"
                #Estudio de los cambios
                if mirna_familia not in nodos.keys():
                    nodos[mirna_familia] = {}
                    cambios[mirna_familia] = {}
                    paralogos[mirna_familia] = {}
                if "motifs" not in nodos[mirna_familia].keys():
                    nodos[mirna_familia]["motifs"] = {}
                    cambios[mirna_familia]["motifs"] = {}
                for motivo,status in inf_2.items():
                    if motivo not in nodos[mirna_familia]["motifs"].keys():
                        nodos[mirna_familia]["motifs"][motivo] = {}
                        cambios[mirna_familia]["motifs"][motivo] = {}
                    if status not in nodos[mirna_familia]["motifs"][motivo].keys():
                        nodos[mirna_familia]["motifs"][motivo][status] = []
                        cambios[mirna_familia]["motifs"][motivo][status] = []
                    cambios[mirna_familia]["motifs"][motivo][status].append(e)
            #Estudio de los motivos export_3p
            for mirna, status in rutas[e]["export_3p"].items():
                if mirna.count("_") >= 1:
                    mirna_familia = mirna.split("_",1)[0]
                elif mirna.endswith("a") or mirna.endswith("b") or mirna.endswith("c"):
                    mirna_familia = mirna[:-1]
                else: 
                    mirna_familia = mirna
                if  mirna_familia.count("-") > 1:
                    mirna_familia = mirna.split("-")[0] + "-" + mirna.split("-")[1]
                if mirna_familia.count("Bantam") == 1:
                    mirna_familia = "Bantam"
                #Estudio de los cambios
                if mirna_familia not in nodos.keys():
                    nodos[mirna_familia] = {}
                    cambios[mirna_familia] = {}
                    paralogos[mirna_familia] = {}
                if "export_3p" not in nodos[mirna_familia].keys():
                    nodos[mirna_familia]["export_3p"] = []
                    cambios[mirna_familia]["export_3p"] = []
                cambios[mirna_familia]["export_3p"].append(e)
                
            #Estudio de los motivos export_5p
            for mirna, status in rutas[e]["export_5p"].items():
                if mirna.count("_") >= 1:
                    mirna_familia = mirna.split("_",1)[0]
                elif mirna.endswith("a") or mirna.endswith("b") or mirna.endswith("c"):
                    mirna_familia = mirna[:-1]
                else: 
                    mirna_familia = mirna
                if  mirna_familia.count("-") > 1:
                    mirna_familia = mirna.split("-")[0] + "-" + mirna.split("-")[1]
                if mirna_familia.count("Bantam") == 1:
                    mirna_familia = "Bantam"
                #Estudio de los cambios
                if mirna_familia not in nodos.keys():
                    nodos[mirna_familia] = {}
                    cambios[mirna_familia] = {}
                    paralogos[mirna_familia] = {}
                if "export_5p" not in nodos[mirna_familia].keys():
                    nodos[mirna_familia]["export_5p"] = []
                    cambios[mirna_familia]["export_5p"] = []
                cambios[mirna_familia]["export_5p"].append(e)
                
            #Estudio de los motivos loop_motifs
            for mirna, inf_2 in rutas[e]["loop_motifs"].items():
                if mirna.count("_") >= 1:
                    mirna_familia = mirna.split("_",1)[0]
                elif mirna.endswith("a") or mirna.endswith("b") or mirna.endswith("c"):
                    mirna_familia = mirna[:-1]
                else: 
                    mirna_familia = mirna
                if  mirna_familia.count("-") > 1:
                    mirna_familia = mirna.split("-")[0] + "-" + mirna.split("-")[1]
                if mirna_familia.count("Bantam") == 1:
                    mirna_familia = "Bantam"
                #Estudio de los cambios
                if mirna_familia not in nodos.keys():
                    nodos[mirna_familia] = {}
                    cambios[mirna_familia] = {}
                    paralogos[mirna_familia] = {}
                if "loop_motifs" not in nodos[mirna_familia].keys():
                    nodos[mirna_familia]["loop_motifs"] = {}
                    cambios[mirna_familia]["loop_motifs"] = {}
                for motivo,status in inf_2.items():
                    if motivo not in nodos[mirna_familia]["loop_motifs"].keys():
                        nodos[mirna_familia]["loop_motifs"][motivo] = {}
                        cambios[mirna_familia]["loop_motifs"][motivo] = {}
                    if status not in nodos[mirna_familia]["loop_motifs"][motivo].keys():
                        nodos[mirna_familia]["loop_motifs"][motivo][status] = []
                        cambios[mirna_familia]["loop_motifs"][motivo][status] = []
                    cambios[mirna_familia]["loop_motifs"][motivo][status].append(e)

        for e,inf in rutas.items():
            mirnas = {}
            if "Parasitiformes" in rutas[e]["linaje"]:
                for mirna in rutas[e]["motifs_gen"].keys():
                    #miRNAs familia
                    if mirna.count("_") >= 1:
                        mirna_familia = mirna.split("_",1)[0]
                    elif mirna.endswith("a") or mirna.endswith("b") or mirna.endswith("c"):
                        mirna_familia = mirna[:-1]
                    else:
                        mirna_familia = mirna
                    if  mirna_familia.count("-") > 1:
                        mirna_familia = mirna.split("-")[0] + "-" + mirna.split("-")[1]
                    if mirna_familia.count("Bantam") == 1:
                        mirna_familia = "Bantam"
                    #Adición datos
                    if mirna_familia not in mirnas:
                        mirnas[mirna_familia] = 0
                    if "parasitiformes" not in paralogos[mirna_familia].keys():
                        paralogos[mirna_familia]["parasitiformes"] = 0
                        paralogos[mirna_familia]["parasitiformes_especies"] = []
                        paralogos[mirna_familia]["parasitiformes_max"] = 0
                    paralogos[mirna_familia]["parasitiformes"] += 1
                    mirnas[mirna_familia] += 1
                    if e not in paralogos[mirna_familia]["parasitiformes_especies"]:
                        paralogos[mirna_familia]["parasitiformes_especies"].append(e)
                #max mirna
                for mirna_fam, inf in mirnas.items():
                    if inf > paralogos[mirna_fam]["parasitiformes_max"]:
                        paralogos[mirna_fam]["parasitiformes_max"] = inf 
            else:
                for mirna in rutas[e]["motifs_gen"].keys():
                    #miRNAs familia
                    if mirna.count("_") >= 1:
                        mirna_familia = mirna.split("_",1)[0]
                    elif mirna.endswith("a") or mirna.endswith("b") or mirna.endswith("c"):
                        mirna_familia = mirna[:-1]
                    else:
                        mirna_familia = mirna
                    if  mirna_familia.count("-") > 1:
                        mirna_familia = mirna.split("-")[0] + "-" + mirna.split("-")[1]
                    if mirna_familia.count("Bantam") == 1:
                        mirna_familia = "Bantam"
                     #Adición datos
                    if mirna_familia not in mirnas:
                        mirnas[mirna_familia] = 0
                    if "resto" not in paralogos[mirna_familia].keys():
                        paralogos[mirna_familia]["resto"] = 0
                        paralogos[mirna_familia]["resto_especies"] = []
                        paralogos[mirna_familia]["resto_max"] = 0
                    paralogos[mirna_familia]["resto"] += 1
                    mirnas[mirna_familia] += 1
                    if e not in paralogos[mirna_familia]["resto_especies"]:
                        paralogos[mirna_familia]["resto_especies"].append(e)
                #max mirna
                for mirna_fam, inf in mirnas.items():
                    if inf > paralogos[mirna_fam]["resto_max"]:
                        paralogos[mirna_fam]["resto_max"] = inf

        for mirna,inf in paralogos.items():
            if "parasitiformes" in inf.keys() and "parasitiformes_especies" in inf.keys():
                paralogos[mirna]["parasitiformes_media"] = paralogos[mirna]["parasitiformes"]/len(paralogos[mirna]["parasitiformes_especies"])
            if "resto" in inf.keys() and "resto_especies" in inf.keys():
                paralogos[mirna]["resto_media"] = paralogos[mirna]["resto"]/len(paralogos[mirna]["resto_especies"])


        for mirna,inf in cambios.items():
            sum = 0
            list = []
            conteo["mirnas_por_especies"][mirna] = {}
            for tipos,inf_2 in inf.items():
                if tipos == "export_3p":
                    nodo_comun = tree.common_ancestor(inf_2)
                    nodos[mirna]["export_3p"] = [nodo_comun.name,len(inf_2)]
                    for i in inf_2:
                        if i not in list:
                            list.append(i)
                            sum += 1
                            niveles = [n.strip() for n in rutas[i]["linaje"].split(";") if n.strip()]
                            for nivel in niveles:
                                if nivel not in conteo["mirnas_por_especies"][mirna].keys():
                                    conteo["mirnas_por_especies"][mirna][nivel] = 0
                                conteo["mirnas_por_especies"][mirna][nivel] += 1
                elif  tipos == "export_5p":
                    nodo_comun = tree.common_ancestor(inf_2)
                    nodos[mirna]["export_5p"] = [nodo_comun.name,len(inf_2)]
                    for i in inf_2:
                        if i not in list:
                            list.append(i)
                            sum += 1
                            niveles = [n.strip() for n in rutas[i]["linaje"].split(";") if n.strip()]
                            for nivel in niveles:
                                if nivel not in conteo["mirnas_por_especies"][mirna].keys():
                                    conteo["mirnas_por_especies"][mirna][nivel] = 0
                                conteo["mirnas_por_especies"][mirna][nivel] += 1
                else:
                    for motivo,inf_3 in inf_2.items():
                        if tipos == "motifs":
                            for status, e in inf_3.items(): 
                                nodo_comun = tree.common_ancestor(e)
                                nodos[mirna]["motifs"][motivo][status] = [nodo_comun.name,len(e)]
                                for i in e:
                                    if i not in list:
                                        list.append(i)
                                        sum += 1
                                        niveles = [n.strip() for n in rutas[i]["linaje"].split(";") if n.strip()]
                                        for nivel in niveles:
                                            if nivel not in conteo["mirnas_por_especies"][mirna].keys():
                                                conteo["mirnas_por_especies"][mirna][nivel] = 0
                                            conteo["mirnas_por_especies"][mirna][nivel] += 1
                        elif tipos == "loop_motifs":
                            for status, e in inf_3.items(): 
                                nodo_comun = tree.common_ancestor(e)
                                nodos[mirna]["loop_motifs"][motivo][status] = [nodo_comun.name,len(e)]
                                for i in e:
                                    if i not in list:
                                        list.append(i)
                                        sum += 1
                                        niveles = [n.strip() for n in rutas[i]["linaje"].split(";") if n.strip()]
                                        for nivel in niveles:
                                            if nivel not in conteo["mirnas_por_especies"][mirna].keys():
                                                conteo["mirnas_por_especies"][mirna][nivel] = 0
                                            conteo["mirnas_por_especies"][mirna][nivel] += 1
            if mirna not in conteo["mirnas"].keys():
                conteo["mirnas"][mirna] = sum
                conteo_esp[mirna] = list

        for e, info in rutas.items():
            for mirna in info["sequence"].keys():
                if mirna.count("_") >= 1:
                    mirna_familia = mirna.split("_",1)[0]
                elif mirna.endswith("a") or mirna.endswith("b") or mirna.endswith("c"):
                    mirna_familia = mirna[:-1]
                else:
                    mirna_familia = mirna
                if  mirna_familia.count("-") > 1:
                    mirna_familia = mirna.split("-")[0] + "-" + mirna.split("-")[1]
                if mirna_familia.count("Bantam") == 1:
                    mirna_familia = "Bantam"

                if mirna_familia not in conteo["mirnas_alin"].keys():
                    conteo["mirnas_alin"][mirna_familia] = 0
                    conteo["mirnas_alin_especies"][mirna_familia] = {}
                conteo["mirnas_alin"][mirna_familia] += 1
                niveles = [n.strip() for n in rutas[e]["linaje"].split(";") if n.strip()]
                for nivel in niveles:
                    if nivel not in conteo["mirnas_alin_especies"][mirna_familia].keys():
                        conteo["mirnas_alin_especies"][mirna_familia][nivel] = 0
                    conteo["mirnas_alin_especies"][mirna_familia][nivel] += 1

    with open(nombre_json, "w") as f:
        json.dump(rutas, f, indent=4)

    with open(nodos_json, "w") as f:
        json.dump(nodos, f, indent=4)
    
    with open(conteo_json, "w") as f:
        json.dump(conteo, f, indent=4)

    with open(conteo_especies_json, "w") as f:
        json.dump(conteo_esp, f, indent=4)

    with open(paralogos_json, "w") as f:
        json.dump(paralogos, f, indent=4)


interpretacion(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])