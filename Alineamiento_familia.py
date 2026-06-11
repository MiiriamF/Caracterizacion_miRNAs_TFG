import subprocess
import sys
import os
import json

def alineamiento(nombre_json, resultados_json):
    ruta = subprocess.run(["pwd"], capture_output=True, text=True, check=True)
    directorio = ruta.stdout.strip() + "/" + "alineamiento"
    
    with open(nombre_json, "r") as f:
        rutas = json.load(f)
        res_json = {}
        distinta = {}
        for e,inf in rutas.items():
            with open(ruta.stdout.strip() + "/consensos.txt","r") as ifile:
                for l in ifile:
                    sec = ""
                    mirna = str(l.strip())
                    consenso = next(ifile).strip()
                    if mirna in rutas[e]["presentes"]:
                        with open (directorio + f"/{mirna}_alineamiento.txt","r") as ifile_2:
                            for l_2 in ifile_2:
                                if l_2.startswith(e): 
                                    sec = l_2.strip().split("\t")[1]
                                    if sec != "":
                                        if mirna not in distinta:
                                            distinta[mirna] = {}
                                        for i in range(0, len(sec)):
                                            if sec[i] != consenso[i]:
                                                if i not in distinta[mirna]:
                                                    if sec[i] == "-":
                                                        distinta[mirna][i] = []
                                                        distinta[mirna][i].append([e,"Delecion"])
                                                    elif consenso[i] == "-":
                                                        distinta[mirna][i] = []
                                                        distinta[mirna][i].append([e,"Insercion"])
                                                    else:
                                                        distinta[mirna][i] = []
                                                        distinta[mirna][i].append([e,"Sustitucion"+"("+consenso[i]+"->"+sec[i]+")"])
                                                else:
                                                    if sec[i] == "-":
                                                        distinta[mirna][i].append([e,"Delecion"])
                                                    elif consenso[i] == "-":
                                                        distinta[mirna][i].append([e,"Insercion"])
                                                    else:
                                                        distinta[mirna][i].append([e,"Sustitucion"+"("+consenso[i]+"->"+sec[i]+")"])
                                        if distinta[mirna] == {}:
                                            del distinta[mirna]

    #añadir número de cambios por miRNA
        for mir,inf in distinta.items():
                res_json[mir] = {}
                for pos,esp in inf.items():
                    if pos not in res_json[mir]:
                        res_json[mir][str(pos + 1)] = {}
                    for x in esp:
                        res_json[mir][str(pos + 1)][x[0]] = x[1]

    with open(nombre_json, "w") as f:
        json.dump(rutas, f, indent=4)

    with open(resultados_json, "w") as f:
        json.dump(res_json, f, indent=4)

alineamiento(sys.argv[1], sys.argv[2])