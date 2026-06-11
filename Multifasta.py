
import subprocess
import sys
import os
import json

def Multifasta(nombre_json):
    ruta = subprocess.run(["pwd"], capture_output=True, text=True, check=True)
    directorio = ruta.stdout.strip() + "/" + "alineamiento"
    subprocess.run(["mkdir", directorio], check=True)

    #Multifasta
    with open(nombre_json, "r") as f:
        rutas = json.load(f)

        for e,inf in rutas.items():
            rutas[e]["presentes"] = []
            for mirna,sec in rutas[e]["sequence"].items():
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
                
                with open(f"{directorio}/{mirna_familia}.fasta","a") as ofile:
                    ofile.write(f">{e}\n{sec}\n")
                    
                if mirna_familia not in rutas[e]["presentes"]:
                    rutas[e]["presentes"].append(mirna_familia)

    with open(nombre_json, "w") as f:
        json.dump(rutas, f, indent=4)

Multifasta(sys.argv[1])