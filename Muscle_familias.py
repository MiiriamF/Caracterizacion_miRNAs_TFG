import subprocess
import sys
import os
import json
from Bio import AlignIO
from Bio import motifs

def MUSCLE_align(nombre_json):
    #MUSCLE
    with open(nombre_json, "r") as f:
        rutas = json.load(f)
        
        ruta = subprocess.run(["pwd"], capture_output=True, text=True, check=True)
        directorio = ruta.stdout.strip() + "/" + "alineamiento"

        for i in os.listdir(directorio):
            if i.endswith(".fasta") and f"{directorio}/{i.split('.')[0]}_muscle.afa" not in os.listdir(directorio):
                subprocess.run(["muscle", "-in", f"{directorio}/{i}", "-out", f"{directorio}/{i.split('.')[0]}_muscle.afa"], cwd=directorio, check=True)

        for i in os.listdir(directorio):
            if i.endswith("_muscle.afa") and f"{directorio}/{i.split('_')[0]}_alineamiento.txt" not in os.listdir(directorio):
                alignment = AlignIO.read(f"{directorio}/{i}", "fasta")
                instancias = []
                for record in alignment:
                    record.id = record.id.split("|")[0]
                    with open(f"{directorio}/{i.split('_')[0]}_alineamiento.txt","a") as ofile:
                        ofile.write(record.id + "\t" + str(record.seq) + "\n")
                    instancias.append(record.seq)
                m = motifs.create(instancias, alphabet="ACGT-")
                consenso = m.consensus
                try:
                    with open(f"{ruta.stdout.strip()}/consensos.txt","r") as ofile_r:
                        presente = False
                        for bus in ofile_r:
                            if i in bus:
                                presente = True
                except FileNotFoundError:
                    presente = False
                    pass
                if not presente:
                    with open(f"{ruta.stdout.strip()}/consensos.txt","a") as ofile:
                            ofile.write(i.split('_')[0] + "\n" + str(consenso) + "\n")

    with open(nombre_json, "w") as f:
        json.dump(rutas, f, indent=4)

MUSCLE_align(sys.argv[1])