import subprocess
import sys
import os
import json
from Bio import Entrez

def descargar(esp_txt, nombre_json, correo, archivos_dir):
    #Inicializar variables
    especies = []
    with open(esp_txt,"r") as ifile:
        for l in ifile:
            if l.strip().startswith("GCF_") or l.strip().startswith("GCA_"): 
                especies.append(l.strip())

    #Obtener linaje
    Entrez.email = str(correo)
    post_handle = Entrez.esearch(db="assembly", term=" OR ".join(especies), usehistory="y")
    post_results = Entrez.read(post_handle)
    post_handle.close()

    webenv = post_results["WebEnv"]
    query_key = post_results["QueryKey"]
    
    fetch_handle = Entrez.esummary(db="assembly", webenv=webenv, query_key=query_key)
    records = Entrez.read(fetch_handle)
    fetch_handle.close()
    resumenes = records["DocumentSummarySet"]["DocumentSummary"]

    taxid = []
    taxid_accesion = {}
    dic_nombre = {}
    for record in resumenes:
        taxid_accesion[record['Taxid']] = record['AssemblyAccession']
        taxid.append(record['Taxid'])
        dic_nombre[record['AssemblyAccession']] = record['SpeciesName']

    tax_post_handle = Entrez.epost(db="taxonomy", id=",".join(taxid))
    tax_post_results = Entrez.read(tax_post_handle)
    tax_webenv = tax_post_results["WebEnv"]
    tax_query_key = tax_post_results["QueryKey"]

    tax_fetch_handle = Entrez.efetch(db="taxonomy", webenv=tax_webenv, query_key=tax_query_key)
    tax_records = Entrez.read(tax_fetch_handle)

    lineage = {}
    for record in tax_records:
        lineage[taxid_accesion[record['TaxId']]] = record['Lineage']

    rutas = { }

    #Análisis mirnas.report
    for e in especies:
        if e in os.listdir(f"/shared/bak/TFG/miriam/{archivos_dir}/expression"):
            if "mirna.json" in os.listdir(f"/shared/bak/TFG/miriam/{archivos_dir}/expression/{e}"):
                with open(f"/shared/bak/TFG/miriam/{archivos_dir}/expression/{e}/mirna.json","r") as ifile:
                    mirnas_report = json.load(ifile)    
                    rutas[e] = {
                        "linaje": lineage[e],
                        "especie": dic_nombre[e],
                        "motifs_gen": {},
                        "motifs": {},
                        "export_3p": {},
                        "export_5p": {},
                        "loop_motifs": {},
                        "sequence": {}
                        }
                    for mirna,inf in mirnas_report.items():
                        rutas[e]["sequence"][mirna.split("-",1)[1]] = mirnas_report[mirna]["hairpin_sequence"]
                        if "motifs" in mirnas_report[mirna].keys():
                            rutas[e]["motifs_gen"][mirna.split("-",1)[1]] = mirnas_report[mirna]["motifs"]
                            rutas[e]["motifs"][mirna.split("-",1)[1]] = {}
                            rutas[e]["loop_motifs"][mirna.split("-",1)[1]] = {}
                            if "processing" in mirnas_report[mirna]["motifs"]:
                                for motivo in mirnas_report[mirna]["motifs"]["processing"].keys():
                                    rutas[e]["motifs"][mirna.split("-",1)[1]][motivo] = mirnas_report[mirna]["motifs"]["processing"][motivo]["status"]
                            if "export_3p" in mirnas_report[mirna]["motifs"]:
                                rutas[e]["export_3p"][mirna.split("-",1)[1]] = [mirnas_report[mirna]["motifs"]["export_3p"][0]["motif"], mirnas_report[mirna]["motifs"]["export_3p"][0]["sequence"]]
                            if "export_5p" in mirnas_report[mirna]["motifs"]:
                                rutas[e]["export_5p"][mirna.split("-",1)[1]] = [mirnas_report[mirna]["motifs"]["export_5p"][0]["motif"], mirnas_report[mirna]["motifs"]["export_5p"][0]["sequence"]]
                            if "loop_motifs" in mirnas_report[mirna]["motifs"]:
                                for motivo in mirnas_report[mirna]["motifs"]["loop_motifs"].keys():
                                    rutas[e]["loop_motifs"][mirna.split("-",1)[1]][motivo] = mirnas_report[mirna]["motifs"]["loop_motifs"][motivo]["status"]
                        else:
                            continue

    with open(nombre_json, "w") as f:
        json.dump(rutas, f, indent=4)

descargar(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])