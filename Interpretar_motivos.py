import sys
import json

def resultados(nodos_json, conteo_json, clean_json):
    with open(nodos_json, "r") as f:
        nodos = json.load(f)
        nodos_resultado = {}
        with open(conteo_json, "r") as f:
            conteo = json.load(f)
            for mir,inf in nodos.items():
                nodos_resultado[mir] = {}
                for tipo,inf_2 in inf.items():
                    nodos_resultado[mir][tipo] = {} 
                    if tipo == "motifs" or tipo == "loop_motifs":
                        for motivo,inf_3 in inf_2.items():
                            nodos_resultado[mir][tipo][motivo] = {}
                            for status,inf_4 in inf_3.items():
                                if inf_4[0].startswith("GCF_") or inf_4[0].startswith("GCA_"):
                                    nodos_resultado[mir][tipo][motivo][status] = [inf_4[0]]
                                else:
                                    nodos_resultado[mir][tipo][motivo][status] = [inf_4[0], inf_4[1]/conteo["mirnas_alin_especies"][mir][inf_4[0]]]
                                    
                    elif tipo == "export_3p" or tipo == "export_5p":
                        if inf_2[0].startswith("GCF_") or inf_2[0].startswith("GCA_"):
                            nodos_resultado[mir][tipo] = [inf_2[0]]
                        else:
                            nodos_resultado[mir][tipo] = [inf_2[0], inf_2[1]/conteo["mirnas_alin_especies"][mir][inf_2[0]]]

    with open(clean_json, "w") as f:
        json.dump(nodos_resultado, f, indent=4)

resultados(sys.argv[1], sys.argv[2], sys.argv[3])
