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
                for pos,esp in inf.items():
                    for tipo,inf_2 in esp.items():
                        if inf_2[1] > 1:
                            if pos not in nodos_resultado[mir]:
                                nodos_resultado[mir][pos] = {}
                            nodos_resultado[mir][pos][tipo] = [inf_2[0], inf_2[1]/conteo["mirnas_alin_especies"][mir][inf_2[0]]]

    with open(clean_json, "w") as f:
        json.dump(nodos_resultado, f, indent=4)

resultados(sys.argv[1], sys.argv[2], sys.argv[3])
