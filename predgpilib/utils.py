
def write_gff_output(acc, sequence, output_file, gpi, cleavage, prob):
    if gpi:
        print(acc, "PredGPI", "GPI-anchor", cleavage, cleavage, prob, ".", ".",
        "Ontology_term=GO:0046658;evidence=ECO:0000256",file = output_file, sep = "\t")
    else:
        print(acc, "PredGPI", "Chain", 1, len(sequence), prob, ".", ".",
        "evidence=ECO:0000256",file = output_file, sep = "\t")

def get_json_output(acc, sequence, gpi, cleavage, prob):
    acc_json = {'accession': acc, 'features': []}
    acc_json['sequence'] = {
                              "length": len(sequence),
                              "sequence": sequence
                           }
    if gpi:
        acc_json['features'].append({
            "type": "LIPID",
            "category": "PTM",
            "description": "GPI-anchor",
            "begin": cleavage,
            "end": cleavage,
            "score": round(float(prob), 2),
            "evidences": [
              {
                "code": "ECO:0000256",
                "source": {
                  "name": "SAM",
                  "id": "PredGPI",
                  "url": "http://gpcr2.biocomp.unibo.it/predgpi"
                }
              }
            ]
        })
    return acc_json
