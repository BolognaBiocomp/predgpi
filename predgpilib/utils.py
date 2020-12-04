
def write_gff_output(acc, sequence, output_file, gpi, cleavage, prob):
    if gpi:
        print(acc, "PredGPI", "GPI-anchor", cleavage, cleavage, prob, ".", ".",
        "Ontology_term=GO:0046658;evidence=ECO:0000256",file = output_file, sep = "\t")
    else:
        print(acc, "PredGPI", "Chain", 1, len(sequence), prob, ".", ".",
        "evidence=ECO:0000256",file = output_file, sep = "\t")
