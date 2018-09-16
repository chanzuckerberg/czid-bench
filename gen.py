#!/usr/bin/env python3
import os

class Genome(object):

    all = dict()

    def __init__(self, category, organism, tax_id, genome_assembly_URL, genbank_URL=None):
        self.category = category
        self.organism = organism
        self.tax_id = tax_id
        self.genome_assembly_URL = genome_assembly_URL
        self.genbank_URL = genbank_URL
        self.key = f"{category}__{organism}__{tax_id}"
        self.filename = f"{self.key}.fasta"
        Genome.all[self.key] = self

    @staticmethod
    def ensure_all_present():
        for g in Genome.all.values():
            assert os.path.isfile(g.filename), f"Missing genome {g.filename}, please download from {g.genome_assembly_URL or f.genbank_URL}"


TOP_6_ID_GENOMES = [
    Genome("bacteria",
           "klebsiella_pneumoniae",
           1125630,
           "https://www.ncbi.nlm.nih.gov/genome/33?genome_assembly_id=22642"),
    Genome("fungi",
           "aspergillus_fumigatus",
           330879,
           "https://www.ncbi.nlm.nih.gov/genome/18?genome_assembly_id=22576"),
    Genome("viruses",
           "rhinovirus_c",
           463676,
           None,
           "https://www.ncbi.nlm.nih.gov/nuccore/MG148341"),
    Genome("viruses",
           "chikungunya",
           37124,
           None,
           "https://www.ncbi.nlm.nih.gov/nuccore/MG049915"),
    Genome("bacteria",
           "staphylococcus_aureus",
           93061,
           "https://www.ncbi.nlm.nih.gov/genome/154?genome_assembly_id=299272"),
    Genome("protista",
           "plasmodium_falciuparum",
           36329,
           "https://www.ncbi.nlm.nih.gov/genome/33?genome_assembly_id=22642")
]
MODELS = ["novaseq", "miseq", "hiseq"]
ABUNDANCES = ["uniform", "log-normal"]


def main():
    print("Generating IDSEQ benchmark data.")
    num_reads = 100 * 1000
    model = MODELS[0]
    abundance = ABUNDANCES[0]
    genome_fastas = " ".join(g.filename for g in TOP_6_ID_GENOMES)
    command = f"iss generate --n_reads {num_reads} --genomes {genome_fastas} --model {model} --gc_bias --output {model}_reads"
    print(command)
    Genome.ensure_all_present()
    os.system(command)


if __name__ == "__main__":
    main()
