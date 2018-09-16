#!/usr/bin/env python3
import os
import traceback

# Pre-requsites:
#
#    pip3 install InSilicoSeq
#    pip3 install ncbi-acc-download
#

QUIET = False

def remove_safely(fn):
    try:
        if os.path.isfile(fn):
            os.remove(fn)
    except:
        print(f"Ignoring exception {traceback.format_exc()}")


def check_call(command):
    if not QUIET:
        print(command)
    exitcode = os.system(command)
    assert exitcode == 0, f"Command failed: {command}"


class Genome(object):

    all = dict()

    def __init__(self, category, organism, tax_id, genome_assembly_URL, versioned_accession_ids=[]):
        self.category = category
        self.organism = organism
        self.tax_id = tax_id
        self.genome_assembly_URL = genome_assembly_URL
        self.versioned_accession_ids = versioned_accession_ids
        self.key = f"{category}__{organism}__{tax_id}"
        self.filename = f"{self.key}.fasta"
        Genome.all[self.key] = self

    @staticmethod
    def fetch_versioned_accession_id(vaccid):  # e.g., "NC_004325.2"
        output_file = f"{vaccid}.fa"
        if os.path.isfile(output_file):
            print(f"{output_file} already exists, nice")
        else:
            try:
                command = f"ncbi-acc-download --format fasta {vaccid} -e all"
                check_call(command)
            except:
                remove_safely(output_file)
                raise
        return output_file

    @staticmethod
    def ensure_all_present():
        for g in Genome.all.values():
            remove_safely(g.filename)
            accession_fas = []
            for f in (g.versioned_accession_ids):
                af = Genome.fetch_versioned_accession_id(f)
                accession_fas.append(af)
            accession_fastas = " ".join(accession_fas)
            command = f"cat {accession_fastas} > {g.filename}"
            check_call(command)
            assert os.path.isfile(g.filename), f"Failed to download genome {g.filename}"


TOP_6_ID_GENOMES = [
    Genome("bacteria",
           "klebsiella_pneumoniae_HS11286",
           1125630,
           None,
           ["NC_016845.1"]),
    Genome("fungi",
           "aspergillus_fumigatus",
           330879,
           "https://www.ncbi.nlm.nih.gov/genome/18?genome_assembly_id=22576",
           ["NC_007194.1", "NC_007195.1", "NC_007196.1", "NC_007197.1", "NC_007198.1", "NC_007199.1", "NC_007200.1", "NC_007201.1"]),
    Genome("viruses",
           "rhinovirus_c",
           463676,
           None,
           ["MG148341"]),
    Genome("viruses",
           "chikungunya",
           37124,
           None,
           ["MG049915"]),
    Genome("bacteria",
           "staphylococcus_aureus",
           93061,
           "https://www.ncbi.nlm.nih.gov/genome/154?genome_assembly_id=299272",
           ["NC_007795.1"]),
    Genome("protista",
           "plasmodium_falciuparum",
           36329,
           "https://www.ncbi.nlm.nih.gov/genome/33?genome_assembly_id=22642",
           ["NC_004325.2", "NC_037280.1", "NC_000521.4", "NC_004318.2", "NC_004326.2", "NC_004327.3", "NC_004328.3", "NC_004329.3", "NC_004330.2", "NC_037281.1", "NC_037282.1", "NC_037284.1", "NC_004331.3", "NC_037283.1", "NC_036769.1"])
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
    Genome.ensure_all_present()
    check_call(command)


if __name__ == "__main__":
    main()
