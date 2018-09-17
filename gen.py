#!/usr/bin/env python3
import os
import traceback
import gzip
from multiprocessing import cpu_count

# Pre-requsites:
#
#    pip3 install InSilicoSeq
#    pip3 install ncbi-acc-download
#

LOGICAL_VERSION = "1"

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


def annotate_reads_work(input_f, output_f, r):
    # Annotate read ID by appending the consecutive read counter, and stripping the _1 or _2
    # paired-end identifier appended by ISS;  this is required to run correctly through STAR.
    assert len(r) == 3
    assert r in ["_1\n", "_2\n"]
    line_number = 1
    try:
        line = input_f.readline()
        while line:
            # The FASTQ format specifies that each read consists of 4 lines,
            # the first of which begins with @ followed by read ID.
            assert line[0] == "@", f"fastq format requires every 4th line to start with @"
            assert line.endswith(r), f"fastq produced by ISS have read IDS ending with _1\\n or _2\\n"
            zero_padded_read_count = "{:010d}".format(line_number//4)
            augmented_header = line[:-3] + f"__{zero_padded_read_count}\n"
            output_f.write(augmented_header.encode('utf-8'))
            for i in range(4):
                line = input_f.readline()
                line_number += 1
                if i < 3:
                    output_f.write(line.encode('utf-8'))
    except Exception as e:
        e.line_number = line_number
        raise


def opener(filename):
    if filename.endswith(".gz"):
        return gzip.open
    else:
        return open


def smart_open(filename, mode):
    return opener(filename)(filename, mode)


def annotate_reads(input_fastq, output_fastq, r):
    try:
        with smart_open(input_fastq, "r") as input_f:
            with smart_open(output_fastq, "w") as output_f:
                annotate_reads_work(input_f, output_f, r)
    except Exception as e:
        print(f"Error parsing line {e.line_number} in {input_fastq}.")
        raise


def main():
    print("Generating IDSEQ benchmark data.")
    num_cpus = cpu_count()
    num_reads = 100 * 1000
    model = MODELS[0]
    abundance = ABUNDANCES[0]
    genome_fastas = " ".join(g.filename for g in TOP_6_ID_GENOMES)
    Genome.ensure_all_present()
    num_organisms = len(TOP_6_ID_GENOMES)
    num_accessions = sum(len(g.versioned_accession_ids) for g in TOP_6_ID_GENOMES)
    output_prefix = f"norg_{num_organisms}__nacc_{num_accessions}__{abundance}_weight_per_accession__{model}_reads__v{LOGICAL_VERSION}_"
    pid = os.getpid()
    tmp_prefix = f"tmp_{pid}"
    tmp_files = [f"{tmp_prefix}_{r}.fastq" for r in ["R1", "R2"]]
    output_files = [f"{output_prefix}_{r}.fastq.gz" for r in ["R1", "R2"]]
    for f in tmp_files + output_files:
        remove_safely(f)
    # TODO:  Currently each chromosome is treated as a separate organism
    # for relative abundance purposes.  Thus, organisms with greater number
    # of chromosomes will have a lot of extra weight in the mix.
    # Probably should fix here?
    remove_safely("top_6_pathogens.fasta")
    command = f"cat {genome_fastas} > top_6_pathogens.fasta"
    check_call(command)
    command = f"iss generate --n_reads {num_reads} --genomes top_6_pathogens.fasta --model {model} --abundance {abundance} --gc_bias --output {tmp_prefix} --cpus {num_cpus}"
    check_call(command)
    for tmp_fastq, output_fastq, r in zip(tmp_files, output_files, ["_1\n", "_2\n"]):
        annotate_reads(tmp_fastq, output_fastq, r)
        remove_safely(tmp_fastq)


if __name__ == "__main__":
    main()
