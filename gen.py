#!/usr/bin/env python3
#
# IDSeq Benchmark Generator.
#
# Prerequsites:
#
#    pip3 install InSilicoSeq
#    pip3 install ncbi-acc-download
#
#
import os
from multiprocessing import cpu_count
from util import remove_safely, check_call, smart_open
from params import MODELS, ABUNDANCES, TOP_6_ID_GENOMES
from genome import Genome


# Increment this as often as you like;  especially if a code change will result
# in different content for the same output filename.
LOGICAL_VERSION = "3"


def augmented_read_header(line, r, line_number):
    assert line.endswith(r), f"fastq produced by ISS have read id's ending with _1\\n or _2\\n"
    assert len(r) == 3
    iss_read_id = line[:-3]
    zero_padded_read_count = "{:010d}".format(line_number // 4)
    serial_number = f"s{zero_padded_read_count}"
    versioned_accession_id = iss_read_id.rsplit("_", 1)[0][1:]
    g = Genome.by_accid[versioned_accession_id]
    benchmark_lineage = f"benchmark_lineage_{g.subspecies_taxid}_{g.species_taxid}_{g.genus_taxid}_{g.family_taxid}"
    return f"{iss_read_id}__{benchmark_lineage}__{serial_number}\n"


def annotate_reads_work(input_f, output_f, r):
    # Annotate read ID by appending the consecutive read counter, and stripping the _1 or _2
    # paired-end identifier appended by ISS;  this is required to run correctly through STAR.
    line_number = 1
    try:
        line = input_f.readline()
        while line:
            # The FASTQ format specifies that each read consists of 4 lines,
            # the first of which begins with @ followed by read ID.
            assert line[0] == "@", f"fastq format requires every 4th line to start with @"
            output_f.write(augmented_read_header(line, r, line_number).encode('utf-8'))
            for i in range(4):
                line = input_f.readline()
                line_number += 1
                if i < 3:
                    output_f.write(line.encode('utf-8'))
    except:
        print(f"Error parsing line {line_number} in {input_fastq}.")
        raise


def annotate_reads(input_fastq, output_fastq, r):
    with smart_open(input_fastq, "r") as input_f:
        with smart_open(output_fastq, "w") as output_f:
            annotate_reads_work(input_f, output_f, r)


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
