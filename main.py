#!/usr/bin/env python3
#
# IDSeq Benchmark Generator.
#
# Prerequsites:
#
#    pip3 install InSilicoSeq
#    pip3 install ncbi-acc-download
#
# Usage:
#
#    Customize params.py, then run main.py, then use idseq-cli to upload output to idseq.
#
import os
import time
from multiprocessing import cpu_count
from util import remove_safely, check_call, smart_open
from params import MODELS, ABUNDANCES, TOP_6_ID_GENOMES
from genome import Genome


# Increment this as often as you like;  especially if a code change will result
# in different content for the same output filename.
LOGICAL_VERSION = "4"


def augment_read_header(line, r, line_number):
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
    # Also annotate with lineage information, that can later be used to score idseq accuracy.
    line_number = 1
    try:
        line = input_f.readline()
        while line:
            # The FASTQ format specifies that each read consists of 4 lines,
            # the first of which begins with @ followed by read ID.
            assert line[0] == "@", f"fastq format requires every 4th line to start with @"
            augmented_read_header = augment_read_header(line, r, line_number)
            output_f.write(augmented_read_header.encode('utf-8'))
            for i in range(4):
                line = input_f.readline()
                line_number += 1
                if i < 3:
                    output_f.write(line.encode('utf-8'))
    except Exception as e:
        e.line_number = line_number
        raise


def annotate_reads(input_fastq, output_fastq, r):
    with smart_open(input_fastq, "r") as input_f:
        with smart_open(output_fastq, "w") as output_f:
            try:
                annotate_reads_work(input_f, output_f, r)
            except Exception as e:
                print(f"Error parsing line {e.line_number} in {input_fastq}.")
                raise


def paired_files(prefix):
    return [f"{prefix}_{r}.fastq.gz" for r in ["R1", "R2"]]


def annotate_and_move(tmp_fastq, output_fastq):
    for tmp_fastq, output_fastq, r in zip(tmp_files, output_files, ["_1\n", "_2\n"]):
        annotate_reads(tmp_fastq, output_fastq, r)
        remove_safely(tmp_fastq)


def run_iss_single_genome(g, num_reads, model, tmp_prefix, num_cpus):
    abundance = ABUNDANCES[0]
    num_organisms = 1
    num_accessions = len(g.versioned_accession_ids)
    tmp_files = paired_files(tmp_prefix)
    output_prefix_single_genome = f"norg_{num_organisms}__nacc_{num_accessions}__{abundance}_weight_per_accession__{model}_reads__{g.key}__v{LOGICAL_VERSION}_"
    print("GENERATING {output_prefix_single_genome} WITH {num_reads} READS.")
    output_fastqs = paired_files(output_prefix_single_genome)
    for f in tmp_files + output_files:
        remove_safely(f)
    iss_command_single_genome = f"iss generate --n_reads {num_reads} --genomes {g.filename} --model {model} --abundance {abundance} --gc_bias --output {tmp_prefix} --cpus {num_cpus}"
    check_call(iss_command_single_genome)
    annotate_and_move(tmp_fastq, output_fastq)


def uniform_abundance_per_organism(genomes, tmp_prefix):
    # Ensure each organism has equal abundance in the mix, despite the varying
    # numbers of chromosomes and accessions across organisms.
    abundance_file = f"{tmp_prefix}_abundance.txt"
    remove_safely(abundance_file)
    with open(abundance_file, "w") as af:
        sum_of_weights = 0
        for g in genomes:
            vaccid_weight = 1.0 / (len(g.versioned_accession_ids) * len(genomes))
            sum_of_weights += vaccid_weight
            af.write(f"{vaccid} {vaccid_weight}\n")
    assert -0.0005 < sum_of_weights - 1.0 < 0.0005, f"{sum_of_weights} != 1.0"
    return abundance_file


def concatenate_fasta(genomes, tmp_prefix):
    genomes_file = f"{tmp_prefix}_genomes.fasta"
    genome_fastas = " ".join(g.filename for g in genomes)
    remove_safely(genomes_file)
    command = f"cat {genome_fastas} > {genomes_file}"
    check_call(command)
    return genomes_file


def run_iss_multiplexed(genomes, num_reads, model, tmp_prefix, num_cpus):
    num_organisms = len(genomes)
    num_accessions = sum(len(g.versioned_accession_ids) for g in genomes)
    output_prefix_multiplexed = f"norg_{num_organisms}__nacc_{num_accessions}__uniform_weight_per_organism__{model}_reads__v{LOGICAL_VERSION}_"
    print("GENERATING {output_prefix_multiplexed} WITH {num_reads} READS.")
    output_files = paired_files(output_prefix_multiplexed)
    tmp_files = paired_files(tmp_prefix)
    for f in tmp_files + output_files:
        remove_safely(f)
    abundance_file = uniform_abundance_per_organism(genomes, tmp_prefix)
    genomes_file = concatenate_fasta(genomes)
    iss_multiplexed_command = f"iss generate --n_reads {num_reads} --genomes {genomes_file} --model {model} --abundance_file {abundance_file} --gc_bias --output {tmp_prefix} --cpus {num_cpus}"
    check_call(iss_multiplexed_command)
    annotate_and_move(tmp_fastq, output_fastq)
    remove_safely(abundance_file)
    remove_safely(genomes_file)


def main():
    print("Generating IDSEQ benchmark data.")
    num_cpus = cpu_count()
    num_reads = 100 * 1000
    Genome.ensure_all_present()
    pid = os.getpid()
    tmp_prefix = f"tmp_{pid}"
    total_reads = num_reads * (1 + len(TOP_6_ID_GENOMES))
    reads_generated_so_far = 0
    t_start = time.time()
    def tick():
        t_elapsed = time.time() - t_start
        t_remaining = (t_elapsed / reads_generated_so_far) * total_reads
        t_eta = t_start + t_remaining
        t_eta_str = time.strftime("%H:%M:%S", time.localtime(t_eta))
        print(f"{3.1f:reads_so_far/total_reads*100:3.1f} percent done, {t_remaining/60:3.1f} minutes remaining, ETA {t_eta_str}")
    for model in MODELS:
        # First, generate a separate benchmark for each genome.
        for g in TOP_6_ID_GENOMES:
            run_iss_single_genome(g, num_reads, model, tmp_prefix, num_cpus)
            reads_generated_so_far += num_reads
            tick()
        # Then generate a multiplexed benchmark.
        run_iss_multiplexed(TOP_6_ID_GENOMES, num_reads, model, tmp_prefix, num_cpus)
        reads_generated_so_far += num_reads
        tick()


if __name__ == "__main__":
    main()
