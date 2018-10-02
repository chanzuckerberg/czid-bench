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
from collections import defaultdict
from util import remove_safely, check_call, smart_open
from params import MODELS, UNIFORM_ABUNDANCE, READ_SIZE, TOP_6_ID_GENOMES, NUM_READS
from genome import Genome


# Increment this as often as you like;  especially if a code change will result
# in different content for the same output filename.
LOGICAL_VERSION = "5"

class ISSRunContext(object):
    "Execution context for each InsilicoSeq run.  Encapsulates all temp/intermediate/output files, but not the logic."

    def __init__(self, tmp_prefix, output_prefix):
        print(f"GENERATING {output_prefix}__[R1, R2]")
        def paired_files(prefix, suffix):
            return [f"{prefix}_{r}.{suffix}" for r in ["R1", "R2"]]
        self.output_files = paired_files(output_prefix, "fastq.gz")
        self.tmp_files = paired_files(tmp_prefix, "fastq")
        self.summary_file = output_prefix[:-1] + "__summary_counts_and_coverage.txt"
        self.abundance_file = f"{tmp_prefix}_abundance.txt"
        self.genomes_file = f"{tmp_prefix}_genomes.fasta"
        self.clean_slate()

    def clean_slate(self):
        for f in self.tmp_files + self.output_files + [self.summary_file, self.abundance_file, self.genomes_file]:
            remove_safely(f)

    def remove_intermediate_and_temp_files(self):
        for f in self.tmp_files + [self.abundance_file, self.genomes_file]:
            remove_safely(f)

    def cleanup(self):
        self.remove_intermediate_and_temp_files()


class ProgressTracker(object):

    def __init__(self, target):
        self.target = target
        self.current = 0
        self.t_start = time.time()

    def advance(self, amount):
        PESSIMISM = 2.0
        self.current += amount
        t_elapsed = time.time() - self.t_start
        t_remaining = (t_elapsed / self.current) * self.target - t_elapsed
        t_remaining *= PESSIMISM
        t_eta = self.t_start + t_elapsed + t_remaining
        t_eta_str = time.strftime("%H:%M:%S", time.localtime(t_eta))
        print(f"*** {self.current/self.target*100:3.1f} percent done, {t_elapsed/60:3.1f} minutes elapsed, {t_remaining/60:3.1f} minutes remaining, ETA {t_eta_str} ***\n")



def benchmark_lineage_tag(g):
    # TODO:  Make double-blind so the tools can't cheat by inspecting these tags? :)
    return f"benchmark_lineage_{g.subspecies_taxid}_{g.species_taxid}_{g.genus_taxid}_{g.family_taxid}"


def augment_and_count_read_header(line, r, line_number, counters):
    assert line.endswith(r), f"fastq produced by ISS have read id's ending with _1\\n or _2\\n"
    assert len(r) == 3
    iss_read_id = line[:-3]
    zero_padded_read_count = "{:010d}".format(line_number // 4)
    serial_number = f"s{zero_padded_read_count}"
    versioned_accession_id = iss_read_id.rsplit("_", 1)[0][1:]
    g = Genome.by_accid[versioned_accession_id]
    benchmark_lineage = benchmark_lineage_tag(g)
    counters[g.key] += 1
    return f"{iss_read_id}__{benchmark_lineage}__{serial_number}\n"


def annotate_and_count_reads(input_fastq, output_fastq, r, counters):
    """Annotate read IDs by appending the consecutive read counter, after stripping
    that _1 or _2  paired-end indicator appended by ISS.  Required to run correctly
    through STAR.  Also annotate with lineage information, that can later be used
    to score idseq accuracy."""
    with smart_open(input_fastq, "r") as input_f, \
         smart_open(output_fastq, "w") as output_f:
        line_number = 1
        try:
            line = input_f.readline()
            while line:
                # The FASTQ format specifies that each read consists of 4 lines,
                # the first of which begins with @ followed by read ID.
                assert line[0] == "@", f"fastq format requires every 4th line to start with @"
                augmented_read_header = augment_and_count_read_header(line, r, line_number, counters)
                output_f.write(augmented_read_header.encode('utf-8'))
                for i in range(4):
                    line = input_f.readline()
                    line_number += 1
                    if i < 3:
                        output_f.write(line.encode('utf-8'))
        except Exception as e:
            print(f"Error parsing line {line_number} in {input_fastq}.")
            raise


def output_summary_counters(summary_file, counters):
    with smart_open(summary_file, "w") as sf:
        print("")
        print("READ_COUNT\tCOVERAGE\tLINEAGE\tGENOME\n")
        for g_key, read_count in sorted(counters.items(), key=lambda pair: pair[1], reverse=True):
            g = Genome.all[g_key]
            benchmark_lineage = benchmark_lineage_tag(g)
            coverage = read_count * READ_SIZE / g.size
            summary_line = f"{read_count}\t{coverage:3.1f}x\t{benchmark_lineage}\t{g.key}\n"
            print(summary_line)
            sf.write(summary_line)


def run_iss(rc, iss_command):
    check_call(iss_command)
    counters = defaultdict(int)
    for tmp_fastq, output_fastq, r in zip(rc.tmp_files, rc.output_files, ["_1\n", "_2\n"]):
        annotate_and_count_reads(tmp_fastq, output_fastq, r, counters)
    output_summary_counters(rc.summary_file, counters)
    rc.cleanup()


def run_iss_single_genome(g, num_reads, model, tmp_prefix, num_cpus):
    num_organisms = 1
    num_accessions = len(g.versioned_accession_ids)
    abundance = UNIFORM_ABUNDANCE
    output_prefix_single_genome = f"norg_{num_organisms}__nacc_{num_accessions}__{abundance}_weight_per_accession__{model}_reads__{g.key}__v{LOGICAL_VERSION}_"
    iss_command_single_genome = f"iss generate --n_reads {num_reads} --genomes {g.filename} --model {model} --abundance {abundance} --gc_bias --output {tmp_prefix} --cpus {num_cpus}"
    rc = ISSRunContext(tmp_prefix, output_prefix_single_genome)
    run_iss(rc, iss_command_single_genome)


def uniform_abundance_per_organism(genomes, abundance_file):
    # Ensure each organism has equal abundance in the mix, despite the varying
    # numbers of chromosomes and accessions across organisms.
    with open(abundance_file, "w") as af:
        sum_of_weights = 0
        for g in genomes:
            vaccid_weight = 1.0 / (len(g.versioned_accession_ids) * len(genomes))
            for vaccid in g.versioned_accession_ids:
                sum_of_weights += vaccid_weight
                af.write(f"{vaccid} {vaccid_weight}\n")
    assert -0.0005 < sum_of_weights - 1.0 < 0.0005, f"{sum_of_weights} != 1.0"


def concatenate_fasta(genomes, genomes_file):
    genome_fastas = " ".join(g.filename for g in genomes)
    command = f"cat {genome_fastas} > {genomes_file}"
    check_call(command)


def run_iss_multiplexed(genomes, num_reads, model, tmp_prefix, num_cpus):
    num_organisms = len(genomes)
    num_accessions = sum(len(g.versioned_accession_ids) for g in genomes)
    output_prefix_multiplexed = f"norg_{num_organisms}__nacc_{num_accessions}__uniform_weight_per_organism__{model}_reads__v{LOGICAL_VERSION}_"
    rc = ISSRunContext(tmp_prefix, output_prefix_multiplexed)
    uniform_abundance_per_organism(genomes, rc.abundance_file)
    concatenate_fasta(genomes, rc.genomes_file)
    iss_multiplexed_command = f"iss generate --n_reads {num_reads} --genomes {rc.genomes_file} --model {model} --abundance_file {rc.abundance_file} --gc_bias --output {tmp_prefix} --cpus {num_cpus}"
    run_iss(rc, iss_multiplexed_command)


def main():
    print("Generating IDSEQ benchmark data.")
    num_cpus = cpu_count()
    num_reads = NUM_READS
    Genome.fetch_all()
    pid = os.getpid()
    tmp_prefix = f"tmp_{pid}"
    ticker = ProgressTracker(target=num_reads * len(MODELS) * (1 + len(TOP_6_ID_GENOMES)))
    for model in MODELS:
        # First, generate a separate benchmark for each genome.
        for g in TOP_6_ID_GENOMES:
            run_iss_single_genome(g, num_reads, model, tmp_prefix, num_cpus)
            ticker.advance(num_reads)
        # Then generate a multiplexed benchmark.
        run_iss_multiplexed(TOP_6_ID_GENOMES, num_reads, model, tmp_prefix, num_cpus)
        ticker.advance(num_reads)


if __name__ == "__main__":
    main()
