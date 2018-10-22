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
import sys
import json
from multiprocessing import cpu_count
from collections import defaultdict
from util import remove_safely, check_call, smart_open, ProgressTracker
from params import MODELS, UNIFORM_ABUNDANCE, TOP_6_ID_GENOMES, NUM_READS
from genome import Genome


# Increment this as often as you like;  especially if a code change will result
# in different content for the same output filename.
LOGICAL_VERSION = "8"


IRREPRODUCIBLE = any("irreproducible" in arg for arg in sys.argv)
if IRREPRODUCIBLE:
    LOGICAL_VERSION += "irreproducible"


class ISSRunContext:
    "Execution context for each InsilicoSeq run.  Encapsulates all temp/intermediate/output files, but not the logic."

    git_commit_hash = None
    iss_version = None

    def __init__(self, tmp_prefix, output_prefix):
        print(f"GENERATING {output_prefix}__[R1, R2]")
        if ISSRunContext.iss_version == None:
            ISSRunContext.iss_version = get_iss_version()
        if not IRREPRODUCIBLE:
            self.is_reproducible = True
            self.idseq_bench_command = " ".join(sys.argv)
            if ISSRunContext.git_commit_hash == None:
                # TODO:  We should actually require a git branch in github and confirm
                # remote matches local, so that whenever parameters are changed there
                # is a permanent record in github.
                ISSRunContext.git_commit_hash = get_git_hash()
        else:
            print("***************************************************************************************************")
            print("*  WARNING:  OUTPUT WILL NOT BE REPRODUCIBLE AND WILL BE REJECTED BY IDSEQ PORTAL PRODUCTION ENV  *")
            print("***************************************************************************************************")
            self.is_reproducible = False
            self.idseq_bench_command = "irreproducible"
            ISSRunContext.git_commit_hash = "irreproducible"
        subdir = "v" + LOGICAL_VERSION + "/" + output_prefix[:-1]
        def paired_files(prefix, suffix):
            return [f"{subdir}/{prefix}_{r}.{suffix}" for r in ["R1", "R2"]]
        self.subdir = subdir
        self.output_files = paired_files(output_prefix, "fastq.gz")
        self.tmp_files = paired_files(tmp_prefix, "fastq")
        self.summary_file_txt = subdir + "/summary.txt"
        self.metadata_file_json = subdir + "/metadata.json"
        self.abundance_file = f"{subdir}/{tmp_prefix}_abundance.txt"
        self.genomes_file = f"{subdir}/{tmp_prefix}_genomes.fasta"
        self.clean_slate()

    def clean_slate(self):
        os.makedirs(self.subdir, exist_ok=True)
        for f in self.tmp_files + self.output_files + [self.summary_file_txt, self.metadata_file_json, self.abundance_file, self.genomes_file]:
            remove_safely(f)

    def remove_intermediate_and_temp_files(self):
        for f in self.tmp_files + [self.abundance_file, self.genomes_file]:
            remove_safely(f)

    def cleanup(self):
        self.remove_intermediate_and_temp_files()


def benchmark_lineage_tag(g):
    # TODO:  Make double-blind so the tools can't cheat by inspecting these tags? :)
    return f"benchmark_lineage_{g.subspecies_taxid}_{g.species_taxid}_{g.genus_taxid}_{g.family_taxid}"


def augment_and_count_read_header(line, r, line_number):
    assert len(r) == 3
    sep = r[0]
    assert line.endswith(r), f"fastq produced by ISS have read id's ending with {sep}1\\n or {sep}2\\n"
    iss_read_id = line[:-3]
    zero_padded_read_count = "{:010d}".format(line_number // 4)
    serial_number = f"s{zero_padded_read_count}"
    versioned_accession_id = iss_read_id.rsplit(sep, 1)[0].rsplit("_", 1)[0][1:]
    g = Genome.by_accid[versioned_accession_id]
    benchmark_lineage = benchmark_lineage_tag(g)
    return f"{iss_read_id}__{benchmark_lineage}__{serial_number}\n", g.key


def annotate_and_count_reads(input_fastq, output_fastq, r, counters, accumulators):
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
                augmented_read_header, g_key = augment_and_count_read_header(line, r, line_number)
                counters[g_key] += 1
                output_f.write(augmented_read_header.encode('utf-8'))
                for i in range(4):
                    line = input_f.readline()
                    line_number += 1
                    if i == 0:
                        accumulators[g_key] += (len(line) - 1)
                    if i < 3:
                        output_f.write(line.encode('utf-8'))
        except Exception as _:
            print(f"Error parsing line {line_number} in {input_fastq}.")
            raise


def output_summary_counters(rc, iss_command, counters, accumulators):
    contents = []
    total_reads = 0
    with smart_open(rc.summary_file_txt, "w") as mft:
        headers = "READ_COUNT\tREAD_SIZE\tCOVERAGE\tLINEAGE\tGENOME\n"
        mft.write(headers)
        print("")
        print(headers)
        for g_key, read_count in sorted(counters.items(), key=lambda pair: pair[1], reverse=True):
            g = Genome.all[g_key]
            benchmark_lineage = benchmark_lineage_tag(g)
            coverage = accumulators[g.key] / g.size
            read_size = accumulators[g.key] / read_count
            summary_line = f"{read_count}\t{read_size}\t{coverage:3.1f}x\t{benchmark_lineage}\t{g.key}\n"
            print(summary_line)
            mft.write(summary_line)
            contents.append({
                'genome': g.key,
                'benchmark_lineage': benchmark_lineage,
                'coverage': coverage,
                'read_count': read_count,
                'read_size' : read_size,
                'genome_size': g.size
            })
            total_reads += read_count
    metadata = {
        "iss_version": rc.iss_version,
        "iss_command": iss_command,
        "idseq_bench_command": rc.idseq_bench_command,
        "idseq_bench_git_hash": ISSRunContext.git_commit_hash,
        "idseq_bench_reproducible": rc.is_reproducible,
        "idseq_bench_version": LOGICAL_VERSION,
        "prefix": os.path.basename(rc.subdir),
        "fastqs": [os.path.basename(f) for f in rc.output_files],
        "verified_total_reads": total_reads,
        "verified_contents": contents
    }
    with smart_open(rc.metadata_file_json, "w") as mdf:
        mdf.write(json.dumps(metadata, indent=4))


def get_iss_version():
    remove_safely("iss_version.txt")
    check_call("iss --version > iss_version.txt")
    with open("iss_version.txt") as ivf:
        iss_version = ivf.readline().strip()
    assert iss_version.startswith("iss version ")
    iss_version = iss_version[len("iss version "):]
    remove_safely("iss_version.txt")
    return iss_version


def get_git_hash():
    no_changes = os.system("git status --porcelain | grep .")
    assert no_changes, "You have uncommitted changes.  Please commit and push to a github branch first.  Or use --irreproducible flag."
    remove_safely("git_status.txt")
    try:
        check_call("git log | head -1 > git_status.txt")
        with open("git_status.txt", "r") as gsf:
            git_commit = gsf.readline().strip()
    finally:
        remove_safely("git_status.txt")
    assert git_commit.startswith("commit ") and len(git_commit) == 47, "could not find git commit hash;  use --irreproducible flag."
    _, commit_hash = git_commit.split()
    return commit_hash


def strictly_above(va, vb):
    for a, b in zip(va.split("."), vb.split(".")):
        if int(a) > int(b):
            return True
        if int(a) < int(b):
            return False
    return False


def run_iss(rc, iss_command):
    check_call(iss_command)
    counters = defaultdict(int)
    accumulators = defaultdict(int)
    if strictly_above(rc.iss_version, "1.2.0"):
        suffixes = ["/1\n", "/2\n"]
    else:
        suffixes = ["_1\n", "_2\n"]
    for tmp_fastq, output_fastq, r in zip(rc.tmp_files, rc.output_files, suffixes):
        annotate_and_count_reads(tmp_fastq, output_fastq, r, counters, accumulators)
    output_summary_counters(rc, iss_command, counters, accumulators)
    rc.cleanup()


def run_iss_single_genome(g, num_reads, model, tmp_prefix, num_cpus):
    num_organisms = 1
    num_accessions = len(g.versioned_accession_ids)
    abundance = UNIFORM_ABUNDANCE
    output_prefix_single_genome = f"norg_{num_organisms}__nacc_{num_accessions}__{abundance}_weight_per_accession__{model}_reads__{g.key}__v{LOGICAL_VERSION}_"
    rc = ISSRunContext(tmp_prefix, output_prefix_single_genome)
    iss_command_single_genome = f"iss generate --n_reads {num_reads} --genomes {g.filename} --model {model} --abundance {abundance} --gc_bias --output {rc.subdir}/{tmp_prefix} --cpus {num_cpus}"
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
    iss_multiplexed_command = f"iss generate --n_reads {num_reads} --genomes {rc.genomes_file} --model {model} --abundance_file {rc.abundance_file} --gc_bias --output {rc.subdir}/{tmp_prefix} --cpus {num_cpus}"
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
