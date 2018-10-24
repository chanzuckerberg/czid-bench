#!/usr/bin/env python3
#
# IDSeq Benchmark Scorer.
#
# After running an idseq-bench sample through the IDSeq Portal,
# score the portal output as follows:
#
#     python3 score.py s3://idseq-samples-prod/samples/16/8848/results/2.8
#
# This computes the QC and recall rate for each sample organism.
#
import sys
import json
import re
from collections import defaultdict
from util import smarter_open, smarter_readline, smart_glob


# TODO: Use NamedTuple instead.
# List ranks in the same order as benchmark_lineage tags.
TAXID_RANKS = ["subspecies", "species", "genus", "family"]

# Benchmark linage tag format.
BENCHMARK_LINEAGE_PATTERN = r'__benchmark_lineage_\d+_\d+_\d+_\d+__'

# IDSeq linage tag format in taxid_annot.fasta.
TAXID = r'[+,-]?\d+'
IDSEQ_LINEAGE_PATTERN = (
    f"family_nr:{TAXID}:family_nt:{TAXID}:" +
    f"genus_nr:{TAXID}:genus_nt:{TAXID}:" +
    f"species_nr:{TAXID}:species_nt:{TAXID}"
)

def glob_sample_data(sample, version):
    "Enumerate pipeline stage data to be counted."
    return {
        "input_fastq":
            smart_glob(f"{sample}/fastqs/*.fastq*", expected=2),
        "post_qc_fasta":
            smart_glob(f"{sample}/results/{version}/gsnap_filter_[1,2].fa*", expected=2),
        "post_alignment_fasta":
            smart_glob(f"{sample}/postprocess/{version}/taxid_annot.fasta", expected=1)
    }


def parse_result_dir(sample_result_dir):
    "Deconstruct '<sample_path>/results/2.8' into ('<sample_path>', '2.8')"
    try:
        assert "/results/" in sample_result_dir
        before, after = sample_result_dir.split("/results/")
        version = after.split("/")[0]
        for version_part in version.split("."):
            assert str(int(version_part)) == version_part
        return before, version
    except:
        print(f"ERROR: Expected something like '.../results/5.1', got {sample_result_dir}")
        raise


def benchmark_lineage_from_header(header_line):
    matches = re.search(BENCHMARK_LINEAGE_PATTERN, header_line)
    assert matches, "Missing or malformed benchmark_lineage tag."
    benchmark_lineage = matches.group(0)[2:-2]
    return benchmark_lineage


def idseq_lineage_from_header(header_line):
    matches = re.search(IDSEQ_LINEAGE_PATTERN, header_line)
    assert matches, f"Missing or malformed idseq lineage tags in '{header_line}'."
    tokens = matches.group(0).split(":")
    return zip(tokens[0::2], tokens[1::2])


def benchmark_lineage_to_taxid_strs(benchmark_lineage):
    return tuple(taxid_str for taxid_str in benchmark_lineage.split("_")[2:])


def accumulators_new():
    "Use result as follows: accumulators[benchmark_lineage][taxid_rank][taxid_str] += 1"
    return defaultdict(lambda: defaultdict(lambda: defaultdict(int)))


def count_fastq(input_fastq):
    assert ".fastq" in input_fastq or ".fq" in input_fastq
    accumulators = accumulators_new()
    with smarter_open(input_fastq, "rb") as input_f:
        line_number = 1
        try:
            line = smarter_readline(input_f)
            while line:
                # The FASTQ format specifies that each read consists of 4 lines,
                # the first of which begins with @ followed by read ID.
                line = line.decode('utf-8')
                assert line[0] == "@", f"fastq format requires every 4th line to start with @"
                benchmark_lineage = benchmark_lineage_from_header(line)
                taxid_strings = benchmark_lineage_to_taxid_strs(benchmark_lineage)
                for taxid_rank, taxid_str in zip(TAXID_RANKS, taxid_strings):
                    accumulators[benchmark_lineage][taxid_rank][taxid_str] += 1
                for _ in range(4):
                    line = smarter_readline(input_f)
                    line_number += 1
        except Exception as _:
            print(f"Error parsing line {line_number} in {input_fastq}.")
            raise
    return accumulators


def count_fasta(fasta_file):
    assert ".fasta" in fasta_file or ".fa" in fasta_file
    accumulators = accumulators_new()
    with smarter_open(fasta_file, "rb") as input_f:
        line_number = 1
        try:
            line = smarter_readline(input_f)
            while line:
                # The FASTA format specifies that each read header starts with ">"
                line = line.decode('utf-8')
                if line[0] == ">":
                    benchmark_lineage = benchmark_lineage_from_header(line)
                    taxid_strings = benchmark_lineage_to_taxid_strs(benchmark_lineage)
                    for taxid_rank, taxid_str in zip(TAXID_RANKS, taxid_strings):
                        accumulators[benchmark_lineage][taxid_rank][taxid_str] += 1
                line = smarter_readline(input_f)
                line_number += 1
        except Exception as _:
            print(f"Error parsing line {line_number} in {fasta_file}.")
            raise
    return accumulators


def increment(accumulators, delta):
    for benchmark_lineage in delta:
        for taxid_rank in delta[benchmark_lineage]:
            for taxid_str in delta[benchmark_lineage][taxid_rank]:
                accumulators[benchmark_lineage][taxid_rank][taxid_str] += delta[benchmark_lineage][taxid_rank][taxid_str]


def pick_from_equal(values):
    result = None
    for v in values:
        if result == None:
            result = v
        assert result == v
    return result


def condense_equal(accumulators):
    condenser = {}
    for benchmark_lineage in accumulators:
        ranks = []
        read_counts = []
        for taxid_rank in accumulators[benchmark_lineage]:
            ranks.append(taxid_rank)
            this_count = pick_from_equal(accumulators[benchmark_lineage][taxid_rank].values())
            read_counts.append(this_count)
        assert ranks == TAXID_RANKS
        condenser[benchmark_lineage] = pick_from_equal(read_counts)
    return condenser


def count_annot_fasta(fasta_file):
    assert ".fasta" in fasta_file or ".fa" in fasta_file
    accumulators = accumulators_new()
    with smarter_open(fasta_file, "rb") as input_f:
        line_number = 1
        try:
            line = smarter_readline(input_f)
            while line:
                # The FASTA format specifies that each read header starts with ">"
                line = line.decode('utf-8')
                if line[0] == ">":
                    benchmark_linage = benchmark_lineage_from_header(line)
                    idseq_lineage = idseq_lineage_from_header(line)
                    # TODO:  Does the annot fasta format exclude duplicates?   If not, ensure we exclude them here.
                    # Perhaps also assert on recall > 1.
                    for taxid_rank, taxid_str in idseq_lineage:
                        accumulators[benchmark_linage][taxid_rank][taxid_str] += 1
                line = smarter_readline(input_f)
                line_number += 1
        except Exception as _:
            print(f"Error parsing line {line_number} in {fasta_file}.")
            raise
    return accumulators


def main(args):
    assert len(args) == 2, "Sample dir argument is required.  See usage."
    sample_result_dir = args[1]
    print(f"Scoring IDSEQ benchmark output {sample_result_dir}")
    sample, version = parse_result_dir(sample_result_dir)
    sample_data = glob_sample_data(sample, version)
    print(json.dumps(sample_data, indent=4))
    r1, r2 = sample_data["input_fastq"]
    fastq_counts = count_fastq(r1)
    increment(fastq_counts, count_fastq(r2))
    r1, r2 = sample_data["post_qc_fasta"]
    fasta_counts = count_fasta(r1)
    increment(fasta_counts, count_fasta(r2))
    r = sample_data["post_alignment_fasta"][0]
    annot_counts = count_annot_fasta(r)
    tally = {
        "input_fastq": condense_equal(fastq_counts),
        "post_qc_fasta": condense_equal(fasta_counts),
        "post_alignment_fasta": annot_counts,
    }
    print(json.dumps(tally, indent=4))



if __name__ == "__main__":
    main(sys.argv)
