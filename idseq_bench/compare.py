#!/usr/bin/env python3
# Compare json output after run the score.py

import json
import sys

def main():
    args = sys.argv
    assert len(args) >= 2, "Usage: python3 compare.py <score1.json>  [<score2.json> <score3.json> ...]"
    bundles = []

    for json_file in args[1:]:
        benchmark_out = json.load(open(json_file, 'r'))
        bundle = (benchmark_out['sample_data']['sample'],
                  benchmark_out['sample_data']['version'],
                  benchmark_out["stats"])
        bundles.append(bundle)

    headers = ['sample', 'version', 'lineage', 'total_reads', 'survived_qc',
               'correct_family_nt', 'correct_family_nr', 'family_best_post_qc',
               'correct_genus_nt', 'correct_genus_nr', 'genus_best_post_qc',
               'correct_species_nt', 'correct_species_nr', 'species_best_post_qc']
    print("\t".join(headers))


    lineages = bundles[0][2].keys()

    for bl in lineages:
        for bundle in bundles:
            sample = bundle[0]
            version = bundle[1]
            stats = bundle[2]
            details = stats[bl]["recalled_correctly"]
            fields = [
                sample, version, bl, stats[bl]['total_reads']['count'], stats[bl]['survived_qc']['count'],
                details["family"]["nt"], details["family"]["nr"], details["family"]["best_post_qc"],
                details["genus"]["nt"], details["genus"]["nr"], details["genus"]["best_post_qc"],
                details["species"]["nt"], details["species"]["nr"], details["species"]["best_post_qc"]
            ]
            fields = list(map(str, fields))
            print("\t".join(fields))

if __name__ == "__main__":
    main()
