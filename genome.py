import os
from util import remove_safely, check_call


class Genome(object):

    all = dict()
    by_accid = dict()

    def __init__(self, category, organism, lineage, versioned_accession_ids, genome_assembly_URL=None):
        self.category = category
        self.organism = organism
        self.lineage = lineage
        assert ["subspecies", "species", "genus", "family"] == [l[0] for l in lineage]
        self.subspecies_taxid, self.species_taxid, self.genus_taxid, self.family_taxid = [(l[1] or 0) for l in lineage]
        self.taxid = self.subspecies_taxid or self.species_taxid
        self.genome_assembly_URL = genome_assembly_URL
        self.versioned_accession_ids = versioned_accession_ids
        self.key = f"{category}__{organism}__{self.taxid}"
        self.filename = f"{self.key}.fasta"
        Genome.all[self.key] = self
        for accid in versioned_accession_ids:
            vaccid = accid
            if "." not in vaccid:
                vaccid = accid + ".1"
                print(f"INFO:  Changing {accid} to {vaccid} to match ISS headers.")
            Genome.by_accid[vaccid] = self

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
