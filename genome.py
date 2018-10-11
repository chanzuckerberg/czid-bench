import os
from util import remove_safely, check_call

DOWNLOADS_SUBDIR = "downloads"

class Genome:

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
        self.versioned_accession_ids = [Genome.ensure_versioned(vaccid) for vaccid in versioned_accession_ids]
        self.key = f"{category}__{organism}__{self.taxid}"
        self.filename = f"{DOWNLOADS_SUBDIR}/{self.key}.fasta"
        # The size (number of bases) is filled in by fetch_all()
        self.size = None
        Genome.all[self.key] = self
        for vaccid in self.versioned_accession_ids:
            assert vaccid not in Genome.by_accid, "Accession ID {vaccid} should not occur in multiple genomes."
            Genome.by_accid[vaccid] = self

    @staticmethod
    def ensure_versioned(vaccid):
        if "." not in vaccid:
            vaccid += ".1"
            print(f"INFO:  Changed {vaccid[:-2]} to {vaccid} in order to match ISS-generated read headers.")
        return vaccid

    @staticmethod
    def fetch_versioned_accession_id(vaccid):  # e.g., "NC_004325.2"
        output_file = f"{DOWNLOADS_SUBDIR}/{vaccid}.fa"
        os.makedirs(DOWNLOADS_SUBDIR, exist_ok=True)
        if os.path.isfile(output_file):
            print(f"{output_file} already exists, nice")
        else:
            try:
                command = f"cd {DOWNLOADS_SUBDIR}; ncbi-acc-download --format fasta {vaccid} -e all"
                check_call(command)
            except:
                remove_safely(output_file)
                raise
        return output_file

    @staticmethod
    def fetch_all():
        for g in Genome.all.values():
            subdir = g.filename.rsplit("/", 1)[0]
            os.makedirs(subdir, exist_ok=True)
            remove_safely(g.filename)
            accession_fas = []
            for f in (g.versioned_accession_ids):
                af = Genome.fetch_versioned_accession_id(f)
                accession_fas.append(af)
            accession_fastas = " ".join(accession_fas)
            command = f"cat {accession_fastas} > {g.filename}"
            check_call(command)
            assert os.path.isfile(g.filename), f"Failed to download genome {g.filename}"
            remove_safely(f"{g.key}.size")
            command = f"grep -v '^>' {g.filename} | tr -d '\n' | wc > {g.key}.size"
            check_call(command)
            with open(f"{g.key}.size") as f:
                line = f.readline().rstrip()
                g.size = int(line.split()[2])
            remove_safely(f"{g.key}.size")
            print(f"Genome {g.key} size {g.size} bases.")
