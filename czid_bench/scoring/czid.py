import json
import re
from collections import defaultdict
import numpy as np
from smart_open import open as smart_open
from czid_bench.util import smart_glob, smart_ls
from czid_bench.parsers import extract_accession_id, extract_fast_file_type_from_path
from .metrics import adjusted_aupr

STORE = 's3://'

HIT_SUMMARY_READ_ID = r"^(?P<read_id>.*?)\t"
BENCHMARK_LINEAGE_PATTERN = r"__benchmark_lineage_(?P<subspecies>\d+)_(?P<species>\d+)_(?P<genus>\d+)_(?P<family>\d+)__"
IDSEQ_LINEAGE_HIT_SUMMARY_PATTERN = r"\t(?P<species>-?\d+)\t(?P<genus>-?\d+)\t(?P<family>-?\d+)(?:\tfrom_assembly)?$"
RANKS = ['species', 'genus', 'family']

FAST_FILE_TYPE = r"\.(?:fast|f)(?P<type>q|a)(?:\.|$)"


class MalformedBenchmarkLineageException(Exception):
  def __init__(self, line):
    super().__init__(f"Missing or malformed benchmark_lineage tag: {line}")


class MalformedHitSummaryLineageException(Exception):
  def __init__(self, line):
    super().__init__(f"Missing or malformed benchmark_lineage tag: {line}")

class MalformedHitSummaryReadIdException(Exception):
  def __init__(self, line):
    super().__init__(f"Missing or malformed id tag: {line}")

class HitCounters:
  def __init__(self):
    self.counters = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

  def __getitem__(self, rank):
    return self.counters[rank]

  def by_rank(self, rank):
    return self.counters[rank]

  def ranks(self):
    return self.counters.keys()

  def increment(self, benchmark_lineage, lineage):
    for rank, tax_id in lineage.items():
      self.counters[rank][benchmark_lineage[rank]][tax_id] += 1

  def __str__(self):
    return json.dumps(self.counters, indent=4)


class IDseqSampleFileManager():
  """Manage download of files from IDseq
  """

  def __init__(self, project_id, sample_id, pipeline_version, env='prod', local_path=None):
    self.project_id = project_id
    self.sample_id = sample_id
    self.pipeline_version = pipeline_version
    self.env = env or 'prod'
    self.store = local_path or STORE
    self.set_directory_vars()

  def set_directory_vars(self):
    # this enables benchmark scoring to work with both SFN-WDL and original filepaths
    env_dir = f"{self.store}czid-samples-{self.env}"
    samples_dir = f"{env_dir}/samples/{self.project_id}/{self.sample_id}"
    self.input_fastq_file_pattern = rf"{samples_dir}/fastqs/.+\.(?:fast|f)q(?:\..+)?"

    results_dir = f"{samples_dir}/results/{self.pipeline_version}"
    post_process_dir = f"{samples_dir}/postprocess/{self.pipeline_version}/assembly"

    if(len(smart_ls(results_dir)) == 0):
      results_dir = f"{samples_dir}/results/czid-{self.env}-main-1/wdl-1/dag-{self.pipeline_version}"
      self.post_assembly_summary_files = {
        'NT': f"{results_dir}/gsnap.hitsummary2.tab",
        'NR': f"{results_dir}/rapsearch2.hitsummary2.tab"
      }
    else:
      self.post_assembly_summary_files = {
        'NT': f"{post_process_dir}/gsnap.hitsummary2.tab",
        'NR': f"{post_process_dir}/rapsearch2.hitsummary2.tab"
      }

    self.post_qc_fasta_file_pattern = rf"{results_dir}/gsnap_filter_[12]\.fa(?:sta)?"

  @staticmethod
  def parse_benchmark_lineage(line):
    matches = re.search(BENCHMARK_LINEAGE_PATTERN, line)
    if not matches:
      raise MalformedBenchmarkLineageException(line)

    return {
      rank: int(matches.group(rank))
      for rank in ['species', 'genus', 'family']
    }

  @staticmethod
  def parse_hit_summary_lineage(line):
    matches = re.search(IDSEQ_LINEAGE_HIT_SUMMARY_PATTERN, line)
    if not matches:
      raise MalformedHitSummaryLineageException(line)

    return {
      rank: int(matches.group(rank))
      for rank in ['species', 'genus', 'family']
    }

  @staticmethod
  def parse_hit_summary_read_id(line):
    matches = re.search(HIT_SUMMARY_READ_ID, line)
    if not matches:
      raise MalformedHitSummaryReadIdException(line)
    return matches.group("read_id")

  def hit_summary_entries(self, summary_file, skip_benchmark_lineage=False):
    try:
      with smart_open(summary_file, 'rb') as input_file:
        for line in input_file:
          line = line.decode('UTF-8')
          entry = {}
          entry['benchmark_lineage'] = None if skip_benchmark_lineage else self.parse_benchmark_lineage(line)
          entry['hit_summary_lineage'] = self.parse_hit_summary_lineage(line)
          entry['read_id'] = self.parse_hit_summary_read_id(line)
          entry['line'] = line

          yield entry
    except GeneratorExit:
      # If the generator is closing, we should just reraise and not output an error
      raise
    except:
      print(f"[ERROR] Parsing file: {summary_file}")
      raise

  def post_assembly_hit_summary_entries(self, db_type, skip_benchmark_lineage=False):
    return self.hit_summary_entries(
      self.post_assembly_summary_files[db_type],
      skip_benchmark_lineage=skip_benchmark_lineage
    )

  def parse_fastx_entry(self, entry):
    parsed_entry = {
      "lineage": self.parse_benchmark_lineage(entry[0]),
      "accession_id": extract_accession_id(entry[0]),
      "read": entry[1].strip()
    }

    if len(entry) == 4:
      parsed_entry['quality'] = entry[3].strip()

    return parsed_entry

  def fastx_iterator(self, fastx_file, file_type="q"):
    file_type = extract_fast_file_type_from_path(fastx_file)

    lines_per_entry = 4 if file_type == "q" else 2
    read_number = 1
    try:
      with smart_open(fastx_file) as input_file:
        entry_first_line = input_file.readline()
        while entry_first_line:
          entry = self.parse_fastx_entry([entry_first_line] + [
            input_file.readline()
            for _ in range(lines_per_entry - 1)
          ])
          yield entry
          read_number += 1
          entry_first_line = input_file.readline()
    except GeneratorExit:
      # If the generator is closing, we should just reraise and not output an error
      raise
    except:
      print(f"[ERROR] Parsing read number {read_number} in {fastx_file}")
      raise

  def input_files(self):
    return smart_glob(self.input_fastq_file_pattern, expected_num_files=[1, 2])

  def post_qc_files(self):
    return smart_glob(self.post_qc_fasta_file_pattern, expected_num_files=[1, 2])

def lineage_key(lineage_dict):
  return "{species}:{genus}:{family}".format(**lineage_dict)

def key_to_lineage(key):
  return {k: int(v) for k, v in zip(["species", "genus", "family"], key.split(":"))}

def hit_summary_counts_per_benchmark_lineage(czid_file_manager, db_type, counters=None):
  counters = counters or HitCounters()
  for entry in czid_file_manager.post_assembly_hit_summary_entries(db_type):
    counters.increment(entry['benchmark_lineage'], entry['hit_summary_lineage'])
  return counters

def hit_summary_counts_per_tax_id(czid_file_manager, db_type):
  print(f" * Counting hits per tax id for {db_type}")
  counters = defaultdict(lambda: defaultdict(int))
  for entry in czid_file_manager.post_assembly_hit_summary_entries(db_type, skip_benchmark_lineage=True):
    for rank, tax_id in entry['hit_summary_lineage'].items():
      counters[rank][tax_id] += 1
  return counters

def hit_summary_concordance(czid_file_manager):
  concordance_counters = defaultdict(int)
  hit_by_read_id = {}
  # Loop through both hit summary files simultaneously to take advantage of
  # similarly sorted entries
  for nt_hit_summary_entry, nr_hit_summary_entry in zip(
      czid_file_manager.post_assembly_hit_summary_entries('NT'),
      czid_file_manager.post_assembly_hit_summary_entries('NR')
    ):
    nt_czid_lineage, nt_read_id = nt_hit_summary_entry['hit_summary_lineage'], nt_hit_summary_entry['read_id']
    nr_czid_lineage, nr_read_id = nr_hit_summary_entry['hit_summary_lineage'], nr_hit_summary_entry['read_id']
    for czid_lineage, read_id in zip([nt_idseq_lineage, nr_idseq_lineage], [nt_read_id, nr_read_id]):
      if read_id in hit_by_read_id:
        for rank, tax_id in czid_lineage.items():
          if hit_by_read_id[read_id][rank] == tax_id:
            concordance_counters[tax_id] += 1
        del hit_by_read_id[read_id]
      else:
        hit_by_read_id[read_id] = czid_lineage
  return concordance_counters


def count_reads_per_benchmark_lineage(czid_file_manager, fastx_files):
  counters = {}
  for fastx_file in fastx_files:
    for entry in czid_file_manager.fastx_iterator(fastx_file):
      for _, tax_id in entry['lineage'].items():
        counters[tax_id] = counters.get(tax_id, 0) + 1
  return counters

def count_hits_per_benchmark_lineage(czid_file_manager):
  counts_nt = hit_summary_counts_per_benchmark_lineage(czid_file_manager, 'NT')
  counts_nr = hit_summary_counts_per_benchmark_lineage(czid_file_manager, 'NR')
  return counts_nt, counts_nr

def count_hits_per_tax_id(czid_file_manager):
  counts_nt = hit_summary_counts_per_tax_id(czid_file_manager, 'NT')
  counts_nr = hit_summary_counts_per_tax_id(czid_file_manager, 'NR')
  return counts_nt, counts_nr

def score_benchmark(project_id, sample_id, pipeline_version, env='prod', local_path=None, force_monotonic=False):
  czid_file_manager = IDseqSampleFileManager(project_id, sample_id, pipeline_version, env=env, local_path=local_path)

  print(" * Counting reads from input files")
  input_reads_by_tax_id = count_reads_per_benchmark_lineage(czid_file_manager, idseq_file_manager.input_files())
  print(" * Counting reads from post qc files")
  post_qc_reads_by_tax_id = count_reads_per_benchmark_lineage(czid_file_manager, idseq_file_manager.post_qc_files())
  print(" * Counting hits per benchmark lineage")
  hit_counters_nt, hit_counters_nr = count_hits_per_benchmark_lineage(czid_file_manager)
  print(" * Counting corcordant hits per taxon id")
  concordance_by_tax_id = hit_summary_concordance(czid_file_manager)
  stats = {
    'per_rank': {}
  }

  ranks = sorted(set(hit_counters_nt.ranks()) | set(hit_counters_nr.ranks()))
  for rank in ranks:
    stats_per_rank = stats['per_rank'].setdefault(rank, {})

    total_reads_per_rank = sum(
      input_reads_by_tax_id[benchmark_tax_id]
      for benchmark_tax_id in hit_counters_nt.by_rank(rank)
      if benchmark_tax_id > 0)
    total_post_qc_reads_per_rank = sum(
      post_qc_reads_by_tax_id[benchmark_tax_id]
      for benchmark_tax_id in hit_counters_nt.by_rank(rank)
      if benchmark_tax_id > 0)

    for db_type, hit_counters in zip(['NT', 'NR'], [hit_counters_nt, hit_counters_nr]):
      stats_per_db_type = stats_per_rank.setdefault(db_type, {})
      benchmark_hits = hit_counters.by_rank(rank)
      total_correct_reads_per_db_type = 0
      for benchmark_tax_id in benchmark_hits.keys():
        stats_by_tax_id = stats_per_db_type.setdefault(benchmark_tax_id, {})
        stats_by_tax_id['total_reads'] = input_reads_by_tax_id[benchmark_tax_id]
        stats_by_tax_id['post_qc_reads'] = post_qc_reads_by_tax_id[benchmark_tax_id]
        stats_by_tax_id['recall_per_read'] = {
          'count': benchmark_hits[benchmark_tax_id][benchmark_tax_id],
          'value': benchmark_hits[benchmark_tax_id][benchmark_tax_id] / post_qc_reads_by_tax_id[benchmark_tax_id],
        }
        total_correct_reads_per_db_type += benchmark_hits[benchmark_tax_id][benchmark_tax_id]

      stats_per_db_type['accuracy'] = {
        'count': total_correct_reads_per_db_type,
        'value': total_correct_reads_per_db_type/total_post_qc_reads_per_rank
      }

      czid_hit_counters = defaultdict(int)
      bench_hit_counters = defaultdict(int)
      for bench_tax_id, czid_hits in benchmark_hits.items():
        for czid_tax_id, counts in idseq_hits.items():
          czid_hit_counters[idseq_tax_id] += counts
          bench_hit_counters[bench_tax_id] += counts

      truth_taxa = [
        {'tax_id': tax_id, 'abs_abundance': counts}
        for tax_id, counts in bench_hit_counters.items()
      ]
      sample_level_metrics = metrics_per_sample(czid_hit_counters, truth_taxa, force_monotonic=force_monotonic)
      stats_per_db_type.update(sample_level_metrics)

    stats_concordance = {}
    benchmark_tax_ids = set(hit_counters_nt.by_rank(rank).keys()) | set(hit_counters_nr.by_rank(rank).keys())
    for benchmark_tax_id in benchmark_tax_ids:
      stats_concordance[benchmark_tax_id] = {
        "count": concordance_by_tax_id[benchmark_tax_id],
        "value": concordance_by_tax_id[benchmark_tax_id] / post_qc_reads_by_tax_id[benchmark_tax_id]
      }
    stats_per_rank['concordance'] = stats_concordance

    stats_per_rank['total_reads'] = total_reads_per_rank
    stats_per_rank['post_qc_reads'] = total_post_qc_reads_per_rank

  return stats

def metrics_per_sample(hit_counters, truth_taxa, force_monotonic=False):
  stats = {}

  total_simulated_taxa = len(truth_taxa)
  total_correctly_identified_taxa = sum(1 for taxon in truth_taxa if taxon['tax_id'] in hit_counters)
  total_identified_taxa = len(hit_counters)

  stats['total_simulated_taxa'] = total_simulated_taxa
  stats['total_correctly_identified_taxa'] = total_correctly_identified_taxa
  stats['total_identified_taxa'] = total_identified_taxa
  recall = total_correctly_identified_taxa / total_simulated_taxa
  precision = total_correctly_identified_taxa / total_identified_taxa
  stats['recall'] = recall
  stats['precision'] = precision
  stats['f1-score'] = 2 * recall * precision / (recall + precision)

  # AUPR - Area Under Precision and Recall curve
  tax_ids = hit_counters.keys()
  benchmark_tax_ids_set = set(taxon['tax_id'] for taxon in truth_taxa)
  missed_tax_ids = [tax_id for tax_id in benchmark_tax_ids_set if tax_id not in tax_ids]
  y_true = [1 if tax_id in benchmark_tax_ids_set else 0 for tax_id in tax_ids] + [1] * len(missed_tax_ids)
  total_sum = sum(hit_counters.values())
  y_score = [hit_counters[tax_id] / total_sum for tax_id in tax_ids] + [0] * len(missed_tax_ids)
  aupr_results = adjusted_aupr(y_true, y_score, force_monotonic=force_monotonic)
  stats['aupr'] = aupr_results["aupr"]

  total_benchmark_reads = sum(taxon['abs_abundance'] for taxon in truth_taxa)
  relative_abundances_diff = [
    taxon['abs_abundance']/total_benchmark_reads - hit_counters[taxon['tax_id']]/total_benchmark_reads
    for taxon in truth_taxa
  ]
  l1_norm = np.linalg.norm(relative_abundances_diff, ord=1)
  l2_norm = np.linalg.norm(relative_abundances_diff, ord=2)
  stats['l1_norm'] = l1_norm
  stats['l2_norm'] = l2_norm

  return stats


def score_sample(project_id, sample_id, pipeline_version, truth_taxa, env='prod', local_path=None, force_monotonic=False):
  czid_file_manager = IDseqSampleFileManager(project_id, sample_id, pipeline_version, env=env, local_path=local_path)
  hit_counters_nt, hit_counters_nr = count_hits_per_tax_id(czid_file_manager)

  stats = {'per_rank': {}}

  ranks = sorted(truth_taxa.keys())
  for rank in ranks:
    stats_per_rank = stats['per_rank'].setdefault(rank, {})
    for db_type, hit_counters in zip(['NT', 'NR'], [hit_counters_nt, hit_counters_nr]):
      stats_per_rank[db_type] = metrics_per_sample(hit_counters[rank], truth_taxa[rank], force_monotonic=force_monotonic)

  return stats
