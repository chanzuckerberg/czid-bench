import json
import re
from collections import defaultdict
from functools import reduce
from idseq_bench.util import smart_glob, smarter_open, smarter_readline
from idseq_bench.parsers import extract_accession_id, extract_fast_file_type_from_path

STORE = 's3://'
ENV_DIR = f"{{store_path}}idseq-samples-{{env}}"
SAMPLES_DIR = f"{ENV_DIR}/samples/{{project_id}}/{{sample_id}}"
RESULTS_DIR = f"{SAMPLES_DIR}/results/{{pipeline_version}}"
POST_PROCESS_DIR = f"{SAMPLES_DIR}/postprocess/{{pipeline_version}}"

INPUT_FASTQ_FILE_PATTERN = f"{SAMPLES_DIR}/fastqs/*.fastq*"
POST_QC_FASTA_FILE_PATTERN = f"{RESULTS_DIR}/gsnap_filter_[1,2].fa*"
POST_ALIGNMENT_FASTA_FILE_PATTERN = f"{POST_PROCESS_DIR}/taxid_annot.fasta"
POST_ASSEMBLY_SUMMARY_FILES = {
  'NT': f"{POST_PROCESS_DIR}/assembly/gsnap.hitsummary2.tab",
  'NR': f"{POST_PROCESS_DIR}/assembly/rapsearch2.hitsummary2.tab"
}

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
    # print(f"Update lineage: {lineage}")
    for rank, tax_id in lineage.items():
      self.counters[rank][benchmark_lineage[rank]][tax_id] += 1

  def __str__(self):
    return json.dumps(self.counters, indent=4)
    # return str(self.counters)



class IDseqSampleFileManager():
  """Manage download of files from IDseq
  """

  def __init__(self, project_id, sample_id, pipeline_version, env='prod', cached_path=None):
    self.project_id = project_id
    self.sample_id = sample_id
    self.pipeline_version = pipeline_version
    self.env = env
    self.store_path = cached_path or STORE

  def apply_context(self, format_str):
    return format_str.format(
      store=self.store_path,
      env=self.env,
      project_id=self.project_id,
      sample_id=self.sample_id,
      pipeline_version=self.pipeline_version
    )

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

  def hit_summary_entries(self, summary_file):
    with smarter_open(summary_file) as input_file:
      line = smarter_readline(input_file)
      while line:
        line = line.decode('UTF-8')
        try:
          entry = {}
          entry['benchmark_lineage'] = self.parse_benchmark_lineage(line)
          entry['hit_summary_lineage'] = self.parse_hit_summary_lineage(line)
          entry['read_id'] = self.parse_hit_summary_read_id(line)
          entry['line'] = line

          yield entry
          line = smarter_readline(input_file)
        except:
          print(f"[ERROR] Parsing file: {summary_file}")
          raise

  def post_assembly_hit_summary_entries(self, db_type):
    return self.hit_summary_entries(
      self.apply_context(POST_ASSEMBLY_SUMMARY_FILES[db_type])
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
      with smarter_open(fastx_file) as input_file:
        entry_first_line = smarter_readline(input_file).decode("utf-8")
        while entry_first_line:
          entry = self.parse_fastx_entry([entry_first_line] + [
            smarter_readline(input_file).decode("utf-8")
            for _ in range(lines_per_entry - 1)
          ])
          yield entry
          read_number += 1
          entry_first_line = smarter_readline(input_file).decode("utf-8")

    except Exception as e:
      # TODO: move to proper assertion
      print(f"[ERROR] Parsing read number {read_number} in {fastx_file}")
      raise e

  def input_files(self):
    return smart_glob(self.apply_context(INPUT_FASTQ_FILE_PATTERN), expected_num_files=[1, 2])

  def post_qc_files(self):
    return smart_glob(self.apply_context(POST_QC_FASTA_FILE_PATTERN), expected_num_files=[1, 2])


def lineage_key(lineage_dict):
  return "{species}:{genus}:{family}".format(**lineage_dict)

def key_to_lineage(key):
  return {k: int(v) for k, v in zip(["species", "genus", "family"], key.split(":"))}

def hit_summary_counts(idseq_file_manager, db_type, counters=None):
  counters = counters or HitCounters()
  for entry in idseq_file_manager.post_assembly_hit_summary_entries(db_type):
    counters.increment(entry['benchmark_lineage'], entry['hit_summary_lineage'])
  return counters

def hit_summary_concordance(idseq_file_manager):
  concordance_counters = defaultdict(int)
  hit_by_read_id = {}
  # Loop through both hit summary files simultaneously to take advantage of
  # similarly sorted entries
  for nt_hit_summary_entry, nr_hit_summary_entry in zip(
      idseq_file_manager.post_assembly_hit_summary_entries('NT'),
      idseq_file_manager.post_assembly_hit_summary_entries('NR')
    ):
    nt_idseq_lineage, nt_read_id = nt_hit_summary_entry['hit_summary_lineage'], nt_hit_summary_entry['read_id']
    nr_idseq_lineage, nr_read_id = nr_hit_summary_entry['hit_summary_lineage'], nr_hit_summary_entry['read_id']
    for idseq_lineage, read_id in zip([nt_idseq_lineage, nr_idseq_lineage], [nt_read_id, nr_read_id]):
      if read_id in hit_by_read_id:
        for rank, tax_id in idseq_lineage.items():
          if hit_by_read_id[read_id][rank] == tax_id:
            concordance_counters[tax_id] += 1
        del hit_by_read_id[read_id]
      else:
        hit_by_read_id[read_id] = idseq_lineage
  return concordance_counters


def count_reads_per_benchmark_lineage(idseq_file_manager, fastx_files):
  counters = {}
  for fastx_file in fastx_files:
    for entry in idseq_file_manager.fastx_iterator(fastx_file):
      for _, tax_id in entry['lineage'].items():
        counters[tax_id] = counters.get(tax_id, 0) + 1
  return counters

def count_hits(idseq_file_manager):
  counts_nt = hit_summary_counts(idseq_file_manager, 'NT')
  counts_nr = hit_summary_counts(idseq_file_manager, 'NR')
  return counts_nt, counts_nr

def score(project_id, sample_id, pipeline_version):
  idseq_file_manager = IDseqSampleFileManager(project_id, sample_id, pipeline_version)


  print(" * Counting read from input files")
  input_reads_by_tax_id = count_reads_per_benchmark_lineage(idseq_file_manager, idseq_file_manager.input_files())
  print(" * Counting read from post qc files")
  post_qc_reads_by_tax_id = count_reads_per_benchmark_lineage(idseq_file_manager, idseq_file_manager.post_qc_files())
  print(" * Counting hits per benchmark lineage")
  hit_counters_nt, hit_counters_nr = count_hits(idseq_file_manager)
  print(" * Counting corcordant hits per taxon id")
  concordance_by_tax_id = hit_summary_concordance(idseq_file_manager)
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

      total_simulated_taxa = sum(1 for benchmark_tax_id in benchmark_hits.keys() if benchmark_tax_id > 0)
      total_correctly_identified_taxa = sum(1 for benchmark_tax_id in benchmark_hits.keys() if benchmark_tax_id > 0 and benchmark_hits[benchmark_tax_id][benchmark_tax_id])
      total_identified_taxa = len(reduce(
        lambda curr_set, identified_taxa: curr_set | set(identified_taxa.keys()),
        benchmark_hits.values(),
        set()
      ))
      stats_per_db_type['total_simulated_taxa'] = total_simulated_taxa
      stats_per_db_type['total_correctly_identified_taxa'] = total_correctly_identified_taxa
      stats_per_db_type['total_identified_taxa'] = total_identified_taxa
      recall = total_correctly_identified_taxa / total_simulated_taxa
      precision = total_correctly_identified_taxa / total_identified_taxa
      stats_per_db_type['recall'] = recall
      stats_per_db_type['precision'] = precision
      stats_per_db_type['f1-score'] = 2 * recall * precision / (recall + precision)

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
