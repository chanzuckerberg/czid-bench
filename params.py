from genome import Genome


NUM_READS_PER_ORGANISM = 10 * 1000


# TODO: train "iseq" model
MODELS = ["hiseq"] # ["novaseq", "miseq", "hiseq"]


# We're committed to uniform.
UNIFORM_ABUNDANCE = "uniform"


# Used in standard benchmarks.
def top_6_id_genomes():
  return [

    Genome("bacteria", "klebsiella_pneumoniae_HS11286",
      [("subspecies", 1125630), ("species", 573), ("genus", 570), ("family", 543)],
      ["NC_016845.1"]),

    Genome("fungi", "aspergillus_fumigatus",
      [("subspecies", 330879), ("species", 746128), ("genus", 5052), ("family", 1131492)],
      ["NC_007194.1", "NC_007195.1", "NC_007196.1", "NC_007197.1", "NC_007198.1", "NC_007199.1",
       "NC_007200.1", "NC_007201.1"],
      "https://www.ncbi.nlm.nih.gov/genome/18?genome_assembly_id=22576"),

    Genome("viruses", "rhinovirus_c",
      [("subspecies", 0), ("species", 463676), ("genus", 12059), ("family", 12058)],
      ["MG148341.1"]),

    Genome("viruses", "chikungunya",
      [("subspecies", 0), ("species", 37124), ("genus", 11019), ("family", 11018)],
      ["MG049915.1"]),

    Genome("bacteria", "staphylococcus_aureus",
      [("subspecies", 93061), ("species", 1280), ("genus", 1279), ("family", 90964)],
      ["NC_007795.1"],
      "https://www.ncbi.nlm.nih.gov/genome/154?genome_assembly_id=299272"),

    Genome("protista", "plasmodium_falciuparum",
      [("subspecies", 36329), ("species", 5833), ("genus", 5820), ("family", 1639119)],
      ["NC_004325.2", "NC_037280.1", "NC_000521.4", "NC_004318.2", "NC_004326.2",
       "NC_004327.3", "NC_004328.3", "NC_004329.3", "NC_004330.2", "NC_037281.1",
       "NC_037282.1", "NC_037284.1", "NC_004331.3", "NC_037283.1", "NC_036769.1"],
      "https://www.ncbi.nlm.nih.gov/genome/33?genome_assembly_id=22642"),

    Genome("eukaryota", "balamuthia_mandrillaris",
      [("subspecies", 0), ("species", 66527), ("genus", 66526), ("family", 555408)],
      ["NC_027736.1"]),           

    Genome("viruses", "rubella_virus", 
      [("subspecies", 0), ("species", 11041), ("genus", 11040), ("family", 2560066)], 
      ["NC_001545.2"]), 

    Genome("bacteria", "Elizabethkingia_anophelis_NUHP1",
      [("subspecies", 1338011 ), ("species", 1117645), ("genus", 308865), ("family", 49546)],
      ["NZ_CP007547.1"]),                          

    Genome("bacteria", "Neisseria_meningitidis_MC58",
      [("subspecies", 122586 ), ("species", 487), ("genus", 482), ("family", 481)],
      ["NC_003112.2"]),  

    Genome("viruses", "Torque_teno_midi_virus_11",
      [("subspecies", 0 ), ("species", 2065052), ("genus", 687333), ("family", 687329)],
      ["NC_038358.1"]),                 

    Genome("viruses", "Human_immunodeficiency_virus_1",
      [("subspecies", 0 ), ("species", 11676), ("genus", 11646), ("family", 11632)],
      ["NC_001802.1"]),                        

    Genome("viruses", "Hubei_mosquito_virus_2",
      [("subspecies", 0 ), ("species", 1922926), ("genus", 0), ("family", 0)],
      ["NC_033305.1", "NC_033306.1"])

    ]

# test a single common microbe with genomic similarity to other microbes
def ecoli_cov():
  return [
    Genome("bacteria", "Escherichia_coli_str_K-12_substr_MG1655",
      [("subspecies",511145), ("species",562), ("genus",561), ("family",543)],
      ["NC_000913.3"],
      "https://www.ncbi.nlm.nih.gov/genome/167?genome_assembly_id=161521")
  ]

# human host removal filters
def human_host_removal():
  return[
    Genome("eukaryota", "homo_sapiens",
      [("subspecies", 0), ("species", 9606), ("genus", 9605), ("family", 9604)],
      ["NC_000001.11", "NC_000002.12", "NC_000003.12","NC_000004.12","NC_000005.10",
      "NC_000006.12","NC_000007.14","NC_000008.11","NC_000009.12","NC_000010.11",
      "NC_000011.10","NC_000012.12","NC_000013.11","NC_000014.9","NC_000015.10",
      "NC_000016.10","NC_000017.11","NC_000018.10","NC_000019.10","NC_000020.11",
      "NC_000021.9","NC_000022.11","NC_000023.11","NC_000024.10"],
      "https://www.ncbi.nlm.nih.gov/genome/51?genome_assembly_id=412355")
  ]

# closely related bacterial species from enterobacteriaceae family
def enterobacteriaceae_family():
  return[

    Genome("bacteria", "Salmonella_enterica_str_CT18",
      [("subspecies", 220341 ), ("species", 28901), ("genus", 590), ("family", 543)],
      ["NC_003198.1","NC_003384.1","NC_003385.1"], 
      "https://www.ncbi.nlm.nih.gov/genome/152?genome_assembly_id=154363"),    

    Genome("bacteria", "Enterobacter_aerogenes",
      [("subspecies",1028307), ("species",548), ("genus",570), ("family",543)],
      ["NC_015663.1"], 
      "https://www.ncbi.nlm.nih.gov/genome/3417?genome_assembly_id=173101"),

    Genome("bacteria", "Enterobacter_cloacae",
      [("subspecies",716541), ("species",550), ("genus",547), ("family",543)],
      ["NC_014121.1","NC_014107.1","NC_014108.1"],
      "https://www.ncbi.nlm.nih.gov/genome/1219?genome_assembly_id=170919"),

    Genome("bacteria", "Klebsiella_oxytoca",
      [("subspecies",0), ("species",571), ("genus",570), ("family",543)],
      ["NZ_CP011636.1","NZ_CP011627.1","NZ_CP011633.1","NZ_CP011628.1","NZ_CP011625.1","NZ_CP011634.1","NZ_CP011629.1","NZ_CP011630.1","NZ_CP011631.1","NZ_CP011626.1","NZ_CP011632.1","NZ_CP011635.1"],
      "https://www.ncbi.nlm.nih.gov/genome/1165?genome_assembly_id=232487"),

    Genome("bacteria", "Klebsiella_pneumoniae_HS11286",
      [("subspecies",1125630), ("species",573), ("genus",570), ("family",543)],
      ["NC_016845.1","NC_016838.1","NC_016846.1","NC_016839.1","NC_016840.1","NC_016847.1","NC_016841.1"],
      "https://www.ncbi.nlm.nih.gov/genome/815?genome_assembly_id=168877"),

    Genome("bacteria", "Citrobacter_freundii_CFNIH1",
      [("subspecies",1333848), ("species",546), ("genus",544), ("family",543)],
      ["NZ_CP007557.1","NZ_CP007558.1"],
      "https://www.ncbi.nlm.nih.gov/genome/2850?genome_assembly_id=172712"),

    Genome("bacteria", "Citrobacter_koseri",
      [("subspecies",0), ("species",545), ("genus",544), ("family",543)],
      ["NZ_CP022073.1","NZ_CP022075.1"],
      "https://www.ncbi.nlm.nih.gov/genome/1169?genome_assembly_id=366944"),

    Genome("bacteria", "Escherichia_coli_str_K-12_substr_MG1655",
      [("subspecies",511145), ("species",562), ("genus",561), ("family",543)],
      ["NC_000913.3"],
      "https://www.ncbi.nlm.nih.gov/genome/167?genome_assembly_id=161521"),

    Genome("bacteria", "Shigella_flexneri_2a_str_301",
      [("subspecies",198214), ("species", 623), ("genus", 620), ("family", 543)],
      ["NC_004337.2", "NC_004851.1"],
      "https://www.ncbi.nlm.nih.gov/genome/182?genome_assembly_id=165033"),

    Genome("bacteria", "Shigella_boydii",
      [("subspecies", 0), ("species", 621), ("genus", 620), ("family", 543)],
      ["NZ_CP011511.1"],
      "https://www.ncbi.nlm.nih.gov/genome/496?genome_assembly_id=232367")

  ]


# Used in mutated virus benchmarks.
def mutated_rhino_c_genome(percent):
    return Genome("viruses", "rhinovirus_c",
                  [("subspecies", 1000000000 + percent), ("species", 463676), ("genus", 12059), ("family", 12058)],
                  ["DERISI_HRC_{:03d}.1".format(percent)])


def mutated_rhinovirus_c_genomes():
    return  [
        mutated_rhino_c_genome(percent)
        for percent in [100, 99, 95, 90, 85, 80, 75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25]
    ]
