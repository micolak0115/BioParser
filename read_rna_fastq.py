# %%
from fileparser import Fastq
import re

# %%
def get_reverse_complement_sequence(seq):
    complement_dna = {"A": "T", "G": "C", "C": "G", "T": "A", "N": "N"}
    
    rev_comp_seq = ""
    for base in seq:
        if not base in complement_dna.keys():
            print(f"Unrecognized Base : {base}")
            raise Exception
    
        complement_base = complement_dna[base]
        rev_comp_seq += complement_base

    return rev_comp_seq

# %%
input_fastq_R1_file = "/BiO/Store/KOGIC/RNASeq/ProcessedData/Korea10K_1/Results_GENCODE_v43/1_fastp/KU10K-02991/KU10K-02991_L01_R1.trimmed.fq.gz"
input_fastq_R2_file = "/BiO/Store/KOGIC/RNASeq/ProcessedData/Korea10K_1/Results_GENCODE_v43/1_fastp/KU10K-02991/KU10K-02991_L01_R2.trimmed.fq.gz"
pract_fastq_R1_file = "/BiO/Research/GeroExpressome/Resources/Scripts/practice_R1.trimmed.fq"
pract_fastq_R2_file = "/BiO/Research/GeroExpressome/Resources/Scripts/practice_R2.trimmed.fq"

fastq_R1 = Fastq(input_fastq_R1_file, sep=None)
fastq_R2 = Fastq(input_fastq_R2_file, sep=None)

fastq_R1.make_practice_file(pract_fastq_R1_file, line_count=100000, file_obj="Fastq")
fastq_R2.make_practice_file(pract_fastq_R2_file, line_count=100000, file_obj="Fastq")

# %%
def check_if_read_set_is_valid(read_set):
    list_is_valid = list(map(lambda x: x != None, read_set))
    is_set_valid = (sum(list_is_valid) == 4)
    is_set_end = (sum(list_is_valid) == 0)

    return is_set_valid, is_set_end

def get_a_single_read_set(fastq):
    header = fastq.get_nextline()
    sequence = fastq.get_nextline()
    optional = fastq.get_nextline()
    quality = fastq.get_nextline()

    read_set = (header, sequence, optional, quality)

    return read_set

def get_header(read_set):
    header = read_set[0]

    return header

def get_sequence(read_set):
    sequence = read_set[1]

    return sequence

def get_optional(read_set):
    optional = read_set[2]

    return optional

def get_quality(read_set):
    quality = read_set[3]

    return quality

def get_read_pair_info_from_a_read_set(read_set):
    header = get_header(read_set)
    assert header.startswith('@'), f"Sequence ID must starts with '@' : {header}"
    
    read_pair = header.split("/")[-1]

    return read_pair

def check_correct_read_pair(read_set1, read_set2):
    read_pair1 = get_read_pair_info_from_a_read_set(read_set1)
    read_pair2 = get_read_pair_info_from_a_read_set(read_set2)
    if (int(read_pair1) == 1 and int(read_pair2) == 2):

        return True
    
    return False

def merge_sequences_paired_end_reads(seq1, seq2, min_overlap=3):
    def find_overlap_sequences(seq1, seq2, min_overlap):
        max_overlap = 0
        overlap_sequence = ''

        for i in range(min_overlap, min(len(seq1), len(seq2)) + 1):
            if seq1[-i:] == seq2[:i]:
                max_overlap = i
                overlap_sequence = seq1[-i:]

        return max_overlap, overlap_sequence

    overlap_length, overlap_seq = find_overlap_sequences(seq1, seq2, min_overlap)

    if overlap_length > 0:
        merged_sequence = seq1 + seq2[overlap_length:]

    else:
        merged_sequence = seq1 + seq2

    return merged_sequence, overlap_length, overlap_seq

def get_platform_name(read_set):
    dict_instrument = fastq_R1.instrument_id
    list_platform = list()
    for instrument_id, platform in dict_instrument.items():
        header = get_header(read_set)
        match_result = re.match(instrument_id, header)
        if match_result != None:
            list_platform.append(platform)
    if len(list_platform) == 1:
        platform_name = list_platform[0]
        
    elif len(list_platform) == 0:
        platform_name = "Unknown"
    
    else:
        print("Multiple Platform IDs Detected")
        raise Exception

    return platform_name

# %%
fastq_R1_prac = Fastq(input_fastq_R1_file, sep=None, is_gzip=False)
fastq_R2_prac = Fastq(input_fastq_R2_file, sep=None, is_gzip=False)

list_overlap_length = list()
while 1:
    read_set_R1 = get_a_single_read_set(fastq_R1_prac)
    read_set_R2 = get_a_single_read_set(fastq_R2_prac)
    is_set_valid_R1, is_set_end_R1 = check_if_read_set_is_valid(read_set_R1)
    is_set_valid_R2, is_set_end_R2 = check_if_read_set_is_valid(read_set_R2)
    if is_set_end_R2:
        break
    if not is_set_valid_R2:
        print(f"File has unrecognizable read. Pass this read instead.", flush = True)
        continue
    seq1 = get_sequence(read_set_R1)
    seq2 = get_sequence(read_set_R2)
    _, overlap_length, overlap_seq = merge_sequences_paired_end_reads(seq1, seq2, min_overlap=0)

    list_overlap_length.append(overlap_length)

# %%
import matplotlib.pyplot as plt
from collections import Counter

cnt_overlap_length = Counter(list_overlap_length)

plt.figure(figsize=(20, 10))
plt.bar(cnt_overlap_length.keys(), cnt_overlap_length.values(), color="gray")
for key, value in cnt_overlap_length.items():
    plt.annotate(text=str(value), xy=(key, value), ha='center', va='bottom', fontsize=8, fontweight="bold")
plt.margins(x=0.01)
plt.yscale("log")
plt.yticks(fontsize=13)
plt.xticks(fontsize=13)
plt.xlabel("Length of Overlaps (Between Paired-end Reads)", fontsize=15)
plt.ylabel("Number of Overlaps", fontsize=15)
plt.show()
plt.close()


# %%
