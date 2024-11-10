# %%
#!/usr/bin/env/python3

import gzip
import re
from collections import Counter
from itertools import chain

class File:
    def __init__(self, path_input_file, sep, is_gzip, header_startswith):
        self.path_input_file = path_input_file
        self.sep = sep
        self.is_gzip = is_gzip
        self.header_startswith = header_startswith
        self.file_reader = None
        self.set_file_reader()

    def set_file_reader(self):
        if self.path_input_file.endswith(".gz"):
            self.is_gzip = True
            self.file_reader = gzip.open(self.path_input_file, mode="rb")
        else:
            self.file_reader = open(self.path_input_file, mode="r")

    def get_nextline(self, num_shorten=None):
        line = next(self.file_reader, None)
        if line == None: return None

        if self.is_gzip:
            line = line.decode()

        if num_shorten:
            line = line[:num_shorten]

        if self.sep:
            record = line.rstrip("\n").split(self.sep)
        
        else:
            record = line.rstrip("\n")
        
        return record
    
    def make_practice_file(self, file_practice, line_count, file_obj):
        with open(file_practice, 'w') as fw:
            if hasattr(file_obj, "is_VCF"):
                fw.write("\t".join(self.header)+"\n")
                for _ in range(line_count):
                    line_practice_file = self.get_nextline()
                    fw.write("\t".join(line_practice_file) + "\n")
            
            else:
                for _ in range(line_count):
                    line_practice_file = self.get_nextline()
                    fw.write(line_practice_file + "\n")

        print("Done writing file!")

    def skip_header(self):
        while True:
            line = self.get_nextline()
            if line is None:
                raise EOFError("Reached end of file without finding header")
            
            if line[0].startswith(self.header_startswith):
                break
        
        header = line

        return header

    def get_index_query_column_name(self, header, colname):
        if colname in header:
            idx_num = header.index(colname)

            return idx_num
    
    def close_reader(self):
        self.file_reader.close()

    def __del__(self):
        if not self.file_reader.closed:
            self.close_reader()
        del self.file_reader

class VCF(File):
    def __init__(self, path_input_file, sep="\t", is_gzip=True, header_startswith="#CHROM"):
        super().__init__(path_input_file, sep, is_gzip, header_startswith)
        self.header = self.skip_header()
        self.format_index = self.get_index_format()
        self.info_index = self.get_index_info()
        self.samples_start_index = self.get_index_samples_start()

    def is_VCF(self):
        pass

    def get_index_format(self, namecolformat="FORMAT"):
        format_index = int(self.header.index(namecolformat))
        
        return format_index
        
    def get_index_info(self, namecolinfo="INFO"):
        info_index = int(self.header.index(namecolinfo))

        return info_index

    def get_index_samples_start(self):
        samples_start_index = int(self.format_index) + 1

        return samples_start_index

    def get_variant_name(self, vcf_record):
        chr = vcf_record[0]
        pos = vcf_record[1] 
        ref = vcf_record[3]
        alt = vcf_record[4]
        variant = "_".join([chr, pos, ref, alt])

        return variant

    def _get_format_vcf(self, vcf_record):
        vcf_format = vcf_record[self.format_index]

        return vcf_format

    def _get_info_vcf(self, vcf_record):
        vcf_info = vcf_record[self.info_index]

        return vcf_info

    def _get_index_format_attribute(self, vcf_record, nameattr):
        format_vcf = self._get_format_vcf(vcf_record)
        attr_index = format_vcf.split(":").index(nameattr)

        return attr_index

    def get_sample_name(self):
        sample_name = self.header[self.samples_start_index:]

        return sample_name

    def _get_sample_level_information(self, vcf_record):
        vcf_sample_info = vcf_record[self.samples_start_index:]

        return vcf_sample_info

    def get_variant_attr(self, vcf_record, nameattr=None):
        vcf_sample_info = self._get_sample_level_information(vcf_record)
        if nameattr is None:
            variant_attr = list(map(lambda x: x.split(":"), vcf_sample_info))
        else:
            attr_index = self._get_index_format_attribute(vcf_record, nameattr)
            variant_attr = list(map(lambda x: x.split(":")[attr_index], vcf_sample_info))
        
        return variant_attr

    def get_variant_info(self, vcf_record):
        info_vals = self._get_info_vcf(vcf_record).split(";")
        dict_info_variant = dict(list((map(lambda x: (x.split("=")[0], x.split("=")[-1]), info_vals))))

        return dict_info_variant

    def change_info_vcf(self, dict_info_variant, new_info_val, info_attr="AC"):
        dict_info_variant[info_attr] = new_info_val

        return dict_info_variant

    def update_variant_info_attr(self, dict_info_variant, info_attr, attr_val):
        dict_info_attr_to_update = {info_attr: attr_val}
        dict_info_variant.update(dict_info_attr_to_update)
        dict_info_variant_update = dict_info_variant

        return dict_info_variant_update

    def get_sample_alleles(self, vcf_variant_gt):
        compile_gt = re.compile("([\.,0-9])[/,\|]([\.,0-9])")
        list_gt_groups = list(map(lambda x: compile_gt.match(x).groups(), vcf_variant_gt))
        
        return list_gt_groups

    def calc_allele_stats(self, vcf_variant_gt):
        list_gt_groups = self.get_sample_alleles(vcf_variant_gt)
        list_all_gt = list(chain(*list_gt_groups))
        dict_cnt_gt = dict(Counter(list_all_gt))

        return dict_cnt_gt

    def calc_allele_count(self, dict_cnt_gt, flag_allele="1", verbose=False):
        try:
            allele_count = dict_cnt_gt[flag_allele]
        except:
            if verbose:
                print(f"no alt alleles: {dict_cnt_gt}")
            allele_count = 0

        return allele_count

    def calc_allele_freq(self, dict_cnt_gt, flag_allele="1", verbose=False):
        try:
            allele_freq = int(dict_cnt_gt[flag_allele])/int(sum(list(dict_cnt_gt.values())))
            allele_freq = round(float(allele_freq), 7)
        except:
            if verbose:
                print(f"no alt alleles: {dict_cnt_gt}")
            allele_freq = 0            

        return allele_freq      

    def calc_allele_number(self, dict_cnt_gt, remove_allele="."):
        dict_cnt_gt.pop(remove_allele, None)
        allele_num = sum(list(dict_cnt_gt.values()))
        
        return allele_num

    def is_homo_ref(self, gt_group):
        if str(list(set([*gt_group]))[0]) == "0":
            return True
        
        return False

    def is_homo_alt(self, gt_group):
        if str(list(set([*gt_group]))[0]) == "1":
            return True
        
        return False

    def is_hetero(self, gt_group):
        num_alleles = len(dict(Counter([*gt_group])).keys())
        if num_alleles > 1:
            return True

        return False

    def calc_allele_balance(self, record, mode="hetero", filt="on"):
        def calc_ab(ad):
            a1, a2 = ad.split(',')
            a1 = int(a1)
            a2 = int(a2)
            tot = a1+a2

            if tot == 0:
                return None
            
            return a2/(a1+a2)

        vcf_variant_ad = self.get_variant_attr(record, nameattr="AD")
        vcf_variant_gt = self.get_variant_attr(record, nameattr="GT")
        vcf_variant_dp = self.get_variant_attr(record, nameattr="DP")
        vcf_variant_gq = self.get_variant_attr(record, nameattr="GQ")
        gt_groups = self.get_sample_alleles(vcf_variant_gt)

        if filt == "on":
            for dp in vcf_variant_dp:
                if str(dp) == "." or int(dp) <= 10: 
                    continue
                    
            for gq in vcf_variant_gq:
                if str(gq) == "." or int(gq) <= 20:
                    continue

        list_ad = list()
        if mode=="homo_ref":
            for ad, gt_group in zip(vcf_variant_ad, gt_groups):
                if gt_group == ("0", "0"):
                    list_ad.append(ad)
                        
        elif mode=="hetero":
            for ad, gt_group in zip(vcf_variant_ad, gt_groups):
                if gt_group == ("1", "0") or gt_group == ("0", "1"):
                    list_ad.append(ad)        

        elif mode =="homo_alt":
            for ad, gt_group in zip(vcf_variant_ad, gt_groups):
                if gt_group == ("1", "1"):
                    list_ad.append(ad)
            
        list_ab_filt = list(set(list_ad))
        list_ab = list(map(lambda x: calc_ab(x), list_ad))
        list_ab_filt = list(filter(lambda x: x is not None, list_ab))

        return list_ab_filt
        
    def calc_call_rate(self, record):
        vcf_variant_dp = self.get_variant_attr(record, nameattr="DP")
        # depth = 0 if depth = "."
        vcf_variant_dp = list(map(lambda x: 0 if x=="." else x, vcf_variant_dp))
        calls = list(filter(lambda x: int(x) != 0, vcf_variant_dp))
        call_rate = round(len(calls)/len(vcf_variant_dp), 3)
        # if verbose:
        #     print(f"call rate = {len(calls)}/{len(vcf_variant_dp)} ({call_rate}%)")

        return float(call_rate)

    def get_samples_variant_attr_after_filt_cond(self, sample_name, variant_attr, filt_cond, filt_num):
        filt_in_samples = list()
        for name, attr in zip(sample_name, variant_attr):
            try:
                attr = int(attr)
            except:
                attr = -1

            if filt_cond == "eq":
                if attr == filt_num:
                    filt_in_samples.append(name)
                else:
                    continue
            elif filt_cond == "ge":
                if attr > filt_num:
                    filt_in_samples.append(name)
                else:
                    continue
            elif filt_cond == "le":
                if attr < filt_num:
                    filt_in_samples.append(name)
                else:
                    continue
            else:
                continue
        
        return filt_in_samples

    def filter_variant_attr(self, vcf_record, sample_name, nameattr="GQ", filt_cond="ge", filt_num=0, filt_mode="out", verbose=True):
        variant_attr = self.get_variant_attr(vcf_record, nameattr)
        filt_in_samples = self.get_samples_variant_attr_after_filt_cond(sample_name, variant_attr, filt_cond, filt_num)

        if verbose:
            print(f"VCF format attribute: {nameattr}")
            dict_ineq_sym = {"ge":">=", "le":"<=", "eq":"="}
            print(f"filter condition: {dict_ineq_sym[filt_cond]}{filt_num}")
            print(f"filter mode: {filt_mode}")
            tot_num = len(sample_name)
            num_filt_in = len(filt_in_samples)
            num_filt_out = tot_num - num_filt_in
            dict_filt = {"filt_in":num_filt_in, "filt_out":num_filt_out}
            print(f"filter stats: {dict_filt}")
        
        if filt_mode == "in":

            return sorted(filt_in_samples)
        
        else:
            filt_out_samples = list(set(sample_name).difference(set(filt_in_samples)))
            
            return sorted(filt_out_samples)

    def designate_missing_genotype(self, sample_name, vcf_variant_gt, samples_miss_gt):
        for idx, (name, gt) in enumerate(zip(sample_name, vcf_variant_gt)):
            if name in set(samples_miss_gt):
                vcf_variant_gt[idx] = "./."

        vcf_variant_gt_adj_miss_gt = vcf_variant_gt

        return vcf_variant_gt_adj_miss_gt

class Fastq(File):
    def __init__(self, path_input_file, sep=None, is_gzip=True, header_startswith=None):
        super().__init__(path_input_file, sep, is_gzip, header_startswith)
        self.instrument_id = {
                    "HWI-M[0-9]+" : "MiSeq",
                    "HWUSI-[A-z,0-9]+" : "Genome Analyzer IIx",
                    "M[0-9]+" : "MiSeq",
                    "HWI-C[0-9]+" : "HiSeq 1500",
                    "C[0-9]+" : "HiSeq 1500",
                    "HWI-D[0-9]+" : "HiSeq 2000/2500",
                    "HWI-ST[0-9]+" : "HiSeq 2000",
                    "SN[0-9+]" : "HiSeq 2500",
                    "D[0-9]+" : "HiSeq 2500",
                    "J[0-9]+" : "HiSeq 3000",
                    "K[0-9]+" : "HiSeq 3000/4000",
                    "E[0-9]+" : "HiSeq X",
                    "NB[0-9]+": "NextSeq",
                    "NS[0-9]+" : "NextSeq",
                    "N[0-9]+" : "NextSeq 500/550",
                    "V[0-9]+" : "NextSeq 2000",
                    "AA[0-9]" : "NextSeq 2000 P1/P2/P3 or Genome Analyzer",
                    "MN[0-9]+" : "MiniSeq",
                    "A[0-9]+" : "NovaSeq 6000",
                    "H[0-9]+" : "NovaSeq S1/S2/S4",
                    "@([A-z,0-9]+)L([0-9]+)C([0-9]{3})R([0-9]+)/([1,2])":"BGI/DNBSeq"
        }
        self.quality_score = {
                    '!': '0',
                    '"': '1',
                    '#': '2',
                    '$': '3',
                    '%': '4',
                    '&': '5',
                    "'": '6',
                    '(': '7',
                    ')': '8',
                    '*': '9',
                    '+': '10',
                    ',': '11',
                    '-': '12',
                    '.': '13',
                    '/': '14',
                    '0': '15',
                    '1': '16',
                    '2': '17',
                    '3': '18',
                    '4': '19',
                    '5': '20',
                    '6': '21',
                    '7': '22',
                    '8': '23',
                    '9': '24',
                    ':': '25',
                    ';': '26',
                    '<': '27',
                    '=': '28',
                    '>': '29',
                    '?': '30',
                    '@': '31',
                    'A': '32',
                    'B': '33',
                    'C': '34',
                    'D': '35',
                    'E': '36',
                    'F': '37',
                    'G': '38',
                    'H': '39',
                    'I': '40'
        }

    def is_Fastq(self):
        pass

    def get_header(self):
        pass
        
    def parse_sequences(self):
        pass

class BAM(File):
    pass

    def is_BAM(self):
        pass

class GTF(File):
    def __init__(self, path_input_file, sep="\t", is_gzip=False, header_startswith="##date"):
        super().__init__(path_input_file, sep, is_gzip, header_startswith)
        _ = self.skip_header()
        self.header = ["chrname", "source", "feature", "start", "end", "score", "strand", "phase", "attribute"]
    
    def is_GTF(self):
        pass

    def get_dictionary_attribute(self, line):
        idx_attr = self.get_index_query_column_name(self.header, "attribute")
        attribute = line[idx_attr]
        list_attributes = list(map(lambda x: x.lstrip().rstrip(), attribute.split(";")))
        list_attr_key = list(map(lambda x: x.split(" ")[0], list_attributes))
        list_attr_val = list(map(lambda x: x.split(" ")[-1], list_attributes))
        dict_attr = dict(zip(list_attr_key, list_attr_val))
        dict_attr = {k:v.lstrip('"').rstrip('"') for k, v in dict_attr.items() if k != ""}

        return dict_attr
    
class GCT(File):
    def __init__(self, path_input_file, sep="\t", is_gzip=False, header_startswith="Name"):
        super().__init__(path_input_file, sep, is_gzip, header_startswith)
        self.header = self.skip_header()

    def is_GCT(self):
        pass

    def get_dictionary_sample_id_tissue(self, path_metadata_file, colsampleid="SAMPID", coltissue="SMTSD"):
        dict_sample_id_tissue = dict()
        with open(path_metadata_file, mode="r") as fr:
            header = fr.readline().rstrip("\n").split("\t")
            idx_sample_id = self.get_index_query_column_name(header, colname=colsampleid)
            idx_tissue = self.get_index_query_column_name(header, colname=coltissue)
            for line in fr:
                record = line.rstrip("\n").split("\t")
                sample_id = record[idx_sample_id]
                tissue = record[idx_tissue]
                dict_sample_id_tissue[sample_id] = tissue

        return dict_sample_id_tissue

    def get_list_sample_id_only_specific_tissue(self, dict_sample_id_tissue, query_tissue="Whole Blood"):
        list_sample_id_specific_tissue = list()
        for sample_id, sample_tissue in dict_sample_id_tissue.items():
            if sample_tissue == query_tissue:
                list_sample_id_specific_tissue.append(sample_id)
        
        return list_sample_id_specific_tissue
    
    def get_list_idx_query_sample_id(self, list_sample_id_query):
        list_idx_sample_id_whole_blood = list(map(lambda colname: self.get_index_query_column_name(self.header, colname) if colname in self.header else None, list_sample_id_query))
        list_idx_line_query = list(filter(lambda idx: idx!=None, self.list_default_idx + list_idx_sample_id_whole_blood))

        return list_idx_line_query

    def filter_line_by_query(self, line, list_idx_line_query):
        line_filtered = list(map(lambda idx: line[idx], list_idx_line_query))

        return line_filtered

    def select_gene_count_table(self, file_gct_filtered, list_sample_id_query):
        list_idx_line_query = self.get_list_idx_query_sample_id(list_sample_id_query)
        with open(file_gct_filtered, mode="w") as fw:
            header_filtered = self.filter_line_by_query(self.header, list_idx_line_query)
            fw.write("\t".join(header_filtered) + "\n")

            while 1:
                line = self.get_nextline()
                if line == None:
                    break

                line_filtered = self.filter_line_by_query(line, list_idx_line_query)
                fw.write("\t".join(line_filtered) + "\n")

class GEO(File):
    def __init__(self, path_input_file, sep="\t", is_gzip=False, header_startswith="!Sample_title"):
        super().__init__(path_input_file, sep, is_gzip, header_startswith)
        self.header = self.skip_header()

    def is_GEO(self):
        pass
    
    def get_sample_id(self):
        header = self.header
        list_samples = list(map(lambda x: x.lstrip('"').rstrip('"'), header[1:]))
        
        return list_samples
    
    def get_sample_size(self):
        list_samples = self.get_sample_id()
        sample_size = len(list_samples)

        return sample_size

    def get_metadata(self):
        import pandas as pd

        prev = []
        dict_meta = dict()
        dict_meta["Project-ID"] = self.get_sample_id()

        while 1:
            line = self.get_nextline()
            if "protocol" in line[0]:
                break
            
            colname = line[0].lstrip('!')
            values = list(map(lambda x: x.lstrip('"').rstrip('"'), line[1:]))
            
            if colname in prev:
                colname = line[1].split(":")[0].lstrip('"')
                values = list(map(lambda x: x.split(":")[-1].lstrip(" ").rstrip('"'), line[1:]))

            prev.append(colname)
            dict_meta[colname] = values

        df_meta = pd.DataFrame(dict_meta)

        return df_meta
    
class HGNC(File):
    
    def __init__(self, path_input_file, sep="\t", is_gzip=False, header_startswith="hgnc_id"):
        super().__init__(path_input_file, sep, is_gzip, header_startswith)
        self.header = self.skip_header()
        self.ensembl_gene = self.get_dictionary_ensembl_to_gene()
        self.entrez_gene = self.get_dictionary_entrez_to_gene()
        self.entrez_ensembl = self.get_dictionary_entrez_to_ensembl()

    def is_GEO(self):
        pass        

    def get_master_table(self):
        import pandas as pd
        df_master_table = pd.read_csv(self.path_input_file, sep="\t")

        return df_master_table

    def get_dictionary_ensembl_to_gene(self, colensembl = "ensembl_gene_id", colgene = "symbol"):
        df_master_table = self.get_master_table()
        dict_ensembl_gene = dict(zip(df_master_table[colensembl], df_master_table[colgene]))

        return dict_ensembl_gene

    def get_dictionary_entrez_to_gene(self, colentrez = "entrez_id", colgene = "symbol"):
        df_master_table = self.get_master_table()
        dict_entrez_gene = dict(zip(df_master_table[colentrez], df_master_table[colgene]))

        return dict_entrez_gene
    
    def get_dictionary_entrez_to_ensembl(self, colentrez = "entrez_id", colensembl = "ensembl_gene_id"):
        df_master_table = self.get_master_table()
        dict_entrez_ensembl = dict(zip(df_master_table[colentrez], df_master_table[colensembl]))

        return dict_entrez_ensembl