# %%
import argparse
import glob
import os
import subprocess

import numpy as np
import pandas as pd
from Bio import Entrez
from joblib import Parallel, delayed
import xml.etree.ElementTree as ET

# %%
def fetch_run_info(sra_ids):
    ids_string = ",".join(sra_ids)
    handle = Entrez.efetch(db="sra", id=ids_string, rettype="runinfo", retmode="text")
    run_info = handle.read()
    handle.close()
    
    return run_info

def query_sra_metadata(query, outmeta, email="micolak0115@gmail.com", retmax=10000):
    if not os.path.exists(outmeta):
        os.makedirs(outmeta,exist_ok=True)
        
    Entrez.email = email
    search_handle = Entrez.esearch(db="sra", term=query, retmax=retmax)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    sra_ids = search_results["IdList"]
    if sra_ids:
        all_run_info = ""
        chunk_size = 500
        for i in range(0, len(sra_ids), chunk_size):
            chunk_ids = sra_ids[i:i+chunk_size]
            run_info = fetch_run_info(chunk_ids)
            all_run_info += run_info.decode("utf-8")

        record = "\t".join(all_run_info.split(","))+"\n"
        
    sra_metdata = f"{outmeta}/sra_metadata_run_info_{query}.txt"
    
    with open(sra_metdata, "w") as fw:
        fw.write(record)

def attach_tissue_info_sra_metadata_from_biosample_id(file_metadata, outdir, query, email="micolak0115@gmail.com"):
    df_sra_metadata = pd.read_csv(file_metadata, sep="\t")

    Entrez.email = email
    df_sra_metadata["BioSample"].to_list()
    list_biosample_id = list(map(lambda x: x.replace("SAMN", ""), df_sra_metadata["BioSample"].to_list()))
    
    with open(os.path.join(outdir, f"tissue_info_{query}.txt"), mode="w") as fw:
        fw.write("BioSample_ID" + "\t" + "Tissue" + "\n")
        for biosample_id in list_biosample_id:
            handle = Entrez.esummary(db="biosample", id=biosample_id)
            record = handle.read()
            handle.close()

            root = ET.fromstring(record)
            for docsum in root.findall(".//DocumentSummary"):
                for item in docsum:
                    if str(item.text).startswith("<BioSample"):
                        second_root = ET.fromstring(item.text)
                        attributes = second_root.find('Attributes')
                        for attribute in attributes:
                            if (attribute.attrib["attribute_name"]) == "tissue":
                                tissue = attribute.text
                                fw.write(str(biosample_id) + "\t" +  str(tissue) + "\n")


def get_sra_metadata_filt_cond(file_metadata, filt_col="LibrarySource", filt_str="TRANSCRIPT"):
    df_sra = pd.read_csv(file_metadata, sep="\t")
    cond_filt = df_sra[filt_col].str.upper().str.startswith(filt_str)
    df_sra_filt = df_sra[cond_filt]
    
    return df_sra_filt

def get_list_sra_run_id(df_sra_filt, colrunid="Run"):
    list_srr_ids = df_sra_filt[colrunid].to_list()
    
    return list_srr_ids

def prefetch_sra_data(srr_id, dir_outsra, dir_sratool):
    dir_sra_data = os.path.join(dir_outsra, srr_id)
    if not os.path.exists(dir_sra_data):
        os.makedirs(dir_sra_data, exist_ok=True)
    cmd = f"{dir_sratool}/prefetch --output-directory {dir_outsra} {srr_id} --max-size 1t >> {dir_outsra}/prefetch.log"
    proc = subprocess.run(cmd, shell=True, check=True)
    if proc.returncode == 0:
        print(f"Successfully prefetched {srr_id}")
    else:
        print(f"Error prefetching {srr_id}")

def validate_integrity_sra_data(srr_id, dir_outsra, dir_sratool):
    dir_sra_data = os.path.join(dir_outsra, srr_id)

    if os.path.exists(f"{dir_sra_data}/{srr_id}.sra"):
        cmd = f"{dir_sratool}/vdb-validate --verbose {dir_sra_data}/{srr_id}.sra" 
    elif os.path.exists(f"{dir_sra_data}/{srr_id}.sralite"):
        cmd = f"{dir_sratool}/vdb-validate --verbose {dir_sra_data}/{srr_id}.sralite" 
    else:
        with open(f"{dir_outsra}/error.log", mode="a") as fw:
            fw.write(f"{srr_id}\n")
        cmd = f"echo {srr_id} skipped..."

    stdout = f"{dir_sra_data}/{srr_id}_vdb-validate.log"
    file_write_stderr = open(stdout, 'w')
    proc = subprocess.run(cmd, shell=True, check=True, stderr=file_write_stderr)

    if proc.returncode == 0:
        print(f"Successfully validated {srr_id}")
    else:
        print(f"Error validating {srr_id}")
    file_write_stderr.close() 

def dump_fastq_from_sra_data(srr_id, dir_outsra, dir_sratool):
    dir_sra_data = os.path.join(dir_outsra, srr_id)
    if os.path.exists(f"{dir_sra_data}/{srr_id}.sra"):
        cmd = f"{dir_sratool}/fasterq-dump --split-3 --skip-technical --outdir {dir_sra_data} {dir_sra_data}/{srr_id}.sra"
    elif os.path.exists(f"{dir_sra_data}/{srr_id}.sralite"):
        cmd = f"{dir_sratool}/fasterq-dump --split-3 --skip-technical --outdir {dir_sra_data} {dir_sra_data}/{srr_id}.sralite"
    else:
        with open(f"{dir_outsra}/error.log", mode="a") as fw:
            fw.write(f"{srr_id}\n")
        cmd = f"echo {srr_id} skipped..."

    proc = subprocess.run(cmd, shell=True, check=True)
    if proc.returncode == 0:
        print(f"Successfully dumped fastq {srr_id}")
    else:
        print(f"Error dumping fastq {srr_id}")

def run_sratool_parallel(func_tool, list_srr_ids, outdir, dir_sratool, n_jobs):
    if n_jobs is None:
        n_jobs = len(list_srr_ids)
    finish_logs = Parallel(n_jobs=n_jobs)(delayed(func_tool)(srr_id, outdir, dir_sratool) for srr_id in list_srr_ids)
    for finish in finish_logs:
        print(finish)

def run_sratool_parallel_chunks(func_tool, list_srr_ids, outdir, dir_sratool, chunk_size, n_jobs=None):
    for i in range(0, len(list_srr_ids), chunk_size):
        list_srr_ids_chunks = list_srr_ids[i: i+chunk_size]
        run_sratool_parallel(func_tool, list_srr_ids_chunks, outdir, dir_sratool, n_jobs)

def parse_validate_log(dir_outsra):
    list_validate_log = glob.glob(f"{dir_outsra}/**/*_vdb-validate.log", recursive=True)
    for validate_log in list_validate_log:
        set_status = set()
        with open(validate_log, mode="r") as fr:
            for line in fr:
                record = line.rstrip("\n").split()
                status = record[-1]
                set_status.add(status)
        if set_status == set():
            with open(validate_log, mode="r") as fr:
                for line in fr:
                    record = line.rstrip("\n").split()
                    print(record)
        else:
            continue

def main(args):
    file_metadata = os.path.join(args.dir_outsra, f"sra_metadata_run_info_{args.query}.txt")
    if not os.path.exists(file_metadata):
        query_sra_metadata(args.query, args.dir_outsra)
    attach_tissue_info_sra_metadata_from_biosample_id(file_metadata, args.dir_outsra, args.query)
    df_sra_filt = get_sra_metadata_filt_cond(file_metadata, filt_col=args.filt_col, filt_str=args.filt_str)

    list_srr_ids = get_list_sra_run_id(df_sra_filt)
    if args.mode=="prefetch":
        run_sratool_parallel_chunks(prefetch_sra_data, list_srr_ids, args.dir_outsra, args.dir_sratool, args.chunk_size)
    elif args.mode=="validate":
        run_sratool_parallel_chunks(validate_integrity_sra_data, list_srr_ids, args.dir_outsra, args.dir_sratool, args.chunk_size)
    elif args.mode=="fastqdump":
        run_sratool_parallel_chunks(dump_fastq_from_sra_data, list_srr_ids, args.dir_outsra, args.dir_sratool, args.chunk_size)
    elif args.mode=="none":
        pass
    else:
        print("Unsupported mode")

# %%
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-q", "--query", type=str, default="SRP214077")
    parser.add_argument("-do", "--dir_outsra", type=str, default="/BiO/Research/GeroExpressome/Resources/Data/SRA/SRP214077")
    parser.add_argument("-ds", "--dir_sratool", type=str, default="/BiO/Research/GeroPathway/Resources/Tools/sratoolkit.3.1.1-ubuntu64/bin")
    parser.add_argument("-fc", "--filt_col", type=str, default="LibrarySource")
    parser.add_argument("-fs", "--filt_str", type=str, default="TRANSCRIPT")
    parser.add_argument("-cz", "--chunk_size", type=int, default=10)
    parser.add_argument("-m", "--mode", type=str, default="none")
    args = parser.parse_args()
    
    main(args)
# %%
