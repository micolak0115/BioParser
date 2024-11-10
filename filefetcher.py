# %%
import os
import pandas as pd
from Bio import Entrez
import xml.etree.ElementTree as ET

# %%
class DB():
    def __init__(self, input_db):
        self.input_db = input_db

    def fetch_run_info(self, fetch_ids):
        ids_string = ",".join(fetch_ids)
        handle = Entrez.efetch(db=self.input_db, id=ids_string, rettype="runinfo", retmode="text")
        run_info = handle.read()
        handle.close()
        
        return run_info

class SRA(DB):
    def __init__(self, input_db="sra"):
        super().__init__(input_db)
        self.input_db = input_db
        self.email="micolak0115@gmail.com"
        self.retmax=10000

    def query_metadata(self, query, dir_outmeta):
        if not os.path.exists(dir_outmeta):
            os.makedirs(dir_outmeta,exist_ok=True)
            
        Entrez.email = self.email
        search_handle = Entrez.esearch(db=self.input_db, term=query, retmax=self.retmax)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        sra_ids = search_results["IdList"]
        if sra_ids:
            all_run_info = ""
            chunk_size = 500
            for i in range(0, len(sra_ids), chunk_size):
                chunk_ids = sra_ids[i:i+chunk_size]
                run_info = self.fetch_run_info(chunk_ids)
                all_run_info += run_info.decode("utf-8")

            record = "\t".join(all_run_info.split(","))+"\n"
        
        query_filename = "_".join(query.split())
        sra_metadata = f"{dir_outmeta}/sra_metadata_run_info_{query_filename}.txt"
        
        with open(sra_metadata, "w") as fw:
            fw.write(record)
        
        return sra_metadata
    
    def attach_tissue_info_sra_metadata_from_biosample_id(self, file_metadata, outdir, query):
        df_sra_metadata = pd.read_csv(file_metadata, sep="\t")

        Entrez.email = self.email
        df_sra_metadata["BioSample"].to_list()
        list_biosample_id = list(map(lambda x: x.replace("SAMN", ""), df_sra_metadata["BioSample"].to_list()))
        
        with open(os.path.join(outdir, f"tissue_info_{query}.txt"), mode="w") as fw:
            fw.write("BioSample_ID" + "\t" + "Tissue" + "\n")
            for biosample_id in list_biosample_id:
                handle = Entrez.esummary(db="biosample", id=biosample_id)
                record = handle.read()
                handle.close()

                # Parse the XML response to extract details
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

# %%
query = "SRP214077"
dir_outmeta = "/BiO/Research/GeroExpressome/Results"
sra = SRA()
file_metadata = sra.query_metadata(query, dir_outmeta)
df_sra_metadata = pd.read_csv(file_metadata, sep="\t")
# %%
