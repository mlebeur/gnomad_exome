import os, pandas, csv, re
import numpy as np
import hashlib
from biothings.utils.dataload import dict_convert, dict_sweep
from biothings import config
logging = config.logger
def load_gnomad_exome(data_folder):
    infile = os.path.abspath("/opt/biothings/GRCh37/gnomAD_exomes/r2.1/GnomadExomes.tsv")
    assert os.path.exists(infile)
    results = {}
    with open(infile) as fp:
        reader = csv.reader(fp)
        header = reader.next()
        for line in fp:
            rec = dict(zip(header,line))
            var = rec["release"] + "_" + str(rec["chromosome"]) + "_" + str(rec["position"]) + "_" + rec["reference"] + "_" + rec["alternative"]       
            _id = hashlib.sha224(var.encode('ascii')).hexdigest()       
            process_key = lambda k: k.replace(" ","_").lower()
            rec = dict_convert(rec,keyfn=process_key)
            rec = dict_sweep(rec,vals=[np.nan])
            results.setdefault(_id,[]).append(rec)
            doc = {"_id": _id, "gnomad_exome" : docs}
            yield doc

        
