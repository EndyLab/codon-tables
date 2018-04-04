import logging
logging.basicConfig(level=logging.INFO)

from datetime import date
from src.parameterizer import *
import pickle
import os

# import environmental variables from .config_replicates.py
from config_competition import *

# define local and remote paths
pickle_path = PATH + 'data/' + S3_HOME_DIR + 'params/' # SAME AS S3 PATH
s3_upload_path = S3_UPLOAD_DIR + FILEPATH + 'params/'

# generate base parameter dictionary
logging.info("Generating Base Parameter Dictionary")
params = genParamDict(
    SIM_NUM, BATCH_NUM, STRAINS, N_POP,
    T_0, T_SIM, DT, T_EXTRA, N_SIMS, MUT_PARAM,
    DATE, CODE, FILEPATH
)

# prepare set of parameter dictionaries
logging.info("Generating Parameter Batch")
paramDicts = batcher(params, NUM_CORES)

# inform user of the load sharing scheme
for params in paramDicts:
    batch_num = params['batch_num']
    n_sim = params['N_sims']
    print('Batch Num: {0}; n_sims = {1}'.format(batch_num, n_sim))

# store parameter files in pickle format
logging.info("Saving Parameter Files Locally to {0}".format(pickle_path))
paramPickler(paramDicts, pickle_path)

# upload files to s3
logging.info("Uploading Parameter Directory {0} to {1}:{2}".format(
    pickle_path, BUCKETNAME, S3_UPLOAD_DIR
    )
)
paramUpload(pickle_path, BUCKETNAME, s3_upload_path, S3_REGION)
success_string = (
    "Success! Check 'https://s3.console.aws.amazon.com/s3/home?region={0}'"
    + " to see parameter files."
).format(S3_REGION)
logging.info(success_string)
