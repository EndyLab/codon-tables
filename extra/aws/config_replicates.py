from numpy import log10
from datetime import date
from src.bacteria import strain
from src.codonTable import codonTable
from src.codonUtils import utils

# absolute path to codon-table directory
PATH='/Users/jonathancalles/Lab/ATD/codon-tables/'

# population variables
N_POP = 1e6
SIM_NUM = 0
BATCH_NUM = 0
PROMISC20 = utils.promiscuity(utils.RED20)
STRAINS = [strain(N_pop=N_POP, table=PROMISC20, code='PROMISC20')]
T_0 = 0
T_SIM = 1000
DT = 0.1
T_EXTRA = 5
N_SIMS = 1000
MUT_PARAM = [1, 2]
DATE = str(date.today())
CODE = 'RED15'

# s3 variables
BUCKETNAME = 'endylab-codon-table-simulations'
S3_UPLOAD_DIR = 'manuscript/N={0}e{1}_b={2}_l={3}/'.format(
    str(N_POP)[0],
    int(log10(N_POP)),
    MUT_PARAM[0],
    MUT_PARAM[1]
)
FILEPATH = 'RED15_{0}/'.format(SIM_NUM)# this is the name for the local directory for this sim
S3_HOME_DIR = S3_UPLOAD_DIR + FILEPATH
S3_REGION = 'us-west-1'
FILENAME = '{0}_{1}_{2}_'.format(DATE, CODE, SIM_NUM, BATCH_NUM) + '{0}_params.pickle'

# AWS Batch variables
NUM_CORES = 256
RAM = 500 # in MB
ATTEMPTS = 1 # per instance
JOBNAME = '{0}-job'.format(CODE)
JOBDEFINITION = 'arn:aws:batch:us-west-1:508804160265:job-definition/codon-sim-replicates:1'
JOBQUEUE = 'arn:aws:batch:us-west-1:508804160265:job-queue/codon-simulation-jobqueue'
DEPENDS_ON  = []
PARAMETERS = {}
