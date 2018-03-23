from numpy import log10
from datetime import date
from src.bacteria import strain
from src.codonTable import codonTable

# absolute path to codon-table directory
PATH='/home/jonathan/Lab/ATD/codon-tables/'

# population variables
N_POP = 5e3
SIM_NUM = 0
BATCH_NUM = 0
STRAINS = [strain(N_pop=N_POP)]
T_0 = 0
T_SIM = 100
DT = 0.1
T_EXTRA = 5
N_SIMS = 10
MUT_PARAM = [1, 2]
DATE = str(date.today())
CODE = 'Standard Code'

# s3 variables
BUCKETNAME = 'endylab-codon-table-simulations'
S3_UPLOAD_DIR = 'manuscript-test/N={0}e{1}_b={2}_l={3}/'.format(
    str(N_POP)[0],
    int(log10(N_POP)),
    MUT_PARAM[0],
    MUT_PARAM[1]
)
FILEPATH = 'manuscript-test-1/'# this is the name for the local directory for this sim
S3_HOME_DIR = S3_UPLOAD_DIR + FILEPATH
S3_REGION = 'us-west-1'
FILENAME = '{0}_{1}_{2}_'.format(DATE, CODE, SIM_NUM, BATCH_NUM) + '{0}_params.pickle'

# AWS Batch variables
NUM_CORES = 32
RAM = 500 # in MB
ATTEMPTS = 1 # per instance
JOBNAME = 'python-generated-job'
JOBDEFINITION = 'arn:aws:batch:us-west-1:508804160265:job-definition/codon-sim-replicates:1'
JOBQUEUE = 'arn:aws:batch:us-west-1:508804160265:job-queue/codon-simulation-jobqueue'
DEPENDS_ON  = []
PARAMETERS = {}
