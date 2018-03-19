import pickle
import os
from copy import deepcopy as copy

def genParamDict(sim_num, batch_num, strains, N_pop, T_0, T_sim, dt, t_extra,
                 N_sims, mut_param, date, code, filepath):
    '''A function that is used to generate a dictionary of parameters, given a
    set of inputs. Used to pass parameters to fit_traces.py and derivative
    scripts via pickle file

    Parameters
    ----------
    - int sim_num: identity of simulation batch (group of replicates of a
        simulation with the same parameter values)
    - int batch_num: identity within a simulation batch (which replicate)
    - list<strain> strains: a list of initial strains for simulation
    - float N_pop: the ptotal population <Up>size of the simulation
    - float T_0: the initial time of the simulation (in generations)
    - float T_sim: the total simulation time to run (in generations)
    - float dt: the time over which to run each epoch (in generations)
    - float t_extra: the amount of time over T_sim to extend the simulation(in
        generations); used to ensure clean interpolation in range [T_0, T_sim]
    - list<float> mut_param: a list of parameters to pass to mutation module
    - str date: the current date (DD-MM-YYYY)
    - str code: denotes genetic code(s) to be simulated (for I/O filenaming)
    - str filepath: the path to output file directory (for storing simulation
        results)

    Returns
    -------
    dict params: packaged version of what was passed in
    '''
    return locals()

def batcher(params, num_cores):
    ''' A function used to generate a batch of simulation parameter
    dictionaries in order to parallelize computation on AWS

    Parameters
    ----------
    - dict params: one dictionary of params to generate copies of
    - num_cores: the number of cores available in the cluster

    Returns
    -------
    list<dict> paramDicts: a list of parameter dictionaries
    '''
    # initialize and poplate list of parameter dictionaries
    paramDicts = []
    # calculate number of replicates per core
    N_sims = params['N_sims']
    N_per = int(N_sims/num_cores)
    # handle case where num_cores > N_sims
    if N_per == 0:
        for batch_num in range(N_sims):
            # perform one replicate per core
            paramCopy = copy(params)
            paramCopy['N_sims'] = 1
            paramCopy['batch_num'] = batch_num
            paramDicts.append(paramCopy)
    # handle case where num_cores < N_sims
    else:
        modulo = N_sims % num_cores
        for batch_num in range(num_cores):
            # adjust n_sim for each batch iteration
            n_sim = N_per
            if (num_cores - batch_num) <= modulo: n_sim += 1
            # package adjusted parameters
            paramCopy = copy(params)
            paramCopy['N_sims'] = n_sim
            paramCopy['batch_num'] = batch_num
            paramDicts.append(paramCopy)

    return paramDicts

def paramPickler(paramDicts, path='./', filenames=None):
    '''A function used to automagically pickle parameter dictionaries.

    Parameters
    ----------
    - list<dict> paramDict: parameter dictionaries to pickle
    - str path: the path to the pickling directory
    - list<str> filename: list of filenames for each paramDict. Generates names
        if none are given

    Returns
    -------
    None
    '''
    # create directories to pickle path as needed
    os.makedirs(pickle_path, exist_ok=True)
    # if only one dictionary is passed, package into a list
    if type(paramDicts) == dict:
        paramDicts = [paramDicts]
    # handle filename autogeneration cases
    if (filenames==None) or not (len(paramDicts)==len(filenames)):
        filenames = []
        for params in paramDicts:
            # extract parameters and format filename for each dict
            date = params['date']
            code = params['code']
            sim_num = params['sim_num']
            batch_num = params['batch_num']
            pickle_file = '{0}_{1}_{2}_{3}_params.pickle'.format(date, code, sim_num, batch_num)
            filenames.append(pickle_file)
    # loop over dictionaries and filenames to pickle
    for (params, filename) in zip(paramDicts, filenames):
        with open(path + filename, 'wb') as handle:
            pickle.dump(params, handle)
    return

def paramUpload(from_path, bucket, s3_upload_dir, s3_region='us-west-1'):
    ''' A function used to recursively upload files from local directory to a specified location on s3.

    Parameters
    ----------
    - str from_path: path to pickle files to upload
    - str bucket: name of AWS bucket to upload to
    - str s3_upload_dir: path to write param files to
    - str s3_region: name of AWS region where bucket lives

    Returns
    -------
    None
    '''
    # prepare boto3 client
    s3 = boto3.resource('s3', region_name=s3_region)
    # recursively walk through files to upload
    for root, dirs, files in os.walk(from_path):
        for filename in files:
            # get local and remote paths
            local_path = os.path.join(root, filename)
            outfile = s3_upload_dir + 'data/' + filepath + 'params/' + pickle_file
