#!/usr/bin/env python
# encoding: utf-8
"""
run.py - control a run on an aws node

Created by Dave Williams on 2010-11-19
"""

import sys
import os
import time
import uuid
import ujson as json
import pickle
import shelve
import copy
import optparse 
import multiprocessing as mp
import boto
import hs
import numpy.random as random
import numpy as np

## Settings
BUCKET_NAME = 'model_results' # For results files
FOLDER_NAME = 'test_disregard/' # Needs trailing /


## Help message
help_mesg = "Around here, we like to call our functions with a little \n \
      thing called options. Here's how we'd go about it: \n \
      usage: %prog [options] arg1 arg2 "

## Helper functions
def log_it(message):
    """Print message to sys.stdout"""
    sys.stdout.write("run.py " + mp.current_process().name + 
            " ## " + message + "\n")
    sys.stdout.flush()

## Run emission
def emit_run(path_local, path_s3, timestep_length, timestep_number, 
             z_line=None, lattice_spacing=None, actin_permissiveness=None,
             comment = None, write = True):
    """Produce a structured JSON file that will be consumed to create a run
    
    emit_run is intended as the counterpoint to the run class and execute_run.
    Import emit_run into an interactive workspace and populate a directory with
    run configurations to be executed by a cluster of machines. 
    
    
    Parameters
    ----------
    path_local: string
        The local (absolute or relative) directory to which we save both 
        emitted files and run output.
    path_s3: string
        The s3 bucket (and optional folder) to save run output to and to which
        the emitted files should be uploaded.
    timestep_length: float
        Length of a timestep in ms
    timestep_number: int
        Number of timesteps run is simulated for
    z_line: float or tuple, optional
        If not given, default distance specified in hs.hs is used. If given as
        float, the z-line distance for the run. If given as a tuple, used as 
        arguments to the current z-line generation function within. The tuple
        format is currently: (offset, amp, period)
    lattice_spacing: float or tuple, optional
        Same as for z-line. The tuple format is currently: (rest_ls, rest_zln)
    actin_permissiveness: float or tuple, optional
        Same as for z-line. The tuple format is currently: (amp, period, on, 
        duration, sharp)
    comment: string, optional
        Space for comment on the purpose or other characteristics of the run
    write: bool, optional
        True (default) writes file to path_local/name.meta.json. Other values 
        don't. In both cases the dictionary describing the run is returned.
    
    Returns
    -------
    rund: dict
        Copy of run dictionary saved to disk as json.
    
    Examples
    --------
    >>> emit_run('./', None, .1, 100, write=False) 
    {'actin_permissiveness': None,
    ...  'actin_permissiveness_func': None,
    ...  'comment': None,
    ...  'lattice_spacing': None,
    ...  'lattice_spacing_func': None,
    ...  'name': ...,
    ...  'path_local': './',
    ...  'path_s3': None,
    ...  'timestep_length': 0.1,
    ...  'timestep_number': 100,
    ...  'z_line': None,
    ...  'z_line_func': None}
    ...
    >>> emit_run('./eggs', 'EGGS', 5, 10, (1250, 25, 25), (16, 1250), .75, comment='Just a test', write=False)  
    {'actin_permissiveness': 0.75,
    ... 'actin_permissiveness_func': None,
    ... 'comment': 'Just a test',
    ... 'lattice_spacing': array([ 16.        ,  15.92445392,  15.95318346,  16.04723114,
    ... 16.07663156,  16.        ,  15.92445392,  15.95318346,
    ... 16.04723114,  16.07663156]),
    ... 'lattice_spacing_func': ("lambda rest_ls, rest_zln: np.sqrt(np.power(rest_ls,2) * rest_zln / rund['z_line'])",
    ... (16, 1250)),
    ... 'name': ...,
    ... 'path_local': './eggs',
    ... 'path_s3': 'EGGS',
    ... 'timestep_length': 5,
    ... 'timestep_number': 10,
    ... 'z_line': array([ 1250.        ,  1261.88820645,  1257.34731565,  1242.65268435,
    ... 1238.11179355,  1250.        ,  1261.88820645,  1257.34731565,
    ... 1242.65268435,  1238.11179355]),
    ... 'z_line_func': ('lambda offset, amp, period, time: offset + 0.5 * amp * np.sin(2*np.pi*time/period)',
    ... (1250, 25, 25))}
    ... 
    """
    rund = {}
    ## Simple run metadata
    name = str(uuid.uuid1())
    rund['name'] = name
    rund['comment'] = comment
    rund['path_local'] = path_local
    rund['path_s3'] = path_s3
    rund['timestep_length'] = timestep_length
    rund['timestep_number'] = timestep_number
    ## Generate and store variable traces
    time = np.arange(0, timestep_number*timestep_length, timestep_length)
    # Define variable sarcomere length/z-line distance
    def variable_z_line(offset, amp, period):
        """Takes offset from zero, peak-to-peak amp, and period in ms"""
        return offset + 0.5 * amp * np.sin(2*np.pi*time/period)
    string_zln = "lambda offset, amp, period, time: offset + 0.5 * amp * np.sin(2*np.pi*time/period)"
    # Parse sarcomere length
    if hasattr(z_line, "__iter__"):
        rund['z_line'] = variable_z_line(*z_line)
        rund['z_line_func'] = (string_zln, z_line)
    else:
        rund['z_line'] = z_line
        rund['z_line_func'] = None
   # Define variable lattice spacing
    def variable_lattice_spacing(rest_ls, rest_zln):
        """Assumes constant volume relation, takes in rest ls and zline"""
        return np.sqrt(np.power(rest_ls,2) * rest_zln / rund['z_line'])
    string_ls = "lambda rest_ls, rest_zln: np.sqrt(np.power(rest_ls,2) * rest_zln / rund['z_line'])"
    # Parse lattice spacing
    if hasattr(lattice_spacing, "__iter__"):
        rund['lattice_spacing'] = variable_lattice_spacing(*lattice_spacing)
        rund['lattice_spacing_func'] = (string_ls, lattice_spacing)
    else:
        rund['lattice_spacing'] = lattice_spacing
        rund['lattice_spacing_func'] = None
    # Define variable actin permissiveness
    def variable_actin_permissiveness(amp, period, on, duration, sharp):
        """Requires amplitude, period of stim cycle, offset of turn-on from
        start of a cycle, duration of on time, and sharpness (lower is more
        gentle) of transition"""
        sig = lambda on, x: \
                0.5*amp*(np.tanh(sharp*(x-on))-np.tanh(sharp*(x-on-duration)))
        cycle_number = int(time[-1])//period+1
        return np.sum([sig(on+period*i, time) for i in range(cycle_number)], 0)
    string_ap = """
        sigmoid = lambda amp, on, duration, sharp, x: 
            0.5*amp*(np.tanh(sharp*(x-on))-np.tanh(sharp*(x-on-duration)))
        bumps = lambda amp, on, duration, sharp, period, x: 
            np.sum([sigmoid(amp, on+period*i, duration, sharp, x) 
            for i in range(int(x[-1])//period+1)],0) """
    # Parse actin permissiveness 
    if hasattr(actin_permissiveness, "__iter__"):
        rund['actin_permissiveness'] = variable_actin_permissiveness(
            *actin_permissiveness)
        rund['actin_permissiveness_func'] = (string_ap, actin_permissiveness)
    else:
        rund['actin_permissiveness'] = actin_permissiveness
        rund['actin_permissiveness_func'] = None
    ## Write out the run description
    if write is True:
        output_filename = os.path.join(path_local, name+'.meta.json')
        with open(output_filename , 'w') as metafile:
            json.dump(rund, metafile, indent=4)
    return rund

## Run status
def run_status(i, total_steps, sec_left, sec_passed, process_name):
    if i%50==0 or i==0:
        sec_left=int(sec_left)
        sys.stdout.write("\n %s has finished %i/%i steps, %ih%im%is left to go"%(process_name, i+1, total_steps, sec_left/60/60, sec_left/60%60, sec_left%60)) 
        sys.stdout.flush()

## Push file to S3
def push_to_s3(s3conn, file_name, 
               bucket_name=BUCKET_NAME, 
               folder_name=FOLDER_NAME):
    ## Fix S3 folder name if needed
    if folder_name[-1] != '/':
        folder_name += '/'
    ## Connect to bucket on S3
    try:
        bucket = s3conn.get_bucket(bucket_name)
    except boto.exception.S3ResponseError:
        log_it('Bucket connection gives error, trying to create bucket')
        bucket = s3conn.create_bucket(bucket_name)
    ## Upload the file
    #os.system('s3cmd put %s s3://%s/%s/%s'%(file_name, bucket_name, 
    #                                        folder_name,
    #                                        file_name.split('/')[-1]))
    key_name = folder_name+file_name.split('/')[-1]
    key = bucket.new_key(key_name)
    key.set_contents_from_filename(file_name)
    key.close()
  
## Data creation functions, the run and the recording callback
def run_callback(sarc, input={}): 
    """Choose what to log, return it"""
    atpase = lambda s: sum(sum(s.last_transitions,[]),[]).count('31')
    appendd = lambda n,v: input[n].append(v)
    appendd('ax',     sarc.axialforce())
    appendd('ra',     sarc.radialforce())
    appendd('rt',     sarc.radialtension())
    appendd('fb',     sarc.get_frac_in_states())
    appendd('atpase', atpase(sarc))
    return input

def run(sarc, timesteps, queue=None, folder_name=FOLDER_NAME, 
        local_path=os.path.expanduser('~/data/')):
    """Create data, save it locally, upload it to S3 and log to SDB"""
    random.seed() #Must do in process to get proper dice rolls
    ## Create data names and open files
    result_name = local_path + time.strftime('%Y%m%d.') + str(uuid.uuid1())
    meta_name = result_name+'.meta.pkl'
    data_name = result_name+'.data.pkl'
    sarc_name = result_name+'.sarc.shelf'
    if os.path.exists(local_path) is False:
        os.makedirs(local_path)
    meta_file = open(meta_name, 'w')
    data_file = open(data_name, 'w')
    sarc_file = shelve.open(sarc_name, protocol=2)
    ## Make the data and log it
    results = {'ax':[], 'ra':[], 'rt':[], 'fb':[], 'atpase':[]}
    tic = time.time()
    for tstep in range(timesteps):
        sarc.timestep()
        results = run_callback(sarc, results)
        sarc_file[str(tstep)]=copy.deepcopy(sarc)
        # Update on how it is going
        toc = int((time.time()-tic) / (tstep+1) * (timesteps-tstep-1))
        run_status(tstep, timesteps, toc, time.time()-tic,
                   mp.current_process().name)
    ## Write to disk
    # Close and compress the large sarcomere copy file
    sarc_file.close()
    try:
        os.system('7za a -mx=3 %s %s'%(sarc_name+'.7z', sarc_name))
        sarc_name += '.7z'
    except:
        print("Couldn't compress, is 7zip installed?")
    # Saving smaller datas
    results['sarc'] = sarc #Nice to have one copy in your hatband
    pickle.dump(results, data_file, 2)
    data_file.close()
    # Save metadatas
    meta_data = {'folder':folder_name,
                 'name':result_name,
                 'run_length': timesteps, 
                 'lattice_spacing': sarc.lattice_spacing,
                 'z_line': sarc.z_line}
    pickle.dump(meta_data, meta_file, 0)
    meta_file.close()
    ## Pass the data and meta filenames to the queue
    if queue is not None:
        queue.put(data_name) 
        queue.put(meta_name) 
        queue.put(sarc_name)
    return (meta_name, data_name, sarc_name)

## Just like calling the main
def like_main(ls=14, sl=1250, ts=10, ub=FOLDER_NAME, multiproc=1, loops=1):
    return main(['--lattice_spacing', str(ls), 
                 '--z_line', str(sl), 
                 '--timesteps', str(ts), 
                 '--folder', str(ub), 
                 '--multiprocessing', str(multiproc), 
                 '--loops', str(loops)])

## Our main man
def main(argv=None):
    ## Get our args from the command line if not passed directly
    if argv is None:
        argv = sys.argv[1:]
    ## Parse arguments into values
    print(argv)
    parser = optparse.OptionParser(help_mesg)
    parser.add_option('-t', '--timesteps', dest="timesteps",
                      default=10, type='int',
                      help='number of timesteps for each run [10]')
    parser.add_option('-l', '--loops', dest="loops",
                      default=1, type='int',
                      help = 'how many runs to do on each core [1]')
    parser.add_option('-s', '--lattice_spacing', dest="ls",
                      default=None, type='float',
                      help='lattice spacing to a number [14.0]')
    parser.add_option('-z', '--z_line', dest="z_line",
                      default=None, type='int',
                      help='z-line where the thin filaments end [1250]')
    parser.add_option('-b', '--bucket', dest="bucket_name",
                      default=BUCKET_NAME, type='string',
                      help='set the bucket to which results will be uploaded')
    parser.add_option('-f', '--folder', dest="folder_name",
                      default=FOLDER_NAME, type='string',
                      help='set the folder to which results will be uploaded')
    parser.add_option('-m', '--multiprocessing', action="store_const",
                      dest="proc_num", default=1, const=mp.cpu_count(),
                      help='run as many copies as there are cores [False]')
    parser.add_option('--halt', action="store_true",
                      dest="the_end_of_the_end", default=False, 
                      help='shutdown computer on completion [False]')
    (options, args) = parser.parse_args(argv)
    ## For each loop we are supposed to run
    for loop in range(options.loops):
        ## Initialize the half-sarcomeres with passed ls and z-lines
        sarcs = [hs.hs(options.ls, options.z_line) for i in
                 range(options.proc_num)]
        ## Set up queue
        queue = mp.Queue()
        ## Set up processes
        process = [mp.Process(target=run, args=(s, options.timesteps, 
                        queue, options.folder_name)) for s in sarcs] 
        ## Start processes, wait for them to finish
        [p.start() for p in process]
        [p.join() for p in process]
        ## Upload the results
        log_it("Uploading results of this loop")
        s3conn = boto.connect_s3()
        while not queue.empty():
            push_to_s3(s3conn, queue.get(), 
                       options.bucket_name, 
                       options.folder_name)
        log_it("Gonna run runs "+str(options.loops-1-loop)+" more times")
    if options.the_end_of_the_end is True:
        os.system("sudo shutdown now -h")
    return 0 #Successful termination
    

if __name__ == "__main__":
    sys.exit(main())
