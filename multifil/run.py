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
import cPickle as pickle
import shelve
import copy
import optparse 
import multiprocessing as mp
import boto
import hs
import numpy.random as random

## Settings
BUCKET_NAME = 'model_results' # For results files
FOLDER_NAME = 'test_disregard/' # Needs trailing /


## Help message
help_mesg = "Around here, we like to call our functions with a little \n \
      thing called options. Here's how we'd go about it: \n \
      usage: %prog [options] arg1 arg2 "

## Logging
def log_it(message):
    """Print message to sys.stdout"""
    sys.stdout.write("run.py " + mp.current_process().name + 
            " ## " + message + "\n")
    sys.stdout.flush()

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
        print "Couldn't compress, is 7zip installed?"
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
    print argv
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
