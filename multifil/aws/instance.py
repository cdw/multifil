#!/usr/bin/env python
# encoding: utf-8
"""
instance.py - manage the behaviour on an individual instance

Created by Dave Williams on 2016-07-05  
"""

import sys
import os
import time
import optparse
import urllib.request 
import multiprocessing as mp
import boto
from . import run

## Reporting to SQS
# This is a bit hacky and I don't like it
try:
    with urllib.request.urlopen(
        'http://169.254.169.254/latest/meta-data/public-ipv4',
        timeout=3) as response:
        ip4 = response.read().decode()
    log_to_sqs=True
    sqs = boto.connect_sqs()
    logging_queue = sqs.get_queue('status-queue')
except urllib.error.URLError:
    log_to_sqs=False


## Helper functions
def log_it(log_message):
    """Write message to console so that it can be viewed from EC2 monitor"""
    identified_message = "instance.py :" + mp.current_process().name +  \
            " ## " + log_message
    print(identified_message)
    with open('/dev/console', 'w') as console:
        console.write(identified_message + '\n')
    if log_to_sqs:
        oclock = time.strftime('%a,%H:%M') + " - "
        msg = oclock + ip4 + " - "+mp.current_process().name+": "+log_message
        logging_queue.write(logging_queue.new_message(msg))
                
def fatal_error(error_log_message, feed_me = "differently", shutdown=False):
    """Log a message likely to have torpedoed the run"""
    log_it("ERROR: " + error_log_message)
    log_it("SHUTTING DOWN: feed me " + feed_me + " next time")
    if shutdown:
        halt_system()

def running_error(exception):
    """An error that occured while running a job"""
    log_it("### An error occurred while running jobs")
    log_it("Exception of type " + str(type(exception)) + 
           ": " + exception.message)
    exc_type, exc_value, exc_traceback = sys.exc_info()
    log_it(repr(traceback.format_exception(exc_type, exc_value, exc_traceback)))

def halt_system():
    """Shut it down"""
    os.system("shutdown now -h")


## Munch tasks off a queue
class queue_eater(object):
    def __init__(self,  sqs_queue_name, id=None, secret=None, shutdown=True):
        """Consume an SQS queue. The queue consists of metafile locations.
        These locations are spun off into run.manage instances and the queue
        messages are deleted after the run.manage returns.

        Parameters
        ----------
        sqs_queue_name: string
            name of the queue to eat
        id: string
            optional AWS id access key
        secret: string
            optional AWS secret key
        shutdown: boolean
            if True, will shutdown on errors
        """
        self.name = sqs_queue_name
        self.id = id
        self.secret = secret
        self.should_shutdown = shutdown
        self._connect_to_queue()
        self.new_meta() # Load first meta msg
        try:
            if self.meta is not None: # in case queue is empty
                self.new_proc()
            while self.meta is not None:
                if not self.proc_alive():
                    if self.process.exitcode==0:
                        self.delete_meta()
                    self.new_meta()
                    if self.meta is not None:
                        self.new_proc()
        except Exception as e:
            running_error(e)
            self.shutdown()
        log_it("Ate all we can, queue eater out")
        return
    
    def _connect_to_queue(self):
        """Connect to our sqs queue"""
        try:
            log_it("Connecting to SQS queue "+self.name)
            sqs = boto.connect_sqs(self.id, self.secret)
            self.queue = sqs.get_queue(self.name)
            if type(self.queue) is type(None):
                raise KeyError("Provided queue name not found")
        except KeyError:
            fatal_error("Given queue non-existent", "a different queue name",
                        self.should_shutdown)
    
    def new_proc(self):
        """Launch a new process from the current meta message"""
        try:
            message_body = self.meta.get_body()
            log_it("Gonna run "+message_body)
            self.process = mp.Process(target = run.manage,
                                      args = (message_body,))
            self.process.start()
        except Exception as e:
            running_error(e)
            self.shutdown()
    
    def proc_alive(self):
        """Is the process alive?"""
        if self.process.is_alive():
            return True
        else:
            self.process.join() # wait until process really terminates
            return False
    
    def new_meta(self):
        """Read the next meta message"""
        self.meta = self.queue.read()
    
    def delete_meta(self):
        """Delete the current meta message"""
        self.queue.delete_message(self.meta)
    
    def shutdown(self):
        """Turn off instance"""
        if self.should_shutdown:
            log_it("Going no further, shutting down now")
            halt_system()


## Many munching mouths
def multi_eaters(sqs_queue_name, num=None, id=None, secret=None, shutdown=True):
    """Launch a number of queue eaters"""
    if num is None:
        num = mp.cpu_count()
    queue_eater_args = (sqs_queue_name, id, secret, shutdown)
    eaters = [mp.Process(target=queue_eater, args = queue_eater_args) \
              for i in range(num)]
    [e.start() for e in eaters]
    time.sleep(0.5)
    while not all([not e.is_alive() for e in eaters]):
        time.sleep(0.5)
    [e.join() for e in eaters]
    if shutdown is True:
        halt_system()


## Our main man
def main(argv=None):
    ## Get our args from the command line if not passed directly
    if argv is None:
        argv = sys.argv[1:]
    ## Parse arguments into values
    print(argv)
    parser = optparse.OptionParser("Investigate options and re-call")
    parser.add_option('-q', '--queue', dest="queue_name",
                      default="job_queue", type='string',
                      help='set the queue to look for commands in')
    parser.add_option('-i', '--id', dest="id",
                      default=None, type='string',
                      help='what aws id credential to use')
    parser.add_option('-s', '--secret', dest="secret",
                      default=None, type='string',
                      help='what aws secret key to use')
    parser.add_option('-m', '--multiprocessing', action="store_const",
                      dest="proc_num", default=1, const=mp.cpu_count(),
                      help='run as many copies as there are cores [False]')
    parser.add_option('--halt', action="store_true",
                      dest="the_end_of_the_end", default=False, 
                      help='shutdown computer on completion or error [False]')
    (options, args) = parser.parse_args(argv)
    multi_eaters(options.queue_name, options.proc_num, options.id,
                 options.secret, options.the_end_of_the_end)
    return 0 #Successful termination
    

if __name__ == "__main__":
    sys.exit(main())
