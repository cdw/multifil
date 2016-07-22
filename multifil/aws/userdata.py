#!/usr/bin/env python3
# encoding: utf-8
"""
sqs_control_userdata_script.py - control an aws instance from a sqs queue

Run this guy on startup as a userdata script and he will connect to 
s3 to download code to a directory, and run commands in it that are 
provided by an SQS queue, one job at a time per core

Processing as a string template, we replace the following keys with their
equivalents:
    - aws_access_key
    - aws_secret_key
    - job_queue_name
    - code_zip_key

Created by Dave Williams on 2011-02-08
"""

## Import present packages
import os
import sys
import time
import traceback
import subprocess as subp
import multiprocessing as mp

## Handle logging and thrown fatal errors
def log_it(log_message):
    print(log_message)
    with open('/dev/console', 'w') as console:
        console.write("USER DATA: "+log_message+'\n')

def fatal_error(error_log_message, feed_me = "differently"):
    log_it("ERROR: " + error_log_message)
    log_it("SHUTTING DOWN: feed me " + feed_me + " next time")
    #os.system("shutdown now -h")

def try_and_log(command, message):
    out = subp.call(command, shell=True)
    log_it(message + str(out))

## Install extra software on the node
log_it("#"*60 + "\n START OF USERDATA SCRIPT\n"*3 + "#"*60)
try_and_log("apt-get -qq update", "Synced package index with result: ")
try_and_log("apt-get -qq install python3-scipy python3-pip unzip > \\dev\\null",
            "Installed scipy, pip, unzip with result: ")
try_and_log("pip3 install boto ujson", "Installed boto, ujson: ")

## Userdata runs as root, but in /, let's move
os.chdir('/root')
HOMEDIR = os.getcwd()+'/'

## Configure control parameters
ACCESS_KEY = '$aws_access_key'
SECRET_KEY = '$aws_secret_key'
JOB_QUEUE = '$job_queue_name'
CODE_ZIP_KEY = '$code_zip_key'

## Write out boto configuration
lines = """[Credentials]
aws_access_key_id = %s 
aws_secret_access_key = %s \n"""%(ACCESS_KEY, SECRET_KEY)
with open('.boto', 'w') as config_file:
    config_file.writelines(lines)

## Connect to aws with boto
try:
    log_it("Connecting to boto")
    import boto # Had to wait until .boto was written
    S3 = boto.connect_s3()
    SQS = boto.connect_sqs()
    SQS.get_all_queues() # Call to test if our keys were accepted
except (boto.exception.NoAuthHandlerFound, boto.exception.SQSError) as e:
    fatal_error("Probably gave bad aws keys", "valid credentials")

## Download files from passed bucket
try:
    log_it("Downloading from code bucket")
    bucket_name = [n for n in CODE_ZIP_KEY.split('/') if len(n)>3][0] #s3:// & /
    key_name = CODE_ZIP_KEY[len(bucket_name)+CODE_ZIP_KEY.index(bucket_name)+1:]
    code_bucket = S3.get_bucket(bucket_name)
    key = code_bucket.get_key(key_name)
    key.get_contents_to_filename(key_name)
    try_and_log("unzip %s"%key_name, 
                "Unzipped local code file %s with result: "%key_name)
    time.sleep(3) # poor man's race condition control!
except boto.exception.S3ResponseError:
    fatal_error("No bucket with given name %s"%(CODE_ZIP_KEY), "a valid bucket")
except IOError:
    fatal_error("Couldn't write code_bucket contents locally")


## Turn control over to the job queue
try:
    log_it(str(dir()))
    log_it("Turning things over to queue eater processes")
    commandment = "python3 -c \"import multifil;\
        multifil.aws.instance.multi_eaters('%s',shutdown=True)\""%JOB_QUEUE 
    try_and_log(commandment, "Called sub-process to manage queue eaters")
    log_it("All done")
except Exception as e:
    log_it("### An error occurred while running jobs")
    log_it("Exception of type " + str(type(e)))
    exc_type, exc_value, exc_traceback = sys.exc_info()
    log_it(repr(traceback.format_exception(exc_type, exc_value, exc_traceback)))
    log_it("Going no further, shutting down now")
finally:
    os.system('shutdown now -h')
    
    
