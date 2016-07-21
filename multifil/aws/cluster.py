#!/usr/bin/env python
# encoding: utf-8
"""
cluster.py - manage the behaviour and start up of a cluster

Created by Dave Williams on 2016-07-19
"""

import os
import sys
import time
import string
import copy
import subprocess as subp
import configparser
import boto


## Defaults
BASE_PATH = os.path.expanduser('~/code/multifil/')
CODE_DIR = 'multifil'
CODE_LOCATION =  BASE_PATH + CODE_DIR
USER_DATA = CODE_LOCATION + '/aws/userdata.py'
CODE_BUCKET = 'model_code'
JOB_QUEUE = 'job_queue'
KEY_FILE = os.path.expanduser('~/.aws/keys/id_gsg-keypair')
KEY_NAME = 'gsg-keypair' 
SECURITY_GROUP_ID = 'sg-1ddae166'
SUBNET_IDS = ('subnet-945f5ca9', 'subnet-c9ce2b80', # each corresponds to
              'subnet-64e6d34e', 'subnet-bd94a0e5') # an availability zone
#AMI = ('ami-2d39803a', 'c4.xlarge') # Ubuntu
AMI = ('ami-2d39803a', 't2.medium') # Testing
SPOT_BID = 0.209 # bidding the on-demand price


## Helper functions, quite substantial
def print_direct(string):
    """Print the given string straight to the stdout"""
    sys.stdout.truncate(0)
    sys.stdout.write(string)
    sys.stdout.flush()
    return

def get_access_keys(filename=os.path.expanduser('~/.aws/credentials'), 
                    section='cluster'):
    """Parse out the access and secret keys"""
    config = configparser.ConfigParser()
    config.read(filename)
    id = config.get(section,'aws_access_key_id')
    secret = config.get(section,'aws_secret_access_key')
    return id, secret

def load_userdata(filename='userdata.py', queue_name=JOB_QUEUE):
    id, secret = get_access_keys()
    user_data_dict = {
        'aws_access_key': id, 
        'aws_secret_key': secret,
        'job_queue_name': queue_name,
        'code_zip_key': "s3://%s/%s.zip"%(CODE_BUCKET, CODE_DIR)}
    with open(filename, 'r') as udfile:
        ud_template = string.Template(udfile.read())
    return ud_template.substitute(user_data_dict)

def update_code_on_s3():
    """Update the code on s3 from our local copy"""
    zipname = CODE_DIR+'.zip' 
    cmds = (
        "cd ../..; zip -roTFS %s %s"%(zipname, CODE_DIR), 
        "cd ../..; aws s3 cp %s s3://%s/"%(zipname, CODE_BUCKET))
    print(os.getcwd())
    print(cmds)
    [print(subp.call(c, shell=True)) for c in cmds]

def launch_on_demand_instances(ec2, num_of, userdata, 
                               ami=AMI[0], inst_type=AMI[1]):
    if len(userdata) > 16*1024:
        print("error: User data file is too big")
        return
    reservation = ec2.run_instances(
        image_id           = ami,
        key_name           = KEY_NAME, 
        security_group_ids = [SECURITY_GROUP_ID],
        user_data          = userdata,
        instance_type      = inst_type,
        min_count          = num_of,
        max_count          = num_of,
        subnet_id          = SUBNET_IDS[1])
    time.sleep(.5) # Give the machines time to register
    nodes = copy.copy(reservation.instances)
    return nodes

def launch_spot_instances(ec2, num_of, userdata, bid=SPOT_BID,
                          ami=AMI[0], inst_type=AMI[1]):
    if len(userdata) > 16*1024:
        print("error: User data file is too big")
        return
    reservation = ec2.request_spot_instances(
        price              = bid, 
        image_id           = ami,
        key_name           = KEY_NAME, 
        security_group_ids = [SECURITY_GROUP_ID],
        user_data          = userdata,
        instance_type      = inst_type,
        min_count          = num_of,
        max_count          = num_of,
        subnet_id          = SUBNET_IDS[1])
    time.sleep(.5) # Give the machines time to register
    return reservation


class cluster(object):
    def __init__(self, 
                 number_of_instances, 
                 queue_name=JOB_QUEUE, 
                 userdata=USER_DATA,
                 use_spot=True):
        """A cluster management object"""
        self.number_of_instances = number_of_instances
        self.queue_name = queue_name
        self.userdata = load_userdata(userdata, queue_name)
        self.use_spot = use_spot
        self.s3 = boto.connect_s3()
        self.ec2 = boto.connect_ec2()

    def launch(self):
        """Get the cluster rolling, manual to make you think a bit"""
        print("Uploading code to S3")
        update_code_on_s3()
        print("Creating reservation")
        # TODO: This next bit could be refactored to be prettier 
        if self.use_spot is True:
            nodes = launch_spot_instances(self.ec2, 
                                          self.number_of_instances,
                                          self.userdata)
            ids = [node.id for node in nodes]
            node_states = lambda nodes: [node.state == 'active' 
                            for node in nodes]
            node_update = lambda : self.ec2.get_all_spot_instance_requests(ids)
        else:
            nodes = launch_on_demand_instances(self.ec2, 
                                               self.number_of_instances,
                                               self.userdata)
            ids = [node.id for node in nodes]
            node_states = lambda nodes: [node.state_code == 16 
                                         for node in nodes]
            node_update = lambda : [inst for res in 
                                    self.ec2.get_all_instances(ids)
                                    for inst in res.instances]
        print("Nodes are starting...")
        while not all(node_states(nodes)):
            nodes = node_update()
            ready = sum(node_states(nodes))
            print_direct("\r%i of %i nodes are ready"%(ready, len(nodes)))
            time.sleep(1)
        print_direct("\nAll nodes ready \n")
        self.nodes = nodes
        return nodes

    def kill_cluster(self):
        """Terminate the cluster nodes"""
        try:
            [node.terminate() for node in self.nodes]
        except:
            [node.cancel() for node in self.nodes]
            ids = [node.instance_id for node in self.nodes]
            [instance.terminate() for reservation in EC2.get_all_instances(ids)
             for instance in reservation.instances]

    def node_ip_addresses(self):
        """Print the ip addresses for each node"""
        [print(instance.ip_address) for instance in self.nodes]
    
    
