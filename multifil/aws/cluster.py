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
import itertools
import subprocess as subp
import configparser
import boto


## Defaults
BASE_PATH = os.path.abspath(os.path.split(__file__)[0]+'/../..')+'/'
CODE_DIR = 'multifil'
CODE_LOCATION =  BASE_PATH + CODE_DIR
USER_DATA = CODE_LOCATION + '/aws/userdata.py'
CODE_BUCKET = 'model-code'
JOB_QUEUE = 'job-queue'
STATUS_QUEUE = 'status-queue'
KEY_FILE = os.path.expanduser('~/.aws/keys/id_aws')
KEY_NAME = 'id_aws' 
SECURITY_GROUP_ID = 'sg-2a31b650'
SUBNET_IDS = {'us-east-1a':'subnet-7653873f', # map an availability zone
              'us-east-1e':'subnet-a5957299', # to the right VPC
              'us-east-1d':'subnet-00ff1b2d', 
              'us-east-1b':'subnet-39a5bf61'} 
AMI = ('ami-2d39803a', 'c4.xlarge') # Ubuntu
HD_SIZE = 200 # primary drive size in GB
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

def get_bdm(ec2=boto.connect_ec2(), ami=AMI[0], size=HD_SIZE):
    bdm = ec2.get_image(ami).block_device_mapping
    bdm['/dev/sda1'].size = size
    bdm['/dev/sda1'].encrypted = None
    return bdm

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
    # Not fragile at all... ha
    cmds = (
        "cd %s; zip -roTFS -imultifil/\* -isetup* %s ./"%(BASE_PATH, zipname), 
        "cd %s; aws s3 cp %s s3://%s/"%(BASE_PATH, zipname, CODE_BUCKET))
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
        subnet_id          = SUBNET_IDS['us-east-1a'], 
        block_device_map   = get_bdm(ec2))
    time.sleep(.5) # Give the machines time to register
    nodes = copy.copy(reservation.instances)
    return nodes

def launch_spot_instances(ec2, num_of, userdata, bid=SPOT_BID,
                          ami=AMI[0], inst_type=AMI[1]):
    if len(userdata) > 16*1024:
        print("error: User data file is too big")
        return
    # Choose cheapest availability zone
    sphs = ec2.get_spot_price_history(filters={'instance_type':inst_type})
    prices = [sph.price for sph in sphs]
    availability_zone = sphs[prices.index(min(prices))].availability_zone
    reservation = ec2.request_spot_instances(
        price              = bid, 
        image_id           = ami,
        key_name           = KEY_NAME, 
        security_group_ids = [SECURITY_GROUP_ID],
        user_data          = userdata.encode('ascii'),
        instance_type      = inst_type,
        count              = num_of,
        placement          = availability_zone,
        subnet_id          = SUBNET_IDS[availability_zone],
        block_device_map   = get_bdm(ec2))
    time.sleep(.5) # Give the machines time to register
    return reservation

def watch_cluster():
    """Give real-time updates on what is happening aboard our cluster"""
    #Make it pretty, or pretty trippy
    range_plus =lambda ri,re,s: [str(i)+s for i in range(ri,re)]
    styles = ["\033["+''.join(style) for style in itertools.product(
        range_plus(0,3,';'), range_plus(30,38,';'), range_plus(40,48,'m'))]
    #Make it work
    sqs = boto.connect_sqs()
    status_queue = sqs.get_queue(STATUS_QUEUE)
    ec2 = boto.connect_ec2()
    print("Starting cluster watch, ^c to stop")
    while True: #quit via ^C
        try:
            # Gather and report messages
            if status_queue.count()>0:
                while True:
                    msg = status_queue.read()
                    body = msg.get_body()
                    last_ip = int(body.split('-')[1].split('.')[-1])
                    style = styles[last_ip%len(styles)]
                    print(style+body)
                    status_queue.delete_message(msg)
            # Make sure some instances are running
            running_instances = ec2.get_all_instances(
                filters=({'instance-state-code':0, 
                          'instance-state-code':16}))
            if len(running_instances) == 0:
                print("\nNo running instances found")
                break
            # Don't hammer the connection
            time.sleep(3)
        except KeyboardInterrupt: #^c pressed
            print("\nMy watch has ended")
            break
        except AttributeError: #no message to read body from
            pass


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
        if self.use_spot:
            nodes = self.ec2.get_only_instances([n.instance_id for n in nodes])
        self.nodes = nodes
        return nodes

    def kill_cluster(self):
        """Terminate the cluster nodes"""
        [node.terminate() for node in self.nodes]

    def node_ip_addresses(self):
        """Print the ip addresses for each node"""
        [print(instance.ip_address) for instance in nodes]
    
    
