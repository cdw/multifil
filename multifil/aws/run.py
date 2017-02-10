#!/usr/bin/env python
# encoding: utf-8
"""
run.py - control a run 

run.emit_meta produces a meta file that describes what we want a run to do: the
values of the z_line, lattice spacing, and actin permissiveness through the run
and where it will be stored after completion.

run.manage controls a run, locally on laptop or locally on an aws node, in
either case without needing to know where it is running. It takes in a meta file
produced by emit_run and uses it to configure a sarcomere 

Created by Dave Williams on 2016-07-02
"""

import sys
import os
import shutil
import subprocess
import time
import uuid
import ujson as json
import zipfile
import multiprocessing as mp
import boto
import numpy as np

from .. import hs

## Configure a run via a saved meta file
def emit_meta(path_local, path_s3, timestep_length, timestep_number, 
              z_line=None, lattice_spacing=None, actin_permissiveness=None,
              comment = None, write = True, name=None):
    """Produce a structured JSON file that will be consumed to create a run
    
    emit_meta is intended as the counterpoint to the manage class and 
    execute_run. Import emit_run into an interactive workspace and populate 
    a directory with run configurations to be executed by a cluster.
    
    
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
        Same as for z-line. The tuple format is currently: (cycle_period, phase, 
        stim_duration, influx_time, half_life)
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
    >>> emit_run('./eggs', 'EGGS', 5, 10, (1250, 25, 25), (16, 1250), .75, comment='Just a test', write=False)  #doctest: +ELLIPSIS
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
    if name is None:
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
        return offset + 0.5 * amp * np.cos(2*np.pi*time/period)
    string_zln = "lambda offset, amp, period, time: offset + 0.5 * amp * np.cos(2*np.pi*time/period)"
    # Parse sarcomere length
    rund['z_line_args'] = str(z_line) # For easier parsing by pandas
    if hasattr(z_line, "__iter__"):
        rund['z_line'] = variable_z_line(*z_line)
        rund['z_line_func'] = string_zln
    else:
        rund['z_line'] = z_line
        rund['z_line_func'] = None
   # Define variable lattice spacing
    def variable_lattice_spacing(rest_ls, rest_zln):
        """Assumes constant volume relation, takes in rest ls and zline"""
        return np.sqrt(np.power(rest_ls,2) * rest_zln / rund['z_line'])
    string_ls = "lambda rest_ls, rest_zln: np.sqrt(np.power(rest_ls,2) * rest_zln / rund['z_line'])"
    # Parse lattice spacing
    rund['lattice_spacing_args'] = str(lattice_spacing) # For pandas
    if hasattr(lattice_spacing, "__iter__"):
        rund['lattice_spacing'] = variable_lattice_spacing(*lattice_spacing)
        rund['lattice_spacing_func'] = string_ls
    else:
        rund['lattice_spacing'] = lattice_spacing
        rund['lattice_spacing_func'] = None
    # Define variable actin permissiveness
    def variable_actin_permissiveness(cycle_period, phase, stim_duration, 
                                      influx_time, half_life):
        """Requires period of stim cycle, phase relative to longest length
        point, duration of on time, time from 10 to 90% influx level, and 
        the half-life of the Ca2+ out-pumping.
        """
        # Things we need to know for the shape of a single cycle
        decay_rate = np.log(1/2)/half_life
        growth_rate = 0.5*influx_time
        max_signal = 1.0
        # Things we need to know for the cyclical nature of the signal 
        run_step_number = timestep_number #for ease of reading
        cycle_step_number = int(cycle_period/timestep_length)
        cycle_time_trace = np.arange(0, cycle_period, timestep_length)
        try:
            steps_before_stim = np.argwhere(
                cycle_time_trace>=(cycle_period*(phase%1)))[0][0]
        except IndexError:
            assert 0 == len(np.argwhere(
                cycle_time_trace>=(cycle_period*(phase%1))))
            steps_before_stim = 0 #b/c phase was 0.999 or similar
        stim_step_number = int(stim_duration/timestep_length)
        no_stim_step_number = cycle_step_number - stim_step_number
        # Things we need to know for smoothing
        sd = 1 #standard deviation of smoothing window in ms
        sw = 3 #smoothing window in ms
        base_normal = np.exp(-np.arange(-sw,sw,timestep_length)**2/(2*sd**2))
        normal = base_normal/sum(base_normal)
        # Step through, generating signal
        out = [0.1]
        for i in range(steps_before_stim):
            out.append(out[-1])
        while len(out)<(4*cycle_step_number+run_step_number):
            for i in range(stim_step_number):
                growth = timestep_length * out[-1] * (growth_rate) *\
                        (1-out[-1]/max_signal)
                out.append(out[-1]+growth)
            for i in range(no_stim_step_number):
                decay = timestep_length * out[-1] * decay_rate
                out.append(out[-1]+decay)
        # Smooth signal 
        out = np.convolve(normal, out)
        return out[2*cycle_step_number:2*cycle_step_number+run_step_number]
    string_ap = """ def variable_actin_permissiveness(cycle_period, phase, stim_duration, 
                                      influx_time, half_life):
        \"\"\"Requires period of stim cycle, phase relative to longest length
        point, duration of on time, time from 10 to 90% influx level, and 
        the half-life of the Ca2+ out-pumping.
        \"\"\"
        # Things we need to know for the shape of a single cycle
        decay_rate = np.log(1/2)/half_life
        growth_rate = 0.5*influx_time
        max_signal = 1.0
        # Things we need to know for the cyclical nature of the signal 
        run_step_number = timestep_number #for ease of reading
        cycle_step_number = int(cycle_period/timestep_length)
        cycle_time_trace = np.arange(0, cycle_period, timestep_length)
        steps_before_stim = np.argwhere(
            cycle_time_trace>=(cycle_period*phase))[0][0]
        stim_step_number = int(stim_duration/timestep_length)
        no_stim_step_number = cycle_step_number - stim_step_number
        # Things we need to know for smoothing
        sd = 1 #standard deviation of smoothing window in ms
        sw = 3 #smoothing window in ms
        base_normal = np.exp(-np.arange(-sw,sw,timestep_length)**2/(2*sd**2))
        normal = base_normal/sum(base_normal)
        # Step through, generating signal
        out = [0.1]
        for i in range(steps_before_stim):
            out.append(out[-1])
        while len(out)<(4*cycle_step_number+run_step_number):
            for i in range(stim_step_number):
                growth = timestep_length * out[-1] * (growth_rate) *\
                        (1-out[-1]/max_signal)
                out.append(out[-1]+growth)
            for i in range(no_stim_step_number):
                decay = timestep_length * out[-1] * decay_rate
                out.append(out[-1]+decay)
        # Smooth signal 
        out = np.convolve(normal, out)
        return out[2*cycle_step_number:2*cycle_step_number+run_step_number]"""
    # Parse actin permissiveness 
    rund['actin_permissiveness_args'] = str(actin_permissiveness) # For pandas
    if hasattr(actin_permissiveness, "__iter__"):
        rund['actin_permissiveness'] = variable_actin_permissiveness(
            *actin_permissiveness)
        rund['actin_permissiveness_func'] = string_ap
    else:
        rund['actin_permissiveness'] = actin_permissiveness
        rund['actin_permissiveness_func'] = None
    ## Write out the run description
    if write is True:
        output_filename = os.path.join(path_local, name+'.meta.json')
        with open(output_filename , 'w') as metafile:
            json.dump(rund, metafile, indent=4)
    return rund


## Manage a local run
class manage(object):
    """Run, now with extra objor flavor"""
    def __init__(self, metafile, unattended=True):
        """Create a managed instance of the sarc, optionally running it

        Parameters
        ----------
        metafile: string
            The location of the metafile describing the run to be worked 
            through. Can be local or on S3. Assumed to be on S3 if the local
            file does not exist.
        unattended: boolean
            Whether to complete the run without further intervention or treat 
            as an interactive session. 
        """
        self.s3 = s3()
        self.uuid = metafile.split('/')[-1].split('.')[0]
        self.working_dir = self._make_working_dir(self.uuid)
        self.metafile = self._parse_metafile_location(metafile)
        self.meta = self.unpack_meta(self.metafile)
        self.sarc = self.unpack_meta_to_sarc(self.meta)
        if unattended:
            try:
                self.run_and_save()
            except:
                mp.current_process().terminate()
    
    @staticmethod
    def _make_working_dir(name):
        """Create a temporary working directory and return the name"""
        wdname = '/tmp/'+name
        os.makedirs(wdname, exist_ok=True)
        return wdname
    
    def _parse_metafile_location(self, metafile):
        """Parse the passed location, downloading the metafile if necessary"""
        if not os.path.exists(metafile):
            return self.s3.pull_from_s3(metafile, self.working_dir)
        else:
            mfn = '/'+metafile.split('/')[-1]
            return shutil.copyfile(metafile, self.working_dir+mfn)
    
    @staticmethod
    def unpack_meta(metafilename):
        """Unpack the local meta file to a dictionary"""
        with open(metafilename, 'r') as metafile:
            meta = json.load(metafile)
        return meta
    
    @staticmethod
    def unpack_meta_to_sarc(meta):
        """Unpack the local meta file and instantiate a sarc as defined 
        in the meta file
        """
        # Prep single values for instantiation of hs
        none_if_list = lambda s: None if type(meta[s]) is list else meta[s]
        z_line = none_if_list('z_line')
        lattice_spacing = none_if_list('lattice_spacing')
        actin_permissiveness = none_if_list('actin_permissiveness')
        # Time dependent values
        time_dep_dict = {}
        for prop in ['lattice_spacing', 'z_line', 'actin_permissiveness']:
            if type(meta[prop]) is list:
                time_dep_dict[prop] = meta[prop]
        # Instantiate sarcomere
        sarc = hs.hs(
            lattice_spacing = lattice_spacing,
            z_line = z_line,
            actin_permissiveness = actin_permissiveness,
            timestep_len = meta['timestep_length'],
            time_dependence = time_dep_dict,
            )
        return sarc
    
    def _copy_file_to_final_location(self, temp_full_fn, final_loc=None):
        """Copy file from the temporary location to the final resting places
        
        Parameters
        ----------
        temp_full_fn: string
            the full location of the temporary file to be copied
        final_loc: string
            an optional string for an extra local directory to copy the file to
        """
        temp_loc = temp_full_fn
        file_name = '/' + temp_loc.split('/')[-1]
        # Upload to S3
        if self.meta['path_s3'] is not None:
            s3_loc = self.meta['path_s3']
            self.s3.push_to_s3(temp_loc, s3_loc)
        # Store in final local path 
        if self.meta['path_local'] is not None:
            local_loc = os.path.abspath(os.path.expanduser(
                self.meta['path_local'])) + file_name
            shutil.copyfile(temp_loc, local_loc)
        # Save to passes local location
        if final_loc is not None:
            local_loc = os.path.abspath(os.path.expanduser(location)) \
                    + file_name
            shutil.copyfile(temp_loc, local_loc)
    
    def run_and_save(self):
        """Complete a run according to the loaded meta configuration and save
        results to meta-specified s3 and local locations"""
        # Initialize data and sarc
        self.sarcfile = sarc_file(self.sarc, self.meta, self.working_dir)
        self.datafile = data_file(self.sarc, self.meta, self.working_dir)
        # Run away
        np.random.seed()
        tic = time.time()
        for timestep in range(self.meta['timestep_number']):
            self.sarc.timestep(timestep)
            self.datafile.append()
            self.sarcfile.append()
            # Update on how it is going
            self._run_status(timestep, tic, 100)
        # Finalize and save files to final locations
        self._log_it("model finished, uploading")
        self._copy_file_to_final_location(self.metafile)
        data_final_name = self.datafile.finalize()
        self._copy_file_to_final_location(data_final_name)
        self.datafile.delete() # clean up temp files
        sarc_final_name = self.sarcfile.finalize()
        self._copy_file_to_final_location(sarc_final_name)
        self.sarcfile.delete() # clean up temp files
        self._log_it("uploading finished, done with this run")
    
    def _run_status(self, timestep, start, every):
        """Report the run status"""
        if timestep%every==0 or timestep==0:
            total_steps = self.meta['timestep_number']
            sec_passed = time.time()-start
            sec_left = int(sec_passed/(timestep+1)*(total_steps-timestep-1))
            self._log_it("finished %i/%i steps, %ih%im%is left"%(
                timestep+1, total_steps, 
                sec_left/60/60, sec_left/60%60, sec_left%60))
    
    @staticmethod
    def _log_it(message):
        """Print message to sys.stdout"""
        sys.stdout.write("run.py " + mp.current_process().name + 
                " ## " + message + "\n")
        sys.stdout.flush()


## File management 
class sarc_file(object):
    def __init__(self, sarc, meta, working_dir):
        """Handles recording a sarcomere dict to disk at each timestep"""
        self.sarc = sarc 
        self.meta = meta
        self.working_directory = working_dir
        sarc_name = '/'+meta['name']+'.sarc.json'
        self.working_filename = self.working_directory + sarc_name
        self.working_file = open(self.working_filename, 'a')
        self.next_write = '[\n'
        self.append(True)
    
    def append(self, first=False):
        """Add the current timestep sarcomere to the sarc file"""
        if not first:
            self.next_write +=',\n'
        self.next_write += json.dumps(self.sarc.to_dict(), sort_keys=True)
        self.working_file.write(self.next_write)
        self.next_write = ''
    
    def finalize(self):
        """Close the current sarcomere file for proper JSON formatting"""
        self.working_file.write('\n]')
        self.working_file.close()
        self.zip_filename = self.working_filename[:-4]+'tar.gz'
        cp = subprocess.run(['tar', 'czf', self.zip_filename, self.working_filename])
        os.remove(self.working_filename)
        return self.zip_filename
    
    def delete(self):
        """Delete the sarc zip file from disk"""
        os.remove(self.zip_filename)


class data_file(object):
    def __init__(self, sarc, meta, working_dir):
        """Generate the dictionary for use with the below data callback"""
        self.sarc = sarc 
        self.meta = meta
        self.working_directory = working_dir
        self.data_dict = {
            'name': self.meta['name'],
            'timestep_length': self.sarc.timestep_len,
            'timestep': [],
            'z_line': [],
            'lattice_spacing': [],
            'axial_force': [], 
            'radial_force_y': [],
            'radial_force_z': [],
            'radial_tension': [],
            'xb_fraction_free': [],
            'xb_fraction_loose': [],
            'xb_fraction_tight': [],
            'xb_trans_12': [],
            'xb_trans_23': [],
            'xb_trans_31': [],
            'xb_trans_21': [],
            'xb_trans_32': [],
            'xb_trans_13': [],
            'xb_trans_static': [],
            'actin_permissiveness': [],
            'thick_displace_mean': [],
            'thick_displace_max': [],
            'thick_displace_min': [],
            'thick_displace_std': [],
            'thin_displace_mean': [],
            'thin_displace_max': [],
            'thin_displace_min': [],
            'thin_displace_std': [],
        }
    
    def append(self):
        """Digest out the non-vector values we want to record for each 
        timestep and append them to the data_dict. This is called at each 
        timestep to build a dict for inclusion in a pandas dataframe.
        """
        ## Lambda helpers
        ad = lambda n,v: self.data_dict[n].append(v)
        ## Calculated components
        radial_force = self.sarc.radialforce()
        xb_fracs = self.sarc.get_frac_in_states()
        xb_trans = sum(sum(self.sarc.last_transitions,[]),[])
        act_perm = np.mean(self.sarc.actin_permissiveness)
        thick_d = np.hstack([t.displacement_per_crown() 
                             for t in self.sarc.thick])
        thin_d = np.hstack([t.displacement_per_node() 
                            for t in self.sarc.thin])
        ## Dictionary work
        ad('timestep', self.sarc.current_timestep)
        ad('z_line', self.sarc.z_line)
        ad('lattice_spacing', self.sarc.lattice_spacing)
        ad('axial_force', self.sarc.axialforce())
        ad('radial_force_y', radial_force[0])
        ad('radial_force_z', radial_force[1])
        ad('radial_tension', self.sarc.radialtension())
        ad('xb_fraction_free', xb_fracs[0])
        ad('xb_fraction_loose', xb_fracs[1])
        ad('xb_fraction_tight', xb_fracs[2])
        ad('xb_trans_12', xb_trans.count('12'))
        ad('xb_trans_23', xb_trans.count('23'))
        ad('xb_trans_31', xb_trans.count('31'))
        ad('xb_trans_21', xb_trans.count('21'))
        ad('xb_trans_32', xb_trans.count('32'))
        ad('xb_trans_13', xb_trans.count('13'))
        ad('xb_trans_static', xb_trans.count(None))
        ad('actin_permissiveness', act_perm)
        ad('thick_displace_mean', np.mean(thick_d))
        ad('thick_displace_max', np.max(thick_d))
        ad('thick_displace_min', np.min(thick_d))
        ad('thick_displace_std', np.std(thick_d))
        ad('thin_displace_mean', np.mean(thin_d))
        ad('thin_displace_max', np.max(thin_d))
        ad('thin_displace_min', np.min(thin_d))
        ad('thin_displace_std', np.std(thin_d))

    def finalize(self):
        """Write the data dict to the temporary file location"""
        data_name = '/'+self.meta['name']+'.data.json'
        self.working_filename = self.working_directory + data_name
        with open(self.working_filename, 'w') as datafile:
            json.dump(self.data_dict, datafile, sort_keys=True)
        return self.working_filename
    
    def delete(self):
        """Delete the data file from disk"""
        try: 
            os.remove(self.working_filename)
        except FileNotFoundError:
            print("File not created yet")


class s3(object):
    def __init__(self):
        """Provide an interface to to S3 that hides some error handling"""
        self._refresh_s3_connection()
    
    def _refresh_s3_connection(self):
        """Reconnect to s3, the connection gets dropped sometimes"""
        self.s3 = boto.connect_s3()
    
    def _get_bucket(self, bucket_name):
        """Return link to a bucket name"""
        try:
            bucket = self.s3.get_bucket(bucket_name)
        except (boto.exception.BotoClientError, 
                boto.exception.BotoClientError) as e:
            print(e)
            print("Trying to reconnect to s3")
            self._refresh_s3_connection()
            bucket = self.s3.get_bucket(bucket_name)
        return bucket
    
    def pull_from_s3(self, name, local='./'):
        """Given a key on S3, download it to a local file 
        
        Parameters
        ----------
        name: string
            bucket/keyname to download. can be prefixed by 's3://', '/', or 
            by nothing
        local: string
            local directory to download key into, defaults to current directory
        
        Returns
        -------
        None
        
        Examples
        --------
        >>>pull_from_s3('s3://just_a_test_bucket/test')
        >>>os.remove('test')
        >>>pull_from_s3('just_a_test_bucket/test')
        >>>os.remove('test')
        """
        # Parse name
        bucket_name = [n for n in name.split('/') if len(n)>3][0] # rm s3:// & /
        key_name = name[len(bucket_name)+name.index(bucket_name):]
        file_name = key_name.split('/')[-1]
        # Connect to bucket
        bucket = self._get_bucket(bucket_name)
        # Connect to key/file
        key = bucket.get_key(key_name)
        # Prep local dirs to receive key
        local = os.path.abspath(os.path.expanduser(local))
        os.makedirs(local, exist_ok=True)
        # Download key
        downloaded_name = local+'/'+file_name
        key.get_contents_to_filename(downloaded_name)
        if key.size != os.stat(downloaded_name).st_size:
            print("Size mismatch, downloading again for %s: "%downloaded_name)
            key.get_contents_to_filename(downloaded_name)
        return downloaded_name
    
    def push_to_s3(self, local, remote):
        """Given a local file, push it to a location on s3
        
        Parameters
        ----------
        local: string
            path to local file to be uploaded to s3
        remote: string
            bucket and (optional) folder to upload local file to, the resulting
            key will be remote+/+local_filename. Can be in formats:
                s3://bucket/optional_folder/key
                /bucket/optional_folder/key
                bucket/optional_folder/key
        
        Returns
        -------
        None
        """
        # Parse names
        file_name = local.split('/')[-1]
        bucket_name = [n for n in remote.split('/') if len(n)>3][0] 
        key_name = remote[len(bucket_name)+remote.index(bucket_name):]
        if len(key_name)==0 or key_name[-1] != '/':
            key_name += '/'
        # Parse bucket and folders
        bucket = self._get_bucket(bucket_name)
        key = bucket.new_key(key_name+file_name)
        key.set_contents_from_filename(local)
        if key.size != os.stat(local).st_size:
            print("Size mismatch, uploading again for %s: "%local)
            key.set_contents_from_filename(local)
        return
