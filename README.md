[![Build Status](https://travis-ci.org/cdw/multifil.svg?branch=master)](https://travis-ci.org/cdw/multifil)

# multifil - a simulation of a half sarcomere

To install, in the top level directory of a local copy of this repository:

- run `pip install -e .` to install a symlinked copy that can be updated with a `git pull`
- or run `pip install .` to install a local copy that isn't symlinked. I recommend the former as this model is in heavy flux and so you'll probably want recent changes directly from the repo. 

## Starting a run

The simplest run would be:

``` python
import multifil
sarc = multifil.hs.hs()
axial_force = sarc.run()
```

More complicated runs (workloops, length-velocity calculations) require the modification of sarcomere parameters such as z-line-to-m-line distance and Calcium activation during the run. This is managed by the `aws.run.py` module through the reading of JSON formatted meta files. The metafiles are specified in `aws.metas.py` module. 

Creating a meta file for a workloop would look like:

``` python
# Imports
import numpy as np 
import multifil as mf
# Set run parameters
z_line_rest, z_line_amp, poisson_ratio = 1250, 50, 0.5 
freq, phase = 12, 0.8
act_time, act_rise, act_fall = 20, 3, 3
local_path, s3_path = './', None
# Create time, length, activation traces
time = np.arange(0,200,0.05) #in ms
z_line = mf.aws.metas.zline_workloop(z_line_rest, z_line_amp, freq, time)
activation = mf.aws.metas.actin_permissiveness_workloop(freq, phase, act_time, act_rise, act_fall, time)
# Emit metafile
meta = mf.aws.metas.emit(local_path, s3_path, time, poisson_ratio, z_line=z_line, actin_permissiveness=activation, comment="Example workloop run", phase=phase, frequency=freq)
```
