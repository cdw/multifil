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

More complicated runs (workloops, length-velocity calculations) require the modification of sarcomere parameters such as z-line-to-m-line distance and Calcium activation during the run. This is managed by the `aws.run.py` module through the reading of JSON formatted meta files. The metafiles are specified in the same. 