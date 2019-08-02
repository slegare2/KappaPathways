# KappaPathways

No setup script is provided yet.

To use KappaPathways, clone or download git repository.
Then, in any python script that would use KappaPathways,
add it to sys.path before importing like:

```python
import sys
sys.path.append('/home/user/path_to/KappaPathways')
import kappapathway
```

KappaPathways requires KaSim and KaFlow.
Graphviz and xdot are also recommended to visualize the results.

An example usage is provided in script run_kappapathways.py.
To try it, create a new working directory.
Copy run_kappapathways.py and ptyr-model2-act.ka to the new directory.
In run_kappapathways.py, change the paths to KappaPathways, KaSim and KaFlow.
Then do

```
python3 run_kappapathways.py
```

This should create a dot file of the pathway to the selected event of interest
and a directory containing the intermediary files used to build the pathway.

Pathways to different event of interest can be obtained by changing the value of
variable `eoi`.

