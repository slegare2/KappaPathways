# KappaPathways

No setup script is provided yet.

To use KappaPathways, clone or download git repository.
Then, in any python script that would use KappaPathways,
add it to sys.path before importing:

```python
import sys
sys.path.append('/home/user/path_to/KappaPathways')
import kappapathways
```

KappaPathways requires KaSim [https://github.com/Kappa-Dev/KaSim] and KaFlow [https://github.com/jonathan-laurent/KaFlow].
To visualize the results, Graphviz [https://www.graphviz.org] and xdot [https://pypi.org/project/xdot] are also recommended.

An example usage is provided in script *run_kappapathways.py*.
To try it, create a new working directory.
Copy *run_kappapathways.py* and *ptyr-model2-act.ka* to the new directory.
Inside *run_kappapathways.py*, change the paths to KappaPathways, KaSim and KaFlow.

```python
# Set path to KaSim and KaFlow.
kasimpath = "/home/user/path_to/KaSim"
kaflowpath = "/home/user/path_to/KaFlow"
```

Then do

```
python3 run_kappapathways.py
```

This should create a dot file of the pathway to the selected event of interest
and a directory containing the intermediary files used to build the pathway.

Pathways to different events of interest can be obtained by changing the value of
variable `eoi` in *run_kappapathways.py*.

