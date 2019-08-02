"""
Find the pathway to a given event of interest from a Kappa simulation
using KappaPathways.
"""

import sys
sys.path.append('/home/user/path_to/KappaPathways')
import kappapathways

# Set path to KaSim and KaFlow.
kasimpath = "/home/user/path_to/KaSim"
kaflowpath = "/home/user/path_to/KaFlow"

# Set event of interest.
eoi = "EGFR(Y1092{p})"

# Other examples of valid event of interest in model ptyr-model2-act.ka.
#eoi = "PXN(Y118{p})"
#eoi = "BCR(Y177{p})"
#eoi = "MET(Y1356{p})"
#eoi = "ZAP70(Y319{p})"
#eoi = "GRB2(sh2[_])"

# Set Kappa Model
kappamodel = "ptyr-model2-act.ka"

# Set list of strings to ignore in rule names.
ignorelist = [" unbinds", " ina", " dephos"]

# Simulation parameters.
simtime = 3600
simseed = None


# Run KappaPathways.
kappapathways.findpathway(eoi, kappamodel, kasimpath, kaflowpath,
                          simtime, simseed, ignorelist, edgelabels=False)

# Optionally, make png files from the dot graphs. Requires Graphviz.
#graphvizpath = "path_to/dot"
#kappapathways.drawpngs(eoi, graphvizpath)

