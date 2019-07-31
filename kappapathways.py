#! /usr/bin/python3

"""
Find the pathway of a given event of interest in a Kappa simulation.

# Example usage:

import kappapathways

# Set event of interest and kappa model.
eoi = "EGFR(Y1092{p})"
kappamodel = "model.ka"

# Set list of rule strings to ignore.
ignorelist = [" unbinds", " ina"]

# Simulation parameters.
simtime = 3600
simseed = 235866

# Set path to KaSim and KaFlow.
kasimpath = "/home/user/.opam/4.07.0/bin/KaSim"
kaflowpath = "/home/user/programfiles/KaFlow/KaFlow"

# Run KappaPathways.
kappapathways.findpathway(eoi, kappamodel, kasimpath, kaflowpath, 
                          simtime, simseed, ignorelist)
"""

import os
import shutil
import subprocess
import warnings
import json


class CausalNode(object):
    """ An individual event node to use in causal graphs. """

    def __init__(self, nodeid, label, rank=None, weight=None):
        """ Initialize class CausalNode. """

        self.nodeid = nodeid
        self.label = label
        self.rank = rank
        self.weight = weight
        self.check_types()


    def check_types(self):
        """ Check that CausalNode attributes have proper types. """

        if not isinstance(self.nodeid, str):
            raise TypeError("nodeid should be a string.")
        if not isinstance(self.label, str):
            raise TypeError("label should be a string.")
        if self.rank != None:
            if not isinstance(self.rank, int):
                raise TypeError("rank should be an integer.")
        if self.weight != None:
            if not isinstance(self.weight, int):
                if not isinstance(self.weight, float):
                    raise TypeError("weight should be an integer or float.")


    def __repr__(self):
        """ Representation of the CausalNode object. """

        res =  "Node "
        res += 'id: "{}",  label: "{}"'.format(self.nodeid, self.label)
        if self.rank != None:
            res += ",  rank: {}".format(self.rank)
        if self.weight != None:
            res += ",  weight: {}".format(self.weight) 

        return res


class CausalEdge(object):
    """ An individual causal relationship to use in causal graphs. """

    def __init__(self, source, target, weight=1):
        """ Initialize class CausalEdge. """

        self.source = source
        self.target = target
        self.weight = weight
        self.check_types()


    def check_types(self):
        """ Check that CausalEdge attributes have proper types. """

        if not isinstance(self.source, CausalNode):
            raise TypeError("source should be a CausalNode.")
        if not isinstance(self.target, CausalNode):
            raise TypeError("target should be a CausalNode.")
        if self.weight != None:
            if not isinstance(self.weight, int):
                if not isinstance(self.weight, float):
                    raise TypeError("weight should be an integer or float.")


    def __repr__(self):
        """ Representation of the CausalEdge object. """

        res =  "Edge "
        res += "source) {},  target) {}".format(self.source, self.target)
        if self.weight != None:
            res += ",  weight: {}".format(self.weight)

        return res


class CausalGraph(object):
    """ General data structure for causal graphs. """

    def __init__(self, filename=None, eoi=None, nodestype="event"):
        """ Initialize class CausalGraph. """

        self.filename = filename
        self.eoi = eoi
        self.nodestype = nodestype # event or species
        self.nodes = []
        self.edges = []
        self.occurrence = 1
        self.max_rank = None
        self.prevcores = None
        if self.filename != None:
            self.read_dot(self.filename)


    def read_dot(self, dotpath):
        """
        Read rules (nodes) and causal links (edges) from input causal core.
        """

        rank = None
        self.label_mapping = {}
        dotfile = open(dotpath, "r").readlines()
        for line in dotfile:
            if "nodestype=" in line:
                type_index = line.index("nodestype")
                quote = line.rfind('"')
                self.nodestype = line[type_index+11:quote]
            if "eoi=" in line:
                eoi_index = line.index("eoi")
                quote = line.rfind('"')
                self.eoi = line[eoi_index+5:quote]
            if "Occurrence" in line:
                occu = line.index("Occurrence")
                quote = line.rfind('"')
                self.occurrence = int(line[occu+12:quote])
            if "rank = same" in line:
                open_quote = line.index('"')
                close_quote = line.rfind('"')
                rank = int(line[open_quote+1:close_quote])
            if "label=" in line and "Occurrence" not in line:
                if "->" not in line:
                    tokens = line.split()
                    node_id = tokens[0]
                    if '"' in node_id:
                        node_id = node_id[1:-1]
                    if "node" not in node_id:
                        node_id = "node{}".format(node_id)
                    label_start = line.index("label=")+7
                    label_end = line.index(",")-1
                    label = "{}".format(line[label_start:label_end])
                    self.nodes.append(CausalNode(node_id, label, rank))
                    self.label_mapping[node_id] = label
        tmp_edges = []
        for line in dotfile:
            if "->" in line and 'style="invis"' not in line:
                tokens = line.split()
                source_id = tokens[0]
                if '"' in source_id:
                    source_id = source_id[1:-1]
                if "node" not in source_id:
                    source_id = "node{}".format(source_id)
                target_id = tokens[2]
                if '"' in target_id:
                    target_id = target_id[1:-1]
                if "node" not in target_id:
                    target_id = "node{}".format(target_id)
                for node in self.nodes:
                    if node.nodeid == source_id:
                        source = node
                    if node.nodeid == target_id:
                        target = node
                if "weight=" in line:
                    weight_start = line.index("weight=")+7
                    weight_end = line.index("]")
                    weight = int(line[weight_start:weight_end])
                else:
                    weight = 1
                tmp_edges.append(CausalEdge(source, target, weight))
        for edge in tmp_edges:
            self.edges.insert(0, edge)
        if rank == None:
            self.find_ranks()
        if self.eoi == None:
            self.get_max_rank()
            for node in self.nodes:
                if node.rank == self.max_rank:
                    self.eoi = node.label
        self.build_dot_file()


    def find_ranks(self):
        """ Find the rank of each node in an initial causal core. """

        self.get_start_nodes()
        allowed_nodes = self.start_nodes.copy()
        nodes_in_the_air = []
        for node in self.nodes:
            if node not in allowed_nodes:
                nodes_in_the_air.append(node)
        while len(nodes_in_the_air) > 0:
            placed_nodes = []
            for node in nodes_in_the_air:
                # Determine if node can be placed.
                can_place = True
                required_nodes = []
                for edge in self.edges:
                    if edge.target == node:
                       required_nodes.append(edge.source)
                for required in required_nodes:
                    if required not in allowed_nodes:
                        can_place = False
                        break
                if can_place == True:
                    # Find the rank of the required nodes and assign highest
                    # rank + 1 to the newly placed node
                    req_ranks = []
                    for required in required_nodes:
                        if required.rank != None:
                            req_ranks.append(required.rank)
                    if len(req_ranks) == 0:
                        node.rank = 2
                        for required in required_nodes:
                            required.rank = 1
                    else:
                        node.rank = max(req_ranks) + 1
                        for required in required_nodes:
                            if required.rank == None:
                                required.rank = max(req_ranks)
                    placed_nodes.append(node)
            for placed in placed_nodes:
                allowed_nodes.append(placed)
                nodes_in_the_air.remove(placed)


    def get_max_rank(self):
        """ Find the highest rank of a node in CausalGraph. """

        all_ranks = []
        for node in self.nodes:
            if node.rank != None:
                all_ranks.append(node.rank)
        if len(all_ranks) > 0:
            self.max_rank = max(all_ranks)


    def get_start_nodes(self):
        """
        Find start nodes. If nodes have ranks assigned, start nodes are nodes
        of rank one. If ranks are not assigned, starting nodes are nodes that
        never appear as targets. May fail if the graph contains loops while 
        having no ranks assigned, which should not occur in current workflow.
        """

        self.start_nodes = []
        for node in self.nodes:
            if node.rank == 1:
                self.start_nodes.append(node)
        if len(self.start_nodes) == 0:
            for node in self.nodes:
                rule_is_target = False
                for edge in self.edges:
                    if node == edge.target:
                        rule_is_target = True
                        break
                if rule_is_target == False:
                    self.start_nodes.append(node)


    def get_end_nodes(self):
        """
        Find end nodes. If nodes have ranks assigned, end nodes are nodes
        with max_rank. If ranks are not assigned,  end nodes are nodes that
        never appear as source. May fail if the graph contains loops while 
        having no ranks assigned, which should not occur in current workflow.
        """

        self.end_nodes = []
        if self.max_rank != None:
            for node in self.nodes:
                if node.rank == self.max_rank:
                    self.end_nodes.append(node)
        if len(self.end_nodes) == 0:
            for node in self.nodes:
                rule_is_source = False
                for edge in self.edges:
                    if node == edge.source:
                        rule_is_source = True
                        break
                if rule_is_source == False:
                    self.end_nodes.append(node)


    def sequentialize_ids(self):
        """ Assign sequential node ids, getting rid of event numbers. """

        node_number = 1
        for current_rank in range(1, self.max_rank+1):
            for node in self.nodes:
                if node.rank == current_rank:
                    node.nodeid = "node{}".format(node_number)
                    node_number += 1
        # Also sort edges.
        sorted_edges = sorted(self.edges, key=lambda x: x.source.rank)
        self.edges = sorted_edges


    def update(self):
        """ Update data after some changes were made on nodes and edges. """

        self.get_max_rank()
        self.get_start_nodes()
        self.get_end_nodes()
        self.sequentialize_ids()


    def build_dot_file(self, edgelabels=False):
        """ build a dot file of the CausalGraph. """

        self.get_max_rank()
        self.get_start_nodes()
        self.get_end_nodes()
        self.sequentialize_ids()
        dot_str = "digraph G{\n"
        dot_str += '  nodestype="{}"\n'.format(self.nodestype)
        if self.eoi != None:
            dot_str += '  eoi="{}"\n'.format(self.eoi)
        if self.occurrence != None:
            dot_str += '  label="Occurrence = {}" '.format(self.occurrence)
            dot_str += "fontsize=28 ;\n"
        if self.prevcores != None:
            dot_str += '  prevcores="{}"\n'.format(self.prevcores)
        dot_str += '  labelloc="t" ;\n'
        if "core" in self.filename:
            dot_str += "  ranksep=0.5 ;\n"
        else:
            dot_str += "  ranksep=1.0 ;\n"
        # Draw nodes.
        for current_rank in range(1, self.max_rank+1):
            rank_str = "{}".format(current_rank)
            dot_str += ('{{ rank = same ; "{}" [shape=plaintext] ;\n'
                        .format(rank_str))
            for node in self.nodes:
                if node.rank == current_rank:
                    node_shape = "invhouse"
                    node_color = "lightblue"
                    if "Intro" in node.label:
                        node_shape = "rectangle"
                        node_color = "white"
                    if node.label == self.eoi:
                        node_shape = "ellipse"
                        node_color = "indianred2"
                    if self.nodestype == "species":
                        node_shape = "ellipse"
                    dot_str += ('"{}" [label="{}", '
                                .format(node.nodeid, node.label))
                    dot_str += "shape={}, style=filled, ".format(node_shape)
                    dot_str += "fillcolor={}] ;\n".format(node_color)
            dot_str += "}\n"
        for rank in range(1, self.max_rank):
            rank_str = "{}".format(rank)
            next_rank = "{}".format(rank+1)
            dot_str += ('"{}" -> "{}" [style="invis"] ;\n'
                        .format(rank_str, next_rank))
        # Draw edges.
        all_weights = []
        for edge in self.edges:
            all_weights.append(edge.weight)
        if len(all_weights) > 0:
            minweight = min(all_weights)
        else:
            minweight = 1
        minpenwidth = 1
        maxpenwidth = 20
        for edge in self.edges:
            dot_str += ('"{}" -> "{}" '
                        .format(edge.source.nodeid, edge.target.nodeid))
            edge_color = "black"
            pensize = edge.weight/minweight * minpenwidth
            if pensize > maxpenwidth:
                pensize = maxpenwidth
            dot_str += "[penwidth={}".format(pensize)
            dot_str += ", color={}".format(edge_color)
            if edgelabels == True:
                dot_str += ', label="  {}"'.format(edge.weight)
            dot_str += ", weight={}] ;\n".format(edge.weight)
        dot_str += "}"
        self.dot_file = dot_str
           

    def __repr__(self):
        """ Representation of the CausalGraph object. """

        res = "CausalGraph"
        if self.filename != None:
            res += " from file {}".format(self.filename)
        res += " ,  Occurrence = {}\n".format(self.occurrence)
        res += "\n"
        #for node in self.nodes:
        #    res+="{}\n".format(node.__repr__())
        for edge in self.edges:
            res+="{}\n".format(edge.__repr__())

        return res

# -------------------- Causal Cores Generation Section ------------------------

def getcausalcores(eoi, kappamodel, kasimpath, kaflowpath, simtime=1000, simseed=None):
    """ 
    Generate initial causal cores of given event of interest by running KaSim 
    and then KaFlow.
    """

    new_model = add_eoi(eoi, kappamodel)
    trace_path = run_kasim(new_model, kasimpath, simtime, simseed)
    run_kaflow(eoi, trace_path, kaflowpath)


def add_eoi(eoi, kappamodel):
    """ Create a new Kappa model where the EOI is added. """

    if not os.path.exists(eoi):
            os.mkdir(eoi)
    last_dot = kappamodel.rfind(".")
    prefix = kappamodel[:last_dot]
    new_path = "{}/{}-eoi.ka".format(eoi, prefix)
    shutil.copyfile(kappamodel, new_path)
    new_file = open(new_path, "a")
    new_file.write("%obs: '{}' |{}|\n".format(eoi, eoi))
    new_file.write("%mod: [true] do $TRACK '{}' [true];\n".format(eoi))
    ## Ask for DIN.
    #din_file = "{}/din.json".format(eoi)
    #new_file.write('\n%mod: [true] do $DIN "{}" [true];\n'.format(din_file))
    new_file.close()
    
    return new_path


def run_kasim(kappa_with_eoi, kasimpath, simtime, simseed):
    """ Run simulation with added EOI to produce trace. """

    last_dot = kappa_with_eoi.rfind(".")
    prefix = kappa_with_eoi[:last_dot]
    output_path = "{}.csv".format(prefix)
    trace_path = "{}.json".format(prefix)
    if simtime <= 100:
        plt_period = 0.1
    else:
        plt_period = 1
    if simseed == None:
        subprocess.run(("{}".format(kasimpath),
                        "-mode", "batch", "--no-log", "-u", "t",
                        "-p", "{}".format(plt_period),
                        "-l", "{}".format(simtime),
                        "-i", "{}".format(kappa_with_eoi),
                        "-o", "{}".format(output_path),
                        "-trace", "{}".format(trace_path)))
    else:
        subprocess.run(("{}".format(kasimpath),
                        "-mode", "batch", "--no-log", "-u", "t",
                        "-p", "{}".format(plt_period),
                        "-l", "{}".format(simtime),
                        "-i", "{}".format(kappa_with_eoi),
                        "-o", "{}".format(output_path),
                        "-trace", "{}".format(trace_path),
                        "-seed", "{}".format(simseed))) 
    
    return trace_path


def run_kaflow(eoi, trace_path, kaflowpath):
   """ Run KaFlow on the trace containing the EOI. """

   subprocess.run(("{}".format(kaflowpath),
                   "-o", "{}/causalcore-".format(eoi),
                   "{}".format(trace_path)))

# ---------------- End of Causal Cores Generation Section  --------------------

# ==================== Causal Cores Merging Section ===========================

def mergecores(eoi, causalgraphs=None, edgelabels=False, writedots=True,
               rmprev=False, printmsg=True):
    """ Merge equivalent causal cores and count occurrence. """

    if causalgraphs == None:
        causal_core_files = get_dot_files(eoi, "causalcore")
        causal_cores = []
        for core_file in causal_core_files:
            core_path = "{}/{}".format(eoi, core_file)
            causal_cores.append(CausalGraph(core_path, eoi))
    else:
       causal_cores = causalgraphs
       causal_core_files = None
    #print("Evaluating possible merges among {} causal graphs."
    #      .format(len(causal_cores)))
    merged_cores = []
    while len(causal_cores) > 0:
        current_core = causal_cores[0]
        equivalent_list = [0]
        for i in range(1, len(causal_cores)):
            same_core, equi_edges = equivalent_graphs(current_core,
                                                      causal_cores[i])
            if same_core == True:
                equivalent_list.append(i)
                current_core.occurrence += causal_cores[i].occurrence
                for j in range(len(current_core.edges)):
                    equi_index = equi_edges[j]
                    w = causal_cores[i].edges[equi_index].weight
                    current_core.edges[j].weight += w
        prevcores = []
        for index in equivalent_list:
            file_name = causal_cores[index].filename
            dash = file_name.rfind("-")
            period = file_name.rfind(".")
            number = int(file_name[dash+1:period])
            prevcores.append(number)
        current_core.prevcores = prevcores
        merged_cores.append(current_core)
        for i in range(len(equivalent_list)-1, -1, -1):
            index = equivalent_list[i]
            del(causal_cores[index])
    sorted_cores = sorted(merged_cores, key=lambda x: x.occurrence,
                          reverse=True)
    for i in range(len(sorted_cores)):
        sorted_cores[i].filename = "core-{}.dot".format(i+1)
    for graph in sorted_cores:
        graph.sequentialize_ids()
        graph.build_dot_file(edgelabels)
    if writedots == True:
        for graph in sorted_cores:
            output_path = "{}/{}".format(eoi, graph.filename)
            outfile = open(output_path, "w")
            outfile.write(graph.dot_file)
            outfile.close()
    if rmprev == True:
        if causal_core_files == None:
            causal_core_files = get_dot_files(eoi, "causalcore")
        for core_file in causal_core_files:
            file_path = "{}/{}".format(eoi, core_file)
            os.remove(file_path)
    if printmsg == True:
        print("Merging equivalent causal cores, {} unique cores obtained."
              .format(len(sorted_cores)))

    return sorted_cores


def get_dot_files(eoi, prefix=None):
    """ Get the number of the first and last stories. """

    tmp_file_list = os.listdir("{}".format(eoi))
    file_list = []
    for file_name in tmp_file_list:
        if "dot" in file_name:
            if prefix == None:
                file_list.append(file_name)
            else:
                dash = file_name.rfind("-")
                if file_name[:dash] == prefix:
                    file_list.append(file_name)
    file_dicts = []
    for file_name in file_list:
        dash = file_name.rfind("-")
        period = file_name.rfind(".")
        number = int(file_name[dash+1:period])
        file_dicts.append({"file": file_name, "num": number})
    sorted_dicts = sorted(file_dicts, key=lambda x: x["num"])
    sorted_list = []
    for d in sorted_dicts:
        sorted_list.append(d["file"])

    return sorted_list


def equivalent_graphs(graph1, graph2):
    """
    Equivalent causal graphs have the same edges between the same rules, but
    the exact time of their events can differ.
    """

    equi_edges = []
    if graph1.max_rank == graph2.max_rank:
        graph2_indexes = list(range(len(graph2.edges)))
        all_edges_found = True
        for edge1 in graph1.edges:
            for i in graph2_indexes:
                edge2 = graph2.edges[i]
                same_edge = equivalent_edges(edge1, edge2)
                if same_edge == True:
                    equi_edges.append(i)
                    graph2_indexes.remove(i)
                    break
            if same_edge == False:
                all_edges_found = False
                break
        # All the edges from graph2 should have been used
        # at this point for both graphs to be equivalent.
        if all_edges_found == True:
            if len(graph2_indexes) > 0:
                equi_graphs = False
            else:
                equi_graphs = True
        else:
            equi_graphs = False
    else:
        equi_graphs = False

    return equi_graphs, equi_edges


def equivalent_edges(edge1, edge2):
    """ Find whether two edges are between the same rules at same rank. """

    same_source = False
    if edge1.source.label == edge2.source.label:
        if edge1.source.rank == edge2.source.rank:
            same_source = True
    same_target = False
    if edge1.target.label == edge2.target.label:
        if edge1.target.rank == edge2.target.rank:
            same_target = True
    if same_source == True and same_target == True:
        equi_edges = True
    else:
        equi_edges = False

    return equi_edges

# ================ End of Causal Cores Merging Section ========================

# +++++++++++++++++++++++ Cores Looping Section +++++++++++++++++++++++++++++++

def loopcores(eoi, causalgraphs=None, ignorelist=None, edgelabels=False,
              writedots=True, rmprev=False, writepremerge=False):
    """ Build looped event paths by merging identical nodes within cores. """

    if causalgraphs == None:
        core_files = get_dot_files(eoi, "core")
        cores = []
        for core_file in core_files:
            core_path = "{}/{}".format(eoi, core_file)
            cores.append(CausalGraph(core_path, eoi))
    else:
       cores = causalgraphs
       core_files = None
    for core in cores:
        remove_intro(core)
        remove_ignored(core, ignorelist)
        merge_same_labels(core)
        fuse_edges(core)
        rerank_nodes(core)
    for core in cores:
        core.build_dot_file(edgelabels)
    # Write looped cores before they are merged for debugging.
    if writepremerge == True:
        for i in range(len(cores)):
            cores[i].filename = "loopedcore-{}.dot".format(i+1)
        for graph in cores:
                output_path = "{}/{}".format(eoi, graph.filename)
                outfile = open(output_path, "w")
                outfile.write(graph.dot_file)
                outfile.close()
    looped_paths = mergecores(eoi, cores, writedots=False,
                              edgelabels=edgelabels, printmsg=False)
    for i in range(len(looped_paths)):
        looped_paths[i].filename = "eventpath-{}.dot".format(i+1)
    if writedots == True:
        for graph in looped_paths:
            output_path = "{}/{}".format(eoi, graph.filename)
            outfile = open(output_path, "w")
            outfile.write(graph.dot_file)
            outfile.close()
    if rmprev == True:
        if core_files == None:
            core_files = get_dot_files(eoi, "core")
        for core_file in core_files:
            file_path = "{}/{}".format(eoi, core_file)
            os.remove(file_path)
    print("Finding loops in cores.")
    print("Merging equivalent looped cores, {} event paths obtained."
          .format(len(looped_paths)))

    return looped_paths


def remove_intro(graph):
    """ Remove all Intro nodes. """

    intro_nodes = []
    for i in range(len(graph.nodes)):
        if "Intro" in graph.nodes[i].label:
            intro_nodes.insert(0, i)
    for i in intro_nodes:
        del(graph.nodes[i])
    intro_edges = []
    for i in range(len(graph.edges)):
        source = graph.edges[i].source.label
        target = graph.edges[i].target.label
        if "Intro" in source or "Intro" in target:
            intro_edges.insert(0, i)
    for i in intro_edges:
        del(graph.edges[i])
    for node in graph.nodes:
        node.rank = node.rank - 1
    graph.update()


def remove_ignored(graph, ignorelist):
    """ Remove nodes whose label contains a string defined in ignorelist. """

    ignored_nodes = []
    for i in range(len(graph.nodes)):
        if any(ignorestr in graph.nodes[i].label for ignorestr in ignorelist):
            ignored_nodes.insert(0, i)
    for i in ignored_nodes:
        del(graph.nodes[i])
    ignored_edges = []
    for i in range(len(graph.edges)):
        source = graph.edges[i].source.label
        target = graph.edges[i].target.label
        if any(ignorestr in source for ignorestr in ignorelist):
            ignored_edges.insert(0, i)
        if any(ignorestr in target for ignorestr in ignorelist):
            if i not in ignored_edges:
                ignored_edges.insert(0, i)
    for i in ignored_edges:
        del(graph.edges[i])
    graph.update()


def merge_same_labels(graph):
    """
    Merge every node with same label within a graph. Merging between two
    nodes is done by deleting the second node and changing any edge which was
    coming from or ging to the second node to now come from or go to the
    second node.
    """

    same_label_remains = True
    while same_label_remains == True:
        same_label_remains = False
        for node1 in graph.nodes:
            same_label_nodes = [node1]
            for node2 in graph.nodes:
                if node2 != node1:
                    if node1.label == node2.label:
                        same_label_nodes.append(node2)
            if len(same_label_nodes) > 1:
                same_label_remains = True
                merge_nodes(same_label_nodes, graph)
                break
    graph.update()
            

def merge_nodes(node_list, graph):
    """ Merge every node from the list onto the node of lowest rank. """

    ranks = []
    for node in node_list:
        ranks.append(node.rank)
    lowest_rank = min(ranks)
    for i in range(len(node_list)):
        if node_list[i].rank == lowest_rank:
            main_node = node_list[i]
            del(node_list[i])
            break
    for edge in graph.edges:
        if edge.source in node_list:
            edge.source = main_node
        if edge.target in node_list:
            edge.target = main_node
    for i in range(len(graph.nodes)-1, -1, -1):
        if graph.nodes[i] in node_list:
            del(graph.nodes[i])
            

def fuse_edges(graph):
    """ Remove duplicate edges between two same nodes but sum weights. """

    unique_edges = []
    for edge1 in graph.edges:
        source1 = edge1.source.nodeid
        target1 = edge1.target.nodeid
        new_edge = True
        for edge2 in unique_edges:
            source2 = edge2.source.nodeid
            target2 = edge2.target.nodeid
            if source1 == source2 and target1 == target2:
                new_edge = False
                break
        if new_edge == True:
            unique_edges.append(edge1)
    for unique_edge in unique_edges:
        unique_source = unique_edge.source.nodeid
        unique_target = unique_edge.target.nodeid
        w = 0
        for edge in graph.edges:
            source = edge.source.nodeid
            target = edge.target.nodeid
            if unique_source == source and unique_target == target:
                w += edge.weight
        unique_edge.weight = w
    graph.edges = unique_edges
    graph.update()


def rerank_nodes(graph):
    """
    Adjust the ranks of every node in a graph. The rank of each node is given
    by the longest loopless path to a node of rank 1.
    """

    for node in graph.nodes:
        paths = climb_up(node, graph.start_nodes, graph)
        lengths = []
        for path in paths:
            lengths.append(len(path))
        node.rank = max(lengths)
    graph.update()


def climb_up(bottom, top, graph):
    """
    Return the list of possible loopless paths from bottom to top nodes.
    Edges are followed in reverse, from target to source.
    """

    if isinstance(top, list) == False:
        top = [top]
    all_paths = [[bottom]]
    top_reached = False
    while top_reached == False:
        top_reached = True
        for i in range(len(all_paths)):
            path = all_paths[i]
            if path[-1] not in top and path[-1] != "root":
                top_reached = False
                next_nodes = []
                for edge in graph.edges:
                    if edge.target == path[-1]:
                        next_nodes.append(edge.source)
                if len(next_nodes) == 0:
                    path.append("root")
                    warnings.warn("Root of graph reached before top "
                                  "nodes were found.")
                else:
                    path_copy = path.copy()
                    path.append(next_nodes[0])
                    for i in range(1, len(next_nodes)):
                        new_path = path_copy.copy()
                        new_path.append(next_nodes[i])
                        all_paths.append(new_path)
        # Remove looping paths.
        for i in range(len(all_paths)-1, -1, -1):
            if len(all_paths[i]) != len(set(all_paths[i])):
                del(all_paths[i])
    for i in range(len(all_paths)-1, -1, -1):
        if all_paths[i][-1] not in top:
            del(all_paths[i])
        
    return all_paths

# ++++++++++++++++++ End of Cores Looping Section +++++++++++++++++++++++++++++

# .................. Event Paths Merging Section ..............................

def mergepaths(eoi, causalgraphs=None, edgelabels=False, writedot=True,
               rmprev=False):
    """ Merge event paths into a single pathway. """

    if causalgraphs == None:
        path_files = get_dot_files(eoi, "eventpath")
        event_paths = []
        for path_file in path_files:
            path_path = "{}/{}".format(eoi, path_file)
            event_paths.append(CausalGraph(path_path, eoi))
    else:
       event_paths = causalgraphs
       path_files = None
    pathway = CausalGraph(eoi=eoi)
    node_number = 1
    seen_labels = []
    for event_path in event_paths:
        for node in event_path.nodes:
            if node.label not in seen_labels:
                seen_labels.append(node.label)
                n_id = "node{}".format(node_number)
                pathway.nodes.append(CausalNode(n_id, node.label, node.rank))
                node_number += 1
    for event_path in event_paths:
        for edge in event_path.edges:
            for node in pathway.nodes:
                if node.label == edge.source.label:
                    source = node
                if node.label == edge.target.label:
                    target = node
            pathway.edges.append(CausalEdge(source, target, edge.weight))
    fuse_edges(pathway)
    rerank_nodes(pathway) 
    pathway.occurrence = None
    pathway.build_dot_file(edgelabels)
    pathway.filename = "eventpathway.dot"
    if writedot == True:
        output_path1 = "{}/{}".format(eoi, pathway.filename)
        outfile1 = open(output_path1, "w")
        outfile1.write(pathway.dot_file)
        outfile1.close()
    if rmprev == True:
        if path_files == None:
            path_files = get_dot_files(eoi, "evpath")
        for path_file in path_files:
            file_path = "{}/{}".format(eoi, path_file)
            os.remove(file_path)
    print("Merging all event paths into one event pathway.")

    return pathway        

# ............... End of Event Paths Merging Section ..........................

# """"""""""""""" Species Pathway Conversion Section """"""""""""""""""""""""""

def speciespathway(eoi, kappamodel, causalgraph=None, edgelabels=False):
    """
    Convert a CausalGraph where node are events to a pathway where nodes are
    species.
    """

    if causalgraph == None:
        graph_path = "{}/eventpathway.dot".format(eoi)
        pathway = CausalGraph(graph_path, eoi)
    elif isinstance(causalgraph, CausalGraph):
        pathway = causalgraph
    else:
        pathway = CausalGraph(causalgraph, eoi)
    kappa_rules = get_kappa_rules(kappamodel)
    kappa_rules["{}".format(eoi)] = "|{}|".format(eoi)
    mod_nodes = get_mod_nodes(pathway, kappa_rules)
    added_nodes = add_first_nodes(pathway, kappa_rules, mod_nodes)
    for node in added_nodes:
        mod_nodes.append(node)
    for node in pathway.nodes:
        if node in mod_nodes:
            node.label = node.species
    rebranch(pathway, mod_nodes)
    merge_same_labels(pathway)
    fuse_edges(pathway)
    rerank_nodes(pathway)
    pathway.occurrence = None
    pathway.filename = "pathway.dot"
    pathway.nodestype = "species"
    pathway.build_dot_file(edgelabels)
    output_path1 = "{}/{}".format(eoi, pathway.filename)
    outfile1 = open(output_path1, "w")
    outfile1.write(pathway.dot_file)
    outfile1.close()
    pathway.filename = "{}-{}.dot".format(pathway.filename[:-4], eoi)
    output_path2 = "{}".format(pathway.filename)
    outfile2 = open(output_path2, "w")
    outfile2.write(pathway.dot_file)
    outfile2.close()
    print("Converting event pathway into species pathway.")
    print("File {} created.".format(pathway.filename))

    return pathway

    
def get_kappa_rules(kappamodel):
    """ Build a dictionary of the rules from the input kappa model. """
    
    kappa_file = open(kappamodel, "r").readlines()
    kappa_rules = {}
    for line in kappa_file:
        if line[0] == "'":
            quote = line.rfind("'")
            rule_name = line[1:quote]
            rule_strt = 0
            for i in range(quote+1, len(line)):
                if line[i] != " ":
                    rule_strt = i
                    break
            rule = line[rule_strt:-1]
            kappa_rules[rule_name] = rule

    return kappa_rules


def parse_rule(rule):
    """
    Create a dict for given rule.
    agent_dict = {"type": X ,
                  "sites": {"site1": {"binding": "1", "state": "p"},
                            "site2": {"binding": "3", "state": "u"}
                           }
                 }
    """

    parsed_agents = []
    if "@" in rule:
        a = rule.index("@")
        rate = rule[a+1:]
        agents_list = rule[:a-1].split(', ')
    elif "|" in rule:
        agents_list = rule[1:-1].split(', ')
    else:
        agents_list = rule.split(', ')
    for agent in agents_list:
        agent_dict = {}
        parenthesis = agent.index("(")
        agent_type = agent[:parenthesis]
        agent_dict["type"] = agent_type
        sites = agent[parenthesis+1:-1].split()
        site_dict = {}
        for site in sites:
            if "[" in site:
                open_bracket = site.index("[")
                close_bracket = site.index("]")
                site_id = site[:open_bracket]
                binding = site[open_bracket+1:close_bracket]
            else:
                binding = None
            if "{" in site:
                open_curl = site.index("{")
                close_curl = site.index("}")
                if "[" not in site:
                    site_id = site[:open_curl]
                state = site[open_curl+1:close_curl]
            else:
                state = None
            site_dict[site_id] = {"binding": binding, "state": state}
        agent_dict["sites"] = site_dict
        parsed_agents.append(agent_dict)

    return parsed_agents


def get_mod_nodes(graph, kappa_rules):
    """
    Create a list of the nodes that correspond to a rule where a state is
    modified. Also change the label of those nodes to the species they produce.
    """

    modification_rules = []
    for node in graph.nodes:
        rule = kappa_rules[node.label]
        rule_agents = parse_rule(rule)
        modified_agents = []
        for agent in rule_agents:
            modified_agent = {}
            modified_sites = {}
            sites = agent["sites"]
            for site in sites.keys():
                state = sites[site]["state"]
                if state != None:
                    if "/" in state:
                        slash = state.index("/")
                        final_state = state[slash+1:]
                        modified_sites[site] = {"state": final_state}
            if len(modified_sites.keys()) > 0:
                modified_agent["type"] = agent["type"]
                modified_agent["sites"] = modified_sites
                modified_agents.append(modified_agent)
        if len(modified_agents) > 0:
            species, kappa_species = label_species(modified_agents)
            node.species = species
            if graph.eoi == kappa_species:
                graph.eoi = species
            #node.species = kappa_species
            modification_rules.append(node)

    return modification_rules


def label_species(agent_list):
    """ Write species string. """

    species = ""
    for i in range(len(agent_list)):
        agent = agent_list[i]
        if i > 0:
            species += ", "
        species += "{}".format(agent["type"])
        sites = agent["sites"]
        for site in sites.keys():
            if "act" not in site:
                species += "-{}".format(site)
                #species += " {}".format(sites[site]["state"])

    kappa_species = ""
    for i in range(len(agent_list)):
        agent = agent_list[i]
        if i > 0:
            kappa_species += ", "
        kappa_species += "{}(".format(agent["type"])
        sites = agent["sites"]
        j = 0
        for site in sites.keys():
            if j > 0:
                species += " "
            kappa_species += "{}".format(site)
            kappa_species += "{{{}}}".format(sites[site]["state"])
            j += 1
        kappa_species += ")"

    return species, kappa_species


def add_first_nodes(graph, kappa_rules, mod_nodes):
    """
    Add the first node to each path. The label of the first node of a given
    path corresponds to all sites that were seen in the nodes before the first
    mod_node, except the modified site itself.
    This part will bug if complex stuff occurs. Will need to improve it a lot.
    """

    added_nodes = []
    all_paths = []
    seen_agents = []
    modified_agents = []
    path_weights = []
    for start_node in graph.start_nodes:
        all_paths.append([start_node])
        seen_agents.append([])
        path_weights.append(0)
    all_complete = False
    while all_complete == False:
        all_complete = True
        for i in range(len(all_paths)):
            current_node = all_paths[i][-1]
            if current_node != "mod_reached":
                all_complete = False
                rule = kappa_rules[current_node.label]
                rule_agents = parse_rule(rule)
                for rule_agent in rule_agents:
                    seen_agents[i].append(rule_agent)
                if current_node in mod_nodes:
                    modified_agents.append(rule_agents)
                    all_paths[i].append("mod_reached")
                else:
                    for edge in graph.edges:
                        if edge.source == current_node:
                            all_paths[i].append(edge.target)
                            path_weights[i] = edge.weight
    # Now sort through the seen_agents to find which one is assumed to
    # be responsible for the modification of the modified_agents.
    new_node_id = 1
    for i in range(len(all_paths)):
        modified_node = all_paths[i][-2]
        modified_agent = modified_agents[i]
        modified_types = []
        for mod_ag in modified_agent:
            modified_types.append(mod_ag["type"])
        all_agents = seen_agents[i]
        required_agents = []
        for agent in all_agents:
            if agent["type"] not in modified_types:             
                if agent["type"] not in required_agents:
                    required_agents.append(agent["type"])
        if len(required_agents) == 0:
            new_edge = CausalEdge(modified_node, modified_node,
                                  path_weights[i])
            graph.edges.append(new_edge)
        else: 
            new_str = "added{}".format(new_node_id)
            new_label = ""
            for i in range(len(required_agents)):
                if i > 0:
                    new_label += ", "
                new_label += "{}".format(required_agents[i])
            new_node = CausalNode(new_str, new_label)
            graph.nodes.append(new_node)
            added_nodes.append(new_node)
            new_edge = CausalEdge(new_node, modified_node, path_weights[i])
            graph.edges.append(new_edge)
        modified_node.rank = 1
            
    return added_nodes
    

def rebranch(graph, mod_nodes):
    """
    Remove non-modification nodes and rebranch their upstream nodes to their
    downstream nodes.
    """

    for j in range(len(graph.nodes)-1, -1, -1):
        node = graph.nodes[j]
        if node not in mod_nodes:
            up_edges = []
            for i in range(len(graph.edges)-1, -1, -1):
                edge = graph.edges[i]
                if edge.target == node:
                    up_edges.append(edge)
                    del(graph.edges[i])
            down_edges = []
            for i in range(len(graph.edges)-1, -1, -1):
                edge = graph.edges[i]
                if edge.source == node:
                    down_edges.append(edge)
                    del(graph.edges[i])
            for up_edge in up_edges:
                for down_edge in down_edges:
                    new_edge = CausalEdge(up_edge.source, down_edge.target,
                                          up_edge.weight)
                    graph.edges.append(new_edge)
            del(graph.nodes[j])
    graph.update()

# """"""""""" End of Species Pathway Conversion Section """""""""""""""""""""""

def findpathway(eoi, kappamodel, kasimpath, kaflowpath, simtime=1000,
                simseed=None, ignorelist=None, edgelabels=False):
    """ Use KappaPathways functions to get the pathway to given EOI. """

    # 1) Get Causal Cores.
    getcausalcores(eoi, kappamodel, kasimpath, kaflowpath, simtime, simseed)

    # 2) Merge Equivalent Causal Cores.
    mergedcores = mergecores(eoi, edgelabels=edgelabels, writedots=True,
                             rmprev=True)

    # 3) Loop causal cores to create event paths.
    loopedcores = loopcores(eoi, causalgraphs=mergedcores,
                            ignorelist=ignorelist, edgelabels=edgelabels,
                            writedots=True, rmprev=False, writepremerge=False)
    #loopedcores = loopcores(eoi,
    #                        ignorelist=ignorelist, edgelabels=edgelabels,
    #                        writedots=True, rmprev=False, writepremerge=False)

    # 4) Merge all event paths into one event pathway.
    eventpath = mergepaths(eoi, causalgraphs=loopedcores,
                           edgelabels=edgelabels, writedot=True, rmprev=False)

    # 5) Convert an event pathway into a species pathway.
    speciespathway(eoi, kappamodel, causalgraph=eventpath,
                   edgelabels=edgelabels)



def drawpngs(eoi, graphvizpath):
    """ Draw a png for every dot file found with given event of interest. """

    dot_files = get_dot_files(eoi)
    print("Drawing {} graphs.".format(len(dot_files)))
    for dot_file in dot_files:
        period = dot_file.rfind(".")
        file_path = "{}/{}".format(eoi, dot_file[:period])
        try:
            subprocess.run(("/usr/bin/dot", "-Tpng",
                            "{}.dot".format(file_path),
                            "-o", "{}.png".format(file_path)))
        except KeyboardInterrupt:
            exit(0)


def togglelabels():
    """
    Add or remove edge labels on all dot files found in current directory
    and all other directories found in current directory.
    """

    current_files = os.listdir(".")
    for current_file in current_files:
        if "dot" in current_file:
            toggle_edge_labels(".", current_file)
    directories = filter(os.path.isdir, os.listdir('.'))
    for directory in directories:
        dot_files = os.listdir("{}".format(directory))
        for dot_file in dot_files:
            if "dot" in dot_file:
                toggle_edge_labels(directory, dot_file)


def toggle_edge_labels(dir_path, dot_file):
    """ Add or remove edge labels in dot file. """

    input_path = "{}/{}".format(dir_path, dot_file)
    input_file = open(input_path, "r").readlines()
    new_file = ""
    for line in input_file:
        if "->" in line and "penwidth" in line:
            if "label" not in line: # Add labels.
                wpos = line.index("weight=")
                bracket = line.index("]")
                weight = line[wpos+7:bracket]
                comma = line.rfind(",")
                new_line = line[:comma+2]
                new_line += 'label="  {}", '.format(weight)
                new_line += line[wpos:]
                new_file += new_line
            else: # Remove labels.
                labelpos = line.index("label=")
                comma = line.rfind(",")
                new_line = line[:labelpos] + line[comma+1:]
                new_file += new_line
        else:
            new_file += line

    period = dot_file.rfind(".")
    prefix = dot_file[:period]
    save_path = "{}/{}-save.dot".format(dir_path, prefix)
    shutil.copyfile(input_path, save_path)
    os.remove(input_path)
    output_file = open(input_path, "w")
    output_file.write(new_file)
    os.remove(save_path)

