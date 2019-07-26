#! /usr/bin/python3

"""
Find the pathway of a given event of interest in a Kappa simulation.

# Example usage:

import kappapathways

# Set kappa model and event of interest.
kappa_model = "signaling-model.ka"
eoi = "EGFR(Y1092{p})"

# Simulation parameters.
simtime = 3600
simseed = 235866

# Run KappaPathways.
kappapathways.findpathway(kappa_model, eoi, simtime, simseed)
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

    def __init__(self, filename=None, eoi=None):
        """ Initialize class CausalGraph. """

        self.filename = filename
        self.eoi = eoi
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


    def build_dot_file(self, edge_labels=False):
        """ build a dot file of the CausalGraph. """

        self.get_max_rank()
        self.get_start_nodes()
        self.get_end_nodes()
        self.sequentialize_ids()
        dot_str = "digraph G{\n"
        if self.occurrence != None:
            dot_str += '  label="Occurrence = {}" '.format(self.occurrence)
            dot_str += "fontsize=28 ;\n"
        if self.prevcores != None:
            dot_str += '  prevcores="{}"\n'.format(self.prevcores)
        dot_str += '  labelloc="t" ;\n'
        dot_str += "  ranksep=.3 ;\n"
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
        minweight = min(all_weights)
        minpenwidth = 1
        maxpenwidth = 20
        for edge in self.edges:
            dot_str += ('"{}" -> "{}" '
                        .format(edge.source.nodeid, edge.target.nodeid))
            edge_color = "black"
            pensize = edge.weight/minweight * minpenwidth
            if pensize > maxpenwidth:
                pensize = maxpenwidth
                edge_color = "red"
            dot_str += "[penwidth={}".format(pensize)
            dot_str += ", color={}".format(edge_color)
            if edge_labels == True:
                dot_str += ', label="   {}"'.format(edge.weight)
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

def getcausalcores(kappamodel, eoi, kasimpath, kaflowpath, simtime=None, simseed=None):
    """ 
    Generate initial causal cores of given event of interest by running KaSim 
    and then KaFlow.
    """

    new_model = add_eoi(kappamodel, eoi)
    trace_path = run_kasim(new_model, kasimpath, simtime, simseed)
    run_kaflow(trace_path, eoi, kaflowpath)


def add_eoi(kappamodel, eoi):
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


def run_kasim(kappa_with_eoi, kasimpath, simtime=1000, simseed=None):
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


def run_kaflow(trace_path, eoi, kaflowpath):
   """ Run KaFlow on the trace containing the EOI. """

   subprocess.run(("{}".format(kaflowpath),
                   "-o", "{}/causalcore-".format(eoi),
                   "{}".format(trace_path)))

# ---------------- End of Causal Cores Generation Section  --------------------

# ==================== Causal Cores Merging Section ===========================

def mergecores(eoi, causalgraphs=None, writedots=True, rmprev=False,
               edge_labels=False):
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
        graph.build_dot_file(edge_labels)
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

def loopcores(eoi, ignorelist=None, causalgraphs=None, writedots=True,
              rmprev=False, writepremerge=False, edge_labels=False):
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
        core.build_dot_file(edge_labels)
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
                              edge_labels=edge_labels)
    for i in range(len(looped_paths)):
        looped_paths[i].filename = "evpath-{}.dot".format(i+1)
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

def mergepaths(eoi, causalgraphs=None, writedot=True, rmprev=False,
               edge_labels=False):
    """ Merge event paths into a single pathway. """

    if causalgraphs == None:
        path_files = get_dot_files(eoi, "evpath")
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

    pathway.update()
    rerank_nodes(pathway)
    pathway.occurrence = None
    pathway.build_dot_file(edge_labels)
    pathway.filename = "event_pathway.dot"
    if writedot == True:
        output_path1 = "{}/{}".format(eoi, pathway.filename)
        outfile1 = open(output_path1, "w")
        outfile1.write(pathway.dot_file)
        outfile1.close()
        output_path2 = "{}-{}.dot".format(pathway.filename[:-4], eoi)
        outfile2 = open(output_path2, "w")
        outfile2.write(pathway.dot_file)
        outfile2.close()
    if rmprev == True:
        if path_files == None:
            path_files = get_dot_files(eoi, "evpath")
        for path_file in path_files:
            file_path = "{}/{}".format(eoi, path_file)
            os.remove(file_path)

    return pathway

# ............... End of Event Paths Merging Section ..........................

# """"""""""""""" Species-Pathway Conversion Section """"""""""""""""""""""""""

def convertpathway(eoi, causalgraph, kappamodel):
    """
    Convert a CausalGraph where node are events to a pathway where nodes are
    species.
    """

    if not isinstance(causalgraph, CausalGraph):
        graph_file = causalgraph
        causalgraph = CausalGraph(graph_file, eoi)
    kappa_rules = get_kappa_rules(kappamodel)
    build_pathway(causalgraph)


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


def build_pathway(graph):
    """
    Put nodes for the agent and sites for every rule
    that switches a state.
    """

    species_pathway = CausalGraph()

    self.species_nodes = []
    self.species_edges = []
    current_nodes = self.start_nodes
    while len(current_nodes) > 0:
        next_nodes = []
        for node in current_nodes:
            rule_name = self.rule_names[node][1:-1]
            rule = ""
            if rule_name in self.kappa_rules.keys():
                rule = self.kappa_rules[rule_name]
            rule_agents = self.parse_rule(rule)
            modified_agents = []
            for agent in rule_agents:
                for site in agent["sites"].keys():
                    state = agent["sites"][site]["state"]
                    if state != None:
                        if "/" in state:
                            modified_agents.append(agent["agent_type"])
            if len(modified_agents) > 0:
                ranked_node = {"node_id": target,
                               "rule_name": self.rule_names[target],
                               "rank": self.node_ranks[target]}
                print("---", modified_agents)
            for i in range(len(rule)):
                if rule[i] == "{":
                    if rule[i+2] == "/":
                        mod_rule = True
                        break
            #if mod_rule == True:
            #    print(rule_name, "     ", rule)
            for edge in self.story_edges:
                if edge["source"] == node:
                    next_nodes.append(edge["target"])
        current_nodes = next_nodes


def parse_rule(self, rule):
    """ Create a dict for given rule. """

    a = rule.index("@")
    rate = rule[a+1:]
    agents_list = rule[:a-1].split(', ')
    parsed_agents = []
    for agent in agents_list:
        agent_dict = {}
        parenthesis = agent.index("(")
        agent_type = agent[:parenthesis]
        agent_dict["agent_type"] = agent_type
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

# """"""""""" End of Species-Pathway Conversion Section """""""""""""""""""""""

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


class KappaPathway:
    """ Produce the pathway to EOI from given Kappa simulation. """

    def __init__(self, eoi, kappa_path, sim_time=1000, out_period=1,
                 seed=None, write_looped=False,
                 node_weight=None, edge_weight="count"):
        """ Initialize KappaPathway class. """

        self.eoi = eoi
        self.kappa_path = kappa_path
        self.sim_time = sim_time
        self.out_period = out_period
        self.seed = seed
        self.write_looped = write_looped
        self.node_weigh = node_weight
        self.edge_weight = edge_weight
        self.kasim_path = "/home/slegare/.opam/4.07.0/bin/KaSim"
        self.kaflow_path = "/home/slegare/programfiles/KaFlow/KaFlow"
        # Run class methods.
        self.add_eoi()
        #self.run_kasim()
        #self.run_kaflow()
        self.get_stories()
        self.extract_linear_paths()
        self.push_path()
        self.path_edges()
        self.count_occur()
        self.sort_by_rank()
        self.build_new_dot()


    def add_eoi(self):
        """ Create a new Kappa model where the EOI is added. """

        last_dot = self.kappa_path.rfind(".")
        self.prefix = self.kappa_path[:last_dot]
        self.kappa_eoi_path = "{}-eoi.ka".format(self.prefix)
        shutil.copyfile(self.kappa_path, self.kappa_eoi_path)
        kappa_eoi_file = open(self.kappa_eoi_path, "a")
        kappa_eoi_file.write("%obs: '{}' |{}|\n".format(self.eoi, self.eoi))
        kappa_eoi_file.write("%mod: [true] do $TRACK '{}' [true];\n".format(self.eoi))
        # Also ask for DIN.
        kappa_eoi_file.write('\n%mod: [true] do $DIN "din.json" [true];\n')
        kappa_eoi_file.close()


    def run_kasim(self):
        """ Run simulation with added EOI. """

        output_path = "{}-eoi.csv".format(self.prefix)
        self.trace_path = "{}-eoi.json".format(self.prefix)
        command_line = ("{} -mode batch --no-log --no-log -u t -p {} -l {} "
                        "-i {} -o {} -trace {}"
                        .format(self.kasim_path, self.out_period,
                                self.sim_time, self.kappa_eoi_path,
                                output_path, self.trace_path))
        if self.seed != None:
            command_line += " -seed {}".format(self.seed)
        os.system(command_line)


    def run_kaflow(self):
        """ Run KaFlow on the trace containing the EOI. """

        command_line = ("{} -o story-'{}'- {}"
                        .format(self.kaflow_path, self.eoi, self.trace_path))
        os.system(command_line)


    def get_stories(self):
        """ 
        Get a full and a simplified version of all stories. The simplified
        stories are used to find the rank of rules in the pathway. The full
        stories are used to get edges and numbers of occurence.
        """

        self.get_stories_range()
        self.full_stories = {}
        self.simple_stories = {}
        for story_num in range(self.first_story, self.last_story+1):
        #for story_num in range(10, 10+1):
            story_name = "story-{}-{}".format(self.eoi, story_num)
            input_path = "{}.dot".format(story_name)
            story_obj = SimplifyStory(input_path)
            simple_dict = {}
            simple_dict["nodes"] = story_obj.cleared_nodes
            simple_dict["edges"] = story_obj.cleared_edges
            end_nodes = self.get_end_nodes(simple_dict)
            simple_dict["end_nodes"] = end_nodes
            simple_dict["rule_names"] = story_obj.rule_names
            simple_name = "simple{}".format(story_name)
            self.simple_stories[simple_name] = simple_dict
            full_dict = {}
            full_dict["nodes"] = story_obj.story_nodes
            full_dict["edges"] = story_obj.story_edges
            full_dict["rule_names"] = story_obj.rule_names
            self.full_stories[story_name] = full_dict
            # Optionally, write simplified stories to dot file for debugging.
            output_path = "{}.dot".format(simple_name)
            output_file = open(output_path, "w")
            output_file.write(story_obj.dot_file)
            output_file.close()


    def get_stories_range(self):
        """ Get the number of the first and last stories. """

        file_list = os.listdir(".")
        num_list = []
        for f in file_list:
            if "story" in f and ".dot" in f:
                if "looped" not in f and "simple" not in f:
                    if self.eoi in f:
                        last_dash = f.rfind("-")
                        period = f.rfind(".")
                        num = int(f[last_dash+1:period])
                        num_list.append(num)
        self.first_story = min(num_list)
        self.last_story = max(num_list)


    def get_end_nodes(self, sto_dict):
        """ Find end nodes as nodes that never appear as sources. """

        end_nodes = []
        for node in sto_dict["nodes"]:
            rule_is_source = False
            for edge in sto_dict["edges"]:
                if node["node_id"] == edge["source"]:
                    rule_is_source = True
                    break
            if rule_is_source == False:
                end_nodes.append(node)
        return end_nodes


    def extract_linear_paths(self):
        """
        Transform every story in possibly many linear paths. That is:

            1 5               1  1  5
            | |               |  |  |
            2 6               2  2  6
           /|/     becomes    |  |  |   
          7 3                 7  3  3
          | |                 |  |  |
          8 4                 8  4  4

        """

        self.linear_paths = []
        for story_name in self.simple_stories.keys():
            story = self.simple_stories[story_name]
            for end_point in story["end_nodes"]:
                paths = [[end_point]]
                top_reached = False
                while top_reached == False:
                    top_reached = True
                    for i in range(len(paths)):
                        if paths[i][-1] != "complete":
                            top_reached = False
                            sources = []
                            for edge in story["edges"]:
                                if edge["target"] == paths[i][-1]["node_id"]:
                                    for node in story["nodes"]:
                                        if node["node_id"] == edge["source"]:
                                            sources.append(node)
                            if len(sources) == 0:
                                paths[i].append("complete")
                            else:
                                for j in range(1, len(sources)):
                                    paths.append(paths[i].copy())
                                paths[i].append(sources[0])
                                for k in range(1, len(sources)):
                                    paths[-k].append(sources[k])
                for path in paths:
                    reverse_path = []
                    for i in range(len(path)-2, -1, -1):
                        reverse_path.append(path[i])
                    self.linear_paths.append(reverse_path)
        #for path in self.linear_paths:
        #    print("---")
        #    print(path)
            

    def push_path(self):
        """
        Increase the rank of rules in stories until every rule has the
        same rank in every story.
        """

        self.pathway_nodes = []
        self.pathway_node_ids = {}
        self.pathway_rule_names = {}
        # Initialize current_rules at rank 1.
        current_rank = 1
        current_rules = []
        for path in self.linear_paths:
            for node in path:
                if node["rank"] == current_rank:
                    if node["rule_name"] not in current_rules:
                        current_rules.append(node["rule_name"])
        # Main path pushing loop.
        node_id = 1
        while len(current_rules) > 0:
            for rule in current_rules:
                # Check if rule appears at higher rank in other paths.
                farther_node = False
                for k in range(len(self.linear_paths)):
                    path = self.linear_paths[k]
                    for node in path:
                        if node["rule_name"] == rule:
                            node_rank = node["rank"]
                            #print(current_rank, node_rank, node["rule_name"], k)
                            break
                    if node_rank > current_rank:
                        farther_node = True
                        break
                # If so, increment the rank of that rule and subsequent ones.
                if farther_node == True:
                    for i in range(len(self.linear_paths)):
                        rule_is_present = False
                        #print(len(self.linear_paths[i]))
                        for j in range(len(self.linear_paths[i])):
                            node = self.linear_paths[i][j]
                            if node["rule_name"] == rule:
                                rule_is_present = True
                            if rule_is_present == True:
                                #print(self.linear_paths[i][j])
                                #self.linear_paths[j][i]["rank"] += 1
                                self.linear_paths[i][j]["rank"] += 1
                # Otherwise, copy rule to pathway_nodes with current_rank.
                elif farther_node == False:
                    self.pathway_nodes.append({"node_id": node_id,
                                               "rule_name": rule,
                                               "rank": current_rank})
                    self.pathway_node_ids[rule] = node_id
                    self.pathway_rule_names[node_id] = rule
                    node_id += 1
            self.end_rank = current_rank
            current_rank += 1
            current_rules = []
            for path in self.linear_paths:
                for node in path:
                    if node["rank"] == current_rank:
                        if node["rule_name"] not in current_rules:
                            current_rules.append(node["rule_name"])
            #print(current_rank, current_rules)


    def path_edges(self):
        """ Get all possible edges from every story. """

        self.pathway_edges = []
        for story_name in self.full_stories.keys():
            story = self.full_stories[story_name]
            for edge in story["edges"]:
                source_rule = story["rule_names"][edge["source"]]
                target_rule = story["rule_names"][edge["target"]]
                if source_rule in self.pathway_node_ids.keys():
                    if target_rule in self.pathway_node_ids.keys():
                        pathway_source = self.pathway_node_ids[source_rule]
                        pathway_target = self.pathway_node_ids[target_rule]
                        new_edge = {"source": pathway_source,
                                    "target": pathway_target}
                        if new_edge not in self.pathway_edges:
                            self.pathway_edges.append(new_edge)


    def count_occur(self):
        """
        Count the number of occurences and uses of nodes and edges from all
        stories. Compute the following values as:
        rule_counts: For a given rule, number of different event_ids that were
        instances of that rule in all stories.
        rule_uses: For a given rule, number of times it occured in all stories.
        edge_counts: For a given edge between two rules, number of different
        pairs of event_ids that composed the edge from all stories.
        edge_uses: For a given edge between two rules, number of times that edge
        occured in all stories.
        rule_counts_tot: Number of times a rule occured in total in simulation
        computed from DIN file.
        """

        rule_occurences = {}
        edge_occurences = {}
        for story_name in self.full_stories.keys():
            story = self.full_stories[story_name]
            for node_id in story["nodes"]:
                node_rule = story["rule_names"][node_id]
                if node_rule not in rule_occurences:
                    rule_occurences[node_rule] = [node_id]
                else:
                    rule_occurences[node_rule].append(node_id)
            for edge in story["edges"]:
                edge_ids = "{}-{}".format(edge["source"], edge["target"])
                source_rule = story["rule_names"][edge["source"]]
                target_rule = story["rule_names"][edge["target"]]
                edge_rules = "{}-{}".format(source_rule, target_rule)
                if edge_rules not in edge_occurences:
                    edge_occurences[edge_rules] = [edge_ids]
                else:
                    edge_occurences[edge_rules].append(edge_ids)
        self.rule_counts = {}
        self.rule_uses = {}
        counts = []
        uses = []
        for node_rule in rule_occurences.keys():
            occurences = rule_occurences[node_rule]
            unique_ids = set(occurences)
            self.rule_counts[node_rule] = len(unique_ids)
            self.rule_uses[node_rule] = len(occurences)
            counts.append(len(unique_ids))
            uses.append(len(occurences))
        self.min_rule_count = min(counts)
        self.max_rule_count = max(counts)
        self.min_rule_use = min(uses)
        self.max_rule_use = max(uses)
        self.edge_counts = {}
        self.edge_uses = {}
        counts = []
        uses = []
        for edge_rules in edge_occurences.keys():
            occurences = edge_occurences[edge_rules]
            unique_ids = set(occurences)
            self.edge_counts[edge_rules] = len(unique_ids)
            self.edge_uses[edge_rules] = len(occurences)
            counts.append(len(unique_ids))
            uses.append(len(occurences))
        self.min_edge_count = min(counts)
        self.max_edge_count = max(counts)
        self.min_edge_use = min(uses)
        self.max_edge_use = max(uses)
        # Compute rule_counts_tot from DIN file.
        din_file = open("din.json")
        din_data = json.load(din_file)
        din_rules = din_data["din_rules"]
        din_hits = din_data["din_hits"]
        self.rule_counts_tot = {}
        counts = []
        for node in self.pathway_nodes:
            rule_name = node["rule_name"][1:-1]
            if rule_name in din_rules:
                rule_ind = din_rules.index(rule_name)
                rule_hits = din_hits[rule_ind]
                self.rule_counts_tot[node["rule_name"]] = int(rule_hits)
                counts.append(int(rule_hits))
        self.min_rule_count_tot = min(counts)
        self.max_rule_count_tot = max(counts)


    def sort_by_rank(self):
        """ Group story nodes by rank. """

        self.sorted_nodes = {}
        self.sorted_edges = []
        for rank in range(1, self.end_rank+1):
            node_list = []
            for node in self.pathway_nodes:
                if node["rank"] == rank:
                    node_list.append(node)
                    for edge in self.pathway_edges:
                        if edge["source"] == node["node_id"]:
                            self.sorted_edges.append(edge)
            rank_str = "{}".format(rank)
            self.sorted_nodes[rank_str] = node_list


    def build_new_dot(self):
        """ Build the dot file that contains the new story. """

        dot_str = "digraph G{\n"
        dot_str += "  ranksep=.3 ;\n"
        for rank in range(1, self.end_rank+1):
            rank_str = "{}".format(rank)
            dot_str += ('{{ rank = same ; "{}" [shape=plaintext] ;\n'
                        .format(rank_str))
            node_list = self.sorted_nodes[rank_str]
            for node in node_list:
                if node["rule_name"] == '"{}"'.format(self.eoi):
                    node_shape = "ellipse"
                    node_color = "indianred2"
                else:
                    node_shape = "invhouse"
                    # Compute node colors based on counts
                    count = self.rule_counts_tot[node["rule_name"]]
                    high = self.max_rule_count_tot
                    frac = count / high
                    #print("{:40} {}".format(node["rule_name"], count))
                    #node_color = '"0.55 {} 1.000"'.format(frac)
                    node_color = '"0.55 0.6 1.000"'.format(frac)
                dot_str += '"node{}" [label={}, '.format(node["node_id"],
                                                   node["rule_name"]) 
                dot_str += 'shape={}, style=filled, '.format(node_shape)
                dot_str += 'fillcolor={}] ;\n'.format(node_color)
            dot_str += "}\n"
        for rank in range(1, self.end_rank):
            rank_str = "{}".format(rank)
            next_rank = "{}".format(rank+1)
            dot_str += ('"{}" -> "{}" [style="invis"];\n'
                        .format(rank_str, next_rank))
        for edge in self.sorted_edges:
            source_rule = self.pathway_rule_names[edge["source"]]
            target_rule = self.pathway_rule_names[edge["target"]]
            edge_rules = "{}-{}".format(source_rule, target_rule)
            #count = self.edge_counts[edge_rules]
            #low = self.min_edge_count
            #high = self.max_edge_count
            count = self.edge_uses[edge_rules]
            low = self.min_edge_use
            high = self.max_edge_use
            minpen = 1
            maxpen = 10
            dpen = maxpen-minpen
            frac = (count-low) / (high-low)
            weight = minpen + dpen * frac
            pensize = "{}".format(weight)
            #pensize = "1"
            dot_str += ('"node{}" -> "node{}" [penwidth={}] ;\n'
                        .format(edge["source"], edge["target"], pensize))
        dot_str += "}"
        self.dot_file = dot_str


class LoopStory:
    """ Reorganize a story into loops, if loops are present. """

    def __init__(self, story_path):
        """ Initialize LoopStory class. """

        self.story_path = story_path
        self.story_file = open(self.story_path, "r").readlines()
        # Run class methods.
        self.ignored_rules()
        self.read_nodes_edges()
        self.get_inout_edges()
        self.find_start_nodes()
        self.order_story()
        self.remove_ignored()
        self.sort_by_rank()
        self.build_new_dot()


    def ignored_rules(self):
        """ Define a list of terms to ignore in stories. """

        self.ignore_terms = [" unbinds", " ina"]


    def read_nodes_edges(self):
        """
        Read rules (nodes) and causal links (edges) from input story.
        Intro nodes are ignored.
        """

        self.story_nodes = []
        self.rule_names = {}
        intro_nodes = []
        edges = []
        for line in self.story_file:
            if "shape=" in line:
                tokens = line.split()
                node_id = tokens[0]
                name_start = line.index("label=")+7
                name_end = line.index(",")-1
                rule_name = '"{}"'.format(line[name_start:name_end])
                if "Intro" not in rule_name:
                    self.story_nodes.append(node_id)
                    self.rule_names[node_id] = rule_name
                else:
                    intro_nodes.append(node_id)
            if "->" in line:
                tokens = line.split()
                source = tokens[0]
                target = tokens[2]
                if source not in intro_nodes:
                    edges.append({"source": source,
                                  "target": target})
        self.story_edges = []
        for edge in edges:
            self.story_edges.insert(0, edge)


    def get_inout_edges(self):
        """ Build a dictionary of the incoming edges of each node. """

        self.incoming_edges = {}
        for target in self.story_nodes:
            in_edges = []
            for edge in self.story_edges:
                if edge["target"] == target:
                    in_edges.append(edge)
            self.incoming_edges[target] = in_edges
        self.outgoing_edges = {}
        for source in self.story_nodes:
            out_edges = []
            for edge in self.story_edges:
                if edge["source"] == source:
                    out_edges.append(edge)
            self.outgoing_edges[source] = out_edges


    def find_start_nodes(self):
        """
        Find starting nodes as nodes that never appear as targets
        and end nodes as nodes that never appear as sources.
        """

        self.start_nodes = []
        for node in self.story_nodes:
            rule_is_target = False
            for edge in self.story_edges:
                if node == edge["target"]:
                    rule_is_target = True
                    break
            if rule_is_target == False:
                self.start_nodes.append(node)
        self.end_nodes = []
        for node in self.story_nodes:
            rule_is_source = False
            for edge in self.story_edges:
                if node == edge["source"]:
                    rule_is_source = True
                    break
            if rule_is_source == False:
                self.end_nodes.append(node)


    def order_story(self):
        """
        Assign rank to each rule and create a loop every time a rule
        appears more than once.
        """

        # Initialize new ordered story.
        self.ordered_nodes = []
        self.ordered_edges = []
        added_nodes = []
        sources = []
        self.node_ranks = {}
        for node in self.start_nodes:
            sources.append(node)
            added_nodes.append(node)
            self.node_ranks[node] = 1
            ranked_node = {"node_id": node,
                           "rule_name": self.rule_names[node],
                           "rank": 1}
            self.ordered_nodes.append(ranked_node)
        # Construct new ordered story.
        original_nodes = {}
        while len(sources) > 0:
            # 1) Get all edges that come out of current sources. 
            current_edges = []
            for source in sources:
                for edge in self.outgoing_edges[source]:
                    current_edges.append(edge)
            # 2) Get all possible targets.
            current_targets = []
            for edge in current_edges:
                if edge["target"] not in current_targets:
                    current_targets.append(edge["target"])
            # 3) Get targets for which all incoming edges are present.
            satisfied_targets = []
            for target in current_targets:
                edges_satisfied = True
                for needed_edge in self.incoming_edges[target]:
                    if needed_edge not in current_edges:
                        edges_satisfied = False
                if edges_satisfied == True:
                    satisfied_targets.append(target)
            # 4) Check if target was not already added.
            targets_to_add = []
            for target in satisfied_targets:
                if target not in added_nodes:
                    targets_to_add.append(target)
            # 5) Get edges that link a source to an added target).
            edges_to_add = []
            for edge in current_edges:
                if edge["target"] in targets_to_add:
                    edges_to_add.append(edge)
            # 6) Get sources for which all outgoing edges have been used.
            sources_to_remove = []
            for source in sources:
                available_edges = self.outgoing_edges[source]
                all_used = True
                for edge in available_edges:
                    if edge not in edges_to_add:
                        all_used = False
                        break
                if all_used == True:
                    sources_to_remove.append(source)
            # 7) Find the rank of the nodes to add and the original node_id
            #    for rules that are repeated.
            for target in targets_to_add:
                target_name = self.rule_names[target]
                already_there = False
                for node in added_nodes:
                    if self.rule_names[node] == target_name:
                        already_there = True
                        self.node_ranks[target] = self.node_ranks[node]
                        original_nodes[target] = node
                        break
                if already_there == False:
                    in_edges = self.incoming_edges[target]
                    source_ranks = []
                    for in_edge in in_edges:
                        prev_node = in_edge["source"]
                        source_ranks.append(self.node_ranks[prev_node])
                    rank = max(source_ranks) + 1
                    self.node_ranks[target] = rank
                    self.end_rank = rank
            # 8) Add new target nodes and edges to ordered story.
            for target in targets_to_add:
                if target not in original_nodes.keys():
                    ranked_node = {"node_id": target,
                                   "rule_name": self.rule_names[target],
                                   "rank": self.node_ranks[target]}
                    self.ordered_nodes.append(ranked_node)
                    added_nodes.append(target)
            for edge in edges_to_add:
                so = edge["source"]
                ta = edge["target"]
                if so in original_nodes.keys():
                    so = original_nodes[so]
                if ta in original_nodes.keys():
                    ta = original_nodes[ta]
                self.ordered_edges.append({"source": so,
                                           "target": ta})
            # 9) Update source nodes.
            new_sources = []
            for source in sources:
                if source not in sources_to_remove:
                    new_sources.append(source)
            for target in targets_to_add:
                new_sources.append(target)
            sources = new_sources.copy()


    def remove_ignored(self):
        """ Remove nodes that contain any ignored_term in their rule_name. """

        self.cleared_nodes = []
        self.cleared_edges = []
        for node in self.ordered_nodes:
            ignored = False
            for term in self.ignore_terms:
                if term in node["rule_name"]:
                    ignored = True
                    break
            if ignored == False:
                self.cleared_nodes.append(node)
                for edge in self.ordered_edges:
                    if edge["source"] == node["node_id"]:
                        ignored_edge = False
                        for term in self.ignore_terms:
                            if term in self.rule_names[edge["target"]]:
                                ignored_edge = True
                                break
                        if ignored_edge == False:
                            self.cleared_edges.append(edge)


    def sort_by_rank(self):
        """ Group story nodes by rank. """

        self.sorted_nodes = {}
        self.sorted_edges = []
        for rank in range(1, self.end_rank+1):
            node_list = []
            for node in self.cleared_nodes:
                if node["rank"] == rank:
                    node_list.append(node)
                    for edge in self.cleared_edges:
                        if edge["source"] == node["node_id"]:
                            self.sorted_edges.append(edge)
            rank_str = "{}".format(rank)
            self.sorted_nodes[rank_str] = node_list


    def build_new_dot(self):
        """ Build the dot file that contains the new story. """

        dot_str = "digraph G{\n"
        dot_str += "  ranksep=.3 ;\n"
        for rank in range(1, self.end_rank+1):
            rank_str = "{}".format(rank)
            dot_str += ('{{ rank = same ; "{}" [shape=plaintext] ;\n'
                        .format(rank_str))
            node_list = self.sorted_nodes[rank_str]
            for node in node_list:
                if node["node_id"] in self.end_nodes:
                    node_shape = "ellipse"
                    node_color = "indianred2"
                else:
                    node_shape = "invhouse"
                    node_color = "lightblue"
                dot_str += '"node{}" [label={}, '.format(node["node_id"],
                                                   node["rule_name"]) 
                dot_str += 'shape={}, style=filled, '.format(node_shape)
                dot_str += 'fillcolor={}] ;\n'.format(node_color)
            dot_str += "}\n"
        for rank in range(1, self.end_rank):
            rank_str = "{}".format(rank)
            next_rank = "{}".format(rank+1)
            dot_str += '"{}" -> "{}" [style="invis"];\n'.format(rank_str,
                                                              next_rank)
        for edge in self.sorted_edges:
            dot_str += '"node{}" -> "node{}" ;\n'.format(edge["source"], edge["target"])
        dot_str += "}"
        self.dot_file = dot_str


class SimplifyStory(LoopStory):
    """
    Create a simplified version of a story. It actually is the same as a
    looped story but where looping edges and duplicate edges were removed.
    """

    def __init__(self, story_path):
        """ Initialize SimplifyStory class. """

        self.story_path = story_path
        self.story_file = open(self.story_path, "r").readlines()
        # Run class methods.
        self.ignored_rules()
        self.read_nodes_edges()
        self.get_inout_edges()
        self.find_start_nodes()
        self.order_story()
        self.remove_ignored()
        self.simplify_story()
        self.sort_by_rank()
        self.build_new_dot()


    def simplify_story(self):
        """ Remove looping edges and duplicate edges. """

        tmp_edges = self.cleared_edges.copy()
        self.cleared_edges = []
        seen_ids = []
        for edge in tmp_edges:
            source_rank = self.node_ranks[edge["source"]]
            target_rank = self.node_ranks[edge["target"]]
            edge_id = "{}-{}".format(edge["source"], edge["target"])
            if target_rank > source_rank and edge_id not in seen_ids:
                self.cleared_edges.append(edge)
                seen_ids.append(edge_id)
        

class SpeciesPathway(LoopStory):
    """ 
    Convert rule-centric story (looped or not) into a species-centric pathway.
    Two different versions can be made, with explicit site or with sites 
    merged into the agent.
    """

    def __init__(self, story_path, kappa_path):
        """ Initialize SpeciesPathway class. """

        self.story_path = story_path
        self.story_file = open(self.story_path, "r").readlines()
        self.kappa_path = kappa_path
        self.kappa_file = open(self.kappa_path, "r").readlines()
        # Run class methods.
        self.ignored_rules()
        self.read_nodes_edges()
        self.get_inout_edges()
        self.find_start_nodes()
        self.get_kappa_rules()
        self.build_pathway()


    def get_kappa_rules(self):
        """ Build a dictionary of the rules from the input kappa model. """

        self.kappa_rules = {}
        for line in self.kappa_file:
            if line[0] == "'":
                quote = line.rfind("'")
                rule_name = line[1:quote]
                rule_strt = 0
                for i in range(quote+1, len(line)):
                    if line[i] != " ":
                        rule_strt = i
                        break
                rule = line[rule_strt:-1]
                self.kappa_rules[rule_name] = rule


    def build_pathway(self):
        """ 
        Put nodes for the agent and sites for every rule 
        that switches a state.
        """

        self.species_nodes = []
        self.species_edges = []
        current_nodes = self.start_nodes
        while len(current_nodes) > 0:
            next_nodes = []
            for node in current_nodes:
                rule_name = self.rule_names[node][1:-1]
                rule = ""
                if rule_name in self.kappa_rules.keys():
                    rule = self.kappa_rules[rule_name]
                rule_agents = self.parse_rule(rule)
                modified_agents = []
                for agent in rule_agents:
                    for site in agent["sites"].keys():
                        state = agent["sites"][site]["state"]
                        if state != None:
                            if "/" in state:
                                modified_agents.append(agent["agent_type"])
                if len(modified_agents) > 0:
                    ranked_node = {"node_id": target,
                                   "rule_name": self.rule_names[target],
                                   "rank": self.node_ranks[target]}
                    print("---", modified_agents)
                for i in range(len(rule)):
                    if rule[i] == "{":
                        if rule[i+2] == "/":
                            mod_rule = True
                            break
                #if mod_rule == True:
                #    print(rule_name, "     ", rule)
                for edge in self.story_edges:
                    if edge["source"] == node:
                        next_nodes.append(edge["target"])
            current_nodes = next_nodes


    def parse_rule(self, rule):
        """ Create a dict for given rule. """

        a = rule.index("@")
        rate = rule[a+1:]
        agents_list = rule[:a-1].split(', ')
        parsed_agents = []
        for agent in agents_list:
            agent_dict = {}
            parenthesis = agent.index("(")
            agent_type = agent[:parenthesis]
            agent_dict["agent_type"] = agent_type
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


## Use class LoopStory.
#for story_num in range(first_story, last_story+1):
#    input_path = "story-{}.dot".format(story_num)
#    loopedstory = LoopStory(input_path)
#    output_path = "loopedstory-{}.dot".format(story_num)
#    output_file = open(output_path, "w")
#    output_file.write(loopedstory.dot_file)


"""

# Output files.
output_file = open("kaflow-pathway.dot", "w")
count_file = open("rule_count.txt", "w")


# Count rules in stories.
rule_counts = {}
rule_names = []
node_ids = {}
for i in range(1, n_stories+1):
    file_name = "story{}.dot".format(i)
    input_file = open(file_name, "r").readlines()
    for line in input_file:
        if "shape=" in line:
            name_start = line.index("label=")+7
            name_end = line.index(",")-1
            rule = '"{}"'.format(line[name_start:name_end])
            tokens = line.split()
            node_id = tokens[0]
            if "Intro" not in rule:
                node_ids[node_id] = rule
                if rule in rule_counts.keys():
                    rule_counts[rule] += 1
                else:
                    rule_counts[rule] = 1
                    rule_names.append(rule)
# Write counted rules to file, sorted by occurence.
rule_list = []
for key in rule_counts.keys():
    rule_list.append([key, rule_counts[key]])
sorted_rules = sorted(rule_list, key=lambda x: x[1], reverse=True)
for rule in sorted_rules:
    count_file.write("{:40}{}\n".format(rule[0], rule[1]))


# Find the first rank of every rule in each KaFlow stories.
# stories_ranks looks like: 
# {'story1' : { 'FYN unbinds sh2 intra' : 1,
#               'FYN part act' : 2,
#               'FYNkin binds BCR-Y177 slow' : 3
#               ...
#             },
#  'story2' : { ...
#             },
#  ...
# }
stories_ranks = {}
#for i in range(1, n_stories+1):
for i in range(1, 11):
    file_name = "story{}.dot".format(i)
    story_name = "story{}".format(i)
    input_file = open(file_name, "r").readlines()
    # Get Intro and rule node_ids.
    intro_nodes = []
    rule_nodes = []
    for line in input_file:
        if "shape=" in line:
            tokens = line.split()
            if "Intro" in line:
                intro_nodes.append(tokens[0])
            else:
                rule_nodes.append(tokens[0])
    # Read edges.
    edges = []
    for line in input_file:
        if "->" in line:
            tokens = line.split()
            source = tokens[0]
            target = tokens[2]
            if source not in intro_nodes:
                edges.append([source, target])
    edges_revert = []
    for edge in edges:
        edges_revert.insert(0, edge)
    # Find starting nodes as nodes that never appear as target.
    sources = []
    for rule in rule_nodes:
        rule_is_target = False
        for edge in edges_revert:
            if rule == edge[1]:
                rule_is_target = True
                break
        if rule_is_target == False:
            sources.append(rule)
    # Rank rules.
    rank = 1
    event_ranks = {}
    for source in sources:
        event_ranks[source] = [rank, node_ids[source]]
    while len(sources) > 0:
        rank += 1
        targets = []
        for source in sources:
            for edge in edges_revert:
                if source == edge[0]:
                    targets.append(edge[1])
        for target in targets:
            event_ranks[target] = [rank, node_ids[target]]
        sources = targets.copy()
    # Get the first occurence of every rule.
    first_occur = {}
    for rule in rule_names:
        all_ranks = []
        for node_id in event_ranks.keys():
            event = event_ranks[node_id]
            rank = event[0]
            rule_name = event[1]
            if rule_name == rule:
                all_ranks.append(int(rank))
        if len(all_ranks) > 0:
           min_rank = min(all_ranks)
           first_occur[rule] = min_rank
    stories_ranks[story_name] = first_occur

# Find the rank of every rule whose occurence is above threshold.
# A rule's rank is given as the highest possible rank where it first occured
# in the story where it appeared the latest.
rule_ranks = {}
for rule in rule_names:
    include = True
    occurence = rule_counts[rule]
    if occurence < threshold:
        include = False
    for term in hide_terms:
        if term in rule:
            include = False
            break
    if include == True:
        first_seen_list = []
        fs_list = []
        for story in stories_ranks.keys():
            first_occur = stories_ranks[story]
            if rule in first_occur.keys():
                first_seen = first_occur[rule]
                first_seen_list.append(first_seen)
                fs_list.append([story, first_seen])
        if len(first_seen_list) > 0:
            max_first_seen = max(first_seen_list)
        else: max_first_seen = "Never"
        print(rule, "rank", max_first_seen)
        for sto in fs_list:
            print(sto)


## Select rules from full influence map.
#output_file.write("digraph G{\n")
## Put nodes.
#placed_nodes = []
#for line in influence_file:
#    if "[" in line and "->" not in line:
#        bracket = line.index("[")
#        node = line[:bracket-1]
#        if node in rule_names:
#            occurence = rule_counts[node]
#            if occurence >= threshold:
#                show = True
#                for term in hide_terms:
#                    if term in node:
#                        show = False
#                        break
#                if show == True:
#                    output_file.write(line)
#                    placed_nodes.append(node)
## Put edges.
#for line in influence_file:
#    if "[" in line and "->" in line:
#        # Remove edge label.
#        bracket = line.index("[")
#        color_ind = line.index("color")
#        after_label = line[color_ind:]
#        before_label = line[:bracket+1]
#        line_nolabel = before_label + after_label
#        # Get nodes.
#        arrow = line_nolabel.index("->")
#        node1 = line_nolabel[:arrow-1]
#        node2 = line_nolabel[arrow+3:bracket-1]
##        agent1 = get_agent(node1)
##        agent2 = get_agent(node2)
#        if node1 in placed_nodes and node2 in placed_nodes:
#            show = True
#            # Hide binding conflicts.
#            if 'arrowhead="tee"' in line_nolabel:
#                if " binds" in node1 and " binds" in node2:
#                    show = False
#            ## Hide influence between phos (unbind) and rebind.
#            #if "phosphorylates" in node1 and " binds" in node2:
#            #    agent1 = get_agent(node1)
#            #    agent2 = get_agent(node2)
#            #    target1 = node1[:-1].split()[2]
#            #    target2 = node2[:-1].split()[2]
#            #    print("++++")
#            #    print(line_nolabel[:-1])
#            #    print(agent1, target1, "|", agent2, target2)
#            #    if target1 == target2:
#            #        show = False
#            if show == True:
#                output_file.write(line_nolabel)
#output_file.write("}")
"""
