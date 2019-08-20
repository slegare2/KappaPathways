#! /usr/bin/python3

"""
Module KappaPathways.
Obtain the pathway to a given event of interest from a Kappa simulation.
"""

import os
import shutil
import subprocess
import warnings
import json
import math
import statistics


class CausalNode(object):
    """ An individual event node to use in causal graphs. """

    def __init__(self, nodeid, label, rank=None, weight=None, intro=False):
        """ Initialize class CausalNode. """

        self.nodeid = nodeid
        self.label = label
        self.rank = rank
        self.weight = weight
        self.intro = intro
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
        res += ",  intro: {}".format(self.intro) 

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
        self.maxrank = None
        self.prevcores = None
        if self.filename != None:
            self.read_dot(self.filename)


    def read_dot(self, dotpath):
        """
        Read rules (nodes) and causal links (edges) from input causal core.
        """

        rank = None
        intros = False
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
                quote = line[occu:].index('"')+occu
                self.occurrence = int(line[occu+12:quote])
            if "maxrank=" in line:
                maxrank_index = line.index("maxrank")
                quote = line.rfind('"')
                self.maxrank = int(line[maxrank_index+9:quote])
            if "rank = same" in line:
                open_quote = line.index('"')
                close_quote = line.rfind('"')
                rank = int(line[open_quote+1:close_quote])
            if "label=" in line and "Occurrence" not in line:
                if "->" not in line:
                    if line[0:2] == "//":
                       read_line = line[2:]
                    else:
                       read_line = line
                    tokens = read_line.split()
                    node_id = tokens[0]
                    if '"' in node_id:
                        node_id = node_id[1:-1]
                    if "node" not in node_id:
                        node_id = "node{}".format(node_id)
                    label_start = read_line.index("label=")+7
                    label_end = read_line[label_start:].index('"')+label_start
                    label = "{}".format(read_line[label_start:label_end])
                    #if "," in label:
                    #    comma = label.index(",")
                    #    label = label[:comma]
                    if "intro=" in read_line:
                        intros = True
                        intro_start = read_line.index("intro=")+6
                        intro_end = read_line[intro_start:].index("e")+intro_start+1
                        intro_str = read_line[intro_start:intro_end]
                        if intro_str == "True":
                            is_intro = True
                        elif intro_str == "False":
                            is_intro = False
                    else:
                        is_intro = False
                    self.nodes.append(CausalNode(node_id, label, rank,
                                                 intro=is_intro))
                    self.label_mapping[node_id] = label
        tmp_edges = []
        for line in dotfile:
            if "->" in line and 'style="invis"' not in line:
                if line[0:2] == "//":
                    read_line = line[2:]
                else:
                    read_line = line
                tokens = read_line.split()
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
                    weight_start = read_line.index("weight=")+7
                    weight_end = read_line.index("]")
                    weight = int(read_line[weight_start:weight_end])
                else:
                    weight = 1
                tmp_edges.append(CausalEdge(source, target, weight))
        for edge in tmp_edges:
            self.edges.insert(0, edge)
        if intros == False:
            for node in self.nodes:
                if "Intro" in node.label:
                    node.intro = True
        if rank == None:
            self.rank_nodes()
        if self.eoi == None:
            self.get_maxrank()
            for node in self.nodes:
                if node.rank == self.maxrank:
                    self.eoi = node.label


    def rank_nodes(self):
        """
        Find the rank of each rule node as the longest acyclic path to an
        initial node (intro node of rank 0). Then, the rank of each intro
        node is the lowest rank of any rule node it points to minus one.
        """

        self.find_init_nodes()
        for node in self.nodes:
            if node.intro == False:
                paths = self.climb_up(node, self.init_nodes)
                lengths = []
                for path in paths:
                    lengths.append(len(path)-1)
                node.rank = max(lengths)
        for node in self.nodes:
            if node.intro == True:
                target_ranks = []
                for edge in self.edges:
                    if edge.source == node:
                        target_ranks.append(edge.target.rank)
                node.rank = min(target_ranks)-1
        self.get_maxrank()
        self.sequentialize_ids()


    def find_init_nodes(self):
        """
        Find the nodes of rank 0 as the intro nodes which point to any rule
        node that requires only intro nodes.
        """

        self.init_nodes = []
        first_rules = []
        for node in self.nodes:
            if node.intro == False:
                input_nodes = []
                for edge in self.edges:
                    if edge.target == node:
                        input_nodes.append(edge.source)
                all_intro = True
                for input_node in input_nodes:
                    if input_node.intro == False:
                        all_intro = False
                        break
                if all_intro == True:
                    first_rules.append(node)
        for node in first_rules:
            for edge in self.edges:
                if edge.target == node:
                    if edge.source not in self.init_nodes:
                        self.init_nodes.append(edge.source)


    def climb_up(self, bottom_node, top_nodes):
        """
        Return the list of possible acyclic paths from bottom to top nodes.
        Edges are followed in reverse, from target to source.
        """
    
        all_paths = [[bottom_node]]
        top_reached = False
        while top_reached == False:
            top_reached = True
            for i in range(len(all_paths)):
                path = all_paths[i]
                if path[-1] not in top_nodes:
                    top_reached = False
                    next_nodes = []
                    for edge in self.edges:
                        if edge.target == path[-1]:
                            next_nodes.append(edge.source)
                    if len(next_nodes) == 0:
                        path.append("root")
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
            if all_paths[i][-1] not in top_nodes:
                del(all_paths[i])
    
        return all_paths


    def slide_down(self, top_node, bottom_nodes):
        """
        Return the list of possible acyclic paths from top to bottom nodes.
        Edges are followed from source to target.
        """

        all_paths = [[top_node]]
        bottom_reached = False
        while bottom_reached == False:
            bottom_reached = True
            for i in range(len(all_paths)):
                path = all_paths[i]
                if path[-1] not in bottom_nodes and path[-1] != "end":
                    bottom_reached = False
                    next_nodes = []
                    for edge in self.edges:
                        if edge.source == path[-1]:
                            next_nodes.append(edge.target)
                    if len(next_nodes) == 0:
                        path.append("end")
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
            if all_paths[i][-1] not in bottom_nodes:
                del(all_paths[i])

        return all_paths


    def get_maxrank(self):
        """ Find the highest rank of a node in CausalGraph. """

        all_ranks = []
        for node in self.nodes:
            if node.rank != None:
                all_ranks.append(node.rank)
        if len(all_ranks) > 0:
            self.maxrank = max(all_ranks)


    def sequentialize_ids(self):
        """ Assign sequential node ids, getting rid of event numbers. """

        node_number = 1
        for current_rank in range(self.maxrank+1):
            for node in self.nodes:
                if node.rank == current_rank:
                    node.nodeid = "node{}".format(node_number)
                    node_number += 1
        # Also sort edges.
        sorted_edges = sorted(self.edges, key=lambda x: x.source.rank)
        self.edges = sorted_edges


    def cleanup(self):
        """
        Remove nodes that do not have a path to an intro node and an eoi node
        and any edge that points to or from these nodes
        """

        self.intro_nodes = []
        self.eoi_nodes = []
        for node in self.nodes:
            if node.intro == True:
                self.intro_nodes.append(node)
            if node.label == self.eoi:
                self.eoi_nodes.append(node)
        nodes_to_clean = []
        for i in range(len(self.nodes)):
            node = self.nodes[i]
            if node.intro == False and node.label != self.eoi:
                paths_up = self.climb_up(node, self.intro_nodes)
                paths_down = self.slide_down(node, self.eoi_nodes)
                if len(paths_up) == 0 or len(paths_down) == 0:
                    nodes_to_clean.insert(0,i)
        edges_to_clean = []
        for j in range(len(self.edges)):
            source = self.edges[j].source
            target = self.edges[j].target
            for i in nodes_to_clean:
                node = self.nodes[i]
                if source == node or target == node:
                    edges_to_clean.insert(0,j)
        for j in edges_to_clean:
            del(self.edges[j])
        for i in nodes_to_clean:
            del(self.nodes[i])
                

    def build_dot_file(self, edgelabels=False, hideintro=False):
        """ build a dot file of the CausalGraph. """

        dot_str = "digraph G{\n"
        dot_str += '  nodestype="{}" ;\n'.format(self.nodestype)
        if self.eoi != None:
            dot_str += '  eoi="{}" ;\n'.format(self.eoi)
        if self.occurrence != None:
            dot_str += '  label="Occurrence = {}" '.format(self.occurrence)
            dot_str += 'fontsize=28 labelloc="t" ;\n'
        if self.maxrank != None:
            dot_str += '  maxrank="{}" ;\n'.format(self.maxrank)
        if self.prevcores != None:
            dot_str += '  prevcores="{}" ;\n'.format(self.prevcores)
        if "core" in self.filename:
            dot_str += "  ranksep=0.5 ;\n"
        else:
            dot_str += "  ranksep=1.25 ;\n"
        # Draw nodes.
        for current_rank in range(self.maxrank+1):
            rank_str = "{}".format(current_rank)
            if hideintro == True and current_rank == 0:
                dot_str += "//"
            dot_str += ('{{ rank = same ; "{}" [shape=plaintext] ;\n'
                        .format(rank_str))
            for node in self.nodes:
                if node.rank == current_rank:
                    node_shape = "invhouse"
                    node_color = "lightblue"
                    if node.intro == True:
                        node_shape = "rectangle"
                        node_color = "white"
                    if node.label == self.eoi:
                        node_shape = "ellipse"
                        node_color = "indianred2"
                    if self.nodestype == "species":
                        node_shape = "ellipse" # "circle"
                    if hideintro == True and node.intro == True:
                        dot_str += "//"
                    dot_str += ('"{}" [label="{}", '
                                .format(node.nodeid, node.label))
                    dot_str += "shape={}, style=filled, ".format(node_shape)
                    dot_str += "fillcolor={}, ".format(node_color)
                    dot_str += "intro={}] ;\n".format(node.intro)
            if hideintro == True and current_rank == 0:
                dot_str += "//"
            dot_str += "}\n"
        for rank in range(self.maxrank):
            if hideintro == True and rank == 0:
                dot_str += "//"
            rank_str = "{}".format(rank)
            next_rank = "{}".format(rank+1)
            dot_str += ('"{}" -> "{}" [style="invis"] ;\n'
                        .format(rank_str, next_rank))
        # Draw edges.
        minpenwidth = 0.5
        medpenwidth = 3
        maxpenwidth = 6.5
        all_weights = []
        for edge in self.edges:
            all_weights.append(edge.weight)
        average_weight = statistics.mean(all_weights)
        for edge in self.edges:
            if hideintro == True and edge.source.intro == True:
                dot_str += "//"
            dot_str += ('"{}" -> "{}" '
                        .format(edge.source.nodeid, edge.target.nodeid))
            edge_color = "black"
            ratio = edge.weight/average_weight
            pensize = math.log(ratio,2) + medpenwidth
            if pensize < minpenwidth:
                pensize = minpenwidth
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

def mergecores(eoi, causalgraphs=None, edgelabels=False, hideintro=False,
               writedots=True, rmprev=False, printmsg=True):
    """ Merge equivalent causal cores and count occurrence. """

    # Reading section.
    if causalgraphs == None:
        causal_core_files = get_dot_files(eoi, "causalcore")
        causal_cores = []
        for core_file in causal_core_files:
            core_path = "{}/{}".format(eoi, core_file)
            causal_cores.append(CausalGraph(core_path, eoi))
    else:
       causal_cores = causalgraphs
       causal_core_files = None
    # Doing the work.
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
        graph.build_dot_file(edgelabels, hideintro)
    # Writing section.
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
    if graph1.maxrank == graph2.maxrank:
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
              hideintro=True, writedots=True, rmprev=False,
              writepremerge=False):
    """ Build looped event paths by merging identical nodes within cores. """

    # Reading section.
    if causalgraphs == None:
        core_files = get_dot_files(eoi, "core")
        cores = []
        for core_file in core_files:
            core_path = "{}/{}".format(eoi, core_file)
            cores.append(CausalGraph(core_path, eoi))
    else:
       cores = causalgraphs
       core_files = None
    # Doing the work. 
    for core in cores:
        remove_ignored(core, ignorelist)
        merge_same_labels(core)
        fuse_edges(core)
        core.rank_nodes()
        core.build_dot_file(edgelabels, hideintro)
    # Writing section.
    if writepremerge == True:
        for i in range(len(cores)):
            cores[i].filename = "loopedcore-{}.dot".format(i+1)
        for graph in cores:
                output_path = "{}/{}".format(eoi, graph.filename)
                outfile = open(output_path, "w")
                outfile.write(graph.dot_file)
                outfile.close()
    looped_paths = mergecores(eoi, cores, writedots=False,
                              edgelabels=edgelabels, hideintro=hideintro,
                              printmsg=False)
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


def merge_same_labels(graph):
    """
    Merge every node with same label within a graph. Merging between two
    nodes is done by deleting the second node and changing any edge which was
    coming from or going to the second node to now come from or go to the
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

# ++++++++++++++++++ End of Cores Looping Section +++++++++++++++++++++++++++++

# .................. Event Paths Merging Section ..............................

def mergepaths(eoi, causalgraphs=None, threshold=0.0, edgelabels=False,
               hideintro=True, writedot=True, rmprev=False):
    """ Merge event paths into a single pathway. """

    # Reading section.
    if causalgraphs == None:
        path_files = get_dot_files(eoi, "eventpath")
        event_paths = []
        for path_file in path_files:
            path_path = "{}/{}".format(eoi, path_file)
            event_paths.append(CausalGraph(path_path, eoi))
    else:
        event_paths = causalgraphs
        path_files = None
    # Doing the work.
    pathway = CausalGraph(eoi=eoi)
    node_number = 1
    seen_labels = []
    for event_path in event_paths:
        for node in event_path.nodes:
            if node.label not in seen_labels:
                seen_labels.append(node.label)
                n_id = "node{}".format(node_number)
                pathway.nodes.append(CausalNode(n_id, node.label, node.rank,
                                                intro=node.intro))
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
    all_weights = []
    for edge in pathway.edges:
        all_weights.append(edge.weight)
    average_weight = statistics.mean(all_weights)
    for i in range(len(pathway.edges)-1, -1, -1):
        if pathway.edges[i].weight < average_weight*threshold:
            del(pathway.edges[i])
    pathway.cleanup()
    pathway.rank_nodes()
    pathway.occurrence = None
    pathway.filename = "eventpathway.dot"
    pathway.build_dot_file(edgelabels, hideintro)
    # Writing section.
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

def speciespathway3(eoi, kappamodel, causalgraph=None, edgelabels=False,
                   hideintro=False):
    """
    Convert a CausalGraph where node are events to a pathway where nodes are
    species.
    """

    # Reading section.
    if causalgraph == None:
        graph_path = "{}/eventpathway.dot".format(eoi)
        event_pathway = CausalGraph(graph_path, eoi)
    elif isinstance(causalgraph, CausalGraph):
        event_pathway = causalgraph
    else:
        event_pathway = CausalGraph(causalgraph, eoi)
    # Doing the work
    get_rules(eoi, kappamodel, event_pathway)
    get_req_res(eoi, event_pathway)
    mod_nodes = get_mod_nodes(eoi, event_pathway)
    complete_req(event_pathway, mod_nodes)
    rebranch(event_pathway, mod_nodes)
    fuse_edges(event_pathway)
    simplify_req_res(event_pathway)
    new_edges = build_species_and_edges(event_pathway)
    species_pathway = build_species_graph(eoi, event_pathway, new_edges)
    species_pathway.rank_nodes()
    species_pathway.build_dot_file(edgelabels, hideintro)
    # Writing section.
    output_path1 = "{}/{}".format(eoi, species_pathway.filename)
    outfile1 = open(output_path1, "w")
    outfile1.write(species_pathway.dot_file)
    outfile1.close()
    species_pathway.filename = ("{}-{}.dot"
                                .format(species_pathway.filename[:-4], eoi))
    output_path2 = "{}".format(species_pathway.filename)
    outfile2 = open(output_path2, "w")
    outfile2.write(species_pathway.dot_file)
    outfile2.close()
    print("Converting event pathway into species pathway.")
    print("File {} created.".format(species_pathway.filename))


def build_species_graph(eoi, graph, edge_list):
    """
    Create the new causal graph with species instead of rules as nodes.
    """

    species_pathway = CausalGraph()
    species_pathway.filename = "pathway.dot"
    species_pathway.eoi = eoi
    species_pathway.nodestype = "species"
    species_pathway.occurrence = None
    for event_node in graph.nodes:
       for species_node in event_node.species_nodes:
           found_as_source = False
           found_as_target = False
           for edge in edge_list:
               if edge.source == species_node:
                   found_as_source = True
               if edge.target == species_node:
                   found_as_target = True
               if species_node.intro == True and found_as_source == True:
                   include_node = True
               elif species_node.label == eoi and found_as_target == True:
                   include_node = True
               elif found_as_source == True and found_as_target == True:
                   include_node = True
               else:
                   include_node = False
               if include_node == True:
                   species_pathway.nodes.append(species_node)
    species_pathway.edges = edge_list

    return species_pathway


def build_species_and_edges(graph):
    """ Build nodes from res species and the edges to connecte them. """

    # Create nodes.
    species_id = 1
    for mod_node in graph.nodes:
        species_nodes = []
        for res in mod_node.res_species:
            node_id = "species{}".format(species_id)
            if res["state"] != None:
                label = "{}({}{{{}}})".format(res["agent"], res["site"],
                                              res["state"])
            else:
                label = "{}({}[{}])".format(res["agent"], res["site"],
                                            res["bound_agent"])
            new_node = CausalNode(node_id, label, rank=mod_node.rank,
                                  intro=mod_node.intro)
            new_node.species = res
            species_nodes.append(new_node)
            species_id += 1
        mod_node.species_nodes = species_nodes
    # Create edges.
    new_edges = []
    for mod_node in graph.nodes:
        for edge in graph.edges:
            if edge.source == mod_node:
                target_mod_node = edge.target
                for src in mod_node.species_nodes:
                    if species_in(src.species, target_mod_node.full_req):
                        for trgt in target_mod_node.species_nodes:
                            new_edges.append(CausalEdge(src, trgt,
                                                        weight=edge.weight))

    return new_edges


def simplify_req_res(graph):
    """
    Remove bindings and unbindings from res. Also remove in req all
    agents found in res, except if there is only one agent type in req.
    """

    for event_node in graph.nodes:
        if event_node.intro == False:
            res_to_remove = []
            for i in range(len(event_node.res_species)):
                if event_node.res_species[i]["state"] == None:
                    res_to_remove.insert(0, i)
            for i in res_to_remove:
                del(event_node.res_species[i])
            req_agents = []
            for req in event_node.full_req:
                if req["agent"] not in req_agents:
                    req_agents.append(req["agent"])
            if len(req_agents) > 1:
                res_agents = []
                for res in event_node.res_species:
                    if res["agent"] not in res_agents:
                        res_agents.append(res["agent"])
                req_to_remove = []
                for i in range(len(event_node.full_req)):
                    if event_node.full_req[i]["agent"] in res_agents:
                        req_to_remove.insert(0, i)
                for i in req_to_remove:
                    del(event_node.full_req[i])
            else:
                req_to_remove = []
                for i in range(len(event_node.full_req)):
                    req_ag = event_node.full_req[i]["agent"]
                    req_site = event_node.full_req[i]["site"]
                    for res in event_node.res_species:
                        if req_ag == res["agent"] and req_site == res["site"]:
                            req_to_remove.insert(0, i)
                            break
                for i in req_to_remove:
                    del(event_node.full_req[i])


def complete_req(graph, mod_nodes):
    """
    Extend req_species of each mod node with the req of upstream nodes.
    """

    mod_no_intro = []
    for mod_node in mod_nodes:
        if mod_node.intro == False:
            mod_no_intro.append(mod_node)
    for mod_node in mod_no_intro:
        top_nodes = []
        for top_node in mod_nodes:
            if top_node != mod_node:
                top_nodes.append(top_node)
        paths = graph.climb_up(mod_node, top_nodes)
        all_reqs = []
        for path in paths:
            path_reqs = []
            for i in range(len(path)-1):
                current_node = path[i]
                for current_req in current_node.req_species:
                    if not species_in(current_req, path_reqs):
                        path_reqs.append(current_req.copy())
                up_node = path[i+1]
                for up_res in up_node.res_species:
                    res_ag = up_res["agent"]
                    res_site = up_res["site"]
                    for path_req in path_reqs:
                        req_ag = path_req["agent"]
                        req_site = path_req["site"]
                        if req_ag == res_ag and req_site == res_site:
                            if up_res["bound_agent"] != None:
                                if path_req["bound_agent"] == "_":
                                    bnd_ag = up_res["bound_agent"]
                                    bnd_site = up_res["bound_site"]
                                    path_req["bound_agent"] = bnd_ag
                                    path_req["bound_site"] = bnd_site
                                    partner = {"agent": bnd_ag,
                                               "site": bnd_site,
                                               "bound_agent": req_ag,
                                               "bound_site": req_site,
                                               "state": None}
                                    path_reqs.append(partner)
            all_reqs.append(path_reqs)
        # Assemble all the species collected from each path into one set.
        comp_req_set = []
        for path_reqs in all_reqs:
            for path_req in path_reqs:
                if not species_in(path_req, comp_req_set):
                    comp_req_set.append(path_req)
        mod_node.full_req = comp_req_set


def get_mod_nodes(eoi, graph):
    """ Get modification nodes based on resulting species with state. """

    mod_nodes = []
    for event_node in graph.nodes:
        if event_node.intro == True or event_node.label == eoi:
            mod_nodes.append(event_node)
        else:
            is_mod = False
            for res in event_node.res_species:
                if res["state"] != None:
                    is_mod = True
                    break
            if is_mod == True:
                mod_nodes.append(event_node)
    ## Remove any species that does not denote a state change
    ## from the res of mod nodes.
    #for mod_node in mod_nodes:
    #    if mod_node.intro == False:
    #        species_to_remove = []
    #        for i in range(len(mod_node.res_species)):
    #            if mod_node.res_species[i]["state"] == None:
    #                species_to_remove.insert(0, i)
    #        for i in species_to_remove:
    #            del(mod_node.res_species[i])

    return mod_nodes


def get_req_res(eoi, graph):
    """
    For each event node, get the requirements and results of
    the associated rule.
    """

    for event_node in graph.nodes:
        req_species = []
        res_species = []
        species_list = build_species(event_node.rule)
        for species in species_list:
            for char in ["binding", "state"]:
                if species[char] != None:
                    if "/" in species[char]:
                        slash = species[char].index("/")
                        before = species[char][:slash]
                        after = species[char][slash+1:]
                        before_species = species.copy()
                        after_species = species.copy()
                        before_species[char] = before
                        after_species[char] = after
                        req_species.append(before_species)
                        res_species.append(after_species)
                    else:
                        req_species.append(species)
        req_type, res_type = type_bonds(req_species, res_species)
        req_set = []
        res_set = []
        for req_entry in req_type:
            is_in = False
            for req_set_entry in req_set:
                are_equal = compare_species(req_entry, req_set_entry)
                if are_equal == True:
                    is_in = True
                    break
            if is_in == False:
                req_set.append(req_entry)
        for res_entry in res_type:
            is_in = False
            for res_set_entry in res_set:
                are_equal = compare_species(res_entry, res_set_entry)
                if are_equal == True:
                    is_in = True
                    break
            if is_in == False:
                res_set.append(res_entry)
        if event_node.intro == False:
            event_node.req_species = req_set
            event_node.res_species = res_set
        elif event_node.intro == True:
            event_node.req_species = []
            event_node.res_species = req_set


def build_species(rule_str):
    """
    Build a list of dictionaries that each represent an individual
    species from a rule.
    """

    species_list = []
    agents_list = []
    ag_tmp = rule_str.split(",")
    for ag in ag_tmp:
        agents_list.append(ag.strip())
    for agent in agents_list:
        par = agent.index("(")
        agent_name = agent[:par]
        sites_list = agent[par+1:-1].split()
        for site in sites_list:
            open_bracket = 10000
            open_curl = 10000
            if "[" in site:
                open_bracket = site.index("[")
                close_bracket = site.index("]")
                binding = site[open_bracket+1:close_bracket]
            else:
                binding = None
            if "{" in site:
                open_curl = site.index("{")
                close_curl = site.index("}")
                state = site[open_curl+1:close_curl]
            else:
                state = None
            if open_bracket < 10000 or open_curl < 10000:
                name_end = min([open_bracket, open_curl])
                site_name = site[:name_end]
            else:
                site_name = site
            if binding != None:
                species = {"agent": agent_name, "site": site_name,
                           "binding": binding, "state": None}
                species_list.append(species)
            if state != None:
                species = {"agent": agent_name, "site": site_name,
                           "binding": None, "state": state}
                species_list.append(species)
    return species_list


def type_bonds(req_species, res_species):
    """ Change link numbers to semi-link with type. """

    bond_numbers = get_bond_numbers(req_species, res_species)
    list_index = 0
    for species_list in [req_species, res_species]:
        new_species = []
        for species in species_list:
            if species["binding"] != None:
                number = species["binding"]
                if number != "." and number != "_":
                    current_site = "{}.{}".format(species["site"],
                                                  species["agent"])
                    new_bond = bond_numbers[number][current_site]
                    period = new_bond.index(".")
                    bnd_agent = new_bond[period+1:]
                    bnd_site = new_bond[:period]
                else:
                    bnd_agent = number
                    bnd_site = number
            else:
                bnd_agent = None
                bnd_site = None
            new_sp = {"agent": species["agent"],
                      "site": species["site"],
                      "bound_agent": bnd_agent,
                      "bound_site": bnd_site,
                      "state": species["state"]}
            new_species.append(new_sp)
        if list_index == 0:
            new_req = new_species
            list_index += 1
        else:
            new_res = new_species

    return new_req, new_res


def get_bond_numbers(req_species, res_species):
    """ Find which agent binds to which agent. """

    bond_numbers_tmp = {}
    for species_list in [req_species, res_species]:
        for species in species_list:
            if species["binding"] != None:
                number = species["binding"]
                if number != "." and number != "_":
                    bond_type = "{}.{}".format(species["site"],
                                               species["agent"])
                    if number not in bond_numbers_tmp.keys():
                        bond_numbers_tmp[number] = [bond_type]
                    else:
                        bond_numbers_tmp[number].append(bond_type)
    bond_numbers = {}
    for number in bond_numbers_tmp.keys():
        partners = bond_numbers_tmp[number]
        bond_numbers[number] = {partners[0]: partners[1],
                                partners[1]: partners[0]}

    return bond_numbers


def species_in(species, species_list):
    """ Check if a species in found inside a species list. """

    is_in = False
    for sp in species_list:
        are_equal = compare_species(species, sp)
        if are_equal == True:
            is_in = True
            break

    return is_in


def compare_species(species1, species2):
    """ Check if two species with type_bonds are identical. """

    are_identical = True
    for char in ["agent", "site", "bound_agent", "bound_site", "state"]:
        if species1[char] != species2[char]:
            are_identical = False
            break

    return are_identical


def speciespathway(eoi, kappamodel, causalgraph=None, edgelabels=False,
                   hideintro=False):
    """
    Convert a CausalGraph where node are events to a pathway where nodes are
    species.
    """

    # Reading section.
    if causalgraph == None:
        graph_path = "{}/eventpathway.dot".format(eoi)
        pathway = CausalGraph(graph_path, eoi)
    elif isinstance(causalgraph, CausalGraph):
        pathway = causalgraph
    else:
        pathway = CausalGraph(causalgraph, eoi)
    # Doing the work.
    gather_rules(eoi, kappamodel, pathway)
    #mod_nodes = get_mod_nodes(eoi, pathway)
    #rebranch(pathway, mod_nodes)
    #fuse_edges(pathway)
    arrow_relationships(eoi, pathway)
    new_links = linkresnodes(pathway)
    # Build new causal graph with the sites as nodes.
    species_pathway = CausalGraph()
    species_pathway.filename = "pathway.dot"
    species_pathway.eoi = eoi
    species_pathway.nodestype = "species"
    species_pathway.occurrence = None
    mod_sites = []
    for rule_node in pathway.nodes:
        for rule_res in rule_node.res:
            found_as_source = False
            found_as_target = False
            for link in new_links:
                if link.source == rule_res:
                    found_as_source = True
                if link.target == rule_res:
                    found_as_target = True
            if rule_res.intro == True and found_as_source == True:
                include_node = True
            elif rule_res.label == eoi and found_as_target == True:
                include_node = True
            elif found_as_source == True and found_as_target == True:
                include_node = True
            else:
                include_node = False
            if include_node == True:
                species_pathway.nodes.append(rule_res)
                if "{" in rule_res.label or rule_node.intro == True:
                    mod_sites.append(rule_res)
    for link in new_links:
        if link.source in species_pathway.nodes:
            if link.target in species_pathway.nodes:
                species_pathway.edges.append(link)
    print(">>>>", mod_sites)
    rebranch(species_pathway, mod_sites)
    merge_same_labels(species_pathway)
    fuse_edges(species_pathway)
    species_pathway.rank_nodes()
    ## Remove intro nodes if the agent in their label is the same as their target.
    #nodes_to_remove = []
    #for i in range(len(species_pathway.nodes)):
    #    node = species_pathway.nodes[i]
    #    if node.intro == True:
    #        targets = []
    #        for edge in species_pathway.edges:
    #            if edge.source == node:
    #                targets.append(edge.target)
    #        target_agents = []
    #        for target in targets:
    #            par = target.label.index("(")
    #            agent = target.label[:par]
    #            target_agents.append(agent)
    #        par = node.label.index("(")
    #        source_agent = node.label[:par]
    #        if source_agent in target_agents:
    #            nodes_to_remove.insert(0, i)
    #edges_to_remove = []
    #for j in range(len(species_pathway.edges)):
    #    edge = species_pathway.edges[j]
    #    for i in nodes_to_remove:
    #        node = species_pathway.nodes[i]
    #        if edge.source == node:
    #            edges_to_remove.insert(0, j)
    #            break
    #for i in edges_to_remove:
    #    del(species_pathway.edges[i])
    #for i in nodes_to_remove:
    #    del(species_pathway.nodes[i])
    # Simplify labels.
    #for node in species_pathway.nodes:
    #    par = node.label.index("(")
    #    agent = node.label[:par]
    #    if "[" in node.label:
    #        bracket = node.label.index("[")
    #        site = node.label[par+1:bracket]
    #        new_label = "{}-{}".format(agent, site)
    #    elif "{" in node.label:
    #        open_brace = node.label.index("{")
    #        close_brace = node.label.index("}")
    #        site = node.label[par+1:open_brace]
    #        state = node.label[open_brace+1:close_brace]
    #        new_label = "{}-{} {}".format(agent, site, state)
    #    if node.label == species_pathway.eoi:
    #        species_pathway.eoi = new_label
    #    node.label = new_label
    species_pathway.build_dot_file(edgelabels)
    # Writing section.
    output_path1 = "{}/{}".format(eoi, species_pathway.filename)
    outfile1 = open(output_path1, "w")
    outfile1.write(species_pathway.dot_file)
    outfile1.close()
    species_pathway.filename = ("{}-{}.dot"
                                .format(species_pathway.filename[:-4], eoi))
    output_path2 = "{}".format(species_pathway.filename)
    outfile2 = open(output_path2, "w")
    outfile2.write(species_pathway.dot_file)
    outfile2.close()
    print("Converting event pathway into species pathway.")
    print("File {} created.".format(pathway.filename))


def get_rules(eoi, kappamodel, graph):
    """ Assign rule to every node. """

    kappa = read_kappa_file(kappamodel)
    kappa["eoi"] = read_eoi(eoi, kappamodel)
    kappa["intros"] = build_creation_rules(kappa)
    for node in graph.nodes:
        if node.intro == True:
            node.rule = kappa["intros"][node.label]
        elif node.label == eoi:
            node.rule = kappa["eoi"][node.label]
        else:
            node.rule = kappa["rules"][node.label]


def read_kappa_file(kappamodel):
    """ Build a dictionary of the rules from the original kappa model. """

    kappa_file = open(kappamodel, "r").readlines()
    kappa = {}
    kappa["agents"] = {}
    kappa["inits"] = {}
    kappa["rules"] = {}
    for line in kappa_file:
        if line[:7] == "%agent:":
            agent_def = line[7:-1].strip()
            par = agent_def.index("(")
            agent_type = agent_def[:par]
            kappa["agents"][agent_type] = agent_def
        if line[:6] == "%init:":
            amount = line[6:-1].strip()
            space = amount.index(" ")
            init_def = amount[space:].strip()
            init_agents = init_def.split(",")
            agent_type = ""
            first_agent = True
            for init_agent in init_agents:
                if first_agent == False:
                    agent_type += ", "
                else:
                    first_agent = False
                agent_str = init_agent.strip()
                par = agent_str.index("(")
                agent_type += init_def[:par]
            kappa["inits"][agent_type] = init_def
        if line[0] == "'":
            quote = line.rfind("'")
            rule_name = line[1:quote]
            a = line.index("@")
            rule = line[quote+1:a].strip()
            kappa["rules"][rule_name] = rule

    return kappa


def read_eoi(eoi, kappamodel):
    """ Read EOI from the kappa file with added event of interest. """

    period = kappamodel.rfind(".")
    prefix = kappamodel[:period]
    kappa_path = "{}/{}-eoi.ka".format(eoi, prefix)
    kappa_file = open(kappa_path, "r").readlines()
    eoi_dict = {}
    for line in kappa_file:
        if line[:5] == "%obs:":
            open_quote = line.index("'")
            close_quote = line.rfind("'")
            obs = line[open_quote+1:close_quote]
            if obs == eoi:
                obs_def = line[close_quote+1:-1].strip()
                if "|" in obs_def:
                    obs_def = obs_def[1:-1]
                eoi_dict[obs] = obs_def

    return eoi_dict


def build_creation_rules(kappa_dict):
    """
    Build a creation rule for each Intro node from the init and agent
    definitions.
    """

    creations = {}
    for intro in kappa_dict["inits"].keys():
        init = kappa_dict["inits"][intro]
        init_agents = init.split(",")
        intro_rule = ""
        first_agent = True
        for init_agent_tmp in init_agents:
            init_agent = init_agent_tmp.strip()
            if first_agent == False:
                intro_rule += ", "
            else:
                first_agent = False
            init_dict = build_site_dict(init_agent)
            agent_name = init_dict["name"]
            def_agent = kappa_dict["agents"][agent_name]
            def_dict = build_site_dict(def_agent)
            intro_rule += "{}(".format(agent_name)
            init_names = init_dict.keys()
            default_bind = "."
            first_site = True
            for init_name in init_names:
                if init_name != "name":
                    if first_site == False:
                        intro_rule += " "
                    else:
                       first_site = False
                    intro_rule += "{}".format(init_name)
                    init_site = init_dict[init_name]
                    def_site = def_dict[init_name]
                    init_bind = init_site["binding"]
                    init_state = init_site["state"]
                    if def_site["state"] == None:
                        if init_bind == None:
                            intro_rule += "[.]"
                        else:
                            intro_rule += "[{}]".format(init_bind)
                    elif def_site["state"] != None:
                        comma = def_site["state"].index(",")
                        default_state = def_site["state"][:comma]
                        if init_bind == None and init_state == None:
                            intro_rule += "[{}]{{{}}}".format(default_bind,
                                                              default_state)
                        elif init_bind != None and init_state == None:
                            intro_rule += "[{}]{{{}}}".format(init_bind,
                                                              default_state)
                        elif init_bind == None and init_state != None:
                            intro_rule += "[{}]{{{}}}".format(default_bind,
                                                              init_state)
                        elif init_bind != None and init_state != None:
                            intro_rule += "[{}]{{{}}}".format(init_bind,
                                                              init_state)
            for def_name in def_dict.keys():
                if def_name not in init_names:
                    if first_site == False:
                        intro_rule += " "
                    else:
                       first_site = False
                    def_site = def_dict[def_name]
                    intro_rule += "{}[{}]".format(def_name, default_bind)
                    if def_site["state"] != None:
                        comma = def_site["state"].index(",")
                        default_state = def_site["state"][:comma]
                        intro_rule += "{{{}}}".format(default_state)
            intro_rule += ")"
        intro_label = "Intro {}" .format(intro)
        creations[intro_label] = intro_rule

    return creations


def build_site_dict(agent_str):
    """
    Build a dictionary of the sites of an agent given a string of that agent.
    """

    site_dict = {}
    par = agent_str.index("(")
    site_dict["name"] = agent_str[:par]
    site_list = agent_str[par+1:-1].split()
    for site in site_list:
        open_bracket = 10000
        open_curl = 10000
        if "[" in site:
            open_bracket = site.index("[")
            close_bracket = site.index("]")
            binding = site[open_bracket+1:close_bracket]
        else:
            binding = None
        if "{" in site:
            open_curl = site.index("{")
            close_curl = site.index("}")
            state = site[open_curl+1:close_curl]
        else:
            state = None
        if open_bracket < 10000 or open_curl < 10000:
            name_end = min([open_bracket, open_curl])
            site_name = site[:name_end]
        else:
            site_name = site
        site_dict[site_name] = {"binding": binding, "state": state}

    return site_dict


def build_site_str(site_dict):
    """
    Build a string of and agent given a dict of its sites.
    This is the reverse of what build_site_dict does.
    """

    agent_str = "{}(".format(site_dict["name"])
    first_site = True
    for site in site_dict.keys():
        if site != "name":
            if first_site == False:
                agent_str += " "
            else:
                first_site = False
            agent_str += "{}".format(site)
            if site_dict[site]["binding"] != None:
                agent_str += "[{}]".format(site_dict[site]["binding"])
            if site_dict[site]["state"] != None:
                agent_str += "{{{}}}".format(site_dict[site]["state"])
    agent_str += ")"

    return agent_str


#def get_mod_nodes(eoi, graph):
#    """ Get modification nodes based on corresponding rule. """
#
#    mod_nodes = []
#    for node in graph.nodes:
#        is_mod = False
#        parts = node.rule.split()
#        for part in parts:
#            if "{" in part:
#                open_brace = part.index("{")
#                close_brace = part.index("}")
#                if "/" in part[open_brace+1:close_brace]:
#                    is_mod = True
#                    break
#        if node.intro == True:
#            is_mod = True
#        if node.label == eoi:
#            is_mod = True
#        if is_mod == True:
#            mod_nodes.append(node)
#
#    return mod_nodes


def arrow_relationships(eoi, graph):
    """
    Add required site nodes and resulting site nodes to each modification node.
    """

    node_index = 1
    for node in graph.nodes:
        req_sites, res_sites = individual_sites(node.rule)
        bond_numbers = get_bond_numbers(req_sites, res_sites)
        req_list = relative_bonds(req_sites, bond_numbers)
        res_list = relative_bonds(res_sites, bond_numbers)
        req = []
        for site in req_list:
            node_id = "site{}".format(node_index)
            req.append(CausalNode(node_id, site, rank=node.rank,
                                  intro=node.intro))
            node_index += 1
        res = []
        for site in res_list:
            node_id = "site{}".format(node_index)
            res.append(CausalNode(node_id, site, rank=node.rank))
            node_index += 1
        if node.intro == True:
            node.res = req
            node.req = []
        else:
            node.req = req
            node.res = res
        if node.label == eoi:
            for req_node in node.req:
                if "[" in req_node.label:
                    node.res.append(req_node)
        print(node.label,"|", node.rule)
        print(node.req)
        print(node.res)
        print("----")


def individual_sites(rule):
    """
    Return lists of individual required and resulting sites (species) from a kappa rule.
    """

    req_sites = []
    res_sites = []
    agents_list = rule.split(",")
    for agent_tmp in agents_list:
        agent = agent_tmp.strip()
        site_dict = build_site_dict(agent)
        for site_key in site_dict.keys():
            site = site_dict[site_key]
            if site_key != "name":
                if site["binding"] != None:
                    if "/" in site["binding"]:
                        slash = site["binding"].index("/")
                        preslash = site["binding"][:slash]
                        #indiv_site = {"name": site_dict["name"],
                        #              site_key: {"binding": preslash}}
                        site_str = "{}({}".format(site_dict["name"], site_key)
                        site_str += "[{}])".format(preslash)
                        req_sites.append(site_str)
                        postslash = site["binding"][slash+1:]
                        #indiv_site = {"name": site_dict["name"],
                        #              site_key: {"binding": postslash}}
                        site_str = "{}({}".format(site_dict["name"], site_key)
                        site_str += "[{}])".format(postslash)
                        res_sites.append(site_str)
                    else:
                        #indiv_site = {"name": site_dict["name"],
                        #              site_key: {"binding": site["binding"]}}
                        site_str = "{}({}".format(site_dict["name"], site_key)
                        site_str += "[{}])".format(site["binding"])
                        req_sites.append(site_str)
                if site["state"] != None:
                    if "/" in site["state"]:
                        slash = site["state"].index("/")
                        preslash = site["state"][:slash]
                        #indiv_site = {"name": site_dict["name"],
                        #              site_key: {"state": preslash}}
                        site_str = "{}({}".format(site_dict["name"], site_key)
                        site_str += "{{{}}})".format(preslash)
                        req_sites.append(site_str)
                        postslash = site["state"][slash+1:]
                        #indiv_site = {"name": site_dict["name"],
                        #              site_key: {"state": postslash}}
                        site_str = "{}({}".format(site_dict["name"], site_key)
                        site_str += "{{{}}})".format(postslash)
                        res_sites.append(site_str)
                    else:
                        #indiv_site = {"name": site_dict["name"],
                        #              site_key: {"state": site["state"]}}
                        site_str = "{}({}".format(site_dict["name"], site_key)
                        site_str += "{{{}}})".format(site["state"])
                        req_sites.append(site_str)

    return req_sites, res_sites


#def get_bond_numbers(req_sites, res_sites):
#    """ Find which agent binds to which agent. """
#
#    bond_numbers_tmp = {}
#    for site_list in [req_sites, res_sites]:
#        for site in site_list:
#            if "[" in site:
#                bracket = site.index("[")
#                number = site[bracket+1:-2]
#                if number != "." and number != "_":
#                    par = site.index("(")
#                    ag = site[:par]
#                    s = site[par+1:bracket]
#                    if number not in bond_numbers_tmp.keys():
#                        bond_numbers_tmp[number] = ["{}.{}".format(s, ag)]
#                    else:
#                        bond_numbers_tmp[number].append("{}.{}".format(s, ag))
#    bond_numbers = {}
#    for number in bond_numbers_tmp.keys():
#        partners = bond_numbers_tmp[number]
#        bond_numbers[number] = {partners[0]: partners[1],
#                                partners[1]: partners[0]}
#
#    return bond_numbers


def type_bonds2(sites, bond_numbers):
    """ Change link numbers to semi-link with type and remove duplicates. """

    new_sites = []
    for site in sites:
        if "[" in site:
            bracket = site.index("[")
            number = site[bracket+1:-2]
            if number != "." and number != "_":
                par = site.index("(")
                ag = site[:par]
                s = site[par+1:bracket]
                current_site = "{}.{}".format(s, ag)
                new_site = "{}({}".format(ag, s)
                new_bond = bond_numbers[number][current_site]
                new_site += "[{}])".format(new_bond)
                new_sites.append(new_site)
            else:
                new_sites.append(site)
        else:
            new_sites.append(site)
    sites_set = list(set(new_sites))
    #site_dicts = []
    #for site in sites_set:
    #    site_dict = build_site_dict(site)
    #    site_dicts.append(site_dict)

    return sites_set
    #return site_dicts


def linkresnodes(graph):
    """
    For each rule, link any of its res nodes with all the res nodes of a
    subsequent rule if the given res node from the first rule has the same
    label as any req node from the subsequent rule.
    """

    links = []
    for node in graph.nodes:
        target_rules = []
        for edge in graph.edges:
            if edge.source == node:
                target_rule = edge.target
                w = edge.weight
                for node_res in node.res:
                    link_res_nodes = False
                    for target_req in target_rule.req:
                        if "[_]" in target_req.label:
                            req_par = target_req.label.index("(")
                            req_bracket = target_req.label.index("[")
                            req_agent = target_req.label[:req_par]
                            req_site = target_req.label[req_par+1:req_bracket]
                            res_par = node_res.label.index("(")
                            res_bracket = node_res.label.index("[")
                            res_agent = node_res.label[:res_par]
                            res_site = node_res.label[res_par+1:res_bracket]
                            if req_agent == res_agent and req_site == res_site:
                                link_res_nodes = True
                                break
                        elif node_res.label == target_req.label:
                            link_res_nodes = True
                            break
                    if link_res_nodes == True:
                        for target_res in target_rule.res:
                            links.append(CausalEdge(node_res, target_res,
                                                    weight=w))

    return links


def oldspeciespathway(eoi, kappamodel, causalgraph=None, edgelabels=False,
                   hideintro=False):
    """
    Convert a CausalGraph where node are events to a pathway where nodes are
    species.
    """

    # Reading section.
    if causalgraph == None:
        graph_path = "{}/eventpathway.dot".format(eoi)
        pathway = CausalGraph(graph_path, eoi)
    elif isinstance(causalgraph, CausalGraph):
        pathway = causalgraph
    else:
        pathway = CausalGraph(causalgraph, eoi)
    # Doing the work.
    kappa_rules = get_kappa_rules(kappamodel)
    kappa_rules["{}".format(eoi)] = "|{}|".format(eoi)
    mod_nodes = get_mod_nodes(eoi, pathway, kappa_rules)
    #added_nodes = add_first_nodes(pathway, kappa_rules, mod_nodes)
    added_nodes = []
    for node in added_nodes:
        mod_nodes.append(node)
    for node in pathway.nodes:
        print(node)
        if node in mod_nodes:
            node.label = node.species
    rebranch(pathway, mod_nodes)
    merge_same_labels(pathway)
    fuse_edges(pathway)
    #pathway.rank_nodes()
    pathway.occurrence = None
    pathway.nodestype = "species"
    pathway.build_dot_file(edgelabels, hideintro)
    # Writing section.
    pathway.filename = "pathway.dot"
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


#def get_mod_nodes(eoi, graph, kappa_rules):
#    """
#    Create a list of the nodes that correspond to a rule where a state is
#    modified. Also change the label of those nodes to the species they produce.
#    """
#
#    modification_rules = []
#    for node in graph.nodes:
#        if node.intro == False:
#            rule = kappa_rules[node.label]
#            rule_agents = parse_rule(rule)
#            modified_agents = []
#            for agent in rule_agents:
#                modified_agent = {}
#                modified_sites = {}
#                sites = agent["sites"]
#                for site in sites.keys():
#                    state = sites[site]["state"]
#                    if state != None:
#                        if "/" in state:
#                            slash = state.index("/")
#                            final_state = state[slash+1:]
#                            modified_sites[site] = {"state": final_state}
#                if len(modified_sites.keys()) > 0:
#                    modified_agent["type"] = agent["type"]
#                    modified_agent["sites"] = modified_sites
#                    modified_agents.append(modified_agent)
#            # Also add the EOI if it is a binding.
#            if node.label == eoi:
#                for agent in rule_agents:
#                    modified_agent = {}
#                    modified_sites = {}
#                    sites = agent["sites"]
#                    for site in sites.keys():
#                        #state = sites[site]["state"]
#                        #if state != None:
#                        #    modified_sites[site] = {"state": state}
#                        binding = sites[site]["binding"]
#                        if binding != None:
#                            modified_sites[site] = {"binding": binding}
#                    if len(modified_sites.keys()) > 0:
#                        modified_agent["type"] = agent["type"]
#                        modified_agent["sites"] = modified_sites
#                        modified_agents.append(modified_agent)
#            if len(modified_agents) > 0:
#                species, kappa_species = label_species(modified_agents)
#                node.species = species
#                if graph.eoi == kappa_species:
#                    graph.eoi = species
#                #node.species = kappa_species
#                modification_rules.append(node)
#
#    return modification_rules


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
            if "binding" in sites[site].keys():
                kappa_species += "[{}]".format(sites[site]["binding"])
            if "state" in sites[site].keys():
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
    intro_nodes = []
    for node in graph.nodes:
        if node.intro == True:
            intro_nodes.append(node)
    for start_node in intro_nodes:
        all_paths.append([start_node])
        seen_agents.append([])
        path_weights.append(0)
    all_complete = False
    while all_complete == False:
        all_complete = True
        for i in range(len(all_paths)):
            current_node = all_paths[i][-1]
            if current_node != "mod_reached":
                if current_node not in intro_nodes:
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
    #graph.update()

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

    # 4) Merge all event paths into one event pathway.
    eventpath = mergepaths(eoi, causalgraphs=loopedcores,
                           edgelabels=edgelabels, writedot=True, rmprev=False)

    # 5) Convert an event pathway into a species pathway.
    speciespathway(eoi, kappamodel, causalgraph=eventpath,
                   edgelabels=edgelabels)

# /////////////////////// Other Useful Functions //////////////////////////////

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


def toggleintros():
    """
    Show or hide intro nodes on all path dot files (excludes cores) found in
    current directory and all other directories found in current directory.
    """

    current_files = os.listdir(".")
    for current_file in current_files:
        if "dot" in current_file and "path" in current_file:
            toggle_intro_nodes(".", current_file)
    directories = filter(os.path.isdir, os.listdir('.'))
    for directory in directories:
        dot_files = os.listdir("{}".format(directory))
        for dot_file in dot_files:
            if "dot" in dot_file and "path" in dot_file:
                toggle_intro_nodes(directory, dot_file)


def toggle_intro_nodes(dir_path, dot_file):
    """ Show or hide intro nodes in dot file. """

    input_path = "{}/{}".format(dir_path, dot_file)
    input_file = open(input_path, "r").readlines()
    new_file = ""
    node_ids = []
    rank0 = False
    for line in input_file:
        if "intro=True" in line:
            new_file += toggle_comment(line)
            tokens = line.split()
            node_ids.append(tokens[0].strip("/"))
        elif 'rank = same ; "0"' in line:
            new_file += toggle_comment(line)
            rank0 = True
        elif rank0 == True and line[-2] == "}":
            new_file += toggle_comment(line)
            rank0 = False
        elif '"0" -> "1" [style="invis"]' in line:
            new_file += toggle_comment(line)
        elif "->" in line:
            tokens = line.split()
            source = tokens[0].strip("/")
            target = tokens[2]
            if source in node_ids or target in node_ids:
                new_file += toggle_comment(line)
            else:
                new_file += line
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


def toggle_comment(line):
    """ Add or remove comment mark. """

    if line[:2] == "//":
        new_line = line[2:]
    else:
        new_line = "//{}".format(line)

    return new_line
