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


class EventNode(object):
    """
    An event node to use in causal graphs. It represents a specific event in
    causal cores and an event type (a rule) in pathways.
    """

    def __init__(self, nodeid, label, rank=None, weight=None, intro=False,
                 first=False):
        """ Initialize class EventNode. """

        self.nodeid = nodeid
        self.label = label
        self.rank = rank
        self.weight = weight
        self.intro = intro
        self.first = first
        self.check_types()


    def check_types(self):
        """ Check that EventNode attributes have proper types. """

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
        res += ",  first: {}".format(self.first)

        return res


class CausalEdge(object):
    """
    A relationship between two event nodes in causal graphs. The relationship
    can be precedence (default), causal or conflict.
    """

    def __init__(self, source, target, weight=1, relationtype="precedence",
                 color="black", underlying=False):
        """ Initialize class CausalEdge. """

        self.source = source
        self.target = target
        self.weight = weight
        self.relationtype = relationtype
        self.color = color
        self.underlying = underlying
        self.check_types()


    def check_types(self):
        """ Check that CausalEdge attributes have proper types. """

        if not isinstance(self.source, EventNode):
            raise TypeError("source should be an EventNode.")
        if not isinstance(self.target, EventNode):
            raise TypeError("target should be an EventNode.")
        if self.weight != None:
            if not isinstance(self.weight, int):
                if not isinstance(self.weight, float):
                    raise TypeError("weight should be an integer or float.")


    def __repr__(self):
        """ Representation of the CausalEdge object. """

        res =  "Edge"
        if self.weight != None:
            res += "  weight = {}".format(self.weight)
        res += "\n"
        res += "source: {}\n".format(self.source)
        res += "target: {}\n".format(self.target)

        return res


class MidNode(object):
    """
    Intermediary node to represent edge groups as CombinedEdge objects.
    MidNodes can be enablings (black by default) or requirements
    (white by default).
    """

    def __init__(self, nodeid, rank=None, midtype="enabling", logic="and"):
        """ Initialize class MidNode. """

        self.nodeid = nodeid
        self.label = ""
        self.rank = rank
        self.midtype = midtype # enabling or requirement
        self.logic = logic
        self.check_types()


    def check_types(self):
        """ Check that MidNode attributes have proper types. """

        if not isinstance(self.nodeid, str):
            raise TypeError("nodeid should be a string.")
        if self.rank != None:
            if not isinstance(self.rank, float):
                raise TypeError("rank should be an float.")


    def __repr__(self):
        """ Representation of the MidNode object. """

        res =  "Node "
        res += 'id: "{}", '.format(self.nodeid)
        if self.rank != None:
            res += ",  rank: {}".format(self.rank)

        return res


class MidEdge(object):
    """
    Intermediary edge to construct combined edges. Source and target can
    be either an EventNode or a MidNode.
    """

    def __init__(self, source, target, weight=1):
        """ Initialize class MidEdge. """

        self.source = source
        self.target = target
        self.weight = weight
        self.check_types()


    def check_types(self):
        """ Check that MidEdge attributes have proper types. """

        if not isinstance(self.source, EventNode):
            if not isinstance(self.source, MidNode):
                raise TypeError("source should be a EventNode or MidNode.")
        if not isinstance(self.target, EventNode):
            if not isinstance(self.target, MidNode):
                raise TypeError("target should be a EventNode or MidNode.")
        if self.weight != None:
            if not isinstance(self.weight, int):
                if not isinstance(self.weight, float):
                    raise TypeError("weight should be an integer or float.")


    def __repr__(self):
        """ Representation of the MidEdge object. """

        res =  "Edge"
        if self.weight != None:
            res += "  weight = {}".format(self.weight)
        res += "\n"
        res += "source: {}\n".format(self.source)
        res += "target: {}\n".format(self.target)

        return res


class CombinedEdge(object):
    """
    Single edge made of a group of edges glued together by intermediary nodes.
    """

    def __init__(self, midnodes, midedges, weight=1):
        """ Initialize class CombinedEdge. """

        self.midnodes = midnodes
        self.midedges = midedges
        self.weight = weight


    def __repr__(self):
        """ Representation of the CombinedEdge object. """

        res =  "MidNodes:\n"
        for midnode in self.midnodes:
            res+="{}\n".format(midnode.__repr__())
        res += "MidEdges:\n"
        for midedge in self.midedges:
            res+="{}\n".format(midedge.__repr__())

        if self.weight != None:
            res += "  weight = {}".format(self.weight)

        return res


class CausalGraph(object):
    """ Data structure for causal graphs. """

    def __init__(self, filename=None, eoi=None, combinedgraph=False,
                 nodestype="event", showintro=True):
        """ Initialize class CausalGraph. """

        # Header variables.
        self.filename = filename
        self.eoi = eoi
        self.combinedgraph = combinedgraph
        self.nodestype = nodestype # event or species
        self.showintro = showintro
        # Main variables.
        self.eventnodes = []
        self.causaledges = []
        self.edgegroups = []
        self.midnodes = []
        self.midedges = []
        self.combinededges = []
        # Post computed variables.
        self.coveredges = []
        self.occurrence = 1
        self.maxrank = None
        self.prevcores = None
        if self.filename != None:
            self.read_dot(self.filename)


    def read_dot(self, dotpath):
        """
        Read event nodes and causal edges from input causal graph.
        """

        rank = None
        self.label_mapping = {}
        dotfile = open(dotpath, "r").readlines()
        for line in dotfile:
            if 'combinedgraph="True"' in line:
                self.combinedgraph = True
            if "nodestype=" in line:
                type_index = line.index("nodestype")
                quote = line.rfind('"')
                self.nodestype = line[type_index+11:quote]
            if 'showintro="False"' in line:
                self.showintro = False
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
                close_quote = line[open_quote+1:].index('"')+open_quote+1
                medrank = float(line[open_quote+1:close_quote])
                rank = int(medrank)
            if "label=" in line and "Occurrence" not in line:
                if "->" not in line and "rank = same" not in line:
                    if 'cover="True"' not in line:
                        if line[0:2] == "//":
                           read_line = line[2:]
                        else:
                           read_line = line
                        tokens = read_line.split()
                        ori_id = tokens[0]
                        if '"' in ori_id:
                            ori_id = ori_id[1:-1]
                        if "node" not in ori_id:
                            node_id = "node{}".format(ori_id)
                        else:
                            node_id = ori_id
                        label_start = read_line.index("label=")+7
                        label_end = read_line[label_start:].index('"')+label_start
                        label = "{}".format(read_line[label_start:label_end])
                        if "intro=True" in read_line:
                            is_intro = True
                        else:
                            is_intro = False
                        if "first=True" in read_line:
                            is_first = True
                        else:
                            is_first = False
                        if "midnode=True" not in line:
                            self.eventnodes.append(EventNode(node_id, label,
                                                             rank,
                                                             intro=is_intro,
                                                             first=is_first))
                            self.label_mapping[node_id] = label
                        else:
                            interm_node = MidNode(ori_id, rank=medrank)
                            self.midnodes.append(interm_node)
        tmp_edges = []
        for line in dotfile:
            if "->" in line and 'style="invis"' not in line:
                if 'cover="True"' not in line and "tee" not in line:
                    if line[0:2] == "//":
                        read_line = line[2:]
                    else:
                        read_line = line
                    tokens = read_line.split()
                    source_id = tokens[0]
                    if '"' in source_id:
                        source_id = source_id[1:-1]
                    if "node" not in source_id and "and" not in source_id:
                        source_id = "node{}".format(source_id)
                    target_id = tokens[2]
                    if '"' in target_id:
                        target_id = target_id[1:-1]
                    if "node" not in target_id and "and" not in target_id:
                        target_id = "node{}".format(target_id)
                    for node in self.eventnodes:
                        if node.nodeid == source_id:
                            source = node
                        if node.nodeid == target_id:
                            target = node
                    for node in self.midnodes:
                        if node.nodeid == source_id:
                            source = node
                        if node.nodeid == target_id:
                            target = node
                    if "weight=" in line:
                        weight_start = read_line.index("weight=")+7
                        rem = read_line[weight_start:]
                        if "," in rem:
                            weight_end = rem.index(",")+weight_start
                        else:
                            weight_end = rem.index("]")+weight_start
                        weight = int(read_line[weight_start:weight_end])
                    else:
                        weight = 1
                    tmp_edges.append(CausalEdge(source, target, weight))
        for edge in tmp_edges:
            self.causaledges.insert(0, edge)
        self.postprocess()


    def postprocess(self):
        """
        Various stuff to do after reading a dot file. Includes the creation
        of hyperedges from the intermediary nodes and edges.
        """

        if self.combinedgraph == False:
            self.create_edgegroups()
            self.create_combinededges()
            for node in self.eventnodes:
                if "Intro" in node.label:
                    node.intro = True
            for c in self.combinededges:
                print("--------")
                for n in c.midnodes:
                    print(n)
                print("========")
                for m in c.midedges:
                    print(m)
            self.find_first_rules()
            self.rank_sequential()
            self.combinedgraph = True
        elif self.hypergraph == True:
           self.read_hyperedges()
        if self.eoi == None:
            self.get_maxrank()
            for node in self.nodes:
                if node.rank == self.maxrank:
                    self.eoi = node.label
        self.build_nointro()


    def create_edgegroups(self):
        """
        Group edges together when there is a path between them using only
        head-to-head and tail-to-tail connections.

        Example: The 3 edges in the following graph form a single group.

                 A  B
                 | /|
                 |/ |
                 C  D
        """

        edgescopy = self.causaledges.copy()
        while len(edgescopy) > 0:
            current_group = [edgescopy[0]]
            del(edgescopy[0])
            new_edge_found = True
            while new_edge_found == True:
                # Find sources and targets:
                sources = []
                targets = []
                for current_edge in current_group:
                    if current_edge.source not in sources:
                        sources.append(current_edge.source)
                    if current_edge.target not in targets:
                        targets.append(current_edge.target)
                # Find other edges with same source or target.
                new_edge_found = False
                copy_to_remove = []
                for i in range(len(edgescopy)):
                    other_source = edgescopy[i].source
                    other_target = edgescopy[i].target
                    if other_source in sources or other_target in targets:
                        new_edge_found = True
                        current_group.append(edgescopy[i])
                        copy_to_remove.insert(0, i)
                for i in copy_to_remove:
                    del(edgescopy[i])
                if new_edge_found == False:
                    self.edgegroups.append(current_group)


    def create_combinededges(self):
        """ Create combined edges from edge groups. """

        midid = 1
        for group in self.edgegroups:
            # Collect all sources and targets.
            sources = []
            targets = []
            for edge in group:
                if edge.source not in sources:
                    sources.append(edge.source)
                if edge.target not in targets:
                    targets.append(edge.target)
            # Create intermediary nodes for event nodes with more than one
            # input or output. Also create an edge between those intermediary
            # nodes and the corresponding event nodes.
            #node_map = {}
            midnodes = []
            midedges = []
            req_sources = []
            for source in sources:
                outgoing = []
                for edge in group:
                    if edge.source == source:
                        outgoing.append(edge)
                if len(outgoing) > 1:
                    req_sources.append(source)
                    midnodes.append(MidNode("mid{}".format(midid),
                                            midtype="requirement"))
                    #node_map[source.nodeid] = midnodes[-1]
                    midedges.append(MidEdge(source, midnodes[-1]))
                    midid += 1
            ena_targets = []
            for target in targets:
                ingoing = []
                for edge in group:
                    if edge.target == target:
                        ingoing.append(edge)
                if len(ingoing) > 1:
                    ena_targets.append(target)
                    midnodes.append(MidNode("mid{}".format(midid),
                                            midtype="enabling"))
                    #node_map[source.nodeid] = midnodes[-1]
                    midedges.append(MidEdge(midnodes[-1], target))
                    midid += 1
            # Add the intermediary edges corresponding to the original edges.
            for ori_edge in group:
                if ori_edge.source in req_sources:
                    for midedge in midedges:
                        if midedge.source == ori_edge.source:
                            s = midedge.target
                else:
                    s = ori_edge.source
                if ori_edge.target in ena_targets:
                    for midedge in midedges:
                        if midedge.target == ori_edge.target:
                            t = midedge.source
                else:
                    t = ori_edge.target
                midedges.append(MidEdge(s, t))

            #print("--------")
            #for e in group:
            #    print(e)
            #print("////////")
            #for n in midnodes:
            #    print(n)
            #print("========")
            #for m in midedges:
            #    print(m)
            w = group[0].weight
            self.combinededges.append(CombinedEdge(midnodes, midedges, w))


    def create_hyperedges(self):
        """
        Create hyperedges from the intermediary edges found in dot file.
        If the dot was not already an hypergraph, create an hyperedge for
        every group of edges that point to the same target.
        If the dot was already an hypergraph, use intermediary nodes to
        build hyperedges.
        """

        intermediary_id = 1
        for node in self.nodes:
            incoming_edges = []
            for edge in self.mededges:
                if edge.target == node:
                    incoming_edges.append(edge)
            if len(incoming_edges) > 0:
                source_list = []
                for in_edge in incoming_edges:
                    source_list.append(in_edge.source)
                sources = NodeGroup(source_list, "and")
                w = incoming_edges[0].weight
                mednode = IntermediaryNode("and{}".format(intermediary_id))
                intermediary_id += 1
                self.hyperedges.append(CausalEdge(sources, node, mednode, w))


    def read_hyperedges(self):
        """
        Rebuild hyperedges from intermediary nodes and edges.
        One hyperedge is created per intermediary node.
        """

        for node in self.mednodes:
            source_list = []
            for edge in self.mededges:
                if edge.target == node:
                    source_list.append(edge.source)
                if edge.source == node:
                    target = edge.target
                    w = edge.weight
            sources = NodeGroup(source_list, "and")
            self.hyperedges.append(CausalEdge(sources, target, node, w))
        intermediary_id = len(self.mednodes) + 1
        # Hyperedges with only one source.
        for edge in self.mededges:
            if "and" not in edge.source.nodeid:
                if "and" not in edge.target.nodeid:
                    w = edge.weight
                    src = NodeGroup([edge.source], "and")
                    mednode = IntermediaryNode("and{}".format(intermediary_id))
                    intermediary_id += 1
                    self.hyperedges.append(CausalEdge(src, edge.target,
                                                      mednode, w))


    def find_first_rules(self):
        """
        Find the rules of rank 1 as any rule which has only intros as
        incoming nodes. This should only be applied to acyclic graphs.
        """

        for node in self.eventnodes:
            if node.intro == False:
                incoming_nodes = []
                for edge in self.causaledges:
                    if edge.target == node:
                        incoming_nodes.append(edge.source)
                all_intro = True
                for incoming_node in incoming_nodes:
                    if incoming_node.intro == False:
                        all_intro = False
                        break
                if all_intro == True:
                    node.first = True


    def rank_sequential(self):
        """
        Find the rank of each node, starting with first nodes and then adding
        the other nodes sequentially as soon as the have one secured incoming
        hyperedge.
        """

        # Initialize ranks.
        current_nodes = []
        for node in self.eventnodes:
            if node.first == True:
                node.rank = 1
                current_nodes.append(node)
            else:
                node.rank = None
        for mednode in self.midnodes:
            mednode.rank = None
        while len(current_nodes) > 0:
            # 1) Gather all hyperedges that come out from a current_node.
            current_hyperedges = []
            for edge in self.combinededges:
                for current_node in current_nodes:
                    if current_node in edge.source.nodelist:
                        if edge not in current_hyperedges:
                            current_hyperedges.append(edge)
            # 2) Check which hyperedges are secured (hyperedges with sources
            #    having defined rank, ignoring intro nodes).
            secured_hyperedges = []
            for edge in current_hyperedges:
                secured = True
                for source in edge.source.nodelist:
                    if source.intro == False and source.rank == None:
                        secured = False
                        break
                if secured == True:
                    secured_hyperedges.append(edge)
            # 3) Set rank of the target of secured hyperedges.
            for edge in secured_hyperedges:
                if edge.target.rank == None:
                    source_ranks = []
                    for source in edge.source.nodelist:
                        if source.intro == False:
                            source_ranks.append(source.rank)
                    edge.target.rank = max(source_ranks)+1
                    current_nodes.append(edge.target)
            # 4) Remove current_nodes that point to nodes that are all ranked.
            next_nodes = []
            for current_node in current_nodes:
                keep_node = False
                node_targets = []
                for edge in self.combinededges:
                    if current_node in edge.source.nodelist:
                        if edge.target not in node_targets:
                            node_targets.append(edge.target)
                for node in node_targets:
                    if node.rank == None:
                        keep_node = True
                        break
                if keep_node == True:
                    next_nodes.append(current_node)
            current_nodes = next_nodes
        # Rank intro nodes.
        for node in self.eventnodes:
            if node.intro == True:
                # 1) Find all hyperedges that have current
                #    intro node as source.
                intro_hyperedges = []
                for edge in self.combinededges:
                    if node in edge.source.nodelist:
                        if edge not in intro_hyperedges:
                            intro_hyperedges.append(edge)
                # 2) Set rank of the intro node according to the lower
                #    rank among targets.
                target_ranks = []
                for edge in intro_hyperedges:
                    target_ranks.append(edge.target.rank)
                node.rank = min(target_ranks)-1
        self.rank_intermediary(self.hyperedges)
        self.get_maxrank()
        self.sequentialize_ids()


    def rank_intermediary(self, edgegroup):
        """  Rank intermediary nodes. """
    
        for edge in edgegroup:
            source_ranks = []
            for node in edge.source.nodelist:
                source_ranks.append(node.rank)
            if edge.target.rank > max(source_ranks):
                edge.mednode.rank = ( edge.target.rank + max(source_ranks) ) / 2.0
            elif edge.target.rank < min(source_ranks):
                edge.mednode.rank = ( edge.target.rank + min(source_ranks) ) / 2.0
            else:
                ranks_set = list(set(source_ranks))
                edge.mednode.rank = statistics.mean(ranks_set)


#    def rank_nodes(self):
#        """
#        Find the rank of each node as its longest upstream acyclic path
#        to an initial node, but ignoring upstream paths that contain any node
#        from a downstream acyclic path.
#        Remaining nodes are iteratively placed one rank lower than the lowest
#        rank of the nodes to which they point to and are already ranked.
#        """
#
#        print(self.filename)
#        #self.find_init_nodes()
#        # Rank rule nodes.
#        for node in self.nodes:
#            if node.intro == False:
#                #print(">>>", node.nodeid, node.label)
#                #upstream_paths = self.follow_edges("up", node, self.init_nodes)
#                upstream_paths = self.follow_edges("up", node)
#                if node.label == "SRCkin binds SRC-Y419 fast":
#                    print("UP")
#                    for path in upstream_paths:
#                        print("---")
#                        for n in path:
#                            print(n)
#                downstream_paths = self.follow_edges("down", node)
#                if node.label == "SRCkin binds SRC-Y419 fast":
#                    print("DOWN")
#                    for path in downstream_paths:
#                        print("---")
#                        for n in path:
#                            print(n)
#                paths_to_remove = []
#                for i in range(len(upstream_paths)):
#                    for down_path in downstream_paths:
#                        for down_node in down_path:
#                            if down_node in upstream_paths[i]:
#                                if i not in paths_to_remove:
#                                    paths_to_remove.insert(0, i)
#                if node.label == "SRCkin binds SRC-Y419 fast":
#                    print(paths_to_remove, "|||",len(upstream_paths))
#                #for i in paths_to_remove:
#                #    del(upstream_paths[i])
#                if node.label == "SRCkin binds SRC-Y419 fast":
#                    print("FINAL")
#                    for path in upstream_paths:
#                        print("---")
#                        for n in path:
#                            print(n)
#                lengths = []
#                for path in upstream_paths:
#                    lengths.append(len(path))
#                node.rank = max(lengths)
#        # Rank intro nodes.
#        for node in self.nodes:
#            if node.intro == True:
#                target_ranks = []
#                for edge in self.edges:
#                    if edge.source == node:
#                        target_ranks.append(edge.target.rank)
#                node.rank = min(target_ranks)-1
#        self.get_maxrank()
#        self.sequentialize_ids()
#
#
#    def old_rank_nodes(self):
#        """
#        Find the rank of each rule node in two steps. First, assign a temporary
#        rank as the longest acyclic path to an initial node (intro node of
#        rank 0). Second, repeat the first step but ignoring any edge for which
#        the temporary rank of the source is higher or equal to the temporary
#        rank of the target. Finally, assign the rank of intro nodes as the
#        lowest rank of any rule node it points to, minus 1.
#        """
#
#        #self.find_init_nodes()
#        self.init_nodes = []
#        for node in self.nodes:
#            if node.intro == True:
#                self.init_nodes.append(node)
#        # First pass.
#        for node in self.nodes:
#            if node.intro == False:
#                paths = self.climb_up(node, self.init_nodes)
#                lengths = []
#                for path in paths:
#                    lengths.append(len(path)-1)
#                node.rank = max(lengths)
#            else:
#                node.rank = 0
#        # Second pass ignoring higher or equal temporary ranks.
#        for node in self.nodes:
#            if node.intro == False:
#                paths = self.climb_up(node, self.init_nodes, ignore_heq=True)
#                lengths = []
#                for path in paths:
#                    lengths.append(len(path)-1)
#                node.rank = max(lengths)
#        # Ranking intro nodes.
#        for node in self.nodes:
#            if node.intro == True:
#                target_ranks = []
#                for edge in self.edges:
#                    if edge.source == node:
#                        target_ranks.append(edge.target.rank)
#                node.rank = min(target_ranks)-1
#        self.get_maxrank()
#        self.sequentialize_ids()


#    def find_init_nodes(self):
#        """
#        Find the nodes of rank 0 as the intro nodes which only point to rule
#        nodes that requires only intro nodes.
#        """
#
#        self.init_nodes = []
#        first_rules = []
#        for node in self.nodes:
#            if node.intro == False:
#                input_nodes = []
#                for edge in self.edges:
#                    if edge.target == node:
#                        input_nodes.append(edge.source)
#                all_intro = True
#                for input_node in input_nodes:
#                    if input_node.intro == False:
#                        all_intro = False
#                        break
#                if all_intro == True:
#                    first_rules.append(node)
#        for node in first_rules:
#            for edge in self.edges:
#                if edge.target == node:
#                    if edge.source not in self.init_nodes:
#                        self.init_nodes.append(edge.source)
#        nodes_to_remove = []
#        for i in range(len(self.init_nodes)):
#            for edge in self.edges:
#                if edge.source == self.init_nodes[i]:
#                    if edge.target not in first_rules:
#                        nodes_to_remove.insert(0, i)
#                        break
#        for i in nodes_to_remove:
#            del(self.init_nodes[i])


    def follow_edges(self, direction, from_node, to_nodes=[]):
        """
        Return a list of all acyclic paths from a given node to the top of the
        graph (using direction="up") or to the bottom (using direction="down").
        If to_nodes are provided, return only the paths that go from from_node
        to any of the to_nodes.
        """
    
        all_paths = [[from_node]]
        ends_reached = False
        while ends_reached == False:
            ends_reached = True
            for i in range(len(all_paths)):
                path = all_paths[i]
                next_nodes = []
                for edge in self.hyperedges:
                    if direction == "up":
                        if edge.target == path[-1]:
                            for src_node in edge.source.nodelist:
                                next_nodes.append(src_node)
                    elif direction == "down":
                        if path[-1] in edge.source.nodelist:
                            next_nodes.append(edge.target)
                if len(next_nodes) > 0 and path[-1] not in to_nodes:
                    ends_reached = False
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
        # Remove paths that do not end with one of the to_nodes if to_nodes
        # was defined.
        if len(to_nodes) > 0:
            for i in range(len(all_paths)-1, -1, -1):
                if all_paths[i][-1] not in to_nodes:
                    del(all_paths[i])
        # Remove the from_node in each path (the first node).
        for i in range(len(all_paths)):
            del(all_paths[i][0])
    
        return all_paths


#    def climb_up(self, bottom_node, top_nodes, ignore_heq=False):
#        """
#        Return the list of possible acyclic paths from bottom to top nodes.
#        Edges are followed in reverse, from target to source.
#        """
#    
#        all_paths = [[bottom_node]]
#        top_reached = False
#        while top_reached == False:
#            top_reached = True
#            for i in range(len(all_paths)):
#                path = all_paths[i]
#                if path[-1] not in top_nodes:
#                    top_reached = False
#                    next_nodes = []
#                    for edge in self.edges:
#                        if edge.target == path[-1]:
#                            if ignore_heq == False:
#                                next_nodes.append(edge.source)
#                            elif ignore_heq == True:
#                                if edge.source.rank < edge.target.rank:
#                                    next_nodes.append(edge.source)
#                    if len(next_nodes) == 0:
#                        path.append("root")
#                    else:
#                        path_copy = path.copy()
#                        path.append(next_nodes[0])
#                        for i in range(1, len(next_nodes)):
#                            new_path = path_copy.copy()
#                            new_path.append(next_nodes[i])
#                            all_paths.append(new_path)
#            # Remove looping paths.
#            for i in range(len(all_paths)-1, -1, -1):
#                if len(all_paths[i]) != len(set(all_paths[i])):
#                    del(all_paths[i])
#        for i in range(len(all_paths)-1, -1, -1):
#            if all_paths[i][-1] not in top_nodes:
#                del(all_paths[i])
#    
#        return all_paths
#
#
#    def slide_down(self, top_node, bottom_nodes):
#        """
#        Return the list of possible acyclic paths from top to bottom nodes.
#        Edges are followed from source to target.
#        """
#
#        all_paths = [[top_node]]
#        bottom_reached = False
#        while bottom_reached == False:
#            bottom_reached = True
#            for i in range(len(all_paths)):
#                path = all_paths[i]
#                if path[-1] not in bottom_nodes and path[-1] != "end":
#                    bottom_reached = False
#                    next_nodes = []
#                    for edge in self.edges:
#                        if edge.source == path[-1]:
#                            next_nodes.append(edge.target)
#                    if len(next_nodes) == 0:
#                        path.append("end")
#                    else:
#                        path_copy = path.copy()
#                        path.append(next_nodes[0])
#                        for i in range(1, len(next_nodes)):
#                            new_path = path_copy.copy()
#                            new_path.append(next_nodes[i])
#                            all_paths.append(new_path)
#            # Remove looping paths.
#            for i in range(len(all_paths)-1, -1, -1):
#                if len(all_paths[i]) != len(set(all_paths[i])):
#                    del(all_paths[i])
#        for i in range(len(all_paths)-1, -1, -1):
#            if all_paths[i][-1] not in bottom_nodes:
#                del(all_paths[i])
#
#        return all_paths


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
        # Also sort hyperedges.
        sorted_edges = sorted(self.hyperedges, key=lambda x: x.mednode.rank)
        self.hyperedges = sorted_edges


    def cleanup(self):
        """
        Remove nodes that do not have a path to an intro node and an eoi node
        and any edge that points to or from these nodes. Then remove intro
        nodes that do not have any targets anymore.
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
                #paths_up = self.climb_up(node, self.intro_nodes)
                paths_up = self.follow_edges("up", node, self.intro_nodes)
                #paths_down = self.slide_down(node, self.eoi_nodes)
                paths_down = self.follow_edges("down", node, self.eoi_nodes)
                if len(paths_up) == 0 or len(paths_down) == 0:
                    nodes_to_clean.insert(0, i)
        edges_to_clean = []
        for j in range(len(self.hyperedges)):
            sources = self.hyperedges[j].source.nodelist
            target = self.hyperedges[j].target
            for i in nodes_to_clean:
                node = self.nodes[i]
                if node in sources or target == node:
                    edges_to_clean.insert(0, j)
        for j in edges_to_clean:
            del(self.hyperedges[j])
        for i in nodes_to_clean:
            del(self.nodes[i])
        intros_to_clean = []
        for i in range(len(self.nodes)):
            node = self.nodes[i]
            if node.intro == True:
                remove_intro = True
                for edge in self.hyperedges:
                    if node in edge.source.nodelist:
                        remove_intro = False
                        break
                if remove_intro == True:
                    intros_to_clean.insert(0, i)
        for i in intros_to_clean:
            del(self.nodes[i])
            

    def build_nointro(self):
        """
        Create new hyperedges for the version of the graph that hides
        intro nodes.
        """

        # Reset information about cover edges. 
        for edge in self.hyperedges:
            edge.underlying = False
        self.coveredges = []
        intermediary_id = self.find_max_med_id()+1
        for edge1 in self.hyperedges:
            edge_group = []
            if edge1.underlying == False:
                has_intro1, src1 = self.check_intro_src(edge1)
                if has_intro1 == True:
                    t1 = edge1.target
                    edge_group.append(edge1)
                    # Check all other hyperedges that have the same sources
                    # and target as edge1, ignoring intro nodes.
                    # They will all be grouped inside a single nointro edge.
                    for edge2 in self.hyperedges:
                        if edge2.color == edge1.color:
                            if edge2 != edge1 and edge2.underlying == False:
                                has_intro2, src2 = self.check_intro_src(edge2)
                                t2 = edge2.target
                                if same_nodes(src1, src2) and t1 == t2:
                                    edge_group.append(edge2)
            if len(edge_group) > 0:
                # Compute weight of nointro edge as the sum of all its
                # underlying edges. Also mark edges that were used as
                # underlying.
                w = 0
                for edge in edge_group:
                    w += edge.weight
                    edge.underlying = True
                # Create a new nointro edge.
                if len(src1) > 0:
                    src = NodeGroup(src1, "and")
                    mednode = IntermediaryNode("and{}".format(intermediary_id))
                    intermediary_id += 1
                    self.coveredges.append(CausalEdge(src, t1, mednode, w,
                                                      color=edge1.color))
        self.rank_intermediary(self.coveredges)
                    

    def check_intro_src(self, chk_edge):
        """
        Check if an edge has intro nodes in its sources and return the
        list of non-intro nodes in its sources.
        """

        has_intro = False
        nodes_without_intro = []
        for node in chk_edge.source.nodelist:
            if node.intro == True:
                has_intro = True
            elif node.intro == False:
                nodes_without_intro.append(node)

        return has_intro, nodes_without_intro


    def find_max_med_id(self):
        """ Find the highest intermediary node id in the graph. """

        intermediary_ids = []
        for edge in self.hyperedges:
            intermediary_ids.append(int(edge.mednode.nodeid[3:]))
        max_med_id = max(intermediary_ids)

        return max_med_id


    def build_dot_file(self, edgelabels=False, showintro=True):
        """ build a dot file of the CausalGraph. """

        self.build_nointro()
        # Write info about graph.
        dot_str = 'digraph G{\n'
        dot_str += '  hypergraph="{}" ;\n'.format(self.hypergraph)
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
            dot_str += '  ranksep=0.25 ;\n'
        else:
            dot_str += '  ranksep=0.5 ;\n'
        # Compute some statistics to assign edge and intermediary node width.
        minpenwidth = 1
        medpenwidth = 3
        maxpenwidth = 6.5
        all_weights = []
        for edge in self.hyperedges:
            all_weights.append(edge.weight)
        average_weight = statistics.mean(all_weights)
        # Draw nodes.
        for int_rank in range((self.maxrank+1)*2):
            current_rank = int_rank/2
            rank_str = "{}".format(current_rank)
            if showintro == False and current_rank < 1:
                dot_str += "//"
            if int_rank%2 == 0:
                dot_str += ('{{ rank = same ; "{}" ['
                            'shape=plaintext];\n'.format(int(current_rank)))
            else:
                dot_str += ('{{ rank = same ; "{}" [label="", '
                            'shape=plaintext];\n'.format(current_rank))
            for node in self.nodes:
                if node.rank == current_rank:
                    node_shape = 'invhouse'
                    node_color = 'lightblue'
                    if node.intro == True:
                        node_shape = 'rectangle'
                        node_color = 'white'
                    if node.label == self.eoi:
                        node_shape = 'ellipse'
                        node_color = 'indianred2'
                    if self.nodestype == 'species':
                        node_shape = 'ellipse'
                    if showintro == False and node.intro == True:
                        dot_str += '//'
                    dot_str += ('"{}" [label="{}", '
                                .format(node.nodeid, node.label))
                    dot_str += 'shape={}, style=filled, '.format(node_shape)
                    dot_str += 'fillcolor={}'.format(node_color)
                    if node.intro == True:
                       dot_str += ', intro={}'.format(node.intro)
                    if node.first == True:
                       dot_str += ', first={}'.format(node.first)
                    dot_str += "] ;\n"
            # Draw intermediary nodes that emulate hyperedges if two
            # sources or more are drawn.
            for edge in self.hyperedges:
                if edge.mednode.rank == current_rank:
                    # Include the mednode no matter what, but comment it
                    # if showintro is False and edge is underlying.
                    if len(edge.source.nodelist) > 1:
                        if showintro == False and edge.underlying == True:
                            dot_str += '//'
                        ratio = edge.weight/average_weight
                        pensize = math.log(ratio,2) + medpenwidth
                        if pensize < minpenwidth:
                            pensize = minpenwidth
                        if pensize > maxpenwidth:
                            pensize = maxpenwidth
                        pensize = math.sqrt(pensize)/12
                        dot_str += '"{}" [label="", '.format(edge.mednode.nodeid)
                        dot_str += 'shape=point, style=filled, '
                        dot_str += 'color={}, '.format(edge.color)
                        dot_str += 'fillcolor={}, '.format(edge.color)
                        dot_str += 'hyperand=True, '
                        dot_str += 'width={}, height={}] ;\n'.format(pensize, pensize)
            # Intermediary nodes from cover edges, same as above but only
            # if showintro is False.
            if showintro == False:
                for cedge in self.coveredges:
                    if cedge.mednode.rank == current_rank:
                        if len(cedge.source.nodelist) > 1:
                            ratio = cedge.weight/average_weight
                            pensize = math.log(ratio,2) + medpenwidth
                            if pensize < minpenwidth:
                                pensize = minpenwidth
                            if pensize > maxpenwidth:
                                pensize = maxpenwidth
                            pensize = math.sqrt(pensize)/12
                            dot_str += '"{}" [label="", '.format(cedge.mednode.nodeid)
                            dot_str += 'shape=point, style=filled, '
                            dot_str += 'color={}, '.format(cedge.color)
                            dot_str += 'fillcolor={}, '.format(cedge.color)
                            dot_str += 'hyperand=True, '
                            dot_str += 'width={}, height={}, '.format(pensize, pensize)
                            dot_str += 'cover="True"] ;\n'
            # Close rank braces.
            if showintro == False and current_rank < 1:
                dot_str += "//"
            dot_str += "}\n"
        # Draw invisible ranking edges.
        for int_rank in range(self.maxrank*2):
            rank = int_rank/2
            if showintro == False and rank < 1:
                dot_str += '//'
            next_rank = rank+0.5
            if int_rank%2 == 0:
                rank_str = '{}'.format(int(rank))
                next_str = '{}'.format(next_rank)
            else:
                rank_str = '{}'.format(rank)
                next_str = '{}'.format(int(next_rank))
            dot_str += ('"{}" -> "{}" [style="invis"] ;\n'
                        .format(rank_str, next_str))
        # Draw edges. If there is only one source, draw a direct edge from the
        # source to the target. If there are several sources, draw edges from
        # each source to the mednode and from the mednode to the target.
        # The weight of each intermediary edge should be the same.
        for edge in self.hyperedges:
            sources = edge.source.nodelist
            edge_color = edge.color
            ratio = edge.weight/average_weight
            pensize = math.log(ratio,2) + medpenwidth
            if pensize < minpenwidth:
                pensize = minpenwidth
            if pensize > maxpenwidth:
                pensize = maxpenwidth
            if len(sources) == 1:
                if showintro == False and edge.underlying == True:
                    dot_str += "//" 
                dot_str += ('"{}" -> "{}" '
                            .format(sources[0].nodeid, edge.target.nodeid))
                dot_str += '[penwidth={}'.format(pensize)
                dot_str += ', color={}'.format(edge_color)
                if edgelabels == True:
                    dot_str += ', label="  {}"'.format(edge.weight)
                dot_str += ', weight={}] ;\n'.format(edge.weight)
            elif len(sources) > 1:
                for node in sources:
                    if showintro == False and edge.underlying == True:
                        dot_str += "//"
                    dot_str += ('"{}" -> "{}" '
                            .format(node.nodeid, edge.mednode.nodeid))
                    dot_str += '[dir=none, penwidth={}'.format(pensize)
                    dot_str += ', color={}'.format(edge_color)
                    dot_str += ', weight={}] ;\n'.format(edge.weight)
                # Intermediary edge to target.
                if showintro == False and edge.underlying == True:
                    dot_str += "//"
                dot_str += ('"{}" -> "{}" '
                            .format(edge.mednode.nodeid, edge.target.nodeid))
                dot_str += '[penwidth={}'.format(pensize)
                dot_str += ', color={}'.format(edge_color)
                if edgelabels == True:
                    dot_str += ', label="  {}"'.format(edge.weight)
                dot_str += ', weight={}] ;\n'.format(edge.weight)
        # Draw cover edges if intro nodes are not shown.
        if showintro == False:
            for cedge in self.coveredges:
                sources = cedge.source.nodelist
                edge_color = cedge.color
                ratio = cedge.weight/average_weight
                pensize = math.log(ratio,2) + medpenwidth
                if pensize < minpenwidth:
                    pensize = minpenwidth
                if pensize > maxpenwidth:
                    pensize = maxpenwidth
                if len(sources) == 1:
                    dot_str += ('"{}" -> "{}" '
                                .format(sources[0].nodeid, cedge.target.nodeid))
                    dot_str += '[penwidth={}'.format(pensize)
                    dot_str += ', color={}'.format(edge_color)
                    if edgelabels == True:
                        dot_str += ', label="  {}"'.format(cedge.weight)
                    dot_str += ', weight={}, '.format(cedge.weight)
                    dot_str += 'cover="True"] ;\n'
                elif len(sources) > 1:
                    for node in sources:
                        dot_str += ('"{}" -> "{}" '
                                .format(node.nodeid, cedge.mednode.nodeid))
                        dot_str += '[dir=none, penwidth={}'.format(pensize)
                        dot_str += ', color={}'.format(edge_color)
                        dot_str += ', weight={}, '.format(cedge.weight)
                        dot_str += 'cover="True"] ;\n'
                    # Intermediary edge to target.
                    dot_str += ('"{}" -> "{}" '
                                .format(cedge.mednode.nodeid, cedge.target.nodeid))
                    dot_str += '[penwidth={}'.format(pensize)
                    dot_str += ', color={}'.format(edge_color)
                    if edgelabels == True:
                        dot_str += ', label="  {}"'.format(cedge.weight)
                    dot_str += ', weight={}, '.format(cedge.weight)
                    dot_str += 'cover="True"] ;\n'
        # Close graph.
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
                   "--precedence-only",
                   "-o", "{}/causalcore-".format(eoi),
                   "{}".format(trace_path)))

# ---------------- End of Causal Cores Generation Section  --------------------

# ==================== Causal Cores Merging Section ===========================

def mergecores(eoi, causalgraphs=None, edgelabels=False, showintro=True,
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
                for j in range(len(current_core.hyperedges)):
                    equi_index = equi_edges[j]
                    w = causal_cores[i].hyperedges[equi_index].weight
                    current_core.hyperedges[j].weight += w
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
        graph.build_dot_file(edgelabels, showintro)
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
        graph2_indexes = list(range(len(graph2.hyperedges)))
        all_edges_found = True
        for edge1 in graph1.hyperedges:
            for i in graph2_indexes:
                edge2 = graph2.hyperedges[i]
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

    # Check sources.
    same_sources = same_nodes(edge1.source.nodelist, edge2.source.nodelist)
    # Check target.
    same_target = False
    if edge1.target.label == edge2.target.label:
        if edge1.target.rank == edge2.target.rank:
            same_target = True
    if same_sources == True and same_target == True:
        equi_edges = True
    else:
        equi_edges = False

    return equi_edges


def same_nodes(node_list1, node_list2):
    """ Find if two lists of nodes contain all the same nodes. """

    group1 = node_list1.copy()
    group2 = node_list2.copy()
    found1 = []
    found2 = []
    for i in range(len(group1)):
        for node2 in group2:
            if group1[i].label == node2.label:
                if group1[i].rank == node2.rank:
                    found1.insert(0, i)
                    break
    for j in range(len(group2)):
        for node1 in group1:
            if group2[j].label == node1.label:
                if group2[j].rank == node1.rank:
                    found2.insert(0, j)
                    break
    for i in found1:
        del(group1[i])
    for j in found2:
        del(group2[j])
    if len(group1) == 0 and len(group2) == 0:
        are_same = True
    else:
        are_same = False

    return are_same


# ================ End of Causal Cores Merging Section ========================

# +++++++++++++++++++++++ Cores Looping Section +++++++++++++++++++++++++++++++

def loopcores(eoi, causalgraphs=None, ignorelist=None, edgelabels=False,
              showintro=False, writedots=True, rmprev=False,
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
    flush_ignored(cores, core_files, ignorelist)
    for core in cores:
        #remove_ignored(core, ignorelist)
        merge_same_labels(core)
        fuse_edges(core)
        core.rank_sequential()
        #core.old_rank_nodes()
        core.build_dot_file(edgelabels, showintro)
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
                              edgelabels=edgelabels, showintro=showintro,
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


def flush_ignored(graph_list, graph_files, ignorelist):
    """
    Temporary fix to completely remove all cores that contain
    an ignored term in any of its node.
    """

    init_len = len(graph_list)
    graphs_to_remove = []
    for i in range(len(graph_list)):
        graph = graph_list[i]
        remove_graph = False
        for node in graph.nodes:
            if any(ignorestr in node.label for ignorestr in ignorelist):
                remove_graph = True
                break
        if remove_graph == True:
            graphs_to_remove.insert(0, i)
    for i in graphs_to_remove:
        if graph_files != None:
           slash = graph_list[i].filename.index("/")
           fname = graph_list[i].filename[slash+1:]
           graph_files.remove(fname)
        del(graph_list[i])
    print("Ignoring {} cores out of {} because they contain reverse rules."
          .format(len(graphs_to_remove), init_len))


def remove_ignored(graph, ignorelist):
    """ Remove nodes whose label contains a string defined in ignorelist. """

    ignored_nodes = []
    for i in range(len(graph.nodes)):
        if any(ignorestr in graph.nodes[i].label for ignorestr in ignorelist):
            ignored_nodes.insert(0, i)
    for i in ignored_nodes:
        del(graph.nodes[i])
    ignored_edges = []
    for i in range(len(graph.hyperedges)):
        source = graph.edges[i].source.label
        target = graph.edges[i].target.label
        if any(ignorestr in source for ignorestr in ignorelist):
            ignored_edges.insert(0, i)
        if any(ignorestr in target for ignorestr in ignorelist):
            if i not in ignored_edges:
                ignored_edges.insert(0, i)
    for i in ignored_edges:
        del(graph.hyperedgesedges[i])
    graph.cleanup()


def merge_same_labels(graph):
    """
    Merge every node with same label within a graph. Merging between two
    nodes is done by deleting the second node and changing any edge which was
    coming from or going to the second node to now come from or go to the
    second node.
    If two nodes with same label appear as sources of a same hyperedge, 
    the intermediary edge from the second node is simply removed.
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


def merge_nodes(nodes_to_merge, graph):
    """ Merge every node from the list onto the node of lowest rank. """

    ranks = []
    for node in nodes_to_merge:
        ranks.append(node.rank)
    lowest_rank = min(ranks)
    for i in range(len(nodes_to_merge)):
        if nodes_to_merge[i].rank == lowest_rank:
            main_node = nodes_to_merge[i]
            del(nodes_to_merge[i])
            break
    for edge in graph.hyperedges:
        # Replace sources.
        sources_found = False
        nodes_to_remove = []
        for i in range(len(edge.source.nodelist)):
            src = edge.source.nodelist[i]
            if src in nodes_to_merge:
                nodes_to_remove.insert(0, i)
                sources_found = True
        if sources_found == True:
            for i in nodes_to_remove:
                del(edge.source.nodelist[i])
            label_already_found = False
            for node in edge.source.nodelist:
                if node.label == main_node.label:
                    label_already_found = True
            if label_already_found == False:
                edge.source.nodelist.append(main_node)
        # Replace target.
        if edge.target in nodes_to_merge:
            edge.target = main_node
    for i in range(len(graph.nodes)-1, -1, -1):
        if graph.nodes[i] in nodes_to_merge:
            del(graph.nodes[i])
            

def fuse_edges(graph):
    """ Remove duplicate edges between two same nodes but sum weights. """

    unique_edges = []
    for edge1 in graph.hyperedges:
        sources1 = edge1.source.nodelist
        target1 = edge1.target.nodeid
        new_edge = True
        for edge2 in unique_edges:
            sources2 = edge2.source.nodelist
            target2 = edge2.target.nodeid
            same_sources = same_nodes(sources1, sources2)
            if same_sources == True and target1 == target2:
                new_edge = False
                break
        if new_edge == True:
            unique_edges.append(edge1)
    for unique_edge in unique_edges:
        unique_sources = unique_edge.source.nodelist
        unique_target = unique_edge.target.nodeid
        w = 0
        for edge in graph.hyperedges:
            sources = edge.source.nodelist
            target = edge.target.nodeid
            same_sources = same_nodes(unique_sources, sources)
            if same_sources == True and unique_target == target:
                w += edge.weight
        unique_edge.weight = w
    graph.hyperedges = unique_edges

# ++++++++++++++++++ End of Cores Looping Section +++++++++++++++++++++++++++++

# .................. Event Paths Merging Section ..............................

def mergepaths(eoi, causalgraphs=None, ignorelist=None, edgelabels=False,
               showintro=False, writedot=True, rmprev=False):
    """ Merge event paths into a single pathway. """

    # Reading section.
    if causalgraphs == None:
        # Using cores or eventpaths both work. But it can suffle the nodes
        # horizontally, yielding a different graph, but with same ranks for
        # all nodes.
        #path_files = get_dot_files(eoi, "eventpath")
        path_files = get_dot_files(eoi, "core")
        event_paths = []
        for path_file in path_files:
            path_path = "{}/{}".format(eoi, path_file)
            event_paths.append(CausalGraph(path_path, eoi))
    else:
        event_paths = causalgraphs
        path_files = None
    # Doing the work.
    flush_ignored(event_paths, path_files, ignorelist)
    pathway = CausalGraph(eoi=eoi, hypergraph=True)
    pathway.occurrence = 0
    node_number = 1
    seen_labels = []
    #lowest_ranks = {}
    for event_path in event_paths:
        pathway.occurrence += event_path.occurrence
        for node in event_path.nodes:
            if node.label not in seen_labels:
                seen_labels.append(node.label)
    #            lowest_ranks[node.label] = node.rank
                n_id = "node{}".format(node_number)
                pathway.nodes.append(CausalNode(n_id, node.label,
                                                node.rank,
                                                intro=node.intro,
                                                first=node.first))
                node_number += 1
    #        else:
    #            if node.rank < lowest_ranks[node.label]:
    #                lowest_ranks[node.label] = node.rank
    #for node in pathway.nodes:
    #    node.rank = lowest_ranks[node.label]
    intermediary_id = 1
    for event_path in event_paths:
        for edge in event_path.hyperedges:
            # For each source and target node of the edge in event_path,
            # find the equivalent node in the pathway.
            source_labels = []
            for source_node in edge.source.nodelist:
                source_labels.append(source_node.label)
            source_list = []
            for node in pathway.nodes:
                if node.label in source_labels:
                    source_list.append(node)
                if node.label == edge.target.label:
                    target = node
            sources = NodeGroup(source_list, "and")
            mednode = IntermediaryNode("and{}".format(intermediary_id))
            intermediary_id += 1
            pathway.hyperedges.append(CausalEdge(sources, target, mednode,
                                            edge.weight))
    fuse_edges(pathway)
    # Uncomment the next 3 lines and comment pathway.rank_sequential()
    # to build unranked version of graph
    #pathway.rank_intermediary()
    #pathway.get_maxrank()
    #pathway.sequentialize_ids()
    pathway.rank_sequential()
    pathway.filename = "eventpathway.dot"
    pathway.build_dot_file(edgelabels, showintro)
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

# ///////////////// Event Paths Simplifying Section ///////////////////////////

def simplifypath(eoi, causalgraphs=None, threshold=0.2, edgelabels=False,
                 showintro=False, writedot=True, rmprev=False):
    """
    Simplify event pathway by ignoring eventpaths that contain edges with
    low occurrence.
    """

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
    pathway = CausalGraph("{}/eventpathway.dot".format(eoi), eoi)
    # Doing the work.
    # Remove eventpaths that contain edges with an occurence that
    # is under average_weight*threshold.
    pathway.build_nointro()
    all_weights = []
    for edge in pathway.hyperedges:
        if edge.underlying == False:
            all_weights.append(edge.weight)
    for cedge in pathway.coveredges:
        all_weights.append(cedge.weight)
    average_weight = statistics.mean(all_weights)
    theshold_str = "Average edge weight: {:.2f} , ".format(average_weight)
    theshold_str += "Treshold: {:.2f} , ".format(threshold)
    theshold_str += "Cutoff: {:.2f}\n".format(average_weight*threshold)
    theshold_str += "Simplifying graph; ignoring story types containing "
    theshold_str += ("edges with weight lower than {:.2f}"
                      .format(average_weight*threshold))
    print(theshold_str)
    normals_to_ignore = []
    for edge in pathway.hyperedges:
        if edge.underlying == False:
            if edge.weight < average_weight*threshold:
                normals_to_ignore.append(edge)
    covers_to_ignore = []
    for cedge in pathway.coveredges:
        if cedge.weight < average_weight*threshold:
            covers_to_ignore.append(cedge)
    # Select only eventpaths that do not contain any of the edges to ignore.
    selected_paths = []
    for event_path in event_paths:
        is_ignored = False
        for edge in event_path.hyperedges:
            is_ignored = ignored_edge(edge, normals_to_ignore)
            if is_ignored == True:
                break
        # If no edge to ignore was found in normal edges, check
        # also cover edges.
        if is_ignored == False:
            event_path.build_nointro()
            for cedge in event_path.coveredges:
                is_ignored = ignored_edge(cedge, covers_to_ignore)
                if is_ignored == True:
                    break
        if is_ignored == False:
            selected_paths.append(event_path)
    if len(selected_paths) == 0:
       raise ValueError("Simplification threshold too high, no story left.") 
    # Build simplified graph.
    simplepathway = CausalGraph(eoi=eoi, hypergraph=True)
    simplepathway.occurrence = 0
    node_number = 1
    seen_labels = []
    for event_path in selected_paths:
        simplepathway.occurrence += event_path.occurrence
        for node in event_path.nodes:
            if node.label not in seen_labels:
                seen_labels.append(node.label)
                n_id = "node{}".format(node_number)
                simplepathway.nodes.append(CausalNode(n_id, node.label,
                                                      node.rank,
                                                      intro=node.intro,
                                                      first=node.first))
                node_number += 1
    intermediary_id = 1
    for event_path in selected_paths:
        for edge in event_path.hyperedges:
            # For each source and target node of the edge in event_path,
            # find the equivalent node in the pathway.
            source_labels = []
            for source_node in edge.source.nodelist:
                source_labels.append(source_node.label)
            source_list = []
            for node in simplepathway.nodes:
                if node.label in source_labels:
                    source_list.append(node)
                if node.label == edge.target.label:
                    target = node
            sources = NodeGroup(source_list, "and")
            mednode = IntermediaryNode("and{}".format(intermediary_id))
            intermediary_id += 1
            simplepathway.hyperedges.append(CausalEdge(sources, target,
                                                       mednode,
                                                       edge.weight))
    fuse_edges(simplepathway)
    simplepathway.rank_sequential()
    simplepathway.filename = "eventpathway-simple.dot"
    simplepathway.build_dot_file(edgelabels, showintro)
    # Writing section.
    if writedot == True:
        output_path1 = "{}/{}".format(eoi, simplepathway.filename)
        outfile1 = open(output_path1, "w")
        outfile1.write(simplepathway.dot_file)
        outfile1.close()
    if rmprev == True:
        if path_files == None:
            path_files = get_dot_files(eoi, "evpath")
        for path_file in path_files:
            file_path = "{}/{}".format(eoi, path_file)
            os.remove(file_path)

    return pathway        


def ignored_edge(edge, ignore_list):
    """ Check if edge is contained in edges to ignore based on labels. """

    is_ignored = False
    trg_lbl = edge.target.label
    src_lbls = []
    for node in edge.source.nodelist:
        src_lbls.append(node.label)
    for to_ignore in ignore_list:
        ign_trg = to_ignore.target.label
        ign_srcs = []
        for node in to_ignore.source.nodelist:
            ign_srcs.append(node.label)
        if sorted(src_lbls) == sorted(ign_srcs) and trg_lbl == ign_trg:
            is_ignored = True
            break

    return is_ignored

# ///////////// End of Event Paths Simplifying Section //////////////////////

# ++++++++++++++ Core Mapping on Event Pathway Section ++++++++++++++++++++++

def mapcores(eoi, causalgraphs=None, ignorelist=None, template=None,
             edgelabels=False, showintro=False, writedot=True, rmprev=False):
    """
    Create a new CausalGraph for each core, where the core is mapped on
    the layout of the event pathway.
    """

    # Reading cores.
    if causalgraphs == None:
        path_files = get_dot_files(eoi, "core")
        event_paths = []
        for path_file in path_files:
            path_path = "{}/{}".format(eoi, path_file)
            event_paths.append(CausalGraph(path_path, eoi))
    else:
        event_paths = causalgraphs
        path_files = None
    # Reading event pathway template.
    if template == None:
        template_path = "{}/eventpathway.dot".format(eoi)
    else:
        template_path = template
    mappedcores = []
    # Doing the work.
    flush_ignored(event_paths, path_files, ignorelist) 
    for event_path in event_paths:
        mappedcore = CausalGraph(template_path, eoi)
        mappedcore.occurrence = event_path.occurrence
        for edge in mappedcore.hyperedges:
            edge.color = "grey90"
            edge.weight = mappedcore.occurrence
        edges_to_add = []
        intermediary_id = mappedcore.find_max_med_id()+1
        for i in range(1, event_path.maxrank+1):
            rank_edges = []
            for edge in event_path.hyperedges:
               if edge.target.rank == i:
                   rank_edges.append(edge)
            # Find the nodes with same labels in the event pathway template.
            for rank_edge in rank_edges:
                source_labels = []
                for node in rank_edge.source.nodelist:
                    source_labels.append(node.label)
                mapped_sources = []
                mapped_target = None
                for node in mappedcore.nodes:
                    if node.label in source_labels:
                        mapped_sources.append(node)
                    if node.label == rank_edge.target.label:
                        mapped_target = node
                # Remove original edge in event pathway.
                edges_to_remove = [] 
                for j in range(len(mappedcore.hyperedges)):
                    edge_found = False
                    path_target = mappedcore.hyperedges[j].target
                    if mapped_target == path_target:
                        edge_found = True
                        path_srcs = mappedcore.hyperedges[j].source.nodelist
                        for mapped_source in mapped_sources:
                            if mapped_source not in path_srcs:
                                edge_found = False
                                break
                        for path_src in path_srcs:
                            if path_src not in mapped_sources:
                                edge_found = False
                                break
                    if edge_found == True:
                        edges_to_remove.insert(0, j)
                for j in edges_to_remove:
                    del(mappedcore.hyperedges[j])
                # Add a new edge from mapped_sources to mapped_target.
                mednode = IntermediaryNode("and{}".format(intermediary_id))
                intermediary_id += 1
                w = event_path.occurrence
                sources = NodeGroup(mapped_sources, "and")
                maxr = float(event_path.maxrank)
                col_val = 0.25+(i/maxr)*0.75
                edge_color = '"{:.3} 1 1"'.format(col_val)
                edges_to_add.append(CausalEdge(sources, mapped_target,
                                               mednode, w, color=edge_color))
        for edge_to_add in edges_to_add:
            mappedcore.hyperedges.append(edge_to_add)
        mappedcore.rank_sequential()
        mappedcores.append(mappedcore)
    for i in range(len(mappedcores)):
        mappedcores[i].filename = "mapped-{}.dot".format(i+1)
    for graph in mappedcores:
        graph.build_dot_file(edgelabels, showintro)
    # Writing section.
    if writedot == True:
        for graph in mappedcores:
            output_path = "{}/{}".format(eoi, graph.filename)
            outfile = open(output_path, "w")
            outfile.write(graph.dot_file)
            outfile.close()
    if rmprev == True:
        if path_files == None:
            path_files = get_dot_files(eoi, "core")
        for path_file in path_files:
            file_path = "{}/{}".format(eoi, path_file)
            os.remove(file_path)


# ++++++++++++ End of Core Mapping on Event Pathway Section ++++++++++++++++++

# """"""""""""""" Species Pathway Conversion Section """"""""""""""""""""""""""

def speciespathway3(eoi, kappamodel, causalgraph=None, edgelabels=False,
                   showintro=True):
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
    init_species = get_init_species(event_pathway)
    new_edges = build_species_and_edges(event_pathway)
    species_pathway = build_species_graph(eoi, event_pathway, new_edges)
    clean_remaining_intros(species_pathway)
    change_intro_species(species_pathway)
    remove_initial_bnd(species_pathway, init_species)
    merge_same_labels(species_pathway)
    fuse_edges(species_pathway)
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


def remove_initial_bnd(graph, init_species):
    """ Remove all binding species that were present as init. """

    nodes_to_remove = []
    for i in range(len(graph.nodes)):
        node = graph.nodes[i]
        if species_in(node.species, init_species):
            nodes_to_remove.insert(0, i)
    edges_to_remove = []
    for j in range(len(graph.edges)):
        source = graph.edges[j].source
        target = graph.edges[j].target
        for i in nodes_to_remove:
            node = graph.nodes[i]
            if source == node or target == node:
                if j not in edges_to_remove:
                    edges_to_remove.insert(0, j)
    for j in edges_to_remove:
        del(graph.edges[j])
    for i in nodes_to_remove:
        del(graph.nodes[i])


def remove_free_events(graph):
    """ Remove all free species. """

    nodes_to_remove = []
    for i in range(len(graph.nodes)):
        node = graph.nodes[i]
        if node.intro == False:
            if node.species["bound_agent"] == ".":
                nodes_to_remove.insert(0, i)
    edges_to_remove = []
    for j in range(len(graph.edges)):
        source = graph.edges[j].source
        target = graph.edges[j].target
        for i in nodes_to_remove:
            node = graph.nodes[i]
            if source == node or target == node:
                if j not in edges_to_remove:
                    edges_to_remove.insert(0, j)
    for j in edges_to_remove:
        del(graph.edges[j])
    for i in nodes_to_remove:
        del(graph.nodes[i])


def clean_remaining_intros(graph):
    """
    Remove intro nodes if they represent a state and their agent is the
    same as that of their target.
    """

    intros_to_remove = []
    for i in range(len(graph.nodes)):
        node = graph.nodes[i]
        if node.intro == True and node.species["state"] != None:
            remove_node = False
            for edge in graph.edges:
                if edge.source == node:
                    target = edge.target
                    if target.species["agent"] == node.species["agent"]:
                        remove_node = True
                        break
            if remove_node == True:
                intros_to_remove.insert(0, i)
    edges_to_remove = []
    for j in range(len(graph.edges)):
        source = graph.edges[j].source
        target = graph.edges[j].target
        for i in intros_to_remove:
            node = graph.nodes[i]
            if source == node or target == node:
                edges_to_remove.insert(0, j)
    for j in edges_to_remove:
        del(graph.edges[j])
    for i in intros_to_remove:
        del(graph.nodes[i])


def change_intro_species(graph):
    """
    Change the label of intro nodes to the req of their target instead
    of their own res.
    """

    for node in graph.nodes:
        if node.intro == True:
            target_reqs = []
            for edge in graph.edges:
                if edge.source == node:
                    for req in edge.target.full_req:
                        target_reqs.append(req)
            n_ag = node.species["agent"]
            n_site = node.species["site"]
            n_bnd = node.species["bound_agent"]
            n_state = node.species["state"]
            for req in target_reqs:
                t_ag = req["agent"]
                t_site = req["site"]
                t_bnd = req["bound_agent"]
                t_state =req["state"]
                if n_ag == t_ag and n_site == t_site:
                    if n_bnd == None and t_bnd == None: # This is a state.
                        new_lbl = "{}({}{{{}}})".format(t_ag, t_site, t_state)
                        node.species["state"] = t_state
                    if n_state == None and t_state == None: # This is a bind.
                        new_lbl = "{}({}".format(t_ag, t_site)
                        if t_bnd == ".":
                            new_lbl += "[.]"
                        new_lbl += ")"
                        node.species["bound_agent"] = t_bnd
                        node.species["bound_site"] = req["bound_site"]
            node.label = new_lbl


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
    final_edges = []
    for edge in edge_list:
        if edge.target in species_pathway.nodes:
            final_edges.append(edge)
    species_pathway.edges = final_edges

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
            if mod_node.intro == False:
                new_node.full_req = mod_node.full_req
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
                    #if species_in(src.species, target_mod_node.full_req):
                    src_ag = src.species["agent"]
                    src_site = src.species["site"]
                    src_bnd = src.species["bound_agent"]
                    src_bndsite = src.species["bound_site"]
                    src_state = src.species["state"]
                    add_edge = False
                    for target_req in target_mod_node.full_req:
                        trg_ag = target_req["agent"]
                        trg_site = target_req["site"]
                        trg_bnd = target_req["bound_agent"]
                        trg_bndsite = target_req["bound_site"]
                        trg_state = target_req["state"]
                        if src_ag == trg_ag and src_site == trg_site:
                            if src_bnd == None and trg_bnd == None:
                                if mod_node.intro == True:
                                    add_edge = True
                                elif src_state == trg_state:
                                    add_edge = True
                            elif src_state == None and trg_state == None:
                                if mod_node.intro == True:
                                    add_edge = True
                                elif src_bnd == trg_bnd:
                                    if src_bndsite == trg_bndsite:
                                        add_edge = True
                    if add_edge == True:
                        for trgt in target_mod_node.species_nodes:
                            new_edges.append(CausalEdge(src, trgt,
                                                        weight=edge.weight))

    return new_edges


def get_init_species(graph):
    """ Create a list of all initial species. """

    init_species = []
    for node in graph.nodes:
        if node.intro == True:
            for res in node.res_species:
                init_species.append(res.copy())

    return init_species


def simplify_req_res(graph):
    """
    (Remove bindings and unbindings from res.** Not removed anymore)
    Remove all binding req where the agent is found in res, except if there is only
    one agent type in req.
    """

    for event_node in graph.nodes:
        if event_node.intro == False:
            ##################################
            res_to_remove = []
            for i in range(len(event_node.res_species)):
                if event_node.res_species[i]["state"] == None:
                    res_to_remove.insert(0, i)
            for i in res_to_remove:
                del(event_node.res_species[i])
            ##################################
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
                        if event_node.full_req[i]["state"] == None:
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
        paths = graph.follow_edges("up", mod_node, top_nodes)
        all_reqs = []
        for path in paths:
            path_reqs = []
            for i in range(len(path)-1):
                current_node = path[i]
                #for current_req in current_node.req_species:
                #    if not species_in(current_req, path_reqs):
                #        path_reqs.append(current_req.copy())
                for current_req in current_node.req_species:
                    cur_ag = current_req["agent"]
                    cur_site = current_req["site"]
                    cur_bnd = current_req["bound_agent"]
                    cur_state = current_req["state"]
                    add_req = True
                    for path_req in path_reqs:
                        path_ag = path_req["agent"]
                        path_site = path_req["site"]
                        path_bnd = path_req["bound_agent"]
                        path_state = path_req["state"]
                        if path_ag == cur_ag and path_site == cur_site:
                            if path_bnd == None and cur_bnd == None:
                                add_req = False
                            if path_state == None and cur_state == None:
                                add_req = False
                    if add_req == True:
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

        comp_req_set = []
        for path_reqs in all_reqs:
            for path_req in path_reqs:
                if not species_in(path_req, comp_req_set):
                    comp_req_set.append(path_req)
        mod_node.full_req = comp_req_set
    #for n in mod_nodes:
    #    if n.intro == False:
    #        print(n.rule)
    #        for req in n.full_req:
    #            print(req)
            


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
            for i in range(len(graph.hyperedges)-1, -1, -1):
                edge = graph.hyperedges[i]
                if edge.target == node:
                    up_edges.append(edge)
                    del(graph.hyperedges[i])
            down_edges = []
            for i in range(len(graph.hyperedges)-1, -1, -1):
                edge = graph.hyperedges[i]
                if edge.source == node:
                    down_edges.append(edge)
                    del(graph.hyperedges[i])
            for up_edge in up_edges:
                for down_edge in down_edges:
                    new_edge = CausalEdge(up_edge.source, down_edge.target,
                                          up_edge.weight)
                    graph.hyperedges.append(new_edge)
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
        if "dot" in current_file:
            #if "path" in current_file or "looped" in current_file:
            toggle_intro_nodes(".", current_file)
    directories = filter(os.path.isdir, os.listdir('.'))
    for directory in directories:
        dot_files = os.listdir("{}".format(directory))
        for dot_file in dot_files:
            if "dot" in dot_file:
                #if "path" in dot_file or "looped" in dot_file:
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
