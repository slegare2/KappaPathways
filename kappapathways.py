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
import random
import copy
import textwrap
import time


class EventNode(object):
    """
    An event node to use in causal graphs. It represents a specific event in
    causal cores and an event type (a rule) in pathways.
    """

    def __init__(self, nodeid, label, rank=None, weight=1, rel_wei=1.0,
                 occurrence=1, rel_occ=1.0, intro=False, first=False,
                 highlighted=False, pos=None, eventid=None):
        """ Initialize class EventNode. """

        self.nodeid = nodeid
        self.label = label
        self.rank = rank
        self.weight = weight # Taken from the stories.
        self.rel_wei = rel_wei # = weight / occurrence_of_EOI (num. of cores)
        self.occurrence = occurrence # Taken from the trace.
        self.rel_occ = rel_occ # = occurrence / occurrence_of_EOI
        self.intro = intro
        self.first = first
        self.highlighted = highlighted
        self.pos = pos
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
        if self.occurrence != None:
            if not isinstance(self.occurrence, int):
                raise TypeError("occurrence should be an integer.")


    def __repr__(self):
        """ Representation of the EventNode object. """

        res =  "Node "
        res += 'id: "{}",  label: "{}"'.format(self.nodeid, self.label)
        if self.rank != None:
            res += ",  rank: {}".format(self.rank)
        if self.occurrence != None:
            res += ",  occurrence: {}".format(self.occurrence)
        res += ",  intro: {}".format(self.intro)
        res += ",  first: {}".format(self.first)

        return res


class StateNode(object):
    """
    A state node to use in causal graphs. It represents the state that
    are changed by an event.
    """

    def __init__(self, nodeid, state, label, rank=None, weight=1, rel_wei=1.0,
                 occurrence=1, rel_occ=1.0, intro=False, first=False,
                 highlighted=False, pos=None, eventid=None):
        """ Initialize class StateNode. """

        self.nodeid = nodeid
        self.state = state
        self.label = label
        self.rank = rank
        self.weight = weight # Taken from the stories.
        self.rel_wei = rel_wei # = weight / occurrence_of_EOI (num. of cores)
        self.occurrence = occurrence # Taken from the trace.
        self.rel_occ = rel_occ # = occurrence / occurrence_of_EOI
        self.intro = intro
        self.first = first
        self.highlighted = highlighted
        self.pos = pos
        self.check_types()


    def check_types(self):
        """ Check that StateNode attributes have proper types. """

        if not isinstance(self.nodeid, str):
            raise TypeError("nodeid should be a string.")
        if not isinstance(self.label, str):
            raise TypeError("label should be a string.")
        #if self.rank != None:
        #    if not isinstance(self.rank, int):
        #        raise TypeError("rank should be an integer.")
        if self.occurrence != None:
            if not isinstance(self.occurrence, int):
                raise TypeError("occurrence should be an integer.")


    def __repr__(self):
        """ Representation of the StateNode object. """

        res =  "Node "
        res += 'id: "{}",  label: "{}"'.format(self.nodeid, self.label)
        if self.rank != None:
            res += ",  rank: {}".format(self.rank)
        if self.occurrence != None:
            res += ",  occurrence: {}".format(self.occurrence)
        res += ",  intro: {}".format(self.intro)
        res += ",  first: {}".format(self.first)

        return res


class CausalEdge(object):
    """
    A relationship between two event nodes in causal graphs. The relationship
    can be precedence (default), causal or conflict.
    """

    def __init__(self, source, target, weight=1, rel_wei=1.0, occurrence=None,
                 rel_occ=None, relationtype="causal", color="black",
                 underlying=False, reverse=False, labelcarrier=True,
                 indicator=False, meshid=None, pos=None, labelpos=None,
                 overridewidth=None, overridelabel=None):
        """ Initialize class CausalEdge. """

        self.source = source
        self.target = target
        self.weight = weight # Taken from the stories.
        self.rel_wei = rel_wei # = weight / occurrence_of_EOI (num. of cores)
        self.occurrence = occurrence # Taken from the trace.
        self.rel_occ = rel_occ # = occurrence / occurrence_of_EOI
        self.relationtype = relationtype
        self.color = color
        self.underlying = underlying
        self.reverse = reverse
        self.labelcarrier = labelcarrier
        self.indicator = indicator
        self.meshid = meshid
        self.pos = pos
        self.labelpos = labelpos
        self.overridewidth = overridewidth
        self.overridelabel = overridelabel
        self.check_types()


    def check_types(self):
        """ Check that CausalEdge attributes have proper types. """

        if not isinstance(self.source, EventNode):
            if not isinstance(self.source, StateNode):
                raise TypeError("source should be an EventNode or StateNode.")
        if not isinstance(self.target, EventNode):
            if not isinstance(self.target, StateNode):
                raise TypeError("target should be an EventNode or StateNode.")
        if self.weight != None:
            if not isinstance(self.weight, int):
                raise TypeError("uses should be an integer.")
        if self.occurrence != None:
            if not isinstance(self.occurrence, int):
                raise TypeError("occurrence should be an integer.")


    def __repr__(self):
        """ Representation of the CausalEdge object. """

        res =  "Edge"
        if self.weight != None:
            res += "  weight = {}".format(self.weight)
        if self.occurrence != None:
            res += "  occurrence = {:.3f}".format(self.occurrence)
        res += "\n"
        res += "source: {}\n".format(self.source)
        res += "target: {}\n".format(self.target)

        return res


class HyperEdge(object):
    """
    Hyperedge implemented as a set of CausalEdges, all with the same
    target EventNode.
    """

    def __init__(self, edgelist, weight=1, underlying=False, color="black",
                 hyperid=None):
        """ Initialize class CausalEdge. """

        self.edgelist = edgelist
        self.weight = weight
        self.underlying = underlying
        self.color = color
        self.hyperid = hyperid
        self.check_types()
        self.update()


    def update(self):
        """ Check that all edges within the hyperedge have the same target. """

        self.target = self.edgelist[0].target
        for edge in self.edgelist:
            if edge.target != self.target:
                raise ValueError("Hyperedge has more than one target.")
        self.sources = []
        for edge in self.edgelist:
            self.sources.append(edge.source)
        #self.weight = self.edgelist[0].weight
        #for edge in self.edgelist:
        #    if edge.weight != self.weight:
        #        raise ValueError("Each edge of an hyperedge should have the"
        #                         "same weight.")


    def addedge(self, edge):
        """ Add one more edge to hyperedge. """

        self.edgelist.append(edge)
        self.update()
        

    def check_types(self):
        """ Check that Hyperedge attributes have proper types. """

        if not isinstance(self.edgelist, list):
            raise TypeError("edgelist should be an list.")


    def __repr__(self):
        """ Representation of the HyperEdge object. """

        res = "HyperEdge:\n"
        for edge in self.edgelist:
            res+="{}\n".format(edge.__repr__())

        return res


#class MidNode(object):
#    """
#    Intermediary node to represent edge groups as Mesh objects.
#    MidNodes can be enablings (black by default) or requirements
#    (white by default).
#    """
#
#    def __init__(self, nodeid, rank=None, midtype="enabling", ghost=False,
#                 logic="and", fillcolor="black", bordercolor="black",
#                 pos=None, overridewidth=None):
#        """ Initialize class MidNode. """
#
#        self.nodeid = nodeid
#        self.label = ""
#        self.rank = rank
#        self.midtype = midtype # enabling or involvement
#        self.ghost = ghost
#        self.logic = logic
#        self.fillcolor = fillcolor
#        self.bordercolor = bordercolor
#        self.pos = pos
#        self.overridewidth = overridewidth
#        if self.midtype == "involvement":
#            self.fillcolor = "white"
#        self.check_types()
#
#
#    def check_types(self):
#        """ Check that MidNode attributes have proper types. """
#
#        if not isinstance(self.nodeid, str):
#            raise TypeError("nodeid should be a string.")
#        if self.rank != None:
#            if not isinstance(self.rank, float):
#                raise TypeError("rank should be an float.")
#
#
#    def __repr__(self):
#        """ Representation of the MidNode object. """
#
#        res =  "Node "
#        res += 'id: "{}", '.format(self.nodeid)
#        res += 'type: "{}", '.format(self.midtype)
#        if self.rank != None:
#            res += ",  rank: {}".format(self.rank)
#
#        return res


#class MidEdge(CausalEdge):
#    """
#    Intermediary edge to construct meshes. Source and target can
#    be either an EventNode or a MidNode.
#    """
#
#    def check_types(self):
#        """ Check that MidEdge attributes have proper types. """
#
#        if not isinstance(self.source, EventNode):
#            if not isinstance(self.source, MidNode):
#                raise TypeError("source should be a EventNode or MidNode.")
#        if not isinstance(self.target, EventNode):
#            if not isinstance(self.target, MidNode):
#                raise TypeError("target should be a EventNode or MidNode.")
#        if self.occurrence != None:
#            if not isinstance(self.occurrence, int):
#                raise TypeError("occurrence should be an integer.")
#
#
#class Mesh(object):
#    """
#    A mesh is made of a group of edges glued together by intermediary nodes.
#    """
#
#    def __init__(self, uses=1, usage=1.0, occurrence=1, rel_occ=1.0,
#                 underlying=False, color="black", meshid=None):
#        """ Initialize class Mesh. """
#
#        self.uses = uses
#        self.usage = usage # = uses / occurrence_of_EOI (number of cores)
#        self.occurrence = occurrence
#        self.rel_occ = rel_occ # = occurrence / occurrence_of_EOI
#        self.weight = self.occurrence
#        self.underlying = underlying
#        self.color = color
#        self.meshid = meshid
#        self.midnodes = []
#        self.midedges = []
#
#
#    def get_events(self):
#        """
#        Return de sources and target nodes that are events from mesh object.
#        """
#
#        sources = []
#        targets = []
#        for midedge in self.midedges:
#            if isinstance(midedge.source, EventNode):
#                if midedge.source not in sources:
#                    sources.append(midedge.source)
#            if isinstance(midedge.target, EventNode):
#                if midedge.target not in targets:
#                    targets.append(midedge.target)
#
#        return sources, targets
#
#
#    def extend_midnodes(self):
#        """ Get event nodes connected to each midnode. """
#
#        neighbors = []
#        for midnode in self.midnodes:
#            event_srcs = []
#            event_trgs = []
#            for midedge in self.midedges:
#                if midedge.target == midnode:
#                    if isinstance(midedge.source, EventNode):
#                        event_srcs.append(midedge.source)
#                    elif isinstance(midedge.source, MidNode):
#                        for midedge2 in self.midedges:
#                            if midedge2.target == midedge.source:
#                                event_srcs.append(midedge2.source)
#                if midedge.source == midnode:
#                    if isinstance(midedge.target, EventNode):
#                        event_trgs.append(midedge.target)
#                    elif isinstance(midedge.target, MidNode):
#                        for midedge2 in self.midedges:
#                            if midedge2.source == midedge.target:
#                                event_trgs.append(midedge2.target)
#            neighbors.append({"srcs": event_srcs, "trgs": event_trgs})
#
#        return neighbors
#
#
#    def extend_midedges(self):
#        """
#        Get event nodes connected to each midedge. Consider only midedges whose
#        source is an involvement or whose target is an enabling (This means
#        that we ignore edges that are from an event node to an involvement or
#        from an enabling to an event node).
#        """
#
#        neighbors = []
#        for midedge in self.midedges:
#            # Filter edges.
#            src_is_inv = False
#            if isinstance(midedge.source, MidNode):
#                if midedge.source.midtype == "involvement":
#                    src_is_inv = True
#            trg_is_ena = False
#            if isinstance(midedge.target, MidNode):
#                if midedge.target.midtype == "enabling":
#                    trg_is_ena = True
#            # Get the neighbors of filtered edges.
#            if src_is_inv == True or trg_is_ena == True:
#                event_srcs = []
#                event_trgs = []
#                if isinstance(midedge.source, EventNode):
#                    event_srcs.append(midedge.source)
#                elif isinstance(midedge.source, MidNode):
#                    for midedge2 in self.midedges:
#                        if midedge2.target == midedge.source:
#                            event_srcs.append(midedge2.source)
#                if isinstance(midedge.target, EventNode):
#                    event_trgs.append(midedge.target)
#                elif isinstance(midedge.target, MidNode):
#                    for midedge2 in self.midedges:
#                        if midedge2.source == midedge.target:
#                            event_trgs.append(midedge2.target)
#                neighbors.append({"reltype": midedge.relationtype,
#                                  "srcs": event_srcs, "trgs": event_trgs})
#
#        return neighbors
#
#
#    def get_sources(self, targetnode):
#        """
#        Get the source nodes that enable a specific target node within a mesh.
#        """
#
#        sources = []
#        enablings = []
#        involvements = []
#        for midedge in self.midedges:
#            if midedge.target == targetnode:
#                if isinstance(midedge.source, EventNode):
#                    sources.append(midedge.source)
#                elif isinstance(midedge.source, MidNode):
#                    enablings.append(midedge.source)
#        for midedge in self.midedges:
#            if midedge.target in enablings:
#                if isinstance(midedge.source, EventNode):
#                    sources.append(midedge.source)
#                elif isinstance(midedge.source, MidNode):
#                    involvements.append(midedge.source)
#        for midedge in self.midedges:
#            if midedge.target in involvements:
#                sources.append(midedge.source)
#
#        return sources
#
#
#    def get_targets(self, sourcenode):
#        """
#        Get the targets nodes that a specific source node is involved in within a mesh.
#        """
#
#        targets = []
#        involvements = []
#        enablings = []
#        for midedge in self.midedges:
#            if midedge.source == sourcenode:
#                if isinstance(midedge.target, EventNode):
#                    targets.append(midedge.target)
#                elif isinstance(midedge.target, MidNode):
#                    involvements.append(midedge.target)
#        for midedge in self.midedges:
#            if midedge.source in involvements:
#                if isinstance(midedge.target, EventNode):
#                    targets.append(midedge.target)
#                elif isinstance(midedge.target, MidNode):
#                    enablings.append(midedge.target)
#        for midedge in self.midedges:
#            if midedge.source in enablings:
#                targets.append(midedge.target)
#
#        return targets
#
#
#    def reverse_midedges(self):
#        """
#        Reverse the direction of intermediary edges if the sources they
#        reach within the mesh have an higher average rank than the targets
#        they reach.
#        """
#
#        if len(self.midedges) > 1:
#            for midedge in self.midedges:
#                sources = self.get_sources(midedge.source)
#                targets = self.get_targets(midedge.target)
#                if isinstance(midedge.source, EventNode):
#                    sources.insert(0, midedge.source)
#                if isinstance(midedge.target, EventNode):
#                    targets.insert(0, midedge.target)
#                src_ranks = []
#                trg_ranks = []
#                for source in sources:
#                    src_ranks.append(source.rank)
#                for target in targets:
#                    trg_ranks.append(target.rank)
#                src_ave = 0
#                if len(src_ranks) > 0:
#                    src_ave = statistics.mean(src_ranks)
#                trg_ave = statistics.mean(trg_ranks)
#                if len(src_ranks) > 0 and src_ave >= trg_ave:
#                    midedge.reverse = True
#                else:
#                    midedge.reverse = False
#
#
#    def check_indicators(self):
#        """ Decide if some midedges of the mesh need direction indicators. """
#
#        # Reset indicator data.
#        for midedge in self.midedges:
#            midedge.indicator = False
#        # Assign indicators.
#        neighbors = self.extend_midnodes()
#        for i in range(len(self.midnodes)):
#            source_ranks = []
#            target_ranks = []
#            for source in neighbors[i]["srcs"]:
#                source_ranks.append(source.rank)
#            for target in neighbors[i]["trgs"]:
#                target_ranks.append(target.rank)
#            if len(source_ranks) == 0:
#                source_ranks = [0]
#            if max(source_ranks) > min(target_ranks):
#                for midedge in self.midedges:
#                    if midedge.source == self.midnodes[i]:
#                        midedge.indicator = True
#
#
#    def assign_label_carrier(self):
#        """
#        Choose which midedge will carry the label in meshes that contain more
#        than one midedge. Take the first midedge that has an enabling as source,
#        or the first involvement if there is no enabling.
#        """
#
#        if len(self.midedges) > 1:
#            for midedge in self.midedges:
#                midedge.labelcarrier = False
#            contains_enablings = False
#            for midnode in self.midnodes:
#                if midnode.midtype == "enabling":
#                    contains_enablings = True
#                    break
#            if contains_enablings == True:
#                for midedge in self.midedges:
#                    if isinstance(midedge.source, MidNode):
#                        if midedge.source.midtype == "enabling":
#                            midedge.labelcarrier = True
#                            break
#            elif contains_enablings == False:
#                for midedge in self.midedges:
#                    if isinstance(midedge.source, MidNode):
#                        midedge.labelcarrier = True
#                        break
#        elif len(self.midedges) == 1:
#            self.midedges[0].labelcarrier = True
#
#
#    def __repr__(self):
#        """ Representation of the Mesh object. """
#
#        res = "Mesh"
#        if self.uses != None:
#            res += "  uses = {}".format(self.occurrence)
#        if self.prob != None:
#            res += "  usage = {:.3f}".format(self.prob)
#        res += "\n\n"
#        res +=  "MidNodes:\n\n"
#        for midnode in self.midnodes:
#            res+="{}\n".format(midnode.__repr__())
#        if len(self.midnodes) == 0:
#            res += "None\n"
#        res += "\n"
#        res += "MidEdges:\n\n"
#        for midedge in self.midedges:
#            res+="{}\n".format(midedge.__repr__())
#
#        return res


class CausalGraph(object):
    """ Data structure for causal graphs. """

    def __init__(self, filename=None, eoi=None, meshedgraph=False,
                 nodestype="event", showintro=True, precedenceonly=False,
                 rankposdict=None):
        """ Initialize class CausalGraph. """

        # Header variables.
        self.filename = filename
        self.eoi = eoi
        self.meshedgraph = meshedgraph
        self.nodestype = nodestype # event or species
        self.showintro = showintro
        self.precedenceonly = precedenceonly
        self.rankposdict = rankposdict
        # Main variables.
        self.eventnodes = []
        self.statenodes = []
        self.causaledges = []
        self.hyperedges = []
        #self.midnodes = []
        #self.midedges = []
        #self.meshes = []
        #self.edgegroups = []
        #self.midnodegroups = []
        # Cover edges and midnodes computed from nointro.
        self.coveredges = []
        self.covermidnodes = []
        self.covermidedges = []
        self.covermeshes = []
        self.covermidnodegroups = []
        # Post-computed variables.
        self.occurrence = 1
        self.maxrank = None
        self.prevcores = None
        if self.filename != None:
            self.read_dot(self.filename)


    def read_dot(self, dotpath):
        """
        Read nodes and edges from input causal graph.
        """

        rank = None
        self.label_mapping = {}
        dotfile = open(dotpath, "r").readlines()
        for line in dotfile:
            if 'precedenceonly="True"' in line:
                self.precedenceonly = True
            if 'precedenceonly="False"' in line:
                self.precedenceonly = False
            if 'meshedgraph="True"' in line:
                self.meshedgraph = True
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
                occu_str = line[occu+12:quote]
                if "/" in occu_str:
                    slash = occu_str.index("/")
                    occu_str = occu_str[:slash-1]
                self.occurrence = int(occu_str)
            if "maxrank=" in line:
                maxrank_index = line.index("maxrank")
                quote = line.rfind('"')
                self.maxrank = int(line[maxrank_index+9:quote])
            if "rank = same" in line:
                open_quote = line.index('"')
                close_quote = line[open_quote+1:].index('"')+open_quote+1
                medrank = float(line[open_quote+1:close_quote])
                rank = int(medrank)
            if line[0] == "}":
                rank = None
            # Read nodes.
            if "label=" in line and "Occurrence" not in line:
                if "->" not in line and "rank = same" not in line:
                    if line[0:2] == "//":
                       read_line = line[2:]
                    else:
                       read_line = line
                    tokens = read_line.split()
                    ori_id = tokens[0]
                    if '"' in ori_id:
                        ori_id = ori_id[1:-1]
                    #if "node" not in ori_id:
                    #    node_id = "node{}".format(ori_id)
                    #else:
                    node_id = ori_id
                    label_start = read_line.index("label=")+7
                    label_end = read_line[label_start:].index('"')+label_start
                    label_str = read_line[label_start:label_end].strip()
                    label = label_str.replace("\\n ", "")
                    if "intro=True" in read_line:
                        is_intro = True
                    else:
                        is_intro = False
                    if "first=True" in read_line:
                        is_first = True
                    else:
                        is_first = False
                    if "midtype" in read_line:
                        mid_start = read_line.index('midtype')+8
                        mid_end = read_line[mid_start:].index(',')+mid_start
                        midtype = read_line[mid_start:mid_end]
                        if "style=dotted" in read_line:
                            ghost = True
                        else:
                            ghost = False
                        fillcolor = get_field("fillcolor=", read_line, "black")
                        bordercolor = get_field(" color=", read_line, "black")
                        new_midnode = MidNode(ori_id, rank, midtype,
                                              ghost=ghost, fillcolor=fillcolor,
                                              bordercolor=bordercolor)
                        if 'cover="True"' not in line:
                            self.midnodes.append(new_midnode)
                        elif 'cover="True"' in line:
                            self.covermidnodes.append(new_midnode)
                    else:
                        self.eventnodes.append(EventNode(node_id, label,
                                                         rank,
                                                         intro=is_intro,
                                                         first=is_first))
                        self.label_mapping[node_id] = label
        # Read edges.
        tmp_edges = []
        tmp_midedges = []
        tmp_cedges = []
        tmp_cmidedges = []
        for line in dotfile:
            if "->" in line and '[style="invis"]' not in line:
                if line[0:2] == "//":
                    read_line = line[2:]
                    underlying = True
                else:
                    read_line = line
                    underlying = False
                tokens = read_line.split()
                source_id = tokens[0]
                if '"' in source_id:
                    source_id = source_id[1:-1]
                #if "node" not in source_id and "mid" not in source_id:
                #    source_id = "node{}".format(source_id)
                target_id = tokens[2]
                if '"' in target_id:
                    target_id = target_id[1:-1]
                #if "node" not in target_id and "mid" not in target_id:
                #    target_id = "node{}".format(target_id)
                source = None
                target = None
                for node in self.eventnodes:
                    if node.nodeid == source_id:
                        source = node
                    if node.nodeid == target_id:
                        target = node
                #for node in self.midnodes:
                #    if node.nodeid == source_id:
                #        source = node
                #    if node.nodeid == target_id:
                #        target = node
                #for node in self.covermidnodes:
                #    if node.nodeid == source_id:
                #        source = node
                #    if node.nodeid == target_id:
                #        target = node
                meshid = get_field("meshid=", read_line, 1)
                meshid = int(meshid)
                weight = get_field("weight=", read_line, 1)
                weight = int(weight)
                color = get_field("color=", read_line, "black")
                if "label=" in line:
                    labelcarrier = True
                else:
                    labelcarrier = False
                if self.precedenceonly == False:
                    if self.meshedgraph == False:
                        if "style=dotted" in line:
                            edgetype = "conflict"
                            source_save = source
                            source = target
                            target = source_save
                        else:
                            edgetype = "causal"
                    elif self.meshedgraph == True:
                        if "style=dashed" in line:
                            edgetype = "conflict"
                        else:
                            edgetype = "causal"
                else:
                    edgetype = "precedence"
                if "rev=True" in line:
                   rev = True
                   source_save = source
                   source = target
                   target = source_save
                else:
                   rev = False
                #source_is_mid = isinstance(source, MidNode)
                #target_is_mid = isinstance(target, MidNode)
                #if source_is_mid or target_is_mid:
                #    new_edge = MidEdge(source, target, uses=uses,
                #                       relationtype=edgetype, reverse=rev,
                #                       meshid=meshid, underlying=underlying,
                #                       color=color, labelcarrier=labelcarrier)
                #    if 'cover="True"' not in line:
                #        tmp_midedges.append(new_edge)
                #    elif 'cover="True"' in line:
                #        tmp_cmidedges.append(new_edge)
                #else:
                new_edge = CausalEdge(source, target, weight=weight,
                                      relationtype=edgetype, meshid=meshid,
                                      underlying=underlying, color=color)
                if 'cover="True"' not in line:
                    tmp_edges.append(new_edge)
                elif 'cover="True"' in line:
                    tmp_cedges.append(new_edge)
        for edge in tmp_edges:
            self.causaledges.insert(0, edge)
        for midedge in tmp_midedges:
            self.midedges.insert(0, midedge)
        for cedge in tmp_cedges:
            self.coveredges.insert(0, cedge)
        for cmidedge in tmp_cmidedges:
            self.covermidedges.insert(0, cmidedge)
        self.postprocess()


    def postprocess(self):
        """
        Various stuff to do after reading a dot file. Includes the creation
        of meshes from the intermediary nodes and edges.
        """

        if self.meshedgraph == False:
            self.create_hyperedges()
            #self.find_edgegroups()
            #self.create_meshes(self.edgegroups, 1)
            for node in self.eventnodes:
                if "Intro" in node.label:
                    node.intro = True
                    node.label = node.label[6:]
                    if node.label == "Lig, Lig":
                        node.label = "Lig"
            self.find_first_rules()
            self.rank_sequentially()
            self.rm_superfluous_causal_edges()
            self.align_vertical()
        elif self.meshedgraph == True:
           self.find_midnodegroups(cover=False)
           self.find_midnodegroups(cover=True)
           self.read_meshes(cover=False)
           self.read_meshes(cover=True)
           sorted_meshes = sorted(self.meshes, key=lambda x: x.meshid)
           sorted_covermeshes = sorted(self.covermeshes, key=lambda x: x.meshid)
           self.meshes = sorted_meshes
           self.covermeshes = sorted_covermeshes
        if self.eoi == None:
            self.get_maxrank()
            for node in self.nodes:
                if node.rank == self.maxrank:
                    self.eoi = node.label


    def create_hyperedges(self):
        """ Create hyperedges by grouping edges with the same target. """

        for edge in self.causaledges:
            edge_found = False
            for hyperedge in self.hyperedges:
                if edge.target == hyperedge.target:
                    hyperedge.addedge(edge)
                    edge_found = True
            if edge_found == False:
                self.hyperedges.append(HyperEdge([edge]))


#    def find_edgegroups(self):
#        """
#        Group edges together when there is a path between them using only
#        head-to-head and tail-to-tail connections.
#
#        Example: The 3 edges in the following graph form a single group.
#
#                 A  B
#                 | /|
#                 |/ |
#                 C  D
#        """
#
#        edgescopy = self.causaledges.copy()
#        midedgescopy = self.midedges.copy()
#        all_edges = edgescopy + midedgescopy
#        while len(all_edges) > 0:
#            current_group = [all_edges[0]]
#            del(all_edges[0])
#            new_edge_found = True
#            while new_edge_found == True:
#                # Find sources and targets:
#                sources = []
#                targets = []
#                for current_edge in current_group:
#                    if current_edge.source not in sources:
#                        sources.append(current_edge.source)
#                    if current_edge.target not in targets:
#                        targets.append(current_edge.target)
#                # Find other edges with same source or target.
#                new_edge_found = False
#                copy_to_remove = []
#                for i in range(len(all_edges)):
#                    other_source = all_edges[i].source
#                    other_target = all_edges[i].target
#                    if other_source in sources or other_target in targets:
#                        new_edge_found = True
#                        current_group.append(all_edges[i])
#                        copy_to_remove.insert(0, i)
#                for i in copy_to_remove:
#                    del(all_edges[i])
#                if new_edge_found == False:
#                    self.edgegroups.append(current_group)


#    def create_meshes(self, edgegroups, startid):
#        """ Create meshes from edge groups. """
#
#        midid = startid
#        for edgegroup in edgegroups:
#            new_mesh = Mesh()
#            # Collect all sources and targets.
#            sources = []
#            targets = []
#            for edge in edgegroup:
#                if edge.source not in sources:
#                    sources.append(edge.source)
#                if edge.target not in targets:
#                    targets.append(edge.target)
#            # Create intermediary nodes for event nodes with more than one
#            # input or output. Also create an edge between those intermediary
#            # nodes and the corresponding event nodes. 
#            inv_sources = []
#            for source in sources:
#                outgoing = []
#                for edge in edgegroup:
#                    if edge.source == source:
#                        outgoing.append(edge)
#                if len(outgoing) > 1:
#                    inv_sources.append(source)
#                    new_midnode = MidNode("mid{}".format(midid),
#                                          midtype="involvement")
#                    new_midedge = MidEdge(source, new_midnode)
#                    new_mesh.midnodes.append(new_midnode)
#                    new_mesh.midedges.append(new_midedge)
#                    midid += 1
#            ena_targets = []
#            for target in targets:
#                ingoing = []
#                for edge in edgegroup:
#                    if edge.target == target:
#                        ingoing.append(edge)
#                if len(ingoing) > 1:
#                    ena_targets.append(target)
#                    new_midnode = MidNode("mid{}".format(midid),
#                                          midtype="enabling")
#                    new_midedge = MidEdge(new_midnode, target)
#                    new_mesh.midnodes.append(new_midnode)
#                    new_mesh.midedges.append(new_midedge)
#                    midid += 1
#            # Add the intermediary edges corresponding to the original edges.
#            for ori_edge in edgegroup:
#                if ori_edge.source in inv_sources:
#                    for midedge in new_mesh.midedges:
#                        if midedge.source == ori_edge.source:
#                            s = midedge.target
#                else:
#                    s = ori_edge.source
#                if ori_edge.target in ena_targets:
#                    for midedge in new_mesh.midedges:
#                        if midedge.target == ori_edge.target:
#                            t = midedge.source
#                else:
#                    t = ori_edge.target
#                reltype = ori_edge.relationtype
#                new_mesh.midedges.append(MidEdge(s, t, relationtype=reltype))
#            new_mesh.uses = edgegroup[0].uses
#            new_mesh.weight = edgegroup[0].uses
#            self.meshes.append(new_mesh)
#        self.meshedgraph = True
#
#
#    def find_midnodegroups(self, cover=False):
#        """ Find groups of midnodes connected together by midedges. """
#
#        if cover == False:
#            worknodes = self.midnodes
#            workedges = self.midedges
#            workgroups = self.midnodegroups
#        elif cover == True:
#            worknodes = self.covermidnodes
#            workedges = self.covermidedges
#            workgroups = self.covermidnodegroups
#        midnodescopy = worknodes.copy()
#        while len(midnodescopy) > 0:
#            current_group = [midnodescopy[0]]
#            del(midnodescopy[0])
#            new_midnode_found = True
#            while new_midnode_found == True:
#                new_midnode_found = False
#                # Find midnodes connected to current midnodes through midedges.
#                midnodes_to_add = []
#                for current_midnode in current_group:
#                    for midedge in workedges:
#                        if midedge.source == current_midnode:
#                            if isinstance(midedge.target, MidNode):
#                                if midedge.target not in current_group:
#                                    if midedge.target not in midnodes_to_add:
#                                        midnodes_to_add.append(midedge.target)
#                                        new_midnode_found = True
#                        if midedge.target == current_midnode:
#                            if isinstance(midedge.source, MidNode):
#                                if midedge.source not in current_group:
#                                    if midedge.source not in midnodes_to_add:
#                                        midnodes_to_add.append(midedge.source)
#                                        new_midnode_found = True
#                for midnode_to_add in midnodes_to_add:
#                    current_group.append(midnode_to_add)
#                # Remove midnodes from midnodescopy.
#                copy_to_remove = []
#                for i in range(len(midnodescopy)):
#                    if midnodescopy[i] in midnodes_to_add:
#                        copy_to_remove.insert(0, i)
#                for i in copy_to_remove:
#                    del(midnodescopy[i])
#                if new_midnode_found == False:
#                    workgroups.append(current_group)


#    def read_meshes(self, cover=False):
#        """ Rebuild meshes based on midnodes and their edges. """
#
#        if cover == False:
#            workgroups = self.midnodegroups
#            workedges = self.midedges
#            workmeshes = self.meshes
#            workcausal = self.causaledges
#        elif cover == True:
#            workgroups = self.covermidnodegroups
#            workedges = self.covermidedges
#            workmeshes = self.covermeshes
#            workcausal = self.coveredges
#        for midnodegroup in workgroups:
#            new_mesh = Mesh()
#            for midnode in midnodegroup:
#                new_mesh.midnodes.append(midnode)
#                for midedge in workedges:
#                    if midedge.source == midnode:
#                        if midedge not in new_mesh.midedges:
#                            new_mesh.midedges.append(midedge)
#                    if midedge.target == midnode:
#                        if midedge not in new_mesh.midedges:
#                            new_mesh.midedges.append(midedge)
#            new_mesh.uses = new_mesh.midedges[0].uses
#            new_mesh.weight = new_mesh.midedges[0].uses
#            new_mesh.meshid = new_mesh.midedges[0].meshid
#            new_mesh.underlying = new_mesh.midedges[0].underlying
#            workmeshes.append(new_mesh)
#        for causaledge in workcausal:
#            new_mesh = Mesh(uses=causaledge.uses)
#            new_mesh.midedges.append(causaledge)
#            new_mesh.meshid = new_mesh.midedges[0].meshid
#            new_mesh.underlying = new_mesh.midedges[0].underlying
#            workmeshes.append(new_mesh)


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


    def rank_sequentially(self, intropos="bot", rulepos="top"):
        """
        Find the rank of each node, starting with first nodes and then adding
        the other nodes sequentially as soon as they have a secured enabling.
        """

        # Initialize ranks.
        current_nodes = []
        for node in self.eventnodes:
            if node.first == True:
                node.rank = 1
                current_nodes.append(node)
            else:
                node.rank = None
        while len(current_nodes) > 0:
            # 1) Gather hyperedges that have a current_node in their sources.
            current_hyperedges = []
            for hyperedge in self.hyperedges:
                for current_node in current_nodes:
                    if current_node in hyperedge.sources:
                        if hyperedge not in current_hyperedges:
                            current_hyperedges.append(hyperedge)
            # 2) Gather candidate nodes as any target of current meshes
            #    that is not ranked yet.
            candidate_nodes = []
            for hyperedge in current_hyperedges:
                if hyperedge.target.rank == None:
                    if hyperedge.target not in candidate_nodes:
                        candidate_nodes.append(hyperedge.target)
            # 3) Set rank of all candidate nodes that are secured: all the
            #    nodes pointing to them (ignoring intro nodes) are already
            #    ranked in at least one edge group.
            for candidate_node in candidate_nodes:
                incoming_hedges = []
                for hyperedge in current_hyperedges:
                    if hyperedge.target == candidate_node:
                        incoming_hedges.append(hyperedge)
                nsecured = 0
                possible_ranks = []
                for incoming_hedge in incoming_hedges:
                    secured = True
                    for source in incoming_hedge.sources:
                        if source.intro == False:
                            if source.rank == None:
                                secured = False
                                break
                    if secured == True:
                        nsecured += 1
                        source_ranks = []
                        for source in incoming_hedge.sources:
                            if source.intro == False:
                                source_ranks.append(source.rank)
                        possible_ranks.append(max(source_ranks)+1)
                if rulepos == "top" and nsecured > 0:
                    candidate_node.rank = min(possible_ranks)
                    current_nodes.append(candidate_node)
                if rulepos == "bot" and nsecured == len(incoming_hedge):
                    candidate_node.rank = max(possible_ranks)
                    current_nodes.append(candidate_node)                    
            # 4) Remove all current_nodes for which all outgoing hyperedges
            #    have their target already ranked.
            next_nodes = []
            for current_node in current_nodes:
                keep_node = False
                target_nodes = []
                for hyperedge in self.hyperedges:
                    if current_node in hyperedge.sources:
                        if hyperedge.target not in target_nodes:
                            target_nodes.append(hyperedge.target)
                for target_node in target_nodes:
                    if target_node.rank == None:
                        keep_node = True
                        break
                if keep_node == True:
                    next_nodes.append(current_node)
            current_nodes = next_nodes
        # Rank intro nodes at 0.
        for node in self.eventnodes:
            if node.intro == True:
                node.rank = 0
        # Optionally, push intro nodes down when possible.
        if intropos == "bot":
            gap_found = True
            while gap_found == True:
                gap_found = False
                for node in self.eventnodes:
                    # Find the rank of all targets of that node
                    # (excluding loop targets).
                    target_ranks = []
                    for hyperedge in self.hyperedges:
                        if node in hyperedge.sources:
                            if hyperedge.target.rank > node.rank:
                                target_ranks.append(hyperedge.target.rank)
                    if len(target_ranks) > 0:
                        new_rank = min(target_ranks) - 1
                        if new_rank > node.rank:
                            node.rank = new_rank
                            gap_found = True
        self.get_maxrank()
        #self.sequentialize_nodeids()


#    def rank_intermediary(self, edgegroup):
#        """
#        Rank intermediary nodes. Not used (nor well defined) for the moment.
#        Intermediary node ranking is left to graphviz for now.
#        """
#
#        for edge in edgegroup:
#            source_ranks = []
#            for node in edge.source.nodelist:
#                source_ranks.append(node.rank)
#            if edge.target.rank > max(source_ranks):
#                edge.mednode.rank = (( edge.target.rank +
#                                       max(source_ranks) ) / 2.0)
#            elif edge.target.rank < min(source_ranks):
#                edge.mednode.rank = (( edge.target.rank +
#                                       min(source_ranks) ) / 2.0)
#            else:
#                ranks_set = list(set(source_ranks))
#                edge.mednode.rank = statistics.mean(ranks_set)


    def rm_superfluous_causal_edges(self):
        """
        Remove causal edges when the source was already used earlier in story.
        """

        new_hyperedges = []
        for hyperedge in self.hyperedges:
            if len(hyperedge.edgelist) == 1:
                new_hyperedges.append(hyperedge)
            else:
                new_edgelist = []
                for edge in hyperedge.edgelist:
                    paths = self.follow_edges("down",
                                              edge.source, [edge.target])
                    #print(edge)
                    #print("-------")
                    #print(paths)
                    #print("=======")
                    if len(paths) == 1:
                        new_edgelist.append(edge)
                if len(new_edgelist) > 0:
                    new_hyperedges.append(HyperEdge(new_edgelist))

        self.hyperedges = new_hyperedges


    def align_vertical(self):
        """
        Adjust edge weights such that the event nodes are aligned vertically.
        """

        for hyperedge in self.hyperedges:
            if len(hyperedge.edgelist) > 1:
                nonintro_present = False
                for edge in hyperedge.edgelist:
                    if edge.source.intro == False:
                        nonintro_present = True
                        break
                if nonintro_present == True:
                    for edge in hyperedge.edgelist:
                        if edge.source.intro == True:
                            edge.weight = 0


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
                for hyperedge in self.hyperedges:
                    if direction == "up":
                        if hyperedge.target == path[-1]:
                            for src_node in hyperedge.sources:
                                next_nodes.append(src_node)
                    elif direction == "down":
                        if path[-1] in hyperedge.sources:
                            next_nodes.append(hyperedge.target)
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


    def get_maxrank(self):
        """ Find the highest rank of a node in CausalGraph. """

        all_ranks = []
        for node in self.eventnodes:
            if node.rank != None:
                all_ranks.append(node.rank)
        if len(all_ranks) > 0:
            self.maxrank = max(all_ranks)


    def sequentialize_nodeids(self):
        """ Assign sequential node ids, getting rid of event numbers. """

        node_number = 1
        for current_rank in range(self.maxrank+1):
            for node in self.eventnodes:
                if node.rank == current_rank:
                    node.nodeid = "node{}".format(node_number)
                    node_number += 1
        # Also sort causal edges.
        # (Not needed, only sort grouped edges before building the dot file).
        #sorted_edges = sorted(self.causaledges, key=lambda x: x.source.rank)
        #self.causaledges = sorted_edges


    def build_nointro(self):
        """
        Create new meshes for the version of the graph that hides
        intro nodes.
        """

        # Reset information about cover edges. 
        for mesh in self.meshes:
            mesh.underlying = False
        self.covermeshes = []
        nointro_groups = []
        midid = self.find_max_midid()+1
        for mesh1 in self.meshes:
            mesh_list = []
            if mesh1.underlying == False:
                has_intro1, all_intro1 = self.check_intro(mesh1)
                if has_intro1 == True:
                    if all_intro1 == True:
                        mesh1.underlying = True
                    else:
                        mesh_list.append(mesh1)
                        noin1 = self.nointro_mesh(mesh1, midid)
                        midid += len(noin1.midnodes)
                        if len(noin1.midedges) > 0:
                            # Check all other meshes that have the same
                            # midedges when ignoring edges with intro nodes as
                            # source. They will all be grouped inside a single
                            # mesh without intro nodes as sources.
                            for mesh2 in self.meshes:
                                if mesh2.color == mesh1.color:
                                    if mesh2 != mesh1:
                                        if mesh2.underlying == False:
                                            noin2 = self.nointro_mesh(mesh2,
                                                                      midid)
                                            if self.equivalent_meshes(noin1,
                                                                      noin2):
                                                mesh_list.append(mesh2)
            if len(mesh_list) > 0:
                # Compute occurence of nointro mesh as the sum of all its
                # underlying meshes. Also mark meshes that were used as
                # underlying.
                uses_summ = 0
                for mesh in mesh_list:
                    uses_summ += mesh.uses
                    mesh.underlying = True
                noin1.uses = uses_summ
                noin1.weight = uses_summ
                if len(noin1.midedges) > 0:
                    self.covermeshes.append(noin1)


    def find_max_midid(self, cover=False):
        """ Find the highest intermediary node id in the graph. """

        midids = []
        for mesh in self.meshes:
            for midnode in mesh.midnodes:
                midids.append(int(midnode.nodeid[3:]))
        if cover == True:
            for mesh in self.covermeshes:
                for midnode in mesh.midnodes:
                    midids.append(int(midnode.nodeid[3:]))
        max_midid = max(midids)

        return max_midid


    def find_max_meshid(self, cover=False):
        """ Find the highest mesh id in the graph. """

        meshids = []
        for mesh in self.meshes:
            meshids.append(mesh.meshid)
        if cover == True:
            for mesh in self.covermeshes:
                meshids.append(mesh.meshid)
        max_meshid = max(meshids)

        return max_meshid


    def check_intro(self, mesh):
        """ Check if a mesh has intro nodes in its sources. """

        has_intro = False
        all_intro = True
        for midedge in mesh.midedges:
            if isinstance(midedge.source, EventNode):
                if midedge.source.intro == True:
                    has_intro = True
                if midedge.source.intro == False:
                    all_intro = False

        return has_intro, all_intro


    def nointro_mesh(self, mesh, midid):
        """
        Build a new mesh without midedges from intro nodes.
        Adjust midnodes to the removal of those midedges.
        If the removal of the intro node disconnects the mesh,
        add an midedge between the enablings.
        """

        new_mesh = Mesh(uses=mesh.uses, color=mesh.color)
        # Add midnodes with temporary ids.
        tmpid = 1
        tmpid_map = {}
        for midnode in mesh.midnodes:
            new_mesh.midnodes.append(MidNode("mid{}".format(tmpid),
                                             midtype=midnode.midtype))
            new_mesh.midnodes[-1].bordercolor = midnode.bordercolor
            if midnode.midtype == "enabling":
                new_mesh.midnodes[-1].fillcolor = midnode.fillcolor
            tmpid_map["mid{}".format(tmpid)] = midnode.nodeid
            tmpid += 1
        # Add the midedges that do not have an intro node as source.
        for midedge in mesh.midedges:
            add_edge = True
            if isinstance(midedge.source, EventNode):
                if midedge.source.intro == True:
                    add_edge = False
            if add_edge == True:
                # Get source.
                if isinstance(midedge.source, EventNode):
                    s = midedge.source
                elif isinstance(midedge.source, MidNode):
                    for new_node in new_mesh.midnodes:
                        if tmpid_map[new_node.nodeid] == midedge.source.nodeid:
                            s = new_node
                # Get target.
                if isinstance(midedge.target, EventNode):
                    t = midedge.target
                elif isinstance(midedge.target, MidNode):
                    for new_node in new_mesh.midnodes:
                        if tmpid_map[new_node.nodeid] == midedge.target.nodeid:
                            t = new_node
                use = midedge.uses
                rel = midedge.relationtype
                col = midedge.color
                new_mesh.midedges.append(MidEdge(s, t, uses=use,
                                                 relationtype=rel, color=col))
        # Treat involvement nodes that have no incoming edge from an event
        # node as ghost nodes.
        for midnode in new_mesh.midnodes:
            if midnode.midtype == "involvement":
                has_incoming = False
                for midedge in new_mesh.midedges:
                    if midedge.target == midnode:
                        if isinstance(midedge.source, EventNode):
                            has_incoming = True
                            break
                if has_incoming == False:
                    midnode.ghost = True
        # Remove enablings if they have only one incoming edge and one
        # outgoing edge once intro nodes were removed. Replace the incoming
        # and outgoing edges by a single edge.
        enas_to_remove = []
        midedges_to_remove = []
        midedges_to_add = []
        for i in range(len(new_mesh.midnodes)):
            midnode = new_mesh.midnodes[i]
            if midnode.midtype == "enabling":
                # Keep midnode if it is found in connectors.
                #connected = False
                #for connect in new_mesh.connectors:
                #    if connect.source == midnode or connect.source == midnode:
                #        connected = True
                #if connected == False or connected == True:
                incoming_count = 0
                outgoing_count = 0
                for midedge in new_mesh.midedges:
                    if midedge.target == midnode:
                        incoming_count += 1
                        in_edge = midedge
                    if midedge.source == midnode:
                        outgoing_count += 1
                        out_edge = midedge
                if incoming_count == 1 and outgoing_count == 1:
                    enas_to_remove.insert(0, i)
                    for j in range(len(new_mesh.midedges)):
                        if new_mesh.midedges[j].source == midnode:
                            if j not in midedges_to_remove:
                                midedges_to_remove.append(j)
                        if new_mesh.midedges[j].target == midnode:
                            if j not in midedges_to_remove:
                                midedges_to_remove.append(j)
                    use = in_edge.uses
                    rel = in_edge.relationtype
                    midedges_to_add.append(MidEdge(in_edge.source,
                                                   out_edge.target,
                                                   uses=use,
                                                   relationtype=rel,
                                                   color=mesh.color))
        for j in sorted(midedges_to_remove, reverse=True):
            del(new_mesh.midedges[j])
        for i in enas_to_remove:
            del(new_mesh.midnodes[i])
        for midedge in midedges_to_add:
            new_mesh.midedges.append(midedge)
        # Reassign intermediary node ids in case some were removed.
        for midnode in new_mesh.midnodes:
            midnode.nodeid = "mid{}".format(midid)
            midid += 1

        return new_mesh


    def equivalent_meshes(self, mesh1, mesh2):
        """
        Find whether two meshes connect the same event nodes
        with equivalent midedges.
        """

        nn1 = len(mesh1.midnodes)
        nn2 = len(mesh2.midnodes)
        ne1 = len(mesh1.midedges)
        ne2 = len(mesh2.midedges)
        if nn1 == nn2 and ne1 == ne2:
            are_same = True
        else:
            are_same = False
        if are_same == True:
            sources1, targets1 = mesh1.get_events()
            sources2, targets2 = mesh2.get_events()
            same_sources = same_objects(sources1, sources2)
            same_targets = same_objects(targets1, targets2)
            if same_sources == True and same_targets == True:
                are_same = True
            else:
                are_same = False
        if are_same == True:
            neighbors1 = mesh1.extend_midedges()
            neighbors2 = mesh2.extend_midedges()
            equi_midedges = self.equivalent_midedges(neighbors1, neighbors2)
            if equi_midedges == True:
                are_same = True
            else:
                are_same = False

        return are_same


    def equivalent_midedges(self, neighbors1, neighbors2):
        """
        Find whether two lists of midedges (described as their respective
        neighbors) connect to the same event nodes.
        """
 
        list1 = neighbors1.copy()
        list2 = neighbors2.copy()
        found1 = []
        found2 = []
        for i in range(len(list1)):
            s1 = list1[i]["srcs"]
            t1 = list1[i]["trgs"]
            for j in range(len(list2)):
                s2 = list2[j]["srcs"]
                t2 = list2[j]["trgs"]
                if list1[i]["reltype"] == list2[j]["reltype"]:
                    if same_objects(s1, s2):
                        if same_objects(t1, t2):
                            found1.insert(0, i)
                            break
        for j in range(len(list2)):
            s2 = list2[j]["srcs"]
            t2 = list2[j]["trgs"]
            for i in range(len(list1)):
                s1 = list1[i]["srcs"]
                t1 = list1[i]["trgs"]
                if list2[j]["reltype"] == list1[i]["reltype"]:
                    if same_objects(s2, s1):
                        if same_objects(t2, t1):
                            found2.insert(0, j)
                            break
        for i in found1:
            del(list1[i])
        for j in found2:
            del(list2[j])
        if len(list1) == 0 and len(list2) == 0:
            are_same = True
        else:
            are_same = False
 
        return are_same


    def color_meshes(self, showintro=True):
        """
        Color meshes to better distinguish them when they have midnodes
        and share at least 1 event node.
        """

        meshes_to_color = []
        if showintro == True:
            for mesh in self.meshes:
                meshes_to_color.append(mesh)
        else:
            for mesh in self.meshes:
                if mesh.underlying == False:
                    meshes_to_color.append(mesh)
            for covermesh in self.covermeshes:
                meshes_to_color.append(covermesh)
        # For each mesh to color, find all other meshes that share
        # at least 1 nodes and overlap in their ranks.
        nodeshares = []
        for i in range(len(meshes_to_color)):
            nodeshares.append([])
        for i in range(len(meshes_to_color)):
            s1, t1 = meshes_to_color[i].get_events()
            nodelist1 = s1 + t1
            for j in range(len(meshes_to_color)):
                if i != j:
                    s2, t2 = meshes_to_color[j].get_events()
                    nodelist2 = s2 + t2
                    n_same = 0
                    for node1 in nodelist1:
                        for node2 in nodelist2:
                            if node1 == node2:
                                n_same += 1
                    if n_same >= 1:
                        overlap = self.rank_overlap(meshes_to_color[i],
                                                    meshes_to_color[j])
                        if overlap == True:
                            nodeshares[i].append(j)
        # Create a list of dictionaries containing data from each mesh.
        dictlist = []
        for i in range(len(meshes_to_color)):
            use = meshes_to_color[i].uses
            nshare = len(nodeshares[i])
            d = {"index": i, "use": use, "s": nshare}
            dictlist.append(d)
        sortlist = sorted(dictlist, key=lambda x: (x["use"], x["s"]),
                          reverse=True)
        # Assign color ids.
        color_ids = []
        for i in range(len(meshes_to_color)):
            color_ids.append(0)
        for d in sortlist:
            mesh_index = d["index"]
            neighbor_indexes = nodeshares[mesh_index]
            used_colors = []
            for neighbor_index in neighbor_indexes:
                used_colors.append(color_ids[neighbor_index])
            # Find the lowest non-zero color index that is absent from
            # used_colors.
            col = 1
            color_found = False
            while color_found == False:
                if col not in used_colors:
                    mesh_color_id = col
                    color_found = True
                col += 1
            color_ids[mesh_index] = mesh_color_id
        # Assign color values.
        color_palette = ["black", "blue2", "chartreuse4", "firebrick3",
                         "darkviolet", "darkorange1", "deepskyblue1",
                         "springgreen", "brown2", "magenta", "orange"]
        for i in range(len(meshes_to_color)):
            palette_index = (color_ids[i]-1)%len(color_palette)
            meshes_to_color[i].color = color_palette[palette_index]
            for midedge in meshes_to_color[i].midedges:
                midedge.color = meshes_to_color[i].color
            for midnode in meshes_to_color[i].midnodes:
                midnode.bordercolor = meshes_to_color[i].color
                if midnode.midtype == "enabling":
                    midnode.fillcolor = meshes_to_color[i].color


    def rank_overlap(self, mesh1, mesh2):
        """
        Find if the ranks of the event nodes connected to 2 meshes overlap.

        Examples of overlapping meshes:

                      |                  |     m1s
         m1s  m2s     |     m1s          |      |
          |    |      |      |           |      |   m2s
          V    V      |      V   m2s     |      |    |
         m1t  m2t     |     m1t   |      |      |    V
                      |           V      |      |   m2t
                      |          m2t     |      V
                      |                  |     m1t

        Examples of non-overlapping meshes

         m1s          |     m1s
          |           |      |
          V           |      V
         m1t          |     m1t
              m2s     |      |
               |      |      V
               V      |     m2t
              m2t     |
        """

        sources1, targets1 = mesh1.get_events()
        nodes1 = sources1 + targets1
        ranks1 = []
        for node1 in nodes1:
            ranks1.append(node1.rank)
        min1 = min(ranks1)
        max1 = max(ranks1)
        sources2, targets2 = mesh2.get_events()
        nodes2 = sources2 + targets2
        ranks2 = []
        for node2 in nodes2:
            ranks2.append(node2.rank)
        min2 = min(ranks2)
        max2 = max(ranks2)
        if max1 > min2 and max1 <= max2:
            between1 = True
        else:
            between1 = False
        if max2 > min1 and max2 <= max1:
            between2 = True
        else:
            between2 = False
        if between1 == True or between2 == True:
            overlap = True
        else:
            overlap = False

        return overlap


    def assign_meshid(self, showintro):
        """ Assign a unique numeral id to every mesh and cover mesh. """

        meshid = 1
        for mesh in self.meshes:
            mesh.meshid = meshid
            meshid += 1
        for covermesh in self.covermeshes:
            covermesh.meshid = meshid
            meshid += 1


    def compute_relstats(self):
        """ Compute relative statistics for every mesh. """

        for hyperedge in self.hyperedges:
            hyperedge.rel_wei = hyperedge.weight/self.occurrence
            hyperedge.rel_occ = hyperedge.occurrence/self.occurrence
            for midedge in mesh.midedges:
                midedge.usage = mesh.usage
                midedge.rel_occ = midedge.occurrence/self.occurrence
        for covermesh in self.covermeshes:
            covermesh.usage = covermesh.uses/self.occurrence
            covermesh.rel_occ = covermesh.occurrence/self.occurrence
            for midedge in covermesh.midedges:
                midedge.usage = covermesh.usage
                midedge.rel_occ = midedge.occurrence/self.occurrence


    def compute_visuals(self, showintro=True, color=True):
        """
        Compute visuals like color and label carriers.
        This erases previous visuals.
        """

        # Meshes without intro nodes.
        if showintro == False:
            self.build_nointro()
        # Reverse edges that point upwards.
        for mesh in self.meshes:
            mesh.reverse_midedges()
        for covermesh in self.covermeshes:
            covermesh.reverse_midedges()
        # Assign mesh ids.
        self.assign_meshid(showintro)
        # Color meshes.
        if color == True:
            self.color_meshes(showintro)
        # Assign which edges will carry labels.
        for mesh in self.meshes:
            mesh.assign_label_carrier()
        for covermesh in self.covermeshes:
            covermesh.assign_label_carrier()


    def build_dot_file(self, showintro=True, addedgelabels=True,
                       showedgelabels=True, edgeid=True, edgeocc=False,
                       edgeuse=True, statstype="rel", weightedges=False):
        """ build a dot file of the CausalGraph. """

        # Write info about graph.
        dot_str = 'digraph G{\n'
        dot_str += '  precedenceonly="{}" ;\n'.format(self.precedenceonly)
        dot_str += '  meshedgraph="{}" ;\n'.format(self.meshedgraph)
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
        #dot_str += '  ranksep=0.5 ;\n'
        #dot_str += '  nodesep=0.2 ;\n' # Default nodesep is 0.25
        dot_str += '  splines=true ;\n'
        dot_str += '  outputorder=nodesfirst ;\n'
        dot_str += '  node [pin=true] ;\n'
        #dot_str += '  edge [fontsize=18] ;\n'
        # Compute some statistics to assign edge and intermediary node width.
        minpenwidth = 1
        medpenwidth = 3
        maxpenwidth = 6.5
        all_weights = []
        for hyperedge in self.hyperedges:
            if hyperedge.underlying == False:
                all_weights.append(hyperedge.weight)
        #for coverhyper in self.coverhypers:
        #    all_uses.append(covermesh.uses)
        average_weight = statistics.mean(all_weights)
        # Draw nodes.
        midranks = 1
        for int_rank in range((self.maxrank+1)*(midranks+1)):
            current_rank = int_rank/(midranks+1)
            rank_str = "{}".format(current_rank)
            if showintro == False and current_rank < 1:
                dot_str += "//"
            if current_rank%1 == 0:
                rank_str = str(int(current_rank))
                dot_str += ('{{ rank = same ; "{}" ['
                            'shape=plaintext'.format(rank_str))
                if self.rankposdict != None:
                    if rank_str in self.rankposdict.keys():
                        rankpos = self.rankposdict[rank_str]
                        dot_str += ', pos={}'.format(rankpos)
                dot_str +='];\n'
            else:
                rank_str = "{:.2f}".format(current_rank)
                dot_str += ('{{ rank = same ; "{}" [label="", '
                            'shape=plaintext'.format(rank_str))
                if self.rankposdict != None:
                    if rank_str in self.rankposdict.keys():
                        rankpos = self.rankposdict[rank_str]
                        dot_str += ', pos={}'.format(rankpos)
                dot_str += '];\n'
            for node in self.eventnodes:
                if node.rank == current_rank:
                    #node_shape = 'invhouse'
                    node_shape = 'rectangle'
                    node_color = 'lightblue'
                    if node.intro == True:
                        node_shape = 'ellipse'
                        node_color = 'white'
                    if node.label == self.eoi:
                        node_shape = 'ellipse'
                        node_color = 'indianred2'
                    if self.nodestype == 'species':
                        node_shape = 'ellipse'
                    if showintro == False and node.intro == True:
                        dot_str += '//'
                    node_lines = textwrap.wrap(node.label, 20,
                                              break_long_words=False)
                    node_str = ""
                    for i in range(len(node_lines)):
                        if i == 0:
                            node_str += " {} ".format(node_lines[i])
                        else:
                            node_str += "\\n {} ".format(node_lines[i])
                    dot_str += ('{} [label="{}"'
                                .format(node.nodeid, node_str))
                    dot_str += ', shape={}, style=filled'.format(node_shape)
                    if node.highlighted == True:
                       dot_str += ', fillcolor=gold, penwidth=2'
                    else:
                       dot_str += ', fillcolor={}'.format(node_color)
                    if node.intro == True:
                        dot_str += ', intro={}'.format(node.intro)
                    if node.first == True:
                        dot_str += ', first={}'.format(node.first)
                    if node.pos != None:
                        dot_str += ', pos={}'.format(node.pos)
                    dot_str += ', penwidth=2'
                    dot_str += "] ;\n"
            for node in self.statenodes:
                if node.rank == current_rank:
                    node_shape = 'ellipse'
                    node_color = 'skyblue2'
                    node_lines = textwrap.wrap(node.label, 20,
                                              break_long_words=False)
                    node_str = ""
                    for i in range(len(node_lines)):
                        if i == 0:
                            node_str += " {} ".format(node_lines[i])
                        else:
                            node_str += "\\n {} ".format(node_lines[i])
                    dot_str += ('{} [label="{}"'
                                .format(node.nodeid, node_str))
                    dot_str += ', shape={}, style=filled'.format(node_shape)
                    if node.highlighted == True:
                       dot_str += ', fillcolor=gold, penwidth=2'
                    else:
                       dot_str += ', fillcolor={}'.format(node_color)
                    if node.intro == True:
                        dot_str += ', intro={}'.format(node.intro)
                    if node.first == True:
                        dot_str += ', first={}'.format(node.first)
                    if node.pos != None:
                        dot_str += ', pos={}'.format(node.pos)
                    dot_str += ', penwidth=2'
                    dot_str += "] ;\n"

            ## Draw intermediary nodes that emulate hyperedges if two
            ## sources or more are drawn.
            #for hyperedge in self.hyperedges:
            #    for midnode in mesh.midnodes:
            #        if midnode.rank == current_rank:
            #            # Include the midnode no matter what, but comment it
            #            # if showintro is False and edge is underlying. 
            #            if showintro == False:
            #                if mesh.underlying == True:
            #                    dot_str += '//'
            #            dot_str += self.write_midnode(mesh, midnode,
            #                average_use, minpenwidth, medpenwidth, maxpenwidth)
            #            dot_str += '] ;\n'
            # Intermediary nodes from cover edges, same as above but only
            # if showintro is False.
            if showintro == False:
                for covermesh in self.covermeshes:
                    for midnode in covermesh.midnodes:
                        if midnode.rank == current_rank:
                            dot_str += self.write_midnode(covermesh, midnode,
                                average_use, minpenwidth, medpenwidth, maxpenwidth)
                            dot_str += ', cover="True"] ;\n'
            # Close rank braces.
            if showintro == False and current_rank < 1:
                dot_str += "//"
            dot_str += "}\n"
        ## Draw unranked midnodes.
        #for mesh in self.meshes:
        #    for midnode in mesh.midnodes:
        #        if midnode.rank == None:
        #            # Include the midnode no matter what, but comment it
        #            # if showintro is False and edge is underlying. 
        #            if showintro == False and mesh.underlying == True:
        #                dot_str += '//'
        #            dot_str += self.write_midnode(mesh, midnode, average_use,
        #                minpenwidth, medpenwidth, maxpenwidth)
        #            dot_str += '] ;\n'
        #if showintro == False:
        #    for covermesh in self.covermeshes:
        #        for midnode in covermesh.midnodes:
        #            if midnode.rank == None:
        #                dot_str += self.write_midnode(covermesh, midnode,
        #                    average_use, minpenwidth, medpenwidth, maxpenwidth)
        #                dot_str += ', cover="True"] ;\n'
        # Draw invisible ranking edges.
        for int_rank in range(self.maxrank*(midranks+1)):
            rank = int_rank/(midranks+1)
            if showintro == False and rank < 1:
                dot_str += '//'
            next_rank = rank+(1.0/(midranks+1))
            if rank%1 == 0:
                rank_str = '{}'.format(int(rank))
            else:
                rank_str = '{:.2f}'.format(rank)
            if next_rank%1 == 0:
                next_str = '{}'.format(int(next_rank))
            else:
                next_str = '{:.2f}'.format(next_rank)
            dot_str += ('"{}" -> "{}" [style="invis"'.format(rank_str,
                                                             next_str))
            if self.rankposdict != None:
                edge_str = "{} -> {}".format(rank_str, next_str)
                if edge_str in self.rankposdict.keys():
                    edgerankpos = self.rankposdict[edge_str]
                    dot_str += ', pos={}'.format(edgerankpos)
            dot_str += '] ;\n'
        # Draw each intermediary edge found in each mesh. Comment if
        # Underlying. The occurrence of each intermediary edge within
        # a mesh should be the same.
        #for mesh in self.meshes:
        #    mesh.check_indicators()
        #    for midedge in mesh.midedges:
        #        if showintro == False and mesh.underlying == True:
        #            dot_str += "//"
        #        dot_str += self.write_midedge(mesh, midedge, average_use,
        #            minpenwidth, medpenwidth, maxpenwidth, addedgelabels,
        #            showedgelabels, edgeid, edgeocc, edgeuse, statstype,
        #            weightedges)
        #        dot_str += '] ;\n'
        for hyperedge in self.hyperedges:
            for edge in hyperedge.edgelist:
                dot_str += self.write_edge(edge)
                dot_str += '] ;\n'
        # Draw cover edges if intro nodes are not shown.
        if showintro == False:
            for covermesh in self.covermeshes:
                covermesh.check_indicators()
                for midedge in covermesh.midedges:
                    dot_str += self.write_midedge(covermesh, midedge,
                        average_use, minpenwidth, medpenwidth, maxpenwidth,
                        addedgelabels, showedgelabels, edgeid, edgeocc,
                        edgeuse, statstype, weightedges)
                    dot_str += ', cover="True"] ;\n'
        # Close graph.
        dot_str += "}"
        self.dot_file = dot_str


    def write_midnode(self, mesh, midnode, average_use, minpenwidth,
                      medpenwidth, maxpenwidth):
        """ Write the line of a dot file for a single midnode."""

        ratio = mesh.uses/average_use
        pensize = math.log(ratio, 2) + medpenwidth
        if pensize < minpenwidth:
            pensize = minpenwidth
        if pensize > maxpenwidth:
            pensize = maxpenwidth
        pensize = math.sqrt(pensize)/12
        mid_str = '"{}" [label=""'.format(midnode.nodeid)
        mid_str += ', shape=circle'
        if midnode.ghost == True:
            mid_str += ', style=dotted'
        else:
            mid_str += ', style=filled'
        mid_str += ', color={}'.format(midnode.bordercolor)
        mid_str += ', fillcolor={}'.format(midnode.fillcolor)
        mid_str += ', midtype={}'.format(midnode.midtype)
        if midnode.overridewidth == None:
            mid_str += ', width={:.4f}'.format(pensize)
            mid_str += ', height={:.4f}'.format(pensize)
        else:
            mid_str += ', width={:.4f}'.format(midnode.overridewidth)
            mid_str += ', height={:.4f}'.format(midnode.overridewidth)
        if midnode.pos != None:
            mid_str += ', pos={}'.format(midnode.pos)

        return mid_str


    def write_edge(self, edge):
        """ Write the line of a dot file for a single edge. """

        edge_str = ('{} -> {} '.format(edge.source.nodeid,
                                           edge.target.nodeid))
        edge_str += "[arrowhead=onormal"
        if edge.relationtype == "conflict":
            edge_str += ", style=dashed"
        edge_str += ', weight={}'.format(edge.weight)
        edge_str += ', penwidth=2'

        return edge_str


    def write_midedge(self, mesh, midedge, average_use, minpenwidth,
                      medpenwidth, maxpenwidth, addedgelabels, showedgelabels,
                      edgeid, edgeocc, edgeuse, statstype, weightedges):
        """ Write the line of a dot file for a single midedge. """

        ratio = mesh.uses/average_use
        pensize = math.log(ratio,2) + medpenwidth
        if pensize < minpenwidth:
            pensize = minpenwidth
        if pensize > maxpenwidth:
            pensize = maxpenwidth
        if midedge.reverse == False:
            mid_str = ('"{}" -> "{}" '.format(midedge.source.nodeid,
                                              midedge.target.nodeid))
        elif midedge.reverse == True:
            mid_str = ('"{}" -> "{}" '.format(midedge.target.nodeid,
                                              midedge.source.nodeid))
        mid_str += '[meshid={}'.format(mesh.meshid)
        if midedge.overridewidth == None:
            mid_str += ', penwidth={}'.format(pensize)
        else:
            mid_str += ', penwidth={}'.format(midedge.overridewidth)
        mid_str += ', color={}'.format(midedge.color)
        if statstype == "abs":
            occ_stat = "{}".format(mesh.occurrence)
            use_stat = "{}".format(mesh.uses)
        elif statstype == "rel":
            occ_stat = "{:.2}".format(mesh.rel_occ)
            use_stat = "{:.2}".format(mesh.usage)
        elif statstype == "both":
            occ_stat = "{}".format(mesh.occurrence)
            occ_stat += " ({:.2})".format(mesh.rel_occ)
            use_stat = "{}".format(mesh.uses)
            use_stat += " ({:.2})".format(mesh.usage)
        if addedgelabels == True:
            if midedge.overridelabel == None:
                if midedge.labelcarrier == True:
                    label_str = ""
                    if edgeid == True:
                        label_str += "  #{}".format(mesh.meshid)
                        if edgeocc == True or edgeuse == True:
                            label_str += "\\n"
                    if edgeocc == True:
                        label_str += "  {}".format(occ_stat)
                        if edgeuse == True:
                           label_str += "\\n"
                    if edgeuse == True:
                        label_str += "  {}".format(use_stat)
                    mid_str += ', label="{}"'.format(label_str)
                    if midedge.labelpos != None:
                        mid_str += ', lp={}'.format(midedge.labelpos)
                if showedgelabels == True:
                    mid_str += ', fontcolor={}'.format(midedge.color)
                elif showedgelabels == False:
                    mid_str += ', fontcolor=transparent'
            else:
                mid_str += ', label="{}"'.format(midedge.overridelabel)
                if midedge.labelpos != None:
                    mid_str += ', lp={}'.format(midedge.labelpos)
        if midedge.indicator == True:
            if isinstance(midedge.source, EventNode):
                mid_str += ", dir=none"
            else:
                mid_str += ", dir=both"
            if midedge.reverse == False:
                if isinstance(midedge.target, MidNode):
                    mid_str += ", arrowhead=none"
                if isinstance(midedge.source, MidNode):
                    #mid_str += ", arrowtail=icurve"
                    mid_str += ", arrowtail=crow"
                    #mid_str += ", arrowtail=inv"
            elif midedge.reverse == True:
                if isinstance(midedge.target, MidNode):
                    mid_str += ", arrowtail=none"
                if isinstance(midedge.source, MidNode):
                    #mid_str += ", arrowhead=icurve"
                    mid_str += ", arrowhead=crow"
                    #mid_str += ", arrowhead=inv"
        elif midedge.indicator == False:
            if isinstance(midedge.target, MidNode):
                mid_str += ", dir=none"
        if midedge.reverse == True:
            mid_str += ", rev=True"
        if midedge.relationtype == "conflict":
            mid_str += ", style=dotted"
        mid_str += ', uses={}'.format(mesh.uses)
        mid_str += ', usage={}'.format(mesh.usage)
        if weightedges == True:
            mid_str += ', weight={}'.format(mesh.weight)
        if midedge.pos != None:
            mid_str += ', pos={}'.format(midedge.pos)

        return mid_str


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


def same_objects(list1, list2):
    """
    Find if two lists of objects contain all the same objects.
    (The objects may be found in a different order in the two lists).
    """

    group1 = list1.copy()
    group2 = list2.copy()
    found1 = []
    found2 = []
    for i in range(len(group1)):
        for object2 in group2:
            if group1[i] == object2:
                found1.insert(0, i)
                break
    for j in range(len(group2)):
        for object1 in group1:
            if group2[j] == object1:
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


def get_field(field, read_str, default):
    """ Extract the value of field in dot file line. """

    if field in read_str:
        field_start = read_str.index(field)+len(field)
        rem = read_str[field_start:]
        if "," in rem:
            field_end = rem.index(",")+field_start
        else:
            field_end = rem.index("]")+field_start
        value = read_str[field_start:field_end]
    else:
        value = default

    return value

# -------------------- Causal Cores Generation Section ------------------------

def getcausalcores(eoi, kappamodel, kasimpath, kaflowpath, simtime=1000,
                   simseed=None, precedenceonly=True):
    """ 
    Generate initial causal cores of given event of interest by running KaSim 
    and then KaFlow.
    """

    new_model = add_eoi(eoi, kappamodel)
    trace_path = run_kasim(new_model, kasimpath, simtime, simseed)
    run_kaflow(eoi, trace_path, kaflowpath, precedenceonly)


def getstories(eoi, kappamodel, kasimpath, kastorpath, simtime=1000,
               simseed=None, compression=None):
    """
    Generate stories using the fill_siphon function of KaStor,
    and optionally weak compression.
    """

    new_model = add_eoi(eoi, kappamodel)
    trace_path = run_kasim(new_model, kasimpath, simtime, simseed)
    run_kastor(eoi, trace_path, kastorpath, compression)


def add_eoi(eoi, kappamodel):
    """ Create a new Kappa model where the EOI is added. """

    if not os.path.exists(eoi):
        os.mkdir(eoi)
    if not os.path.exists("{}/tmp".format(eoi)):
        os.mkdir("{}/tmp".format(eoi))
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


def run_kaflow(eoi, trace_path, kaflowpath, precedenceonly):
    """ Run KaFlow on the trace containing the EOI. """

    if precedenceonly == True:
        subprocess.run(("{}".format(kaflowpath),
                        "--precedence-only",
                        "-o", "{}/tmp/precedencecore-".format(eoi),
                        "{}".format(trace_path)))
        # Add a line to indicate that precedenceonly is True in core files.
        core_files = get_dot_files("{}/tmp".format(eoi), "precedencecore")
        for core_file in core_files:
            input_file = open("{}/tmp/{}".format(eoi, core_file), "r")
            content = input_file.readlines()
            input_file.close()
            content.insert(1, '  precedenceonly="True"\n')
            output_file = open("{}/tmp/{}".format(eoi, core_file), "w")
            output_file.writelines(content)
            output_file.close()
    elif precedenceonly == False:
        subprocess.run(("{}".format(kaflowpath),
                        "-o", "{}/tmp/causalcore-".format(eoi),
                        "{}".format(trace_path)))
        # Add a line to indicate that precedenceonly is True in core files.
        core_files = get_dot_files("{}/tmp".format(eoi), "causalcore")
        for core_file in core_files:
            input_file = open("{}/tmp/{}".format(eoi, core_file), "r")
            content = input_file.readlines()
            input_file.close()
            content.insert(1, '  precedenceonly="False"\n')
            output_file = open("{}/tmp/{}".format(eoi, core_file), "w")
            output_file.writelines(content)
            output_file.close()
        

def run_kastor(eoi, trace_path, kastorpath, compression):
    """ Run KaStor with fill_siphon and optionally weak compression. """

    if compression == None:
        subprocess.run(("{}".format(kastorpath), "--none",
                        "{}".format(trace_path)))
    else: # This takes a really long time.
        subprocess.run(("{}".format(kastorpath), "--{}".format(compression),
                        "{}".format(trace_path)))
    
    

# ---------------- End of Causal Cores Generation Section  --------------------

# ........................ Siphon Filling Section .............................

def fillsiphon(eoi, kappamodel, kastorpath):
    """
    Add spurious init event to break causal dependences using KaStor's
    fill_siphon function.
    """

    # Check time.
    time_start = time.perf_counter()
    # Open original trace file.
    period = kappamodel.index(".")
    modelprefix = kappamodel[:period]
    tracefile = open("{}/{}-eoi.json".format(eoi, modelprefix))
    sim = json.load(tracefile)

    causal_core_files = get_dot_files("{}/tmp".format(eoi), "causalcore")
    for causal_core_file in causal_core_files:
        # Create trace from causal core.
        core_path = "{}/tmp/{}".format(eoi, causal_core_file)
        corefile = open(core_path, "r").readlines()
        dash = causal_core_file.index("-")
        period = causal_core_file.index(".")
        filenumber = causal_core_file[dash+1:period]
        trace_path = "{}/tmp/causalcore-{}.json".format(eoi, filenumber)
        tracefile = open(trace_path, 'w')
        # Find the index of the events that are part of the causal core.
        event_indexes = []
        for line in corefile:
            if "style=filled" in line:
                tokens = line.split()
                event_index = int(tokens[0])
                event_indexes.append(event_index)
        # Build core trace with only events from the causal core.
        trace = sim["trace"]
        new_trace = []
        for i in range(len(trace)):
            if i in event_indexes:
                new_trace.append(trace[i])
        new_sim = sim.copy()
        new_sim["trace"] = new_trace
        # Write core trace to file.
        print("\nWriting temporary core trace {}".format(trace_path))
        json.dump(new_sim, tracefile)
        tracefile.close()
        # Get story with fill siphon from core trace using KaStor.
        subprocess.run(("{}".format(kastorpath), "--none", "{}".format(trace_path)))
        os.rename("cflow_none_0.dot", "{}/tmp/siphon-{}.dot".format(eoi, filenumber))
        os.remove("cflow_none_Summary.dat")
        os.remove(trace_path)
    # Clean up and check calculation time.
    os.remove("compression_status.txt")
    os.remove("profiling.html")
    time_stop = time.perf_counter()
    time_diff = time_stop - time_start
    print("\nCalculation time : {:.2f}s\n".format(time_diff))

# ..................... End of Siphon Filling Section .........................

def tweakstories(eoi, showintro=True, addedgelabels=False,
                 showedgelabels=False, edgeid=True, edgeocc=False,
                 edgeuse=False, statstype="abs", writedot=True,
                 weightedges=True):
    """ Tweak stories for deemed better readability. """

    # Reading section.
    story_files = get_dot_files("{}/tmp".format(eoi), "siphon")
    stories = []
    for story_file in story_files:
        story_path = "{}/tmp/{}".format(eoi, story_file)
        stories.append(CausalGraph(story_path, eoi))
    # Tweak each story.
    for i in range(len(stories)):
        stories[i].filename = "story-{}.dot".format(i+1)
    for story in stories:
        #story.compute_relstats()
        #story.compute_visuals(showintro, color=False)
        story.build_dot_file(showintro, addedgelabels, showedgelabels,
                             edgeid, edgeocc, edgeuse, statstype, weightedges)
        output_path = "{}/{}".format(eoi, story.filename)
        outfile = open(output_path, "w")
        outfile.write(story.dot_file)
        outfile.close()

#def getcustomstories(eoi, kappamodel, kasimpath, kaflowpath, kastorpath,
#                     simtime=1000, simseed=None):
#    """
#    Build custom stories using a combination of causal cores from KaFlow and
#    compressed stories from KaStor.
#    """
#
#    new_model = add_eoi(eoi, kappamodel)
#    trace_path = run_kasim(new_model, kasimpath, simtime, simseed)
#    run_kaflow(eoi, trace_path, kaflowpath, precedenceonly=True)
#    run_kaflow(eoi, trace_path, kaflowpath, precedenceonly=False)
#    fillsiphon(eoi, kappamodel, kastorpath)

def showedits(eoi, kappamodel, showintro=True, addedgelabels=False,
              showedgelabels=False, edgeid=True, edgeocc=False,
              edgeuse=False, statstype="abs", writedot=True,
              weightedges=True):
    """ Add dual nodes showing the states edited by every event. """

    # Reading section.
    story_files = get_dot_files("{}".format(eoi), "story")
    stories = []
    for story_file in story_files:
        story_path = "{}/{}".format(eoi, story_file)
        stories.append(CausalGraph(story_path, eoi))
    # Read trace.
    period = kappamodel.rfind(".")
    modelname = kappamodel[:period]
    tracefile = open("{}/{}-eoi.json".format(eoi, modelname))
    sim = json.load(tracefile)
    signatures = sim["model"]["update"]["signatures"]
    #rules = sim["model"]["ast_rules"]
    steps = sim["trace"]
    # Write stories with state edits.
    for story in stories:
        # Get actions for each event node.
        for node in story.eventnodes:
            node.states = []
            step = steps[int(node.nodeid)]
            bnd_num = 1
            if step[0] == 1: # Rule
                actions = step[2][1]
                tmp_states = []
                for action in actions:
                    state, bnd_num = state_from_action(signatures, action,
                                                       bnd_num)
                    tmp_states.append(state)
                node.states = tmp_states
            if step[0] == 3: # Init
                actions = step[1]
                tmp_states = []
                for action in actions:
                    state, bnd_num = state_from_action(signatures, action,
                                                       bnd_num)
                    # Check if last state is the same Bind_to as a previous
                    # one but reversed.
                    add_state = True
                    if action[0] == 3:
                        ag1 = state[0]["agent"]
                        ag2 = state[1]["agent"]
                        id1 = state[0]["agentid"]
                        id2 = state[1]["agentid"]
                        for prev_state in tmp_states:
                            if prev_state[0]["action"] == 3:
                                pa1 = prev_state[0]["agent"]
                                pa2 = prev_state[1]["agent"]
                                pi1 = prev_state[0]["agentid"]
                                pi2 = prev_state[1]["agentid"]
                                if pa1 == ag2 and pi1 == id2:
                                    if pa2 == ag1 and pi2 == id1:
                                        add_state = False
                    if add_state == True:
                        tmp_states.append(state)
                # Remove creations if the agents are present in other actions.
                for tmp_state in tmp_states:
                    keep_state = True
                    if tmp_state[0]["action"] == 0:
                        ag1 = state[0]["agent"]
                        id1 = state[0]["agentid"]
                        for tmp_state2 in tmp_states:
                            if tmp_state2[0]["action"] != 0:
                                for elem in tmp_state2:
                                    ag2 = elem["agent"]
                                    id2 = elem["agentid"]
                                    if ag1 == ag2 and id1 == id2:
                                        keep_state = False
                    if keep_state == True:    
                        node.states.append(tmp_state)
        # Add a StateNode for each state of EventNodes.
        state_id = 1
        for node in story.eventnodes:
            for state in node.states:
                node_id = "state{}".format(state_id)
                rank = node.rank + 0.5
                state_str = ""
                for i in range(len(state)):
                    agent_str = write_kappa_expression(state[i], "num")
                #    state_str += "{}".format(state[i]["agent"])
                #    state_str += ":{}".format(state[i]["agentid"])
                #    state_str += "({}".format(state[i]["site"])
                #    if state[i]["bond"] != None and state[i]["bond"] != ".":
                #        state_str += "[{}]".format(state[i]["bond"]["num"])
                #    if state[i]["value"] != None:
                #        state_str += "{{{}}}".format(state[i]["value"])
                #    state_str += ")"
                    state_str += agent_str
                    if i < len(state)-1:
                        state_str += ", "
                label = state_str
                new_state_node = StateNode(node_id, state, label, rank)
                story.statenodes.append(new_state_node)
                new_edge = CausalEdge(node, new_state_node)
                story.causaledges.append(new_edge)
                story.hyperedges.append(HyperEdge([new_edge]))
                state_id += 1
        # Get tests for each event node.
        for node in story.eventnodes:
            node.tests = []
            step = steps[int(node.nodeid)]
            bnd_num = 1
            tests = None
            if step[0] == 1: # Rule
                tests_tmp = step[2][0]
                tests = []
                for sublist in tests_tmp:
                    for test in sublist:
                        tests.append(test)
            if step[0] == 4: # Obs
                tests = step[1][0]
            tmp_states = []
            if tests != None:
                for test in tests:
                    state, bnd_num = state_from_test(signatures, test,
                                                     bnd_num)
                    print("xxx", state)
                    if state[0]["test"] != 0:
                        tmp_states.append(state)
                node.tests = tmp_states
        # Add edges showing for which events each state node is required.
        for statenode in story.statenodes:
            for eventnode in story.eventnodes:
                editused = False
                for eventtest in eventnode.tests:
                    are_same = compare_state_test(statenode.state, eventtest, statenode.label, eventnode.label)
                    if are_same == True:
                        editused = True
                        break
                if editused == True:
                    new_edge = CausalEdge(statenode, eventnode)
                    story.causaledges.append(new_edge)
                    story.hyperedges.append(HyperEdge([new_edge]))
        # Remove all edges between two event nodes.
        edges_to_remove = []
        for i in range(len(story.causaledges)):
            edge = story.causaledges[i]
            eventsrc = isinstance(edge.source, EventNode)
            eventtrg = isinstance(edge.target, EventNode)
            if eventsrc == True and eventtrg == True:
                edges_to_remove.insert(0, i)
        for i in edges_to_remove:
            del(story.causaledges[i])
        for hyperedge in story.hyperedges:
            edges_to_remove = []
            for i in range(len(hyperedge.edgelist)):
                edge = hyperedge.edgelist[i]
                eventsrc = isinstance(edge.source, EventNode)
                eventtrg = isinstance(edge.target, EventNode)
                if eventsrc == True and eventtrg == True:
                    edges_to_remove.insert(0, i)
            for i in edges_to_remove:
                del(hyperedge.edgelist[i])
    # Write stories with edited states.
    for i in range(len(stories)):
        stories[i].filename = "edits-{}.dot".format(i+1)
    for story in stories:
        story.build_dot_file(showintro, addedgelabels, showedgelabels,
                             edgeid, edgeocc, edgeuse, statstype, weightedges)
        output_path = "{}/{}".format(eoi, story.filename)
        outfile = open(output_path, "w")
        outfile.write(story.dot_file)
        outfile.close()


#def write_kappa_expression(state, bond="num"):
#    """
#    Write a state as a string using Kappa language.
#    The value of bond can be either 'num' or 'partner'.
#    """
#
#    state_str = ""
#    for i in range(len(state)):
#        state_str += "{}".format(state[i]["agent"])
#        state_str += ":{}".format(state[i]["agentid"])
#        state_str += "({}".format(state[i]["site"])
#        if state[i]["bond"] != None and state[i]["bond"] != ".":
#            if bond == "num":
#                state_str += "[{}]".format(state[i]["bond"]["num"])
#            elif bond == "partner":
#                partner = state[i]["bond"]["partner"]
#                state_str += "[{}.{}:{}]".format(partner["site"],
#                                                 partner["agent"],
#                                                 partner["agentid"])
#        if state[i]["value"] != None:
#            state_str += "{{{}}}".format(state[i]["value"])
#        state_str += ")"
#        if i < len(state)-1:
#            state_str += ", "
#
#    return state_str

def write_kappa_expression(agent, bond="num"):
    """
    Write an agent as a string using Kappa language.
    The value of bond can be either 'num' or 'partner'.
    """

    agent_str = ""
    agent_str += "{}".format(agent["agent"])
    agent_str += ":{}".format(agent["agentid"])
    agent_str += "({}".format(agent["site"])
    if agent["bond"] != None:
        if agent["bond"] == ".":
            agent_str += "[.]"
        else:
            if bond == "num":
                agent_str += "[{}]".format(agent["bond"]["num"])
            elif bond == "partner":
                partner = agent["bond"]["partner"]
                agent_str += "[{}.{}:{}]".format(partner["site"],
                                                 partner["agent"],
                                                 partner["agentid"])
    if agent["value"] != None:
        agent_str += "{{{}}}".format(agent["value"])
    agent_str += ")"

    return agent_str


def compare_state_test(state, test, lab1, lab2):
    """ Determine if a state and test represent the same species. """

#--> [{'agent': 'Rec', 'agentid': 3313, 'site': 'g', 'bond': {'num': 1, 'partner': {'agent': 'Syk', 'agentid': 4093, 'site': 'tSH2'}}, 'value': None, 'action': 2},
#     {'agent': 'Syk', 'agentid': 4093, 'site': 'tSH2', 'bond': {'num': 1, 'partner': {'agent': 'Rec', 'agentid': 3313, 'site': 'g'}}, 'value': None, 'action': 2}]
#
#==> [{'agent': 'Syk', 'agentid': 4093, 'site': 'tSH2', 'bond': {'num': 5, 'partner': {'agent': 'Rec', 'agentid': 3313, 'site': 'g'}}, 'value': None, 'test': 5},
#     {'agent': 'Rec', 'agentid': 3313, 'site': 'g', 'bond': {'num': 5, 'partner': {'agent': 'Syk', 'agentid': 4093, 'site': 'tSH2'}}, 'value': None, 'test': 5}]

    #are_same = False
    #if lab1 == "Rec:3313(g[1]), Syk:4093(tSH2[1])" and lab2 == "Lyn-U phosphorylates Syk-l":
    #    print("-->", state)
    #    print("==>", test)
    list1 = state.copy()
    list2 = test.copy()
    found1 = []
    found2 = []
    for i in range(len(list1)):
        agent1_str = write_kappa_expression(list1[i], bond="partner")
        for agent2 in list2:
            agent2_str = write_kappa_expression(agent2, bond="partner")
            if agent1_str == agent2_str:
                found1.insert(0, i)
                break
    for i in range(len(list2)):
        agent2_str = write_kappa_expression(list2[i], bond="partner")
        for agent1 in list1:
            agent1_str = write_kappa_expression(agent1, bond="partner")
            if agent2_str == agent1_str:
                found2.insert(0, i)
                break
    for i in found1:
        del(list1[i])
    for j in found2:
        del(list2[j])
    if len(list1) == 0 and len(list2) == 0:
        are_same = True
    else:
        are_same = False

    return are_same


# same_objects
#    list1 = nodelist1.copy()
#    list2 = nodelist2.copy()
#    found1 = []
#    found2 = []
#    for i in range(len(list1)):
#        for node2 in list2:
#            if list1[i].label == node2.label:
#                if enforcerank == False:
#                    found1.insert(0, i)
#                    break
#                elif enforcerank == True:
#                    if list1[i].rank == node2.rank:
#                        found1.insert(0, i)
#                        break
#    for j in range(len(list2)):
#        for node1 in list1:
#            if list2[j].label == node1.label:
#                if enforcerank == False:
#                    found2.insert(0, j)
#                    break
#                elif enforcerank == True:
#                    if list2[j].rank == node1.rank:
#                        found2.insert(0, j)
#                        break
#    for i in found1:
#        del(list1[i])
#    for j in found2:
#        del(list2[j])
#    if len(list1) == 0 and len(list2) == 0:
#        are_equi = True
#    else:
#        are_equi = False

def state_from_action(signatures, action, bnd_num):
    """ Find the resulting state of an action from the trace file. """

    state = []
    if action[0] == 0: # Create (I only look at the agent, not the sites).
        ag_n = action[1][1]
        agid_n = action[1][0]
        entry = signatures[ag_n]
        agent = entry["name"]
        state.append({"agent":agent, "agentid": agid_n, "site":None,
                      "bond":None, "value":None, "action":0})
    if action[0] == 1: # Mod_internal
        ag_n = action[1][0][1]
        agid_n = action[1][0][0]
        site_n = action[1][1]
        val_n = action[2]
        entry = signatures[ag_n]
        agent = entry["name"]
        site = entry["decl"][site_n]["name"]
        value = entry["decl"][site_n]["decl"][0][val_n]["name"]
        state.append({"agent":agent, "agentid": agid_n, "site":site,
                      "bond":None, "value":value, "action":1})
    if action[0] == 2 or action[0] == 3: # Bind or Bind_to
        ag1_n = action[1][0][1]
        agid1_n = action[1][0][0]
        site1_n = action[1][1]
        entry1 = signatures[ag1_n]
        agent1 = entry1["name"]
        site1 = entry1["decl"][site1_n]["name"]
        ag2_n = action[2][0][1]
        agid2_n = action[2][0][0]
        site2_n = action[2][1]
        entry2 = signatures[ag2_n]
        agent2 = entry2["name"]
        site2 = entry2["decl"][site2_n]["name"]
        partner1 = {"agent":agent1, "agentid": agid1_n, "site":site1}
        partner2 = {"agent":agent2, "agentid": agid2_n, "site":site2}
        state.append({"agent":agent1, "agentid": agid1_n, "site":site1,
                      "bond":{"num": bnd_num, "partner":partner2},
                      "value":None, "action":action[0]})
        state.append({"agent":agent2, "agentid": agid2_n, "site":site2,
                      "bond":{"num": bnd_num, "partner":partner1},
                      "value":None, "action":action[0]})
        bnd_num += 1
    if action[0] == 4: # Free
        ag_n = action[1][0][1]
        agid_n = action[1][0][0]
        site_n = action[1][1]
        entry = signatures[ag_n]
        agent = entry["name"]
        site = entry["decl"][site_n]["name"]
        state.append({"agent":agent, "agentid": agid_n, "site":site,
                      "bond":".", "value":None, "action":4})
    #if action[0] == 5: # Remove (I still do not have any example).

    return state, bnd_num


def state_from_test(signatures, test, bnd_num):
    """ Find the required state of a test from the trace file. """

    state = []
    if test[0] == 0: # Is_Here
        ag_n = test[1][1]
        agid_n = test[1][0]
        entry = signatures[ag_n]
        agent = entry["name"]
        state.append({"agent":agent, "agentid": agid_n, "site":None,
                      "bond":None, "value":None, "test":0})
    if test[0] == 1: # Has_Internal
        ag_n = test[1][0][1]
        agid_n = test[1][0][0]
        site_n = test[1][1]
        val_n = test[2]
        entry = signatures[ag_n]
        agent = entry["name"]
        site = entry["decl"][site_n]["name"]
        value = entry["decl"][site_n]["decl"][0][val_n]["name"]
        state.append({"agent":agent, "agentid": agid_n, "site":site,
                      "bond":None, "value":value, "test":1})
    if test[0] == 2: # Is_Free
        ag_n = test[1][0][1]
        agid_n = test[1][0][0]
        site_n = test[1][1]
        entry = signatures[ag_n]
        agent = entry["name"]
        site = entry["decl"][site_n]["name"]
        state.append({"agent":agent, "agentid": agid_n, "site":site,
                      "bond":".", "value":None, "test":2})
    #if test[0] == 3: # Is_Bound (No example yet).
    #if test[0] == 4: # Has_Binding_type (No example yet).
    if test[0] == 5: # Is_Bound_to
        ag1_n = test[1][0][1]
        agid1_n = test[1][0][0]
        site1_n = test[1][1]
        entry1 = signatures[ag1_n]
        agent1 = entry1["name"]
        site1 = entry1["decl"][site1_n]["name"]
        ag2_n = test[2][0][1]
        agid2_n = test[2][0][0]
        site2_n = test[2][1]
        entry2 = signatures[ag2_n]
        agent2 = entry2["name"]
        site2 = entry2["decl"][site2_n]["name"]
        partner1 = {"agent":agent1, "agentid": agid1_n, "site":site1}
        partner2 = {"agent":agent2, "agentid": agid2_n, "site":site2}
        state.append({"agent":agent1, "agentid": agid1_n, "site":site1,
                      "bond":{"num": bnd_num, "partner":partner2},
                      "value":None, "test":test[0]})
        state.append({"agent":agent2, "agentid": agid2_n, "site":site2,
                      "bond":{"num": bnd_num, "partner":partner1},
                      "value":None, "test":test[0]})
        bnd_num += 1
        print("---", test)
        print("===", state)

    return state, bnd_num


# ==================== Causal Cores Merging Section ===========================

def mergecores(eoi, causalgraphs=None, siphon=False, showintro=True,
               addedgelabels=False, showedgelabels=False, edgeid=True,
               edgeocc=False, edgeuse=False, statstype="abs",
               weightedges=False, color=True, writedot=True, rmprev=False,
               msg=True):
    """
    Merge equivalent causal cores and count occurrence.
    Write the final cores as meshed graphs.
    """

    # Reading section.
    if causalgraphs == None:
        if siphon == False:
            causal_core_files = get_dot_files(eoi, "causalcore")
        elif siphon == True:
            causal_core_files = get_dot_files(eoi, "siphon")
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
        analogous_list = [0]
        for i in range(1, len(causal_cores)):
            same_core, equi_meshes = analogous_graphs(current_core,
                                                      causal_cores[i])
            if same_core == True:
                analogous_list.append(i)
                current_core.occurrence += causal_cores[i].occurrence
                for j in range(len(current_core.meshes)):
                    equi_index = equi_meshes[j]
                    uses = causal_cores[i].meshes[equi_index].uses
                    current_core.meshes[j].uses += uses
        prevcores = []
        for index in analogous_list:
            file_name = causal_cores[index].filename
            dash = file_name.rfind("-")
            period = file_name.rfind(".")
            if "_node" in file_name:
                underscore = file_name.index("_node")
                previd = file_name[:period]
            else:
                previd = file_name[dash+1:period]
            prevcores.append(previd)
        current_core.prevcores = prevcores
        merged_cores.append(current_core)
        for i in range(len(analogous_list)-1, -1, -1):
            index = analogous_list[i]
            del(causal_cores[index])
    sorted_cores = sorted(merged_cores, key=lambda x: x.occurrence,
                          reverse=True)
    for i in range(len(sorted_cores)):
        sorted_cores[i].filename = "meshedcore-{}.dot".format(i+1)
    for graph in sorted_cores:
        graph.compute_relstats()
        graph.compute_visuals(showintro, color)
        graph.build_dot_file(showintro, addedgelabels, showedgelabels,
                             edgeid, edgeocc, edgeuse, statstype, weightedges)
    # Writing section.
    if writedot == True:
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
    if msg == True:
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


def analogous_graphs(graph1, graph2):
    """
    Analogous causal graphs have analogous meshes. That is, all their meshes
    are between nodes with same labels at same ranks.
    """

    equi_meshes = []
    if graph1.maxrank == graph2.maxrank:
        graph2_indexes = list(range(len(graph2.meshes)))
        all_edges_found = True
        for mesh1 in graph1.meshes:
            for i in graph2_indexes:
                mesh2 = graph2.meshes[i]
                are_equi = analogous_meshes(mesh1, mesh2, enforcerank=True)
                if are_equi == True:
                    equi_meshes.append(i)
                    graph2_indexes.remove(i)
                    break
            if are_equi == False:
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

    return equi_graphs, equi_meshes


def analogous_meshes(mesh1, mesh2, enforcerank=True):
    """
    Find whether two meshes connect to event nodes with same labels
    with analogous midedges.
    Optionally, nodes may be also required to be at same ranks.
    """

    nn1 = len(mesh1.midnodes)
    nn2 = len(mesh2.midnodes)
    ne1 = len(mesh1.midedges)
    ne2 = len(mesh2.midedges)
    if nn1 == nn2 and ne1 == ne2:
        are_equi = True
    else:
        are_equi = False
    if are_equi == True:
        sources1, targets1 = mesh1.get_events()
        sources2, targets2 = mesh2.get_events()
        equi_sources = analogous_nodes(sources1, sources2, enforcerank)
        equi_targets = analogous_nodes(targets1, targets2, enforcerank)
        if equi_sources == True and equi_targets == True:
            are_equi = True
        else:
            are_equi = False
    if are_equi == True:
        neighbors1 = mesh1.extend_midedges()
        neighbors2 = mesh2.extend_midedges()
        equi_midedges = analogous_midedges(neighbors1, neighbors2,
                                           enforcerank)
        if equi_midedges == True:
            are_equi = True
        else:
            are_equi = False

    return are_equi


def analogous_nodes(nodelist1, nodelist2, enforcerank=True):
    """
    Find whether two lists of nodes contain nodes with
    same labels.
    Optionally, nodes may be also required to be at same ranks.
    (This is comparable to the function "same_objects" used in
    method "equivalent_meshes".)
    """

    list1 = nodelist1.copy()
    list2 = nodelist2.copy()
    found1 = []
    found2 = []
    for i in range(len(list1)):
        for node2 in list2:
            if list1[i].label == node2.label:
                if enforcerank == False:
                    found1.insert(0, i)
                    break
                elif enforcerank == True:
                    if list1[i].rank == node2.rank:
                        found1.insert(0, i)
                        break
    for j in range(len(list2)):
        for node1 in list1:
            if list2[j].label == node1.label:
                if enforcerank == False:
                    found2.insert(0, j)
                    break
                elif enforcerank == True:
                    if list2[j].rank == node1.rank:
                        found2.insert(0, j)
                        break
    for i in found1:
        del(list1[i])
    for j in found2:
        del(list2[j])
    if len(list1) == 0 and len(list2) == 0:
        are_equi = True
    else:
        are_equi = False

    return are_equi


def analogous_midedges(neighbors1, neighbors2, enforcerank=True):
    """
    Find whether two lists of midedges (described as their respective
    neighbors) connect to event nodes with same labels.
    Optionally, nodes may be also required to be at same ranks.
    (This is comparable to the function "equivalent_midedges" used in
    method "equivalent_meshes".)
    """

    list1 = neighbors1.copy()
    list2 = neighbors2.copy()
    found1 = []
    found2 = []
    for i in range(len(list1)):
        s1 = list1[i]["srcs"]
        t1 = list1[i]["trgs"]
        for j in range(len(list2)):
            s2 = list2[j]["srcs"]
            t2 = list2[j]["trgs"]
            if list1[i]["reltype"] == list2[j]["reltype"]:
                if analogous_nodes(s1, s2, enforcerank):
                    if analogous_nodes(t1, t2, enforcerank):
                        found1.insert(0, i)
                        break
    for j in range(len(list2)):
        s2 = list2[j]["srcs"]
        t2 = list2[j]["trgs"]
        for i in range(len(list1)):
            s1 = list1[i]["srcs"]
            t1 = list1[i]["trgs"]
            if list2[j]["reltype"] == list1[i]["reltype"]:
                if analogous_nodes(s2, s1, enforcerank):
                    if analogous_nodes(t2, t1, enforcerank):
                        found2.insert(0, j)
                        break
    for i in found1:
        del(list1[i])
    for j in found2:
        del(list2[j])
    if len(list1) == 0 and len(list2) == 0:
        are_equi = True
    else:
        are_equi = False

    return are_equi

# ================ End of Causal Cores Merging Section ========================

# .................. Event Paths Merging Section ..............................

def foldcores(eoi, causalgraphs=None, siphon=False, ignorelist=[],
              showintro=False, addedgelabels=True, showedgelabels=True,
              edgeid=True, edgeocc=False, edgeuse=True, statstype="rel",
              weightedges=True, color=True, writedot=True, rmprev=False):
    """ Fold meshed cores into a single event pathway. """

    # Reading section.
    if causalgraphs == None:
        # Using cores or eventpaths both work. But it can suffle the nodes
        # horizontally, yielding a different graph, but with same ranks for
        # all nodes.
        core_files = get_dot_files(eoi, "meshedcore")
        meshedcores = []
        for core_file in core_files:
            core_path = "{}/{}".format(eoi, core_file)
            meshedcores.append(CausalGraph(core_path, eoi))
    else:
        meshedcores = causalgraphs
        core_files = None
    # Doing the work.
    flush_ignored(meshedcores, core_files, ignorelist)
    pathway = CausalGraph(eoi=eoi, meshedgraph=True)
    pathway.occurrence = 0
    node_number = 1
    seen_labels = []
    midid = 1
    for meshedcore in meshedcores:
        # Add nodes.
        pathway.occurrence += meshedcore.occurrence
        for node in meshedcore.eventnodes:
            if node.label not in seen_labels:
                seen_labels.append(node.label)
                n_id = "node{}".format(node_number)
                pathway.eventnodes.append(EventNode(n_id, node.label,
                                                    node.rank,
                                                    intro=node.intro,
                                                    first=node.first))
                node_number += 1
        # Add meshes (edges).
        for mesh in meshedcore.meshes:
            mesh_found = False
            for pathwaymesh in pathway.meshes:
                if analogous_meshes(mesh, pathwaymesh, enforcerank=False):
                    mesh_found = True
                    pathwaymesh.uses += mesh.uses
                    pathwaymesh.weight += mesh.uses
                    break
            if mesh_found == False:
                add_mesh(pathway, mesh, midid)
                midid += len(mesh.midnodes)

    # Uncomment the next 3 lines and comment pathway.rank_sequentially()
    # to build unranked version of graph
    #pathway.rank_intermediary()
    #pathway.get_maxrank()
    #pathway.sequentialize_ids()
    pathway.rank_sequentially()
    pathway.filename = "eventpathway.dot"
    compute_mesh_occurrence(eoi, pathway)
    pathway.compute_visuals(showintro, color)
    pathway.compute_relstats()
    pathway.build_dot_file(showintro, addedgelabels, showedgelabels, edgeid,
                           edgeocc, edgeuse, statstype, weightedges)
    # Writing section.
    if writedot == True:
        output_path1 = "{}/{}".format(eoi, pathway.filename)
        outfile1 = open(output_path1, "w")
        outfile1.write(pathway.dot_file)
        outfile1.close()
    if rmprev == True:
        if path_files == None:
            path_files = get_dot_files(eoi, "meshedcore")
        for path_file in path_files:
            file_path = "{}/{}".format(eoi, path_file)
            os.remove(file_path)
    print("Merging all event paths into one event pathway.")

    return pathway        


def flush_ignored(graph_list, graph_files, ignorelist, msg=True):
    """
    Remove cores that contain any ignored term in any of its nodes.
    """

    init_len = len(graph_list)
    graphs_to_remove = []
    for i in range(len(graph_list)):
        graph = graph_list[i]
        remove_graph = False
        for node in graph.eventnodes:
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
    if msg == True:
        print("Ignoring {} cores out of {} because they contain reverse rules."
              .format(len(graphs_to_remove), init_len))


def add_mesh(graph, mesh, startid, insertpos=None):
    """
    Add a new mesh to a graph based on the labels of the event nodes found in
    the graph and in the source or target of midedges from the mesh.
    """

    new_mesh = Mesh(uses=mesh.uses, color=mesh.color)
    midid = startid
    midid_map = {}
    for midnode in mesh.midnodes:
        new_midnode = MidNode("mid{}".format(midid), midtype=midnode.midtype)
        new_midnode.fillcolor = midnode.fillcolor
        new_midnode.bordercolor = midnode.bordercolor
        new_mesh.midnodes.append(new_midnode)
        midid_map[midnode.nodeid] = "mid{}".format(midid)
        midid += 1
    for midedge in mesh.midedges:
        if isinstance(midedge.source, EventNode):
            for node in graph.eventnodes:
                if node.label == midedge.source.label:
                    source = node
                    break
        elif isinstance(midedge.source, MidNode):
            for new_mid in new_mesh.midnodes:
                if new_mid.nodeid == midid_map[midedge.source.nodeid]:
                    source = new_mid
                    break
        if isinstance(midedge.target, EventNode):
            for node in graph.eventnodes:
                if node.label == midedge.target.label:
                    target = node
                    break
        elif isinstance(midedge.target, MidNode):
            for new_mid in new_mesh.midnodes:
                if new_mid.nodeid == midid_map[midedge.target.nodeid]:
                    target = new_mid
                    break
        new_midedge = MidEdge(source, target, uses=mesh.uses,
                              relationtype=midedge.relationtype,
                              color=midedge.color)
        new_mesh.midedges.append(new_midedge)
    if insertpos == None:
        graph.meshes.append(new_mesh)
    else:
        graph.meshes.insert(insertpos, new_mesh)


def compute_mesh_occurrence(eoi, graph):
    """ Compute the occurrence of every midedge of each mesh of a pathway. """

    # I need the original cores (or the fill_siphon)
    causal_core_files = get_dot_files(eoi, "causalcore")
    causal_cores = []
    for core_file in causal_core_files:
        core_path = "{}/{}".format(eoi, core_file)
        causal_cores.append(CausalGraph(core_path, eoi))
    # For every pathway mesh, gather all analogous causal core mesh.
    core_mesh_lists = []
    for i in range(len(graph.meshes)):
        core_mesh_lists.append([])
    for causal_core in causal_cores:
        for core_mesh in causal_core.meshes:
            for i in range(len(graph.meshes)):
                if analogous_meshes(core_mesh, graph.meshes[i]):
                    core_mesh_lists[i].append(core_mesh)
    for i in range(len(graph.meshes)):
        print(len(core_mesh_lists[i]))

    # Compute occurrence for each central midedge of each pathway mesh.
    for i in range(len(graph.meshes)):
        mesh = graph.meshes[i]
        core_meshes = core_mesh_lists[i]
        # Get midedge touching a midnode (central edges).
        central_edges = []
        for midedge in mesh.midedges:
            s = midedge.source
            t = midedge.target
            if isinstance(s, MidNode) or isinstance(t, MidNode):
                central_edges.append(midedge)
        central_srcs, central_trgs = get_cen_srcs_trgs(central_edges, mesh)
        # For each central edge, get the list of corresponding edges from the
        # causal cores, based on src_trg_dict.
        core_edges = []
        for j in range(len(central_edges)):
            core_edges.append([])
        for core_mesh in core_meshes:
            central_core_edges = []
            for midedge in core_mesh.midedges:
                s = midedge.source
                t = midedge.target
                if isinstance(s, MidNode) or isinstance(t, MidNode):
                    central_core_edges.append(midedge)
            core_srcs, core_trgs = get_cen_srcs_trgs(central_core_edges,
                                                     core_mesh)
            for j in range(len(central_srcs)):
                seen_core_edges = []
                for k in range(len(core_srcs)):
                    if k not in seen_core_edges:
                        same_srcs = same_objects(central_srcs[j], core_srcs[k])
                        same_trgs = same_objects(central_trgs[j], core_trgs[k])
                        if same_srcs and same_trgs:
                            core_edges[j].append(central_core_edges[k])
                            seen_core_edges.append(k)
        # Here, each list inside core_edges should contain all the midedges
        # from causal cores that corresponds to each central midedges of the
        # pathway
            
            
def get_cen_srcs_trgs(edge_list, mesh):
    """ Get the source events and target events of each central edge. """

    src_list = []
    trg_list = []
    for edge in edge_list:
        srcs = []
        if isinstance(edge.source, EventNode):
            srcs.append(edge.source)
        elif isinstance(edge.source, MidNode):
            for midedge in mesh.midedges:
                if midedge.target == edge.source:
                    srcs.append(midedge.source)
        src_list.append(srcs)
        trgs = []
        if isinstance(edge.target, EventNode):
            trgs.append(edge.target)
        elif isinstance(edge.target, MidNode):
            for midedge in mesh.midedges:           
                if midedge.source == edge.target:
                    trgs.append(midedge.target)
        srcs_trgs.append(src_trg_dict)
        trg_list.append(trgs)

    return src_list, trg_list

# ............... End of Event Paths Merging Section ..........................

# ///////////////// Event Paths Simplifying Section ///////////////////////////

def simplifypathway(eoi, causalgraphs=None, threshold=0.2, edgelabels=False,
                    showintro=False, color=True, writedot=True, rmprev=False,
                    weightedges=False):
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
    # is under average_prob*threshold.
    pathway.build_nointro()
    all_probs = []
    for edge in pathway.hyperedges:
        if edge.underlying == False:
            all_probs.append(edge.prob)
    for cedge in pathway.coveredges:
        all_probs.append(cedge.prob)
    average_prob = statistics.mean(all_probs)
    theshold_str = "Average edge prob: {:.2f} , ".format(average_prob)
    theshold_str += "Treshold: {:.2f} , ".format(threshold)
    theshold_str += "Cutoff: {:.2f}\n".format(average_prob*threshold)
    theshold_str += "Simplifying graph; ignoring story types containing "
    theshold_str += ("edges with probablity lower than {:.2f}"
                      .format(average_prob*threshold))
    print(theshold_str)
    normals_to_ignore = []
    for edge in pathway.hyperedges:
        if edge.underlying == False:
            if edge.prob < average_prob*threshold:
                normals_to_ignore.append(edge)
    covers_to_ignore = []
    for cedge in pathway.coveredges:
        if cedge.prob < average_prob*threshold:
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
                                                       edge.prob))
    fuse_edges(simplepathway)
    simplepathway.rank_sequentially()
    simplepathway.filename = "eventpathway-simple.dot"
    simplepathway.build_dot_file(showintro, addedgelabels, showedgelabels,
                                 edgeid, edgeocc, edgeprob, weightedges)
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

def mapcores(eoi, causalgraphs=None, template=None, ignorelist=[],
             showintro=False, addedgelabels=True, showedgelabels=True,
             edgeid=True, edgeocc=False, edgeprob=True, weightedges=True,
             writedot=True, rmprev=False, msg=True, transdist=18):
    """
    Create a new CausalGraph for each core, where the core is mapped on
    the layout of the event pathway.
    """

    # Reading cores.
    if causalgraphs == None:
        core_files = get_dot_files(eoi, "meshedcore")
        meshedcores = []
        for core_file in core_files:
            core_path = "{}/{}".format(eoi, core_file)
            meshedcores.append(CausalGraph(core_path, eoi))
    else:
        meshedcores = causalgraphs
        core_files = None
    flush_ignored(meshedcores, core_files, ignorelist, msg=msg)
    # Reading event pathway template.
    if template == None:
        template_path = "{}/eventpathway.dot".format(eoi)
    else:
        template_path = template
    templategraph = CausalGraph(template_path, eoi)
    produce_layout(eoi, templategraph)
    occtot = templategraph.occurrence
    ##read_layout(eoi, templategraph, showintro)
    ##templategraph.filename = "mapped.dot"
    ##templategraph.build_dot_file(showintro, addedgelabels, showedgelabels,
    ##                                  edgeid, edgeocc, edgeprob, weightedges,
    ##                                  color=True)
    ##output_path = "{}/{}".format(eoi, templategraph.filename)
    ##outfile = open(output_path, "w")
    ##outfile.write(templategraph.dot_file)
    ##outfile.close()
    # Doing the work.
    mappedcores = []
    #for meshedcore in [meshedcores[12]]:
    for meshedcore in meshedcores:
        if len(meshedcore.covermeshes) == 0:
            meshedcore.build_nointro()
        mappedcore = CausalGraph(template_path, eoi)
        read_layout(eoi, mappedcore, showintro)
        mappedcore.compute_relstats()
        mappedcore.occurrence = meshedcore.occurrence
        mappedcore.prevcores = meshedcore.prevcores
        edges_to_add = []
        midid = mappedcore.find_max_midid(cover=True)+1
        meshid = mappedcore.find_max_meshid(cover=True)+1
        maxr = meshedcore.maxrank-1
        # Count the occurrence of every mapped core mesh in the meshed core
        # and translate the mapped core meshes accordingly.
        for mappedmesh in mappedcore.meshes:
            mappedmesh.totcount = 0
            for coremesh in meshedcore.meshes:
                if analogous_meshes(coremesh, mappedmesh, enforcerank=False):
                    mappedmesh.totcount += 1
            if mappedmesh.totcount > 1:
                offset = - (mappedmesh.totcount-1) * transdist / 2
                translate_mesh(mappedmesh, 1, offset)
        for cmappedmesh in mappedcore.covermeshes:
            cmappedmesh.totcount = 0
            for ccoremesh in meshedcore.covermeshes:
                if analogous_meshes(ccoremesh, cmappedmesh, enforcerank=False):
                    cmappedmesh.totcount += 1
            if cmappedmesh.totcount > 1:
                offset = - (cmappedmesh.totcount-1) * transdist / 2
                translate_mesh(cmappedmesh, 1, offset)
        # Find mesh ranks in core.
        for meshlist in [meshedcore.meshes, meshedcore.covermeshes]:
            for mesh in meshlist:
                sources, targets = mesh.get_events()
                source_ranks = []
                for source in sources:
                    source_ranks.append(source.rank)
                if len(source_ranks) == 0:
                    source_ranks = [0]
                mesh.rank = max(source_ranks)
        # Set default color and mesh use counter in mapped core template.
        for meshlist in [mappedcore.meshes, mappedcore.covermeshes]:
            for mesh in meshlist:
                mesh.count = 0
                mesh.color = "grey80"
                for midnode in mesh.midnodes:
                    midnode.bordercolor = "grey80"
                    if midnode.midtype == "enabling":
                        midnode.fillcolor = "grey80"
                for midedge in mesh.midedges:
                    midedge.color = "grey80"
        # Map meshes on mapped core template.
        additional_meshes = []
        additional_covermeshes = []
        for current_rank in range(maxr+1):
            midid, meshid = mapmesh(meshedcore.meshes,
                                    mappedcore.meshes,
                                    additional_meshes,
                                    current_rank, midid, meshid, maxr,
                                    transdist)
            midid, meshid = mapmesh(meshedcore.covermeshes,
                                    mappedcore.covermeshes,
                                    additional_covermeshes,
                                    current_rank, midid, meshid, maxr,
                                    transdist)
        for addmesh in additional_meshes:
            mappedcore.meshes.append(addmesh)
        for addcovermesh in additional_covermeshes:
            mappedcore.covermeshes.append(addcovermesh)
        mappedcores.append(mappedcore)
    for i in range(len(mappedcores)):
        mappedcores[i].filename = "mapped-{}.dot".format(i+1)
    for mappedcore in mappedcores:
        mappedcore.occurrence = "{} / {}".format(mappedcore.occurrence, occtot)
        mappedcore.build_dot_file(showintro, addedgelabels, showedgelabels,
                                  edgeid, edgeocc, edgeprob, weightedges)
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

    return mappedcores


def produce_layout(eoi, graph):
    """ Produce dot file with explicit positions. """

    subprocess.run(("/usr/bin/dot", "{}".format(graph.filename), "-o", 
                    "{}/layout.dot".format(eoi)))


def read_layout(eoi, graph, si):
    """ Read the layout produced by dot. """

    layout_file = open("{}/layout.dot".format(eoi), "r").readlines()
    posdict = {}
    labelposdict = {}
    rankposdict = {}
    for i in range(len(layout_file)):
        if ";" in layout_file[i]:
            start_line = i+1
            break
    for i in range(start_line, len(layout_file)):
        line = layout_file[i]
        tokens = line.split()
        if tokens[0][:4] == "node" or tokens[0][:3] == "mid":
            if "->" not in line and "midtype" not in line:
                idstr = tokens[0]
        if "rank=same" in line:
            tokens2 = layout_file[i+1].split()
            idstr = tokens2[0]
        if "->" in line:
            idstr = "{} -> {}".format(tokens[0], tokens[2])
        if "pos=" in line:
            read_start = line.index("pos=") + 4
            posstr = line[read_start:-2]
            line2 = line
            j = 1
            # The position of long edges may span several lines.
            while line2[-2] == "\\":
                line2 = layout_file[i+j]
                posstr += line2[:-2]
                j += 1
            if "node" in idstr or "mid" in idstr:
                posdict[idstr] = posstr
            else:
                rankposdict[idstr] = posstr
        if tokens[0][:3] == "lp=":
            read_start = line.index("lp=") + 3
            labelposdict[idstr] = line[read_start:-2]
    for eventnode in graph.eventnodes:
        if si == True or (si == False and eventnode.intro == False):
            eventnode.pos = posdict[eventnode.nodeid]
    for meshlist in [graph.meshes, graph.covermeshes]:
        for mesh in meshlist:
            if si == True or (si == False and mesh.underlying == False):
                for midnode in mesh.midnodes:
                    midnode.pos = posdict[midnode.nodeid]
                for midedge in mesh.midedges:
                    if midedge.reverse == False:
                        idstr = "{} -> {}".format(midedge.source.nodeid,
                                                  midedge.target.nodeid)
                    elif midedge.reverse == True:
                        idstr = "{} -> {}".format(midedge.target.nodeid,
                                                  midedge.source.nodeid)
                    midedge.pos = posdict[idstr]
                    if midedge.labelcarrier == True:
                        midedge.labelpos = labelposdict[idstr]
    graph.rankposdict = rankposdict


def mapmesh(coremeshes, mapmeshes, add_list, current_rank, midid, meshid,
            maxr, transdist):
    """ Map a mesh from meshed core onto the mapped core template. """

    for coremesh in coremeshes:
       if coremesh.rank == current_rank:
           # Find analoguous mesh in mapped core template.
           analog_meshes = []
           for i in range(len(mapmeshes)):
               mappedmesh = mapmeshes[i]
               if analogous_meshes(coremesh, mappedmesh, enforcerank=False):
                   analog_meshes.append(i)
           if len(analog_meshes) != 1:
               raise ValueError("Exactly 1 analogous mesh "
                                "expected in template.")
           else:
               mindex = analog_meshes[0]
               if mapmeshes[mindex].count == 0:
                   # Change color and pensize of mapped core mesh.
                   set_core_colors(mapmeshes[mindex], coremesh, maxr)
                   for midnode in mapmeshes[mindex].midnodes:
                       midnode.overridewidth = 0.144 # math.sqrt(3)/12
                   for midedge in mapmeshes[mindex].midedges:
                       midedge.overridewidth = 3
                       midedge.labelcarrier = False
               elif mapmeshes[mindex].count > 0:
                   # Duplicate mesh.
                   dupl_mesh = duplicate_mesh(mapmeshes[mindex],
                                              midid, meshid)
                   midid += len(mapmeshes[mindex].midnodes)
                   meshid += 1
                   set_core_colors(dupl_mesh, coremesh, maxr)
                   for midnode in dupl_mesh.midnodes:
                       midnode.overridewidth = 0.144 # math.sqrt(3)/12
                   for midedge in dupl_mesh.midedges:
                       midedge.overridewidth = 3
                   use_count = mapmeshes[mindex].count
                   translate_mesh(dupl_mesh, use_count, transdist)
                   add_list.append(dupl_mesh)
               mapmeshes[mindex].count += 1

    return midid, meshid


def set_core_colors(mappedmesh, coremesh, max_rank):
    """ Assign colors of mapped mesh based on ranks of core mesh. """

    # Get core source ranks
    msources, mtargets = mappedmesh.get_events()
    csources, ctargets = coremesh.get_events()
    for msource in msources:
        ranks = []
        for csource in csources:
            if csource.label == msource.label:
                ranks.append(csource.rank)
        if len(ranks) < 1:
            raise ValueError("Event node {} not found in core mesh"
                              .format(msource.label))
        else:
            msource.corerank = max(ranks)
    mappedmesh.corerank = coremesh.rank
    # Assign colors based on core ranks.
    # For each midnode and midedge, assign color based on highest
    # rank among reachable sources.
    col_val = 0.35+(mappedmesh.corerank/max_rank)*0.65
    mappedmesh.color = '"{:.3f} 1 1"'.format(col_val)
    for midnode in mappedmesh.midnodes:
        sources = mappedmesh.get_sources(midnode)
        source_ranks = []
        for source in sources:
            source_ranks.append(source.corerank)
        if len(source_ranks) == 0:
            source_ranks = [0]
        rank = max(source_ranks)
        col_val = 0.35+(rank/max_rank)*0.65
        midnode.bordercolor = '"{:.3f} 1 1"'.format(col_val)
        if midnode.midtype == "enabling":
            midnode.fillcolor = '"{:.3f} 1 1"'.format(col_val)
    for midedge in mappedmesh.midedges:
        if isinstance(midedge.source, MidNode):
            midedge.color = midedge.source.bordercolor
        elif isinstance(midedge.source, EventNode):
            col_val = 0.35+(midedge.source.corerank/max_rank)*0.65
            midedge.color = '"{:.3f} 1 1"'.format(col_val)


def duplicate_mesh(mapmesh, startid, meshid):
    """
    Build a copy of given mesh with new midids and a new meshid.
    """

    duplicate_mesh = Mesh(uses=mapmesh.uses, usage=mapmesh.usage,
                          underlying=mapmesh.underlying, color=mapmesh.color,
                          meshid=meshid)
    midid = startid
    midid_map = {}
    for midnode in mapmesh.midnodes:
        new_midnode = MidNode("mid{}".format(midid), midtype=midnode.midtype,
                              ghost=midnode.ghost, pos=midnode.pos)
        new_midnode.fillcolor = midnode.fillcolor
        new_midnode.bordercolor = midnode.bordercolor
        duplicate_mesh.midnodes.append(new_midnode)
        midid_map[midnode.nodeid] = "mid{}".format(midid)
        midid += 1
    for midedge in mapmesh.midedges:
        if isinstance(midedge.source, EventNode):
            source = midedge.source
        elif isinstance(midedge.source, MidNode):
            for new_mid in duplicate_mesh.midnodes:
                if new_mid.nodeid == midid_map[midedge.source.nodeid]:
                    source = new_mid
                    break
        if isinstance(midedge.target, EventNode):
            target = midedge.target
        elif isinstance(midedge.target, MidNode):
            for new_mid in duplicate_mesh.midnodes:
                if new_mid.nodeid == midid_map[midedge.target.nodeid]:
                    target = new_mid
                    break
        new_midedge = MidEdge(source, target, uses=mapmesh.uses,
            relationtype=midedge.relationtype, color=midedge.color,
            underlying=midedge.underlying, reverse=midedge.reverse,
            labelcarrier=midedge.labelcarrier, indicator=midedge.indicator,
            meshid=meshid, pos=midedge.pos, labelpos=midedge.labelpos)
        duplicate_mesh.midedges.append(new_midedge)

    return duplicate_mesh


def translate_mesh(mesh, use_count, dist):
    """
    Translate a mesh depending on how many duplicates
    already exist in the graph.
    """

    for midnode in mesh.midnodes:
        if midnode.pos != None:
            new_pos = change_coordinates(midnode.pos, dist, use_count)
            midnode.pos = new_pos
    for midedge in mesh.midedges:
        if midedge.pos != None:
            new_pos = change_coordinates(midedge.pos, dist, use_count)
            midedge.pos = new_pos
        if midedge.labelpos != None:
            new_pos = change_coordinates(midedge.labelpos, dist, use_count)
            midedge.labelpos = new_pos


def change_coordinates(pos_str, dist, use_count):
    """ Change coordinates to reflect a translation. """

    tokens = pos_str[1:-1].split()
    new_coords = []
    for token in tokens:
        if token[0] == "e" or token[0] == "s":
            coord = token[2:]
        else:
            coord = token
        comma = coord.index(",")
        xpos = int(float(coord[:comma]))
        new_xpos = xpos + (dist * use_count)
        new_coord = "{}{}".format(new_xpos, coord[comma:])
        if token[0] == "e":
            new_coord = "e,{}".format(new_coord)
        if token[0] == "s":
            new_coord = "s,{}".format(new_coord)
        new_coords.append(new_coord)
    new_pos = '"'
    for i in range(len(new_coords)):
        if i > 0:
            new_pos += ' '
        new_pos += '{}'.format(new_coords[i])
    new_pos += '"'

    return new_pos

# ++++++++++++ End of Core Mapping on Event Pathway Section ++++++++++++++++++

# ::::::::::::::::::::::: Node Highlight Section ::::::::::::::::::::::::::::::

def highlightnodes(eoi, nodelabels=None, causalgraphs=None, template=None,
                   ignorelist=[], showintro=False, addedgelabels=True,
                   showedgelabels=True, edgeid=True, edgeocc=False,
                   edgeprob=True, weightedges=True, writedot=True,
                   rmprev=False, transdist=18):
    """
    For each nodelabel, create a new CausalGraph for each path that can lead
    to a node with that label in the cores. Then give the probability of each
    exit mesh from that label on each new CausalGraph.
    If nodelabels is None, do for each node labels found in template pathway.
    """

    # Reading cores.
    if causalgraphs == None:
        core_files = get_dot_files(eoi, "meshedcore")
        meshedcores = []
        for core_file in core_files:
            core_path = "{}/{}".format(eoi, core_file)
            meshedcores.append(CausalGraph(core_path, eoi))
    else:
        meshedcores = causalgraphs
        core_files = None
    flush_ignored(meshedcores, core_files, ignorelist)
    # Reading event pathway template.
    if template == None:
        template_path = "{}/eventpathway.dot".format(eoi)
    else:
        template_path = template
    # Check node labels to highlight.
    if nodelabels == None:
        nodelabels = []
        template_graph = CausalGraph(template_path, eoi)
        for eventnode in template_graph.eventnodes:
            nodelabels.append(eventnode.label)
    elif not isinstance(nodelabels, list):
        nodelabels = [nodelabels]
    for nodelabel in nodelabels:
        graphs_path = "{}/{}".format(eoi, nodelabel)
        if not os.path.exists(graphs_path):
            os.mkdir(graphs_path)
    # Loop on node labels.
    for nodelabel in nodelabels:
        # Extract path to nodes carrying label.
        meshedcorescopy = copy.deepcopy(meshedcores)
        extractedpaths = extractpaths(eoi, nodelabel, meshedcorescopy)
        #dir_path = "{}/{}".format(eoi, nodelabel)
        #for extractedpath in extractedpaths:
        #    extractedpath.build_dot_file(edgeprob=False)
        #    output_path = "{}/{}".format(dir_path, extractedpath.filename)
        #    outfile = open(output_path, "w")
        #    outfile.write(extractedpath.dot_file)
        #    outfile.close()

        # Merge the paths keeping track of the filenames in prevcores.
        # mergecores destroys the causalgraphs used, so use the copy.
        mergedpaths = mergecores(eoi, causalgraphs=extractedpaths,
                                 showintro=True, color=False,
                                 writedot=False, msg=False)
        #for i in range(len(mergedpaths)):
        #    mergedpaths[i].filename = "path-{}.dot".format(i+1)
        #dir_path = "{}/{}".format(eoi, nodelabel)
        #for mergedpath in mergedpaths:
        #    output_path = "{}/{}".format(dir_path, mergedpath.filename)
        #    outfile = open(output_path, "w")
        #    outfile.write(mergedpath.dot_file)
        #    outfile.close()

        # Map merged paths on mapped core template.
        mappedpaths = mapcores(eoi, causalgraphs=mergedpaths, template=None,
                               ignorelist=ignorelist, showintro=showintro,
                               addedgelabels=addedgelabels,
                               showedgelabels=showedgelabels, edgeid=edgeid,
                               edgeocc=edgeocc, edgeprob=edgeprob,
                               weightedges=weightedges, msg=False,
                               transdist=transdist,
                               writedot=False)
        # If writedot=False in the previous command, the file are written
        # in the eoi directory rather than in the nodelabel directory.

        # Add the outgoing edges from the highlighted node.
        #for mappedpath in [mappedpaths[5]]:
        for mappedpath in mappedpaths:
            for eventnode in mappedpath.eventnodes:
                if eventnode.label == nodelabel:
                    eventnode.highlighted = True
            # Get event nodes and meshes from prevcores.
            midid = mappedpath.find_max_midid(cover=True)+1
            meshid = mappedpath.find_max_meshid(cover=True)+1
            for pathmesh in mappedpath.meshes:
                pathmesh.skip = False
                pathmesh.outgoing = False
            for cpathmesh in mappedpath.covermeshes:
                cpathmesh.skip = False
                cpathmesh.outgoing = False
            for prevcore in mappedpath.prevcores:
                underscore = prevcore.index("_")
                originalfile = "{}/{}.dot".format(eoi, prevcore[:underscore])
                highnodeid = prevcore[underscore+1:]
                origaph = None
                for meshedcore in meshedcores:
                    if meshedcore.filename == originalfile:
                        origaph = meshedcore
                        break
                highnode = None
                for eventnode in origaph.eventnodes:
                    if eventnode.nodeid == highnodeid:
                        highnode = eventnode
                        break
                # Compute cover meshes if they are not present in meshed core.
                if len(origaph.covermeshes) == 0:
                    origaph.build_nointro()
                # Add missing event nodes and gather meshes to add.
                meshes_to_add = []
                for mesh in origaph.meshes:
                    sources, targets = mesh.get_events()
                    if highnode in sources:
                        meshes_to_add.append(mesh)
                covermeshes_to_add = []
                for cmesh in origaph.covermeshes:
                    sources, targets = cmesh.get_events()
                    if highnode in sources:
                        covermeshes_to_add.append(cmesh)
                midid, meshid = draw_outgoing(meshes_to_add,
                                              mappedpath.meshes,
                                              midid, meshid, transdist)
                midid, meshid = draw_outgoing(covermeshes_to_add,
                                              mappedpath.covermeshes,
                                              midid, meshid, transdist)
            # Compute outgoing probability (sum should be equals to 1).
            # Also compute average to assign edge width.
            totout = 0
            out_uses = []
            for pathmesh in mappedpath.meshes:
                if pathmesh.outgoing == True:
                    if pathmesh.underlying == False:
                        totout += pathmesh.uses
                        out_uses.append(pathmesh.uses)
            for cpathmesh in mappedpath.covermeshes:
                if cpathmesh.outgoing == True:
                    totout += cpathmesh.uses
                    out_uses.append(cpathmesh.uses)
            aveout = statistics.mean(out_uses)
            # Set edge width and labels for outgoing edges.
            for pathmesh in mappedpath.meshes:
                if pathmesh.outgoing == True:
                    outgoing_mesh(pathmesh, totout, aveout)
            for cpathmesh in mappedpath.covermeshes:
                if cpathmesh.outgoing == True:
                    outgoing_mesh(cpathmesh, totout, aveout)
        for i in range(len(mappedpaths)):
            mappedpaths[i].filename = "highlight-{}.dot".format(i+1)
        dir_path = "{}/{}".format(eoi, nodelabel)
        #for mappedpath in [mappedpaths[5]]:
        for mappedpath in mappedpaths:
            mappedpath.build_dot_file(showintro, addedgelabels, showedgelabels,
                                      edgeid, edgeocc, edgeprob, weightedges)
            output_path = "{}/{}".format(dir_path, mappedpath.filename)
            outfile = open(output_path, "w")
            outfile.write(mappedpath.dot_file)
            outfile.close()


def outgoing_mesh(mesh, totout, aveout):
    """ Set edge width and labels of outgoing meshes. """

    mesh.usage = mesh.uses/totout
    minpenwidth = 1
    medpenwidth = 3
    maxpenwidth = 6.5
    ratio = mesh.uses/aveout
    pensize = math.log(ratio, 2) + medpenwidth
    if pensize < minpenwidth:
        pensize = minpenwidth
    if pensize > maxpenwidth:
        pensize = maxpenwidth
    midnodesize = math.sqrt(pensize)/12
    for midnode in mesh.midnodes:
        midnode.overridewidth = midnodesize
    for midedge in mesh.midedges:
        midedge.overridewidth = pensize
    mesh.assign_label_carrier()
    for midedge in mesh.midedges:
        if midedge.labelcarrier == True:
            midedge.overridelabel = " {:.3f}".format(mesh.usage) 


def draw_outgoing(meshes_to_add, mappedmeshes, midid, meshid, transdist):
    """ Add the outgoing meshes of highlighted node."""

    for mesh_to_add in meshes_to_add:
        analog_meshes = []
        for pathmesh in mappedmeshes:
            if pathmesh.skip == False:
                if analogous_meshes(mesh_to_add, pathmesh,
                                    enforcerank=False):
                    analog_meshes.append(pathmesh)
        if len(analog_meshes) < 1:
            raise ValueError("Mesh not found.")
        elif len(analog_meshes) == 1:
            if analog_meshes[0].color == "grey80":
                paint_mesh_black(analog_meshes[0])
                analog_meshes[0].uses = mesh_to_add.uses
                analog_meshes[0].outgoing = True
            elif analog_meshes[0].color == "black":
                analog_meshes[0].uses += mesh_to_add.uses
            else: # Mesh already used to represent path.
                # Add a copy of that mesh with skip = False
                # and with occurrence = mesh_to_add.occurrence.
                # The copy is translated by transdist.
                dupl_mesh = duplicate_mesh(analog_meshes[0], midid, meshid)
                dupl_mesh.skip = False
                dupl_mesh.outgoing = True
                midid += len(analog_meshes[0].midnodes)
                meshid += 1
                translate_mesh(dupl_mesh, 1, transdist)
                paint_mesh_black(dupl_mesh)
                mappedmeshes.append(dupl_mesh)
                analog_meshes[0].skip = True
        elif len(analog_meshes) > 1:
            # Add a copy of the last of those meshes with skip = False
            # and with uses = mesh_to_add.uses.
            # The copy is translated by transdist.
            dupl_mesh = duplicate_mesh(analog_meshes[-1], midid, meshid)
            dupl_mesh.skip = False
            dupl_mesh.outgoing = True
            midid += len(analog_meshes[-1].midnodes)
            meshid += 1
            translate_mesh(dupl_mesh, 1, transdist)
            paint_mesh_black(dupl_mesh)
            mappedmeshes.append(dupl_mesh)
            for analog_mesh in analog_meshes:
                analog_mesh.skip = True

    return midid, meshid


def paint_mesh_black(mesh):
    """ Set all colors to black in a mesh. """

    mesh.color = "black"
    for midedge in mesh.midedges:
        midedge.color = "black"
    for midnode in mesh.midnodes:
        midnode.bordercolor = "black"
        if midnode.midtype == "enabling":
            midnode.fillcolor = "black"


def extractpaths(eoi, nodelabel, meshedcores):
    """
    For the given node label, create a new CausalGraph for each path
    that can lead to a node with that label in the meshedcores.
    """

    extracted_paths = []
    for meshedcore in meshedcores:
        event_instances = []
        for eventnode in meshedcore.eventnodes:
            if eventnode.label == nodelabel:
                event_instances.append(eventnode)
        for event_instance in event_instances:
            extracted_path = CausalGraph(eoi=eoi, meshedgraph=True)
            extracted_path.occurrence = meshedcore.occurrence
            slash = meshedcore.filename.rfind("/")
            fname = "{}_{}.dot".format(meshedcore.filename[slash+1:-4],
                                       event_instance.nodeid)
            extracted_path.filename = fname
            extracted_path.eventnodes = [event_instance]
            extracted_path.meshes = []
            # Reconstruct the graph upstream of event_instance.
            current_nodes = [event_instance]
            current_meshes = []
            for mesh in meshedcore.meshes:
                sources, targets = mesh.get_events()
                for current_node in current_nodes:
                    if current_node in targets:
                        if mesh not in current_meshes:
                            current_meshes.append(mesh)
            for current_mesh in current_meshes:
                if current_mesh not in extracted_path.meshes:
                    extracted_path.meshes.append(current_mesh)
            while len(current_meshes) > 0:
                current_nodes = []
                for current_mesh in current_meshes:
                    sources, targets = current_mesh.get_events()
                    for source in sources:
                        if source not in current_nodes:
                            current_nodes.append(source)
                current_meshes = []
                for mesh in meshedcore.meshes:
                    sources, targets = mesh.get_events()
                    for current_node in current_nodes:
                        if current_node in targets:
                            if mesh not in current_meshes:
                                current_meshes.append(mesh)
                for current_node in current_nodes:
                    if current_node not in extracted_path.eventnodes:
                        extracted_path.eventnodes.append(current_node)
                for current_mesh in current_meshes:
                    if current_mesh not in extracted_path.meshes:
                        extracted_path.meshes.append(current_mesh)
            for mesh in extracted_path.meshes:
                sources, targets = mesh.get_events()
                for target in targets:
                    if target not in extracted_path.eventnodes:
                        extracted_path.eventnodes.append(target)
            ## Highlight event instance without highlighting instances
            ## of the same evetn type within a given path.
            #instancecopy = EventNode(event_instance.nodeid,
            #    event_instance.label, event_instance.rank,
            #    event_instance.occurrence, event_instance.prob,
            #    event_instance.intro, event_instance.first,
            #    highlighted=True, pos=event_instance.pos)
            extracted_path.get_maxrank()
            extracted_paths.append(extracted_path)

    return extracted_paths


# ::::::::::::::::::: End of Node Highlight Section :::::::::::::::::::::::::::

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
                            occ = edge.occurrence
                            new_edges.append(CausalEdge(src, trgt,
                                                        occurrence=occ))

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
                occ = edge.occurrence
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
                                                    occurrence=occ))

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
    path_probs = []
    intro_nodes = []
    for node in graph.nodes:
        if node.intro == True:
            intro_nodes.append(node)
    for start_node in intro_nodes:
        all_paths.append([start_node])
        seen_agents.append([])
        path_probs.append(0)
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
                                path_probs[i] = edge.prob
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
                                  path_probs[i])
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
            new_edge = CausalEdge(new_node, modified_node, path_probs[i])
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
                                          up_edge.prob)
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
    mergedcores = mergecores(eoi, edgelabels=edgelabels, writedot=True,
                             rmprev=True)

    # 3) Loop causal cores to create event paths.
    loopedcores = loopcores(eoi, causalgraphs=mergedcores,
                            ignorelist=ignorelist, edgelabels=edgelabels,
                            writedot=True, rmprev=False, writepremerge=False)

    # 4) Merge all event paths into one event pathway.
    eventpath = foldcores(eoi, causalgraphs=loopedcores,
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
                ppos = line.index("prob=")
                bracket = line.index("]")
                prob = line[ppos+7:bracket]
                comma = line.rfind(",")
                new_line = line[:comma+2]
                new_line += 'label="  {}", '.format(prob)
                new_line += line[ppos:]
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
