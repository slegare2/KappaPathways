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
import collections
import itertools
import string


class EventNode(object):
    """
    An event node to use in causal graphs. It represents a specific event in
    causal cores and an event type (a rule) in pathways.
    """

    def __init__(self, nodeid, label, rank=None, weight=1, rel_wei=1.0,
                 occurrence=1, rel_occ=1.0, intro=False, first=False,
                 highlighted=False, pos=None, eventid=None, shrink=False,
                 incoming=[], outgoing=[], pdh=False):
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
        self.eventid = eventid
        self.shrink = shrink
        self.incoming = incoming
        self.outgoing = outgoing
        self.pdh = pdh
        self.check_types()
        self.edits = []


    def check_types(self):
        """ Check that EventNode attributes have proper types. """

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
        """ Representation of the EventNode object. """

        res =  "{Node "
        res += 'id: "{}",  label: "{}"'.format(self.nodeid, self.label)
        if self.rank != None:
            res += ",  rank: {}".format(self.rank)
        if self.occurrence != None:
            res += ",  occurrence: {}".format(self.occurrence)
        res += ",  intro: {}".format(self.intro)
        res += ",  first: {}}}".format(self.first)

        return res


class StateNode(object):
    """
    A state node to use in causal graphs. It represents the state that
    are changed by an event.
    """

    def __init__(self, nodeid, label, rank=None, weight=1, rel_wei=1.0,
                 occurrence=1, rel_occ=1.0, intro=False, first=False,
                 highlighted=False, pos=None, eventid=None, edit=None,
                 context=[], shrink=False, incoming=[], outgoing=[],
                 pdh=False, stdedit=None):
        """ Initialize class StateNode. """

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
        self.eventid = eventid
        self.edit = edit
        self.context = context
        self.shrink = shrink
        self.incoming = incoming
        self.outgoing = outgoing
        self.pdh = pdh
        self.stdedit = stdedit
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

        res =  "{Node "
        res += 'id: "{}",  label: "{}"'.format(self.nodeid, self.label)
        if self.rank != None:
            res += ",  rank: {}".format(self.rank)
        if self.occurrence != None:
            res += ",  occurrence: {}".format(self.occurrence)
        res += ",  intro: {}".format(self.intro)
        res += ",  first: {}}}".format(self.first)

        return res


class CausalEdge(object):
    """
    A relationship between event or state nodes in causal graphs. The relationship
    can be causal or conflict.
    """

    def __init__(self, source, target, weight=1, layout_weight=1, rel_wei=1.0,
                 occurrence=1, rel_occ=1.0, number=1, rel_num=1.0,
                 relationtype="causal",
                 color="black", secondary=False, underlying=False,
                 reverse=False, labelcarrier=True, indicator=False,
                 meshid=None, pos=None, labelpos=None, overridewidth=None,
                 overridelabel=None, essential=False):
        """ Initialize class CausalEdge. """

        self.source = source
        self.target = target
        self.weight = weight # Appears as w in dot file
        self.layout_weight = layout_weight # Appears as weight in dot file
        self.rel_wei = rel_wei # = weight / occurrence_of_EOI (num. of cores)
        self.occurrence = occurrence # Taken from the trace.
        self.rel_occ = rel_occ # = occurrence / occurrence_of_EOI
        self.number = number # unique instances in stories.
        self.rel_num = rel_num
        self.relationtype = relationtype
        self.color = color
        self.secondary = secondary
        self.underlying = underlying
        self.reverse = reverse
        self.labelcarrier = labelcarrier
        self.indicator = indicator
        self.meshid = meshid
        self.pos = pos
        self.labelpos = labelpos
        self.overridewidth = overridewidth
        self.overridelabel = overridelabel
        self.essential = essential
        self.check_types()


    def check_types(self):
        """ Check that CausalEdge attributes have proper types. """

        if not isinstance(self.source, EventNode):
            if not isinstance(self.source, StateNode):
                if not isinstance(self.source, MidNode):
                    raise TypeError("Source should be an EventNode, StateNode "
                                    "or MidNode.")
        if not isinstance(self.target, EventNode):
            if not isinstance(self.target, StateNode):
                if not isinstance(self.target, MidNode):
                    raise TypeError("Target should be an EventNode, StateNode "
                                    "or MidNode.")
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

    def __init__(self, edgelist, weight=1, layout_weight=1, rel_wei=1.0, 
                 occurrence=1, rel_occ=1.0, number=1, rel_num=1.0,
                 relationtype="causal", color="black", midcolor="black",
                 secondary=False, underlying=False, reverse=False,
                 labelcarrier=True, indicator=False, hyperid=None, pos=None,
                 labelpos=None, overridewidth=None, overridelabel=None,
                 essential=False, cover=False):
        """ Initialize class CausalEdge. """

        self.edgelist = edgelist
        self.weight = weight
        self.layout_weight = layout_weight
        self.rel_wei = rel_wei
        self.occurrence = occurrence
        self.rel_occ = rel_occ
        self.number = number
        self.rel_num = rel_num
        self.relationtype = relationtype
        self.color = color
        self.midcolor = midcolor
        self.secondary = secondary
        self.underlying = underlying
        self.reverse = reverse
        self.labelcarrier = labelcarrier
        self.indicator = indicator
        self.hyperid = hyperid
        self.pos = pos
        self.labelpos = labelpos
        self.overridewidth = overridewidth
        self.overridelabel = overridelabel
        self.essential = essential
        self.cover = cover
        self.check_types()
        self.update()


    def update(self):
        """ Check that all edges within the hyperedge have the same target. """

        self.target = self.edgelist[0].target
        self.sources = []
        all_weights = []
        all_numbers = []
        all_conflicts = True
        for subedge in self.edgelist:
            if subedge.target != self.target:
                raise ValueError("Hyperedge has more than one target.")
            self.sources.append(subedge.source)
            all_weights.append(subedge.weight)
            all_numbers.append(subedge.number)
            if subedge.relationtype != "conflict":
                all_conflicts = False
        self.weight = min(all_weights)
        self.number = min(all_numbers)
        if all_conflicts == True:
            self.relationtype = "conflict"
        # I do not enforce equal weight among all subedges because I want
        # to allow zero weight on edges with intro source for better layout.
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


class MidNode(object):
    """ Intermediary node to build hyperedges. """
    
    def __init__(self, nodeid):
        """ Initialize class MidNode. """

        self.nodeid = nodeid
        self.check_types()


    def check_types(self):
        """ Check that MidNode attributes have proper types. """

        if not isinstance(self.nodeid, str):
            raise TypeError("nodeid should be a string.")


    def __repr__(self):
        """ Representation of the MidNode object. """

        res =  "MidNode "
        res += 'id: "{}", '.format(self.nodeid)

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


    #def assign_label_carrier(self):
    #    """ Choose which hyperedge and subedge will carry labels. """
    #
    #    for hyperedge in self.hyperedges+self.coverhyperedges: 
    #        print(hyperedge.labelcarrier)


#    def assign_label_carrier_for_mesh(self):
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

    def __init__(self, filename=None, eoi=None, hypergraph=False,
                 producedby="KappaPathways", showintro=True,
                 precedenceonly=False, rankposdict=None):
        """ Initialize class CausalGraph. """

        # Header variables.
        self.filename = filename
        self.eoi = eoi
        self.hypergraph = hypergraph
        self.producedby = producedby
        self.showintro = showintro
        self.precedenceonly = precedenceonly
        self.rankposdict = rankposdict
        # Main variables.
        self.eventnodes = []
        self.statenodes = []
        self.midnodes = []
        self.causaledges = []
        self.hyperedges = []
        #self.midnodes = []
        #self.midedges = []
        #self.meshes = []
        #self.edgegroups = []
        #self.midnodegroups = []
        # Cover edges and midnodes computed from nointro.
        #self.coveredges = []
        #self.covermidnodes = []
        #self.covermidedges = []
        #self.covermeshes = []
        #self.covermidnodegroups = []
        self.coverhyperedges = []
        # Post-computed variables.
        self.occurrence = 1
        self.maxrank = None
        self.minrank = None
        self.prevcores = None
        if self.filename != None:
            self.read_dot(self.filename)


    #def read_dualstory(self, dotpath):
    #    """ Read nodes and edges from input dual story. """


    def read_dot(self, dotpath):
        """ Read nodes and edges from input causal graph. """

        rank = None
        self.label_mapping = {}
        dotfile = open(dotpath, "r").readlines()
        for line in dotfile:
            if 'precedenceonly="True"' in line:
                self.precedenceonly = True
            if 'precedenceonly="False"' in line:
                self.precedenceonly = False
            if 'hypergraph="True"' in line:
                self.hypergraph = True
            #if "nodestype=" in line:
            #    type_index = line.index("nodestype")
            #    quote = line.rfind('"')
            #    self.nodestype = line[type_index+11:quote]
            if 'producedby=' in line:
                prod_index = line.index("producedby")
                quote = line.rfind('"')
                self.producedby = line[prod_index+12:quote]
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
                maxr_str = line[maxrank_index+9:quote]
                if "." in maxr_str:
                    self.maxrank = float(line[maxrank_index+9:quote])
                else:
                    self.maxrank = int(line[maxrank_index+9:quote])
            if "minrank=" in line:
                minrank_index = line.index("minrank")
                quote = line.rfind('"')
                minr_str = line[minrank_index+9:quote]
                if "." in minr_str:
                    self.minrank = float(line[minrank_index+9:quote])
                else:
                    self.minrank = int(line[minrank_index+9:quote])
            if "rank = same" in line:
                open_quote = line.index('"')
                close_quote = line[open_quote+1:].index('"')+open_quote+1
                rank_str = line[open_quote+1:close_quote]
                if "." in rank_str:
                    rank = float(rank_str)
                else:
                    rank = int(rank_str)
                #medrank = float(line[open_quote+1:close_quote])
                #rank = int(medrank)
            if line[0] == "}":
                rank = None
            # Read nodes.
            read_it = False
            if "label=" in line and "Occurrence" not in line:
                if "->" not in line and "rank = same" not in line:
                    if "cover=True" not in line:
                        read_it = True
            if read_it == True:
                if line[0:2] == "//":
                   read_line = line[2:]
                else:
                   read_line = line
                tokens = read_line.split()
                ori_id = tokens[0]
                if '"' in ori_id:
                    ori_id = ori_id[1:-1]
                if any(s in ori_id for s in ["ev", "state", "mid"]):
                    node_id = ori_id
                else:
                    node_id = "ev{}".format(ori_id)
                lbl_start = read_line.index("label=")+7
                stded_start = -1
                if "stded=" in read_line:
                    stded_start = read_line.index("stded=")
                shrk = False
                if "hlabel=" in read_line:
                    lbl_start = read_line.index("hlabel=")+8
                    shrk = True
                if ">" in read_line:
                    lbl_end = (read_line[lbl_start:stded_start]
                                  .rfind('>')+lbl_start)
                else:
                    lbl_end = (read_line[lbl_start:stded_start]
                                  .index('"')+lbl_start)
                label_str = read_line[lbl_start:lbl_end].strip()
                label = label_str.replace("\\n ", "")
                label = label.replace("<br/>", " ")
                if "intro=True" in read_line:
                    is_intro = True
                else:
                    is_intro = False
                if "first=True" in read_line:
                    is_first = True
                else:
                    is_first = False
                stdedit = get_stded(read_line)
                #if "midtype" in read_line:
                #    mid_start = read_line.index('midtype')+8
                #    mid_end = read_line[mid_start:].index(',')+mid_start
                #    midtype = read_line[mid_start:mid_end]
                #    if "style=dotted" in read_line:
                #        ghost = True
                #    else:
                #        ghost = False
                #    fillcolor = get_field("fillcolor=", read_line, "black")
                #    bordercolor = get_field(" color=", read_line, "black")
                #    new_midnode = MidNode(ori_id, rank, midtype,
                #                          ghost=ghost, fillcolor=fillcolor,
                #                          bordercolor=bordercolor)
                #    if 'cover=True' not in line:
                #        self.midnodes.append(new_midnode)
                #    elif 'cover=True' in line:
                #        self.covermidnodes.append(new_midnode)
                if "ev" in node_id:
                    eventid = node_id[2:]
                    self.eventnodes.append(EventNode(node_id, label,
                                                     rank,
                                                     intro=is_intro,
                                                     first=is_first,
                                                     shrink=shrk,
                                                     eventid=eventid))
                    self.label_mapping[node_id] = label

                elif "state" in node_id:
                    eventid = get_field("ev=", read_line, None)
                    self.statenodes.append(StateNode(node_id, label,
                                                     rank,
                                                     intro=is_intro,
                                                     first=is_first,
                                                     eventid=eventid,
                                                     stdedit=stdedit))
                    self.label_mapping[node_id] = label
                elif "mid" in node_id:
                    self.midnodes.append(MidNode(node_id))
        # Read edges.
        tmp_edges = []
        #tmp_midedges = []
        #tmp_cedges = []
        #tmp_cmidedges = []
        for line in dotfile:
            read_it = False
            if "->" in line and '[style="invis"]' not in line:
                if "cover=True" not in line:
                    read_it = True
            if read_it == True:
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
                if "state" not in source_id and "mid" not in source_id:
                    if "ev" not in source_id:
                        source_id = "ev{}".format(source_id)
                target_id = tokens[2]
                if '"' in target_id:
                    target_id = target_id[1:-1]
                if "state" not in target_id and "mid" not in target_id:
                    if "ev" not in target_id:
                        target_id = "ev{}".format(target_id)
                source = None
                target = None
                for eventnode in self.eventnodes:
                    if eventnode.nodeid == source_id:
                        source = eventnode
                    if eventnode.nodeid == target_id:
                        target = eventnode
                for statenode in self.statenodes:
                    if statenode.nodeid == source_id:
                        source = statenode
                    if statenode.nodeid == target_id:
                        target = statenode
                for node in self.midnodes:
                    if node.nodeid == source_id:
                        source = node
                    if node.nodeid == target_id:
                        target = node
                #for node in self.covermidnodes:
                #    if node.nodeid == source_id:
                #        source = node
                #    if node.nodeid == target_id:
                #        target = node
                #meshid = get_field("meshid=", read_line, 1)
                #meshid = int(meshid)
                weight = get_field("w=", read_line, 1)
                weight = int(weight)
                layout_weight = get_field("weight=", read_line, 1)
                layout_weight = int(layout_weight)
                color = get_field("color=", read_line, "black")
                if "label=" in line:
                    labelcarrier = True
                else:
                    labelcarrier = False
                if self.precedenceonly == False:
                    if self.producedby == "KaFlow":
                        if color == "grey":
                            edgetype = "conflict"
                            color = "black"
                        else:
                            edgetype = "causal"
                    elif self.producedby == "KaStor":
                        if "style=dotted" in line:
                            edgetype = "conflict"
                            source_save = source
                            source = target
                            target = source_save
                        else:
                            edgetype = "causal"
                    elif self.producedby == "KappaPathways":
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
                ess = False
                if "ess=True" in line:
                    ess = True
                new_edge = CausalEdge(source, target, weight=weight,
                                      number=weight, relationtype=edgetype,
                                      underlying=underlying,
                                      color=color, essential=ess)
                tmp_edges.append(new_edge)
                #elif 'cover=True' in line:
                #    tmp_cedges.append(new_edge)
        for edge in tmp_edges:
            self.causaledges.insert(0, edge)
        #for midedge in tmp_midedges:
        #    self.midedges.insert(0, midedge)
        #for cedge in tmp_cedges:
        #    self.coveredges.insert(0, cedge)
        #for cmidedge in tmp_cmidedges:
        #    self.covermidedges.insert(0, cmidedge)
        self.postprocess()


    def postprocess(self):
        """
        Various stuff to do after reading a dot file. Includes the creation
        of meshes from the intermediary nodes and edges.
        """

        for node in self.eventnodes:
            if "Intro" in node.label:
                node.intro = True
                node.label = node.label[6:]
                if node.label == "Lig, Lig":
                    node.label = "Lig"
        if self.hypergraph == False:
            self.create_hyperedges()
            self.read_states_from_file()
        elif self.hypergraph == True:
            self.read_hyperedges()
        if self.producedby != "KappaPathways":
            self.find_first_rules()
            self.rank_sequentially()
            self.rm_superfluous_causal_edges()
            self.align_vertical()
        if self.eoi == None:
            self.get_maxrank()
            for node in self.nodes:
                if node.rank == self.maxrank:
                    self.eoi = node.label


    def read_states_from_file(self):
        """ Read node states from separate json file. """

        if len(self.statenodes) > 0:
            dash = self.filename.rfind("-")
            period = self.filename.rfind(".")
            slash = self.filename.rfind("/")
            num = self.filename[dash+1:period]
            prefix = self.filename[:slash]
            json_path = "{}/statefile-{}.json".format(prefix, num)
            statefile = open(json_path, "r")
            statedict = json.load(statefile)
            for statenode in self.statenodes:
                statenode.state = statedict["states"][statenode.nodeid]
                statenode.edit = statedict["edits"][statenode.nodeid]


    def create_hyperedges(self):
        """ Create hyperedges by grouping edges with the same target. """

        self.hyperedges = []
        for edge in self.causaledges:
            edge_found = False
            for hyperedge in self.hyperedges:
                if edge.target == hyperedge.target:
                    hyperedge.addedge(edge)
                    edge_found = True
            if edge_found == False:
                self.hyperedges.append(HyperEdge([edge]))


    def read_hyperedges(self):
        """
        Create hyperedges by grouping edges pointing to a same midnode.
        """

        self.hyperedges = []
        hyperdict = {}
        targetdict = {}
        # Find the target of each midnode.
        for edge in self.causaledges:
            if isinstance(edge.source, MidNode):
                for midnode in self.midnodes:
                    if edge.source == midnode:
                        targetdict[midnode.nodeid] = edge.target
        # Create a hyperedge with a single source for causal edges that
        # directly link events and states without passing thougth midnodes.
        for edge in self.causaledges:
            if not isinstance(edge.source, MidNode):
                if not isinstance(edge.target, MidNode):
                    self.hyperedges.append(HyperEdge([edge]))
        # Create a hyperedge with many sources for each midnode.
        for edge in self.causaledges:
            if isinstance(edge.target, MidNode):
                for midnode in self.midnodes:
                    if edge.target == midnode:
                        new_target = targetdict[midnode.nodeid]
                        edge.target = new_target
                        if midnode.nodeid not in hyperdict.keys():
                            hyperdict[midnode.nodeid] = HyperEdge([edge])
                        else:
                            hyperdict[midnode.nodeid].addedge(edge)
        for midid in hyperdict.keys():
            self.hyperedges.append(hyperdict[midid])

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
        # I now realize that this will not work if I have some already active
        # proteins at the beginning of the simulation (Example some active
        # ABL1 in the pYnet model). This will always be the case if I restart
        # a simulation or start it using protein abundance data.
        # !!! Need to find a way to keep the ranking in these cases !!!

        # Initialize ranks.
        current_nodes = []
        for node in self.eventnodes+self.statenodes:
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
                for hyperedge in self.hyperedges:
                    if hyperedge.target == candidate_node:
                        incoming_hedges.append(hyperedge)
                secured_hedges = []
                potential_hedges = []
                possible_ranks = []
                for incoming_hedge in incoming_hedges:
                    secured = True
                    for source in incoming_hedge.sources:
                        if source.intro == False:
                            if source.rank == None:
                                secured = False
                                break
                    if secured == True:
                        secured_hedges.append(incoming_hedge)
                        potential_hedges.append(incoming_hedge)
                        source_ranks = []
                        all_intro = True
                        for source in incoming_hedge.sources:
                            if source.intro == False:
                                all_intro = False
                                shrk = False
                                if isinstance(source, EventNode):
                                    if source.shrink == True:
                                        shrk = True
                                if shrk == False:
                                    source_ranks.append(source.rank)
                                else:
                                    # If node is shrunk, take the max rank
                                    # among its sources instead of own rank.
                                    subranks = []
                                    for hyperedge2 in self.hyperedges:
                                        if hyperedge2.target == source:
                                            for subsource in hyperedge2.sources:
                                                subranks.append(subsource.rank)
                                    source_ranks.append(max(subranks))
                        #if all_intro == True:
                        #    source_ranks.append(0)
                        possible_ranks.append(max(source_ranks)+1)
                    elif secured == False:
                        # Hyperedges that are not secured still count as
                        # potential hyperedges if they do not loop back to
                        # the candidate node.
                        if rulepos == "bot":
                            # I need to compute this only of rulepos == "bot"
                            # and it is very long to compute if the graph is
                            # large.
                            looping_hedge = False
                            for source in incoming_hedge.sources:
                                paths = self.follow_hyperedges("up", source,
                                                               [candidate_node])
                                if len(paths) > 0:
                                    looping_hedge = True
                            if looping_hedge == False:
                                potential_hedges.append(incoming_hedge)
                nsecured = len(secured_hedges)
                npotential = len(potential_hedges)
                if rulepos == "top" and nsecured > 0:
                    candidate_node.rank = min(possible_ranks)
                    current_nodes.append(candidate_node)
                if rulepos == "bot" and nsecured == npotential:
                    candidate_node.rank = max(possible_ranks)
                    #if 1 in possible_ranks:
                    #    candidate_node.rank = 1
                    #else:
                    #    candidate_node.rank = max(possible_ranks)
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
        # Rank intro nodes at 0 if top is selected.
        for node in self.eventnodes:
            if node.intro == True:
                node.rank = 0
        # Rank intro nodes at the lowest rank among its targets, minus 1.
        if intropos == "bot":
            for node in self.eventnodes:
                if node.intro == True:
                    target_ranks = []
                    for hyperedge in self.hyperedges:
                        if node in hyperedge.sources:
                            if hyperedge.target.shrink == False:
                                target_ranks.append(hyperedge.target.rank)
                            else:
                                for h2 in self.hyperedges:
                                    if hyperedge.target in h2.sources:
                                        target_ranks.append(h2.target.rank)
                    node.rank = min(target_ranks) - 1
        # Optionally, push targets of intro nodes down when possible.
        # This way of doing it does not work very well. Takes a long time on
        # large graphs.
        if intropos == "bot2":
            # Root nodes are first nodes with the longest path to the EOI.
            for eventnode in self.eventnodes:
                if eventnode.label == self.eoi:
                    eoi_node = eventnode
            path_lengths = []
            for node in self.eventnodes+self.statenodes:
                if node.first == True:
                    paths = self.follow_hyperedges("down", node, [eoi_node])
                    for path in paths:
                        path_lengths.append(len(path))
            longest_path = max(path_lengths)
            root_nodes = []
            for node in self.eventnodes+self.statenodes:
                if node.first == True:
                    paths = self.follow_hyperedges("down", node, [eoi_node])
                    path_lengths = []
                    for path in paths:
                        path_lengths.append(len(path))
                    if max(path_lengths) == longest_path:
                        root_nodes.append(node)

            # Fix in place any node that has a path up to a root node.
            fixed_nodes = []
            for node in self.eventnodes+self.statenodes:
               paths = self.follow_hyperedges("up", node, root_nodes)
               if len(paths) > 0:
                   fixed_nodes.append(node)
            # Move down nodes that are not fixed when possible.
            gap_found = True
            while gap_found == True:
                gap_found = False
                for node in self.eventnodes+self.statenodes:
                    if node not in fixed_nodes:
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
                    paths = self.follow_hyperedges("down", edge.source,
                                                   [edge.target])
                    if len(paths) == 1:
                        new_edgelist.append(edge)
                if len(new_edgelist) > 0:
                    new_hyperedges.append(HyperEdge(new_edgelist))

        self.hyperedges = new_hyperedges


    def align_vertical(self):
        """
        Adjust edge layout weights such that the event nodes are aligned vertically.
        """

        for hyperedge in self.hyperedges:
            hyperedge.layout_weight = hyperedge.weight
            for edge in hyperedge.edgelist:
                #edge.layout_weight = hyperedge.layout_weight
                edge.layout_weight = edge.weight
            if len(hyperedge.edgelist) > 1:
                nonintro_present = False
                for edge in hyperedge.edgelist:
                    if edge.source.intro == False:
                        nonintro_present = True
                        break
                if nonintro_present == True:
                    for edge in hyperedge.edgelist:
                        if edge.source.intro == True:
                            edge.layout_weight = 0


    def get_all_reachables(self):
        """
        Assign a list of all reachable nodes to each rule_output node. This
        algorithm seeks to optimize computation time by starting at the bottom
        of the graph. While gathering the reachable nodes of a current
        rule_output node, if a previous rule_output node is found, the
        reachables of that previous rule_output are added to the current
        rule_output's reachables instead of continuing the search.

        The order of the nodes within the lists of reachable nodes cannot be
        guaranteed.

        This algorithm does not check for loops, use only on stories, not
        quotiented pathways.
        """

        # Initialize lists of reachable nodes.
        for node in self.statenodes + self.eventnodes:
            node.reachable = []
        # Find the last rule_outputs (the ones pointing directly to the EOI).
        outputs_fringe = []
        for edge in self.eoi_node.incoming:
            outputs_fringe.append(edge.source)
        # Read graph upstream.
        seen_nodes = []
        while len(outputs_fringe) > 0:
            # Find the reachable nodes of each output_fringe node here,
            # checking previous node.reachable.
            for output_node in outputs_fringe:
                self.get_reachables(output_node)
                #if output_node not in seen_nodes:
                #    seen_nodes.append(output_node)
            # Go up once.
            up_next = []
            for output_node in outputs_fringe:
                if len(output_node.incoming) > 0:
                    for edge in output_node.incoming:
                        up_next.append(edge.source)
            outputs_fringe = up_next
            # Keep going up until all fringe nodes are rule_outputs
            # (or have no incoming edge).
            outputs_reached = False
            while outputs_reached == False:
                outputs_reached = True
                up_next = []
                for up_node in outputs_fringe:
                    # Keep up_node if it is a rule.
                    if isinstance(up_node, EventNode) and up_node.intro == False:
                        up_next.append(up_node)
                    # or if it is in rule_outputs.
                    elif up_node in self.rule_outputs:
                        up_next.append(up_node)
                    elif len(up_node.incoming) > 0:
                        for edge in up_node.incoming:
                            if edge.source not in seen_nodes:
                                up_next.append(edge.source)
                                seen_nodes.append(edge.source)
                outputs_fringe = up_next
                # Check if all fringe nodes are rule_outputs or rule.
                for up_node in outputs_fringe:
                    is_rule = False
                    if isinstance(up_node, EventNode):
                        if up_node.intro == False:
                            is_rule = True
                    if is_rule == False and up_node not in self.rule_outputs:
                        outputs_reached = False
                        break


    def get_reachables(self, from_node):
        """
        Return all downstream reachable nodes from given node. If a node is
        found with already defined reachable list, stop search and add its
        list of reachables to current node's reachable.
        """

        # Initialize fringe nodes as the immediate targets of from_node.
        fringe = []
        for edge in from_node.outgoing:
            fringe.append(edge.target)
        list_of_reachables = []
        while len(fringe) > 0:
            # Add fringe nodes to from_node reachables.
            for node in fringe:
                if node not in list_of_reachables:
                    list_of_reachables.append(node)
            # Also add fringe node's own reachables if it has some.
            next_fringe = []
            for node in fringe:
                if len(node.reachable) > 0:
                    for rnode in node.reachable:
                        if rnode not in list_of_reachables:
                            list_of_reachables.append(rnode)
                # If the fringe node does not have reachables, put its
                # immediate target in the next fringe round.
                else:
                    for edge in node.outgoing:
                        if edge.target not in list_of_reachables:
                            next_fringe.append(edge.target)
            fringe = next_fringe
        from_node.reachable = list_of_reachables


    def reachability_with_block(self, from_node, to_nodes, block):
        """
        Tell if at least one of the to_nodes is reachable from the from_node
        without passing through the block node.
        The to_nodes should be the reachable nodes of the blocking node.
        """

        reachable = False
        # Initialize fringe nodes as the immediate targets of from_node.
        fringe = []
        for edge in from_node.outgoing:
            fringe.append(edge.target)
        list_of_reachables = []
        while len(fringe) > 0:
            # Check if one of the to_nodes is in the fringe.
            for node in fringe:
                if node in to_nodes:
                    reachable = True
            if reachable == True:
                break
            # Add fringe nodes to from_node's reachables.
            for node in fringe:
                if node not in list_of_reachables:
                    list_of_reachables.append(node)
            # Follow edges downstream to find next fringe round. Do not add
            # node if it is a state node that modifies a site found in
            # from_node's edit.
            next_fringe = []
            for node in fringe:
                if isinstance(node, StateNode):
                    is_mod = self.modifies_state(from_node.edit, node.edit)
                else:
                    is_mod = False
                if is_mod == False:
                    if node != block:
                        for edge in node.outgoing:
                            if edge.target not in list_of_reachables:
                                next_fringe.append(edge.target)
            fringe = next_fringe

        return reachable


    def modifies_state(self, state1, state2):
        """
        Check if state2 contains at least one site which is a modification of at
        least one site of state1.
        """

        is_modification = False
        for agent1 in state1:
            n1 = agent1["name"]
            id1 = agent1["id"]
            s1 = agent1["sites"]
            for agent2 in state2:
                n2 = agent2["name"]
                id2 = agent2["id"]
                if n1 == n2 and id1 == id2:
                    s2 = agent2["sites"]
                    for site1 in s1:
                        for site2 in s2:
                            if site1["name"] == site2["name"]:
                                is_modification = True
                                break
                        if is_modification == True:
                            break
                if is_modification == True:
                    break
            if is_modification == True:
                break

        return is_modification


#    def oldreachability_with_block(self, from_node, to_nodes, block):
#        """
#        Tell if at least one of the to_nodes is reachable from the from_node
#        without passing through the block node.
#        The to_nodes should be the reachable nodes of the blocking node.
#        """
#
#        reachable = False
#        # Initialize fringe nodes as the immediate targets of from_node.
#        fringe = []
#        for edge in from_node.outgoing:
#            fringe.append(edge.target)
#        list_of_reachables = []
#        while len(fringe) > 0:
#            # Check if one of the to_nodes is in the fringe.
#            for node in fringe:
#                if node in to_nodes:
#                    reachable = True
#            if reachable == True:
#                break
#            # Add fringe nodes to from_node reachables if they are not the
#            # blocking node.
#            for node in fringe:
#                if node != block:
#                    if node not in list_of_reachables:
#                        list_of_reachables.append(node)
#            # Follow edges downstream to find next fringe round if fringe node
#            # is not blocking.
#            next_fringe = []
#            for node in fringe:
#                if node != block:
#                    for edge in node.outgoing:
#                        if edge.target not in list_of_reachables:
#                            next_fringe.append(edge.target)
#            fringe = next_fringe
#
#        return reachable


#    def reachable(self, from_node, block_nodes=[])
#        """
#        Return all downstream reachable nodes from given node. Optionally
#        provide a list of blocking nodes. Reachability search stops if a
#        blocking node is met.
#
#        This is faster than finding all paths to the EOI in large graphs
#        because it does not produce a combinatorial explosion. However, the
#        order of the nodes within the list of reachable nodes cannot be
#        guaranteed.
#
#        This algorithm does not check for loops, use only on stories, not
#        quotiented pathways.
#        """


    def follow_edges(self, direction, from_node, to_nodes=[], block=None,
                     ignore_conflict=False, stop_at_first=False):
        """
        Return a list of all acyclic paths from a given node to the top of the
        graph (using direction="up") or to the bottom (using direction="down").
        If to_nodes are provided, return only the paths that go from from_node
        to any of the to_nodes.
        """
        # This method takes a lot of time on large graphs and can most probably
        # be improved to speed up calculation.
        
    
        all_paths = [[from_node]]
        ends_reached = False
        while ends_reached == False:
            ends_reached = True
            for i in range(len(all_paths)):
                path = all_paths[i]
                next_nodes = []
                for edge in self.causaledges:
                    skip = False
                    if ignore_conflict == True:
                        if edge.relationtype == "conflict":
                            skip = True
                    if skip == False:
                        if direction == "up":
                            if edge.target == path[-1]:
                                next_nodes.append(edge.source)
                        elif direction == "down":
                            if edge.source == path[-1]:
                                next_nodes.append(edge.target)
                if len(next_nodes) > 0 and path[-1] not in to_nodes:
                    ends_reached = False
                    if len(next_nodes) > 0:
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
            # Remove paths that end with blocking node if defined.
            if block != None:
                for i in range(len(all_paths)-1, -1, -1):
                    if all_paths[i][-1] == block:
                        del(all_paths[i])
            # Exit prematurely if only one path is sufficient.
            if stop_at_first == True:
                for path in all_paths:
                    if path[-1] in to_nodes:
                        ends_reached = True
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


    def follow_hyperedges(self, direction, from_node, to_nodes=[]):
        """
        Return a list of all acyclic paths from a given node to the top of the
        graph (using direction="up") or to the bottom (using direction="down").
        If to_nodes are provided, return only the paths that go from from_node
        to any of the to_nodes.
        """
        # Just as follow_edges, can probably be improved to speed up.
    
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
        for node in self.eventnodes + self.statenodes:
            if node.rank != None:
                all_ranks.append(node.rank)
        if len(all_ranks) > 0:
            self.maxrank = max(all_ranks)
            self.minrank = min(all_ranks)
        self.midranks = 1
        for r in all_ranks:
            if isinstance(r, float):
                if r.is_integer() == False:
                    self.midranks = 3
                    break


    def reverse_subedges(self):
        """
        Reverse the direction of subedges if their source has a higher rank
        than their target. Reverse hyperedge if all its subedges are reversed.
        """

        all_hyperedges = self.hyperedges + self.coverhyperedges
        # Reset reverse edge information.
        for hyperedge in all_hyperedges:
            hyperedge.reverse = False
            for subedge in hyperedge.edgelist:
                subedge.reverse = False
        # Compute new edge reversion.
        for hyperedge in all_hyperedges:
            # Check if source or target is shrunk.
            shrunk_src, shrunk_trg = False, False
            if isinstance(hyperedge.sources[0], EventNode):
                shrunk_src = hyperedge.sources[0].shrink
            if isinstance(hyperedge.target, EventNode):
                shrunk_trg = hyperedge.target.shrink
            # Normal case, no shrunk node as source or target.
            if shrunk_src == False and shrunk_trg == False:
                all_subedges_reversed = True
                for subedge in hyperedge.edgelist:
                   src_rank = subedge.source.rank
                   if src_rank != None and src_rank > subedge.target.rank:
                       subedge.reverse = True
                   else:
                       all_subedges_reversed = False
                if all_subedges_reversed == True:
                    hyperedge.reverse = True
                #all_sources_higher_rank = True
                #for subedge in hyperedge.edgelist:
                #    if subedge.source.rank < subedge.target.rank:
                #        all_sources_higher_rank = False
                #    else:
                #        subedge.reverse = True
                #if all_sources_higher_rank == True:
                #    hyperedge.reverse = True
            # If the source is shrunk, the hyperedge has only one subedge.
            # Reverse edge only if all sources are of higher rank.
            elif shrunk_src == True:
                src_list = []
                for hyperedge2 in all_hyperedges:
                    if hyperedge2.target == hyperedge.sources[0]:
                        src_list += hyperedge2.sources
                all_sources_higher_rank = True
                for source in src_list:
                    src_rank = source.rank
                    if src_rank != None and src_rank < hyperedge.target.rank:
                        all_sources_higher_rank = False
                        break
                if all_sources_higher_rank == True:
                    hyperedge.edgelist[0].reverse = True
                    hyperedge.reverse = True
            # Case where the target is shrunk.
            # Reverse subedge only if all its targets are of lower rank.
            # Reverse the hyperedge if all subedges are reversed.
            elif shrunk_trg == True:
                trg_list = []
                for hyperedge2 in all_hyperedges:
                    if hyperedge2.sources[0] == hyperedge.target:
                        trg_list.append(hyperedge2.target)
                all_subedges_reversed = True
                for subedge in hyperedge.edgelist:
                    all_targets_lower = True
                    src_rank = subedge.source.rank
                    for target in trg_list:
                        if src_rank != None and src_rank < target.rank:
                            all_targets_lower = False
                            break
                    if all_targets_lower == True:
                        subedge.reverse = True
                    else:
                        all_subedges_reversed = False
                if all_subedges_reversed == True:
                    hyperedge.reverse = True


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
        Create new hyperedges for the version of the graph where intro nodes
        are hidden.
        """

        # Reset adjacency lists to avoid infinite loop when doing deepcopy.
        for node in self.eventnodes + self.statenodes:
            node.incoming = []
            node.outgoing = []
            node.reachable = []
        for hyperedge in self.hyperedges:
            hyperedge.underlying = False
        self.coverhyperedges = []
        for hyperedge in self.hyperedges:
            underlying_list = []
            if hyperedge.underlying == False:
                has_intro, all_intro = self.check_intro(hyperedge)
                if has_intro == True:
                    if all_intro == True:
                        hyperedge.underlying = True
                    else:
                        a, b, c = self.find_underlying(hyperedge)
                        nointro_hedge = a
                        underlying_list = b
                        correspondances = c
                        #hedge_list.append(hyperedge)
                        #nointro = self.build_nointro_hyperedge(hyperedge)
                        #if len(nointro.edgelist) > 0:
                        #    # Check all other hyperedges that have the same
                        #    # subedges when ignoring subedges with intro nodes
                        #    # as source. They will all be grouped inside a
                        #    # single hyperedge without any intro source.
                        #    for hyperedge2 in self.hyperedges:
                        #        if hyperedge2.color == hyperedge.color:
                        #            if hyperedge2 != hyperedge:
                        #                if hyperedge2.underlying == False:
                        #                    nointro2 = self.build_nointro_hyperedge(hyperedge2)
                        #                    are_equi, corr = equivalent_hyperedges(nointro, nointro2, return_correspondances=True)
                        #                    if are_equi == True:
                        #                        hedge_list.append(hyperedge2)
                        #                        corrs.append(corr)
            if len(underlying_list) > 0:
                # Compute weight of nointro hyperedge as the sum the weight
                # of all its underlying hyperedges. Also mark meshes that were used as
                # underlying.
                weight_summ = 0
                for hedge in underlying_list:
                    weight_summ += hedge.weight
                    hedge.underlying = True
                nointro_hedge.weight = weight_summ
                if len(nointro_hedge.edgelist) > 0:
                    nointro_hedge.cover = True
                    self.coverhyperedges.append(nointro_hedge)


    def find_underlying(self, hyperedge):
        """
        Find the list of underlying hyperedges to form a cover hyperedge.
        """

        underlying_list = [hyperedge]
        nointro_hedge = self.build_nointro_hyperedge(hyperedge)
        correspondances = []
        if len(nointro_hedge.edgelist) > 0:
            # Check all other hyperedges that have the same subedges when
            # ignoring subedges with intro nodes as source. They will all
            # be grouped inside a single hyperedge without any intro source.
            for hyperedge2 in self.hyperedges:
                if hyperedge2.color == hyperedge.color:
                    if hyperedge2 != hyperedge:
                        if hyperedge2.underlying == False:
                            nointro_hedge2 = (self.
                                build_nointro_hyperedge(hyperedge2))
                            are_equi, corr = equivalent_hyperedges(nointro_hedge,
                                nointro_hedge2, return_correspondances=True)
                            if are_equi == True:
                                underlying_list.append(hyperedge2)
                                correspondances.append(corr)

        return nointro_hedge, underlying_list, correspondances


#    def build_nointro_for_meshes(self):
#        """
#        Create new meshes for the version of the graph that hides
#        intro nodes.
#        """
#
#        # Reset information about cover edges. 
#        for mesh in self.meshes:
#            mesh.underlying = False
#        self.covermeshes = []
#        nointro_groups = []
#        midid = self.find_max_midid()+1
#        for mesh1 in self.meshes:
#            mesh_list = []
#            if mesh1.underlying == False:
#                has_intro1, all_intro1 = self.check_intro(mesh1)
#                if has_intro1 == True:
#                    if all_intro1 == True:
#                        mesh1.underlying = True
#                    else:
#                        mesh_list.append(mesh1)
#                        noin1 = self.nointro_mesh(mesh1, midid)
#                        midid += len(noin1.midnodes)
#                        if len(noin1.midedges) > 0:
#                            # Check all other meshes that have the same
#                            # midedges when ignoring edges with intro nodes as
#                            # source. They will all be grouped inside a single
#                            # mesh without intro nodes as sources.
#                            for mesh2 in self.meshes:
#                                if mesh2.color == mesh1.color:
#                                    if mesh2 != mesh1:
#                                        if mesh2.underlying == False:
#                                            noin2 = self.nointro_mesh(mesh2,
#                                                                      midid)
#                                            if self.equivalent_meshes(noin1,
#                                                                      noin2):
#                                                mesh_list.append(mesh2)
#            if len(mesh_list) > 0:
#                # Compute occurence of nointro mesh as the sum of all its
#                # underlying meshes. Also mark meshes that were used as
#                # underlying.
#                uses_summ = 0
#                for mesh in mesh_list:
#                    uses_summ += mesh.uses
#                    mesh.underlying = True
#                noin1.uses = uses_summ
#                noin1.weight = uses_summ
#                if len(noin1.midedges) > 0:
#                    self.covermeshes.append(noin1)


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


#    def find_max_meshid(self, cover=False):
#        """ Find the highest mesh id in the graph. """
#
#        meshids = []
#        for mesh in self.meshes:
#            meshids.append(mesh.meshid)
#        if cover == True:
#            for mesh in self.covermeshes:
#                meshids.append(mesh.meshid)
#        max_meshid = max(meshids)
#
#        return max_meshid


    def check_intro(self, hyperedge):
        """ Check if a hyperedge has intro nodes in its sources. """

        has_intro = False
        all_intro = True
        for subedge in hyperedge.edgelist:
            if subedge.essential == False:
                #if isinstance(subedge.source, EventNode):
                if subedge.source.intro == True:
                    has_intro = True
                if subedge.source.intro == False:
                    all_intro = False

        return has_intro, all_intro


    #def check_intro_for_meshes(self, mesh):
    #    """ Check if a mesh has intro nodes in its sources. """

    #    has_intro = False
    #    all_intro = True
    #    for midedge in mesh.midedges:
    #        if isinstance(midedge.source, EventNode):
    #            if midedge.source.intro == True:
    #                has_intro = True
    #            if midedge.source.intro == False:
    #                all_intro = False

    #    return has_intro, all_intro


    def build_nointro_hyperedge(self, hyperedge):
        """ Build a new hyperedge without subedges from intro nodes. """

        new_hyperedge = copy.deepcopy(hyperedge)
        new_edgelist = []
        for subedge in hyperedge.edgelist:
            if isinstance(subedge.source, EventNode):
                if subedge.source.intro == False:
                    new_edgelist.append(subedge)
            else:
                new_edgelist.append(subedge)
        if len(new_edgelist) > 0:
            new_hyperedge.edgelist = new_edgelist
            new_hyperedge.update()

        return new_hyperedge
                

    def assign_label_carriers(self, s):
        """ Choose which hyperedge and subedge will carry labels. """

        for hyperedge in self.hyperedges+self.coverhyperedges:
            if (s == False and hyperedge.underlying == False) or s == True:
                if hyperedge.sources[0].shrink == False:
                    if len(hyperedge.edgelist) > 1:
                        for subedge in hyperedge.edgelist:
                            subedge.labelcarrier = False
                elif hyperedge.sources[0].shrink == True:
                    # Show labels of edges from shrunk nodes only if their
                    # weight is different from that of the upstream hyperedge.
                    for hyperedge2 in self.hyperedges+self.coverhyperedges:
                        if hyperedge2.target == hyperedge.sources[0]:
                            up_hedge = hyperedge2
                            break
                    if up_hedge.weight == hyperedge.weight:
                        hyperedge.edgelist[0].labelcarrier = False
                        hyperedge.labelcarrier = False


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


    def build_adjacency(self, hyper=False):
        """
        For each node, build the lists of incoming and outgoing hyperedges.
        """
        
        for node in self.eventnodes + self.statenodes:
            node.incoming = []
            node.outgoing = []
            if hyper == False:
                for edge in self.causaledges:
                    if node == edge.target:
                        node.incoming.append(edge)
                    if node == edge.source:
                        node.outgoing.append(edge)
            elif hyper == True:
                for hyperedge in self.hyperedges:
                    if node == hyperedge.target:
                        node.incoming.append(hyperedge)
                    if node in hyperedge.sources:
                        node.outgoing.append(hyperedge)


    def build_dot_file(self, showintro=True, addedgelabels=True,
                       showedgelabels=True, edgeid=True, edgeocc=False,
                       edgeuse=True, statstype="rel", weightedges=False,
                       edgewidthscale=1):
        """ build a dot file of the CausalGraph. """

        # Write info about graph.
        dot_str = 'digraph G{\n'
        dot_str += '  producedby="{}" ;\n'.format(self.producedby)
        dot_str += '  precedenceonly="{}" ;\n'.format(self.precedenceonly)
        dot_str += '  hypergraph="{}" ;\n'.format(self.hypergraph)
        #dot_str += '  nodestype="{}" ;\n'.format(self.nodestype)
        if self.eoi != None:
            dot_str += '  eoi="{}" ;\n'.format(self.eoi)
        if self.occurrence != None:
            dot_str += '  label="Occurrence = {}" '.format(self.occurrence)
            dot_str += 'fontsize=28 labelloc="t" ;\n'
        if self.maxrank != None:
            dot_str += '  maxrank="{}" ;\n'.format(self.maxrank)
        if self.minrank != None:
            dot_str += '  minrank="{}" ;\n'.format(self.minrank)
        if self.prevcores != None:
            dot_str += '  prevcores="{}" ;\n'.format(self.prevcores)
        #dot_str += '  ranksep=1 ;\n'
        #dot_str += '  nodesep=0.2 ;\n' # Default nodesep is 0.25
        dot_str += '  splines=true ;\n'
        dot_str += '  outputorder=nodesfirst ;\n'
        dot_str += '  node [pin=true] ;\n'
        #dot_str += '  edge [fontsize=18] ;\n'
        # Compute some statistics to assign edge and intermediary node width.
        minpenwidth = 1 * edgewidthscale
        medpenwidth = 3 * edgewidthscale
        maxpenwidth = 6.5 * edgewidthscale
        all_weights = []
        all_numbers = []
        for hyperedge in self.hyperedges:
            if hyperedge.underlying == False:
                all_weights.append(hyperedge.weight)
                all_numbers.append(hyperedge.number)
        #for coverhyper in self.coverhypers:
        #    all_uses.append(covermesh.uses)
        average_weight = statistics.mean(all_weights)
        average_number = statistics.mean(all_numbers)
        # Build drawing parameters dict.
        params = {"average_weight": average_weight,
                  "average_number": average_number,
                  "minpenwidth": minpenwidth,
                  "medpenwidth": medpenwidth,
                  "maxpenwidth": maxpenwidth,
                  "addedgelabels": addedgelabels,
                  "showedgelabels": showedgelabels,
                  "edgeid": edgeid,
                  "edgeocc": edgeocc,
                  "edgeuse": edgeuse,
                  "statstype": statstype,
                  "weightedges": weightedges,
                  "edgewidthscale": edgewidthscale}
        # Draw nodes.
        midranks = self.midranks
        for int_rank in range(int((self.minrank)*(midranks+1)),
                              int((self.maxrank+1)*(midranks+1))):
            current_rank = int_rank/(midranks+1)
            rank_str = "{}".format(current_rank)
            if showintro == False and current_rank < 1:
                dot_str += "//"
            if current_rank%1 == 0 and current_rank <= self.maxrank:
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
                if node.rank == current_rank and node.shrink == False:
                    #node_shape = 'invhouse'
                    node_shape = 'rectangle'
                    node_color = 'lightblue'
                    if node.intro == True:
                        node_shape = 'ellipse'
                        node_color = 'white'
                    if node.label == self.eoi:
                        node_shape = 'ellipse'
                        node_color = 'indianred2'
                    if showintro == False and node.intro == True:
                        dot_str += '//'
                    node_lines = textwrap.wrap(node.label, 20,
                                               break_long_words=False)
                    dot_str += '{} '.format(node.nodeid)
                    node_str = ""
                    for i in range(len(node_lines)):
                        if i == 0:
                            node_str += "{}".format(node_lines[i])
                        else:
                            node_str += "<br/>{}".format(node_lines[i])
                    # Add PDH information if not already present.
                    prefix_num = ""
                    if node.pdh != False and ":" not in node_str:
                        prefix_num = " : {}".format(node.pdh)
                    dot_str += ('[label=<{}{}>'
                                .format(node_str, prefix_num))
                    dot_str += ', shape={}, style="filled'.format(node_shape)
                    if node.pdh == False:
                        dot_str += '"'
                    else:
                        dot_str += ',dashed"'
                    #if node.highlighted == True:
                    #   dot_str += ', fillcolor=gold, penwidth=2'
                    #else:
                    dot_str += ', fillcolor={}'.format(node_color)
                    if node.highlighted == True:
                        dot_str += ', penwidth=4'
                    if node.intro == True:
                        dot_str += ', intro={}'.format(node.intro)
                    if node.first == True:
                        dot_str += ', first={}'.format(node.first)
                    if node.pos != None:
                        dot_str += ', pos={}'.format(node.pos)
                    #dot_str += ', penwidth=2'
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
                            node_str += "{}".format(node_lines[i])
                        else:
                            node_str += "<br/>{}".format(node_lines[i])
                    prefix_num = ""
                    if node.pdh != False and ":" not in node_str:
                        prefix_num = " : {}".format(node.pdh)
                    dot_str += ('{} [label=<{}{}>'
                                .format(node.nodeid, node_str, prefix_num))
                    dot_str += ', shape={}, style="filled'.format(node_shape)
                    if node.pdh == False:
                        dot_str += '"'
                    else:
                        dot_str += ',dashed"'
                    #if node.highlighted == True:
                    #   dot_str += ', fillcolor=gold, penwidth=2'
                    #else:
                    dot_str += ', fillcolor={}'.format(node_color)
                    if node.highlighted == True:
                        dot_str += ', penwidth=4'
                    if node.intro == True:
                        dot_str += ', intro={}'.format(node.intro)
                    if node.first == True:
                        dot_str += ', first={}'.format(node.first)
                    if node.pos != None:
                        dot_str += ', pos={}'.format(node.pos)
                    if node.stdedit != None:
                        dot_str += ', stded="{}"'.format(node.stdedit)
                    dot_str += ', ev={}'.format(node.eventid)
                    #dot_str += ', penwidth=2'
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
            #if showintro == False:
            #    for covermesh in self.covermeshes:
            #        for midnode in covermesh.midnodes:
            #            if midnode.rank == current_rank:
            #                dot_str += self.write_midnode(covermesh, midnode,
            #                    average_use, minpenwidth, medpenwidth, maxpenwidth)
            #                dot_str += ', cover="True"] ;\n'
            # Close rank braces.
            if showintro == False and current_rank < 1:
                dot_str += "//"
            dot_str += "}\n"
        # Draw unranked event nodes and shrank nodes.
        for node in self.eventnodes:
            if node.rank == None or node.shrink == True:
                node_lines = textwrap.wrap(node.label, 20,
                                           break_long_words=False)
                node_str = ""
                for i in range(len(node_lines)):
                    if i == 0:
                        node_str += "{}".format(node_lines[i])
                    else:
                        node_str += "<br/>{}".format(node_lines[i])
                if node.shrink == False:
                    node_shape = 'ellipse'
                    node_color = 'white'
                    if showintro == False and node.intro == True:
                        dot_str += '//'
                    dot_str += '{} '.format(node.nodeid)
                    dot_str += '[label=<{}>'.format(node_str)
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
                    dot_str += "] ;\n"
                elif node.shrink == True:
                    dot_str += '{} '.format(node.nodeid)
                    dot_str += '[label="", hlabel=<{}>'.format(node_str)
                    dot_str += ', shape=circle'
                    dot_str += ', style=filled'
                    dot_str += ', fillcolor=white'
                    dot_str += ', width=0.1'
                    dot_str += ', height=0.1'
                    dot_str += ', penwidth=2'
                    dot_str += "] ;\n"
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
        for int_rank in range(int((self.minrank)*(midranks+1)),
                              int(self.maxrank*(midranks+1))):
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
        # If showintro is True, write underlying edges and do not write cover
        # edges.
        # If showintro is False, write underlying edges as comments and write
        # cover edges.
        # The method read_dot reads all underlying edges, even if commented,
        # and does not read cover edges.
        if showintro == True:
            hyperedges_to_write = self.hyperedges
        elif showintro == False:
            hyperedges_to_write = self.hyperedges + self.coverhyperedges
        edges_str = ""
        midid = 1
        for hyperedge in hyperedges_to_write:
            if self.hypergraph == False:
                under = hyperedge.underlying
                for subedge in hyperedge.edgelist:
                    edges_str += self.write_edge(subedge, params,
                                                 underlying=under,
                                                 cover=hyperedge.cover,
                                                 showintro=showintro)
            elif self.hypergraph == True:
                edges_str += self.write_hyperedge(hyperedge, midid, params,
                                                  showintro=showintro)
                midid += 1
        dot_str += edges_str
        # Close graph.
        dot_str += "}"
        self.dot_file = dot_str


    def write_midnode(self, midid, color, scale, underlying=False, cover=False,
                      showintro=True):
        """ Write the line of a dot file for a single midnode."""

        mid_str = ""
        if showintro == False and underlying == True:
            mid_str += "//"
        mid_str += '"{}" [label=""'.format(midid)
        mid_str += ', shape=circle'
        mid_str += ', style=filled'
        mid_str += ', fillcolor={}, color={}'.format(color, color)
        mid_str += ', width={:.2}'.format(0.1*math.sqrt(scale))
        mid_str += ', height={:.2}'.format(0.1*math.sqrt(scale))
        if cover == True:
            mid_str += ', cover=True'
        mid_str += '] ;\n'

        return mid_str


    def write_edge(self, edge, params, arrow=True, custom_src_id=None,
                   custom_trg_id=None, underlying=False, cover=False,
                   showintro=True):
        """ Write the line of a dot file for a single edge. """

        edge_str = ""
        if showintro == False and underlying == True:
            edge_str += "//"
        # Check if a custom source or target is used.
        source_id = custom_src_id
        if custom_src_id == None:
            source_id = edge.source.nodeid
        target_id = custom_trg_id
        if custom_trg_id == None:
            target_id = edge.target.nodeid
        # Write source and target node ids.
        if edge.reverse == False:
            edge_str += ('{} -> {} ['.format(source_id, target_id))
        elif edge.reverse == True:
            edge_str += ('{} -> {} ['.format(target_id, source_id))
            edge_str += 'rev=True, '
        # Write edge color.
        edge_str += 'color={}'.format(edge.color)
        # Overwrite arrow if target is shrunk.
        if isinstance(edge.target, EventNode) and edge.target.shrink == True:
            arrow = False
        if arrow == False:
            edge_str += ', dir=none'
        elif arrow == True and edge.reverse == True:
            edge_str += ', dir=back'
        # Indicate conflict by a dashed line.
        if edge.relationtype == "conflict":
            edge_str += ", style=dashed"
        # Write statistics.
        edge_str += ', w={}'.format(edge.weight)
        if params["weightedges"] == True:
            edge_str += ', weight={}'.format(edge.layout_weight)
        # Compute penwidth.
        #ratio = edge.weight/params["average_weight"]
        ratio = edge.number/params["average_number"]
        pensize = math.log(ratio,2) + params["medpenwidth"]
        if pensize < params["minpenwidth"]:
            pensize = params["minpenwidth"]
        if pensize > params["maxpenwidth"]:
            pensize = params["maxpenwidth"]
        edge_str += ', penwidth={:.2}'.format(pensize)
        # Write labels.
        if params["addedgelabels"] == True:
            if edge.labelcarrier == True:
                edge_str += ', label=" {}\\n'.format(edge.number)
                edge_str += ' {}"'.format(edge.weight)
            if params["showedgelabels"] == True:
                edge_str += ', fontcolor={}'.format(edge.color)
            elif params["showedgelabels"] == False:
                edge_str += ', fontcolor=transparent'
        # Write whether edge is essential.
        if edge.essential == True:
            edge_str += ', ess=True'
        if cover == True:
            edge_str += ', cover=True'
        #if addedgelabels == True:
        #    if midedge.overridelabel == None:
        #        if midedge.labelcarrier == True:
        #            label_str = ""
        #            if edgeid == True:
        #                label_str += "  #{}".format(mesh.meshid)
        #                if edgeocc == True or edgeuse == True:
        #                    label_str += "\\n"
        #            if edgeocc == True:
        #                label_str += "  {}".format(occ_stat)
        #                if edgeuse == True:
        #                   label_str += "\\n"
        #            if edgeuse == True:
        #                label_str += "  {}".format(use_stat)
        #            mid_str += ', label="{}"'.format(label_str)
        #            if midedge.labelpos != None:
        #                mid_str += ', lp={}'.format(midedge.labelpos)
        #        if showedgelabels == True:
        #            mid_str += ', fontcolor={}'.format(midedge.color)
        #        elif showedgelabels == False:
        #            mid_str += ', fontcolor=transparent'
        #    else:
        #        mid_str += ', label="{}"'.format(midedge.overridelabel)
        #        if midedge.labelpos != None:
        #            mid_str += ', lp={}'.format(midedge.labelpos)
        edge_str += '] ;\n'

        return edge_str


    def write_hyperedge(self, hyperedge, midid, params, showintro=True):
        """
        Write the dot lines for the subedges an midnode of an hyperedge.
        """

        if len(hyperedge.edgelist) == 1:
            hyper_str = self.write_edge(hyperedge.edgelist[0], params,
                                        underlying=hyperedge.underlying,
                                        cover=hyperedge.cover,
                                        showintro=showintro)
        elif len(hyperedge.edgelist) > 1:
            # Write mid node.
            midid_str = "mid{}".format(midid)
            hyper_str = self.write_midnode(midid_str, hyperedge.midcolor,
                                           params["edgewidthscale"],
                                           underlying=hyperedge.underlying,
                                           cover=hyperedge.cover,
                                           showintro=showintro)
            # Write subedges without arrow.
            all_conflict = True
            for subedge in hyperedge.edgelist:
                hyper_str += self.write_edge(subedge, params, arrow=False,
                                             custom_trg_id=midid_str,
                                             underlying=hyperedge.underlying,
                                             cover=hyperedge.cover,
                                             showintro=showintro)
                if subedge.relationtype != "conflict":
                    all_conflict = False
            # Write final edge with arrow (if target not shrunk).
            hyper_str += self.write_edge(hyperedge, params,
                                         custom_src_id=midid_str,
                                         underlying=hyperedge.underlying,
                                         cover=hyperedge.cover,
                                         showintro=showintro)

        return hyper_str


    #def write_midnode(self, mesh, midnode, average_use, minpenwidth,
    #                  medpenwidth, maxpenwidth):
    #    """ Write the line of a dot file for a single midnode."""
    #
    #    ratio = mesh.uses/average_use
    #    pensize = math.log(ratio, 2) + medpenwidth
    #    if pensize < minpenwidth:
    #        pensize = minpenwidth
    #    if pensize > maxpenwidth:
    #        pensize = maxpenwidth
    #    pensize = math.sqrt(pensize)/12
    #    mid_str = '"{}" [label=""'.format(midnode.nodeid)
    #    mid_str += ', shape=circle'
    #    if midnode.ghost == True:
    #        mid_str += ', style=dotted'
    #    else:
    #        mid_str += ', style=filled'
    #    mid_str += ', color={}'.format(midnode.bordercolor)
    #    mid_str += ', fillcolor={}'.format(midnode.fillcolor)
    #    mid_str += ', midtype={}'.format(midnode.midtype)
    #    if midnode.overridewidth == None:
    #        mid_str += ', width={:.4f}'.format(pensize)
    #        mid_str += ', height={:.4f}'.format(pensize)
    #    else:
    #        mid_str += ', width={:.4f}'.format(midnode.overridewidth)
    #        mid_str += ', height={:.4f}'.format(midnode.overridewidth)
    #    if midnode.pos != None:
    #        mid_str += ', pos={}'.format(midnode.pos)
    #
    #    return mid_str
    #
    #
    #def write_midedge(self, mesh, midedge, average_use, minpenwidth,
    #                  medpenwidth, maxpenwidth, addedgelabels, showedgelabels,
    #                  edgeid, edgeocc, edgeuse, statstype, weightedges):
    #    """ Write the line of a dot file for a single midedge. """
    #
    #    ratio = mesh.uses/average_use
    #    pensize = math.log(ratio,2) + medpenwidth
    #    if pensize < minpenwidth:
    #        pensize = minpenwidth
    #    if pensize > maxpenwidth:
    #        pensize = maxpenwidth
    #    if midedge.reverse == False:
    #        mid_str = ('"{}" -> "{}" '.format(midedge.source.nodeid,
    #                                          midedge.target.nodeid))
    #    elif midedge.reverse == True:
    #        mid_str = ('"{}" -> "{}" '.format(midedge.target.nodeid,
    #                                          midedge.source.nodeid))
    #    mid_str += '[meshid={}'.format(mesh.meshid)
    #    if midedge.overridewidth == None:
    #        mid_str += ', penwidth={}'.format(pensize)
    #    else:
    #        mid_str += ', penwidth={}'.format(midedge.overridewidth)
    #    mid_str += ', color={}'.format(midedge.color)
    #    if statstype == "abs":
    #        occ_stat = "{}".format(mesh.occurrence)
    #        use_stat = "{}".format(mesh.uses)
    #    elif statstype == "rel":
    #        occ_stat = "{:.2}".format(mesh.rel_occ)
    #        use_stat = "{:.2}".format(mesh.usage)
    #    elif statstype == "both":
    #        occ_stat = "{}".format(mesh.occurrence)
    #        occ_stat += " ({:.2})".format(mesh.rel_occ)
    #        use_stat = "{}".format(mesh.uses)
    #        use_stat += " ({:.2})".format(mesh.usage)
    #    if addedgelabels == True:
    #        if midedge.overridelabel == None:
    #            if midedge.labelcarrier == True:
    #                label_str = ""
    #                if edgeid == True:
    #                    label_str += "  #{}".format(mesh.meshid)
    #                    if edgeocc == True or edgeuse == True:
    #                        label_str += "\\n"
    #                if edgeocc == True:
    #                    label_str += "  {}".format(occ_stat)
    #                    if edgeuse == True:
    #                       label_str += "\\n"
    #                if edgeuse == True:
    #                    label_str += "  {}".format(use_stat)
    #                mid_str += ', label="{}"'.format(label_str)
    #                if midedge.labelpos != None:
    #                    mid_str += ', lp={}'.format(midedge.labelpos)
    #            if showedgelabels == True:
    #                mid_str += ', fontcolor={}'.format(midedge.color)
    #            elif showedgelabels == False:
    #                mid_str += ', fontcolor=transparent'
    #        else:
    #            mid_str += ', label="{}"'.format(midedge.overridelabel)
    #            if midedge.labelpos != None:
    #                mid_str += ', lp={}'.format(midedge.labelpos)
    #    if midedge.indicator == True:
    #        if isinstance(midedge.source, EventNode):
    #            mid_str += ", dir=none"
    #        else:
    #            mid_str += ", dir=both"
    #        if midedge.reverse == False:
    #            if isinstance(midedge.target, MidNode):
    #                mid_str += ", arrowhead=none"
    #            if isinstance(midedge.source, MidNode):
    #                #mid_str += ", arrowtail=icurve"
    #                mid_str += ", arrowtail=crow"
    #                #mid_str += ", arrowtail=inv"
    #        elif midedge.reverse == True:
    #            if isinstance(midedge.target, MidNode):
    #                mid_str += ", arrowtail=none"
    #            if isinstance(midedge.source, MidNode):
    #                #mid_str += ", arrowhead=icurve"
    #                mid_str += ", arrowhead=crow"
    #                #mid_str += ", arrowhead=inv"
    #    elif midedge.indicator == False:
    #        if isinstance(midedge.target, MidNode):
    #            mid_str += ", dir=none"
    #    if midedge.reverse == True:
    #        mid_str += ", rev=True"
    #    if midedge.relationtype == "conflict":
    #        mid_str += ", style=dotted"
    #    mid_str += ', uses={}'.format(mesh.uses)
    #    mid_str += ', usage={}'.format(mesh.usage)
    #    if weightedges == True:
    #        mid_str += ', weight={}'.format(mesh.weight)
    #    if midedge.pos != None:
    #        mid_str += ', pos={}'.format(midedge.pos)
    #
    #    return mid_str


    def clear_adjacency(self):
        """
        Empty adjacency lists to avoid infinit recursions while using deepcopy.
        """

        for node in self.eventnodes + self.statenodes:
            node.incoming = []
            node.outgoing = []


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


def get_stded(read_str):
    """ Extract the value of field in dot file line. """

    if "stded=" in read_str:
        field_start = read_str.index("stded=")+6+1
        rem = read_str[field_start:]
        field_end = rem.index('"')+field_start
        value = read_str[field_start:field_end]
    else:
        value = None

    return value


# Converting alphabet strings to base 10 numbers and vice versa.
# I copied this from an anonymous post (username 301_Moved_Permanently)
# on StackExchange without trying to understand how it works.
A_LOWERCASE = ord('a')
ALPHABET_SIZE = 26
def _decompose(number):
    """Generate digits from `number` in base alphabet, least significants
    bits first.

    Since A is 1 rather than 0 in base alphabet, we are dealing with
    `number - 1` at each iteration to be able to extract the proper digits.
    """

    while number:
        number, remainder = divmod(number - 1, ALPHABET_SIZE)
        yield remainder

def base_10_to_alphabet(number):
    """Convert a decimal number to its base alphabet representation"""

    return ''.join(
            chr(A_LOWERCASE + part)
            for part in _decompose(number)
    )[::-1]

def base_alphabet_to_10(letters):
    """Convert an alphabet number to its decimal representation"""

    return sum(
            (ord(letter) - A_LOWERCASE + 1) * ALPHABET_SIZE**i
            for i, letter in enumerate(reversed(letters.lower()))
    )

# -------------------- Causal Cores Generation Section ------------------------

def getcausalcores(eoi, kappamodel, kasimpath, kaflowpath,
                   simtime=1000, simseed=None, precedenceonly=True,
                   ignorelist=[], eoi_def=None):
    """ 
    Generate initial causal cores of given event of interest by running KaSim 
    and then KaFlow.
    """

    new_model = add_eoi(eoi, kappamodel, eoi_def)
    trace_path = run_kasim(new_model, kasimpath, simtime, simseed)
    run_kaflow(eoi, trace_path, kaflowpath, precedenceonly, ignorelist)


def getstories(eoi, kappamodel, kasimpath, kastorpath,
               simtime=1000, simseed=None, compression=None, eoi_def=None):
    """
    Generate stories using the fill_siphon function of KaStor,
    and optionally weak compression.
    """

    new_model = add_eoi(eoi, kappamodel, eoi_def)
    trace_path = run_kasim(new_model, kasimpath, simtime, simseed)
    run_kastor(eoi, trace_path, kastorpath, compression)


def add_eoi(eoi, kappamodel, eoi_def=None):
    """ Create a new Kappa model where the EOI is added. """

    if not os.path.exists(eoi):
        os.mkdir(eoi)
    if not os.path.exists("{}/tmp".format(eoi)):
        os.mkdir("{}/tmp".format(eoi))
    if not os.path.exists("{}/unique".format(eoi)):
        os.mkdir("{}/unique".format(eoi))
    last_dot = kappamodel.rfind(".")
    prefix = kappamodel[:last_dot]
    new_path = "{}/{}-eoi.ka".format(eoi, prefix)
    shutil.copyfile(kappamodel, new_path)
    new_file = open(new_path, "a")
    if eoi_def != None:
        new_file.write("%obs: '{}' {}\n".format(eoi, eoi_def))
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


def run_kaflow(eoi, trace_path, kaflowpath, precedenceonly, ignorelist):
    """ Run KaFlow on the trace containing the EOI. """

    if precedenceonly == True:
        subprocess.run(("{}".format(kaflowpath),
                        "--precedence-only",
                        "-o", "{}/tmp/causalcore-".format(eoi),
                        "{}".format(trace_path)))
    elif precedenceonly == False:
        subprocess.run(("{}".format(kaflowpath),
                        "-o", "{}/tmp/causalcore-".format(eoi),
                        "{}".format(trace_path)))
    core_files = get_dot_files("{}/tmp".format(eoi), "causalcore")
    filenum = 1
    nignored = 0
    for core_file in core_files:
        input_path = "{}/tmp/{}".format(eoi, core_file)
        input_file = open(input_path, "r")
        content = input_file.readlines()
        contains_ignored = check_ignored(eoi, input_path, ignorelist)
        input_file.close()
        os.remove(input_path)
        if contains_ignored == True:
            nignored += 1
        elif contains_ignored == False:
            content.insert(1, '  precedenceonly="{}"\n'.format(precedenceonly))
            content.insert(1, '  producedby="KaFlow"\n')
            output_file = open("{}/tmp/causalcore-{}.dot".format(eoi, filenum), "w")
            output_file.writelines(content)
            output_file.close()
            filenum += 1
    print("Ignoring {} cores out of {} because they contain undo rules."
        .format(nignored, len(core_files)))


def check_ignored(eoi, input_path, ignorelist):
    """
    Check if core contains any ignored term in any of its nodes.
    """

    contains_ignored = False
    core = CausalGraph(input_path, eoi)
    for eventnode in core.eventnodes:
        if any(ignorestr in eventnode.label for ignorestr in ignorelist):
            contains_ignored = True

    return contains_ignored


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
        # Using weak compression here will sometimes remove loop events.
        #subprocess.run(("{}".format(kastorpath), "--weak", "{}".format(trace_path)))
        #os.rename("cflow_Weakly_0.dot", "{}/tmp/siphonweak-{}.dot".format(eoi, filenumber))
        #os.remove("cflow_Weakly_Summary.dat")
        os.remove(trace_path)
    # Add a line to indicate that the file was produced by KaStor.
    story_files = get_dot_files("{}/tmp".format(eoi), "siphon")
    #story_files = get_dot_files("{}/tmp".format(eoi), "siphonweak")
    for story_file in story_files:
        input_file = open("{}/tmp/{}".format(eoi, story_file), "r")
        content = input_file.readlines()
        input_file.close()
        content.insert(3, ' producedby="KaStor"\n')
        output_file = open("{}/tmp/{}".format(eoi, story_file), "w")
        output_file.writelines(content)
        output_file.close()
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
    story_files = get_dot_files("{}/tmp".format(eoi), "causalcore")
    # I cannot use KaStor fillsiphon for now because I lose event ids in the
    # process. I need to be able the recover those ids.
    #story_files = get_dot_files("{}/tmp".format(eoi), "siphon")
    stories = []
    for story_file in story_files:
        story_path = "{}/tmp/{}".format(eoi, story_file)
        stories.append(CausalGraph(story_path, eoi))
    # Tweak each story.
    for i in range(len(stories)):
        stories[i].filename = "story-{}.dot".format(i+1)
        stories[i].producedby = "KappaPathways"

    ## Optionally keep only short stories for faster calculation.
    #num = 1
    #stories_tmp = []
    #for story in stories:
    #    if story.maxrank < 20:
    #        story.filename = "story-{}.dot".format(num)
    #        story.producedby = "KappaPathways"
    #        stories_tmp.append(story)
    #        num += 1
    #stories = stories_tmp

    for story in stories:
        #story.compute_relstats()
        #story.compute_visuals(showintro, color=False)
        #story.create_hyperedges()
        story.rank_sequentially(intropos="bot2")
        story.occurrence = None
        story.build_dot_file(showintro, addedgelabels, showedgelabels,
                             edgeid, edgeocc, edgeuse, statstype, weightedges)
        output_path = "{}/{}".format(eoi, story.filename)
        outfile = open(output_path, "w")
        outfile.write(story.dot_file)
        outfile.close()

# ;;;;;;;;;;;;;;;;;;;;;;;;;;;; Context Section ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

# Get StateNodes:
# A state is a list of sites.
# A site is a dict with agent id, site name, binding and internal value, etc.
# Sites also have a type, which can be edit or context.
# EventNodes contain a list of states, each state coming from one action
# from the trace.
# All sites of every states of an EventNode are of type edit.
# Each state of an EventNode gives rise to one StateNode.
# Each StateNode contains one state, which is a list of sites, some of which
# are of type edit and some of type context.
# The context sites of a StateNode are given by the cumulative relevant
# context of the story.

# 2nd try:
# A state is a list of agents.
# An agent is a list of sites.
# A site is a dict with agent id, site name, binding and internal value, etc.
# Sites also have a type, which can be edit or context.
# EventNodes contain a list of states, each state coming from one action
# from the trace.
# All sites of every states of an EventNode are of type edit.
# Each state of an EventNode gives rise to one StateNode.
# Each StateNode contains one state, which is a list of sites, some of which
# are of type edit and some of type context.
# The context sites of a StateNode are given by the cumulative relevant
# context of the story.

# state = [agents]
# agent = {"name": ,"id": ,"sites": [sites], "type": (None if sites exist),
#          "action": (None if sites exist)}
# site = {"name": ,"bond": ,"value", "agentname": ,"agentid": ,
#         "type": (edit, context or parallel), ("action": or "test:")}
# bond = {"num": ,"partner": }
# partner = {"agentname": ,"agentid": ,"sitename":}


def getdualstories(eoi, kappamodel, showintro=True, addedgelabels=False,
                   showedgelabels=False, edgeid=True, edgeocc=False,
                   edgeprob=False, statstype="abs", writedot=True,
                   weightedges=True):
    """
    Add state nodes showing what was edited by rules along with
    the relevant context.
    """

    # Reading section.
    story_files = get_dot_files("{}".format(eoi), "story")
    stories = []
    for story_file in story_files:
        story_path = "{}/{}".format(eoi, story_file)
        stories.append(CausalGraph(story_path, eoi))
    # Keep a copy of the stories.
    #storiescopy = copy.deepcopy(stories)
    # Read trace.
    period = kappamodel.rfind(".")
    modelname = kappamodel[:period]
    tracefile = open("{}/{}-eoi.json".format(eoi, modelname))
    sim = json.load(tracefile)
    signatures = sim["model"]["update"]["signatures"]
    #rules = sim["model"]["ast_rules"]
    steps = sim["trace"]

    # Find edits produced by each event.
    #tmp_stories = stories[21:22]
    #stories = tmp_stories
    for story in stories:
        print(story.filename)
        # Get actions for each event node.
        for eventnode in story.eventnodes:
            step = steps[int(eventnode.nodeid[2:])]
            bnd_num = 1
            editlist = [] # List of states containing only sites of type edit.
            if step[0] == 1: # Rule
                actions = step[2][1]
                for action in actions:
                    edit, bnd_num = state_from_action(signatures, action,
                                                      bnd_num)
                    editlist.append(edit)
                eventnode.edits = editlist
            if step[0] == 3: # Init
                actions = step[1]
                for action in actions:
                    edit, bnd_num = state_from_action(signatures, action,
                                                       bnd_num)
                    # Check if last state is the same Bind_to as a previous
                    # one but reversed.
                    add_edit = True
                    if action[0] == 3:
                        ag1 = edit[0]["name"]
                        ag2 = edit[1]["name"]
                        id1 = edit[0]["id"]
                        id2 = edit[1]["id"]
                        for prev_edit in editlist:
                            if prev_edit[0]["sites"] != None:
                                if prev_edit[0]["sites"][0]["action"] == 3:
                                    pa1 = prev_edit[0]["name"]
                                    pa2 = prev_edit[1]["name"]
                                    pi1 = prev_edit[0]["id"]
                                    pi2 = prev_edit[1]["id"]
                                    if pa1 == ag2 and pi1 == id2:
                                        if pa2 == ag1 and pi2 == id1:
                                            add_edit = False
                    if add_edit == True:
                        editlist.append(edit)
                # Remove creations if the agents are present in other actions.
                for edit in editlist:
                    keep_edit = True
                    if edit[0]["action"] == 0:
                        ag1 = edit[0]["name"]
                        id1 = edit[0]["id"]
                        for edit2 in editlist:
                            if edit2[0]["action"] != 0:
                                for elem in edit2:
                                    ag2 = elem["name"]
                                    id2 = elem["id"]
                                    if ag1 == ag2 and id1 == id2:
                                        keep_edit = False
                    if keep_edit == True:    
                        eventnode.edits.append(edit)
            # Group sites by agents for each edit.
            grouped_edits = []
            for edit in eventnode.edits:
                new_edit = group_sites_by_agent(edit)
                grouped_edits.append(new_edit)
            eventnode.edits = grouped_edits
        # Add a StateNode for each state of EventNodes.
        state_id = 1
        for eventnode in story.eventnodes:
            for edit in eventnode.edits:
                node_id = "state{}".format(state_id)
                rank = eventnode.rank + 0.5
                lbl = write_context_expression(edit)
                new_state_node = StateNode(node_id, lbl, rank, edit=edit,
                                           eventid=eventnode.eventid)
                story.statenodes.append(new_state_node)
                new_edge = CausalEdge(eventnode, new_state_node)
                story.causaledges.append(new_edge)
                state_id += 1
        # Get tests for each event node.
        for eventnode in story.eventnodes:
            eventnode.tests = []
            step = steps[int(eventnode.nodeid[2:])]
            bnd_num = 1
            tests = None
            if step[0] == 1: # Rule
                tests_tmp = step[2][0]
                tests = []
                for sublist in tests_tmp:
                    for test in sublist:
                        tests.append(test)
            if step[0] == 4: # Obs
                tests = step[2][0]
            tmp_states = []
            if tests != None:
                for test in tests:
                    state, bnd_num = state_from_test(signatures, test, bnd_num)
                    if state[0]["test"] != 0:
                        tmp_states.append(state)
                eventnode.tests = tmp_states
            # Group sites by agents for each test.
            grouped_tests = []
            for test in eventnode.tests:
                new_test = group_sites_by_agent(test)
                grouped_tests.append(new_test)
            eventnode.tests = grouped_tests
        # Add edges between state nodes and the event nodes that require them.
        # A given test must have only one edit that satisfy it.
        for eventnode in story.eventnodes:
            for eventtest in eventnode.tests:
                satisfying_edits = []
                for statenode in story.statenodes:
                    if statenode.rank < eventnode.rank:
                        are_rel = edit_vs_test(statenode.edit, eventtest)
                        if are_rel == True:
                            satisfying_edits.append(statenode)
                if len(satisfying_edits) > 0:
                    sorted_edits = sorted(satisfying_edits,
                                          key=lambda x: x.rank, reverse=True)
                    new_edge = CausalEdge(sorted_edits[0], eventnode)
                    story.causaledges.append(new_edge)
        # Add conflict edges between states and target events.
        for edge in story.causaledges:
            if edge.relationtype == "conflict":
                out_state_nodes = []
                for edge2 in story.causaledges:
                    if edge2.source == edge.source:
                        if isinstance(edge2.target, StateNode):
                            out_state_nodes.append(edge2.target)
                for out_state_node in out_state_nodes:
                    new_edge = CausalEdge(out_state_node, edge.target,
                                          relationtype="conflict")
                    story.causaledges.append(new_edge)
        # Find secondary edges.
        for statenode in story.statenodes:
            # Get source events.
            source_events = []
            for edge in story.causaledges:
                if edge.target == statenode:
                    source_events.append(edge.source)
            # Get direct target events.
            direct_target_events = []
            for edge in story.causaledges:
                if edge.source in source_events:
                    if isinstance(edge.target, EventNode):
                        direct_target_events.append(edge.target)
            # Get outgoing edges.
            outgoing_edges = []
            for edge in story.causaledges:
                if edge.source == statenode:
                    outgoing_edges.append(edge)
            # Assign secondary edges.
            for edge in outgoing_edges:
                if edge.target not in direct_target_events:
                    edge.secondary = True
                    edge.color = "gray60"
        # Revert secondary for rule outputs that have only secondary edges.
        story.rule_outputs = []
        for edge in story.causaledges:
            if isinstance(edge.source, EventNode):
                if edge.source.intro == False:
                    if isinstance(edge.target, StateNode):
                        story.rule_outputs.append(edge.target)
        for statenode in story.rule_outputs:
            outgoing_edges = []
            target_ranks = []
            has_normal = False
            for i in range(len(story.causaledges)):
                edge = story.causaledges[i]
                if edge.source == statenode:
                    outgoing_edges.append(i)
                    target_ranks.append(edge.target.rank)
                    if edge.secondary == False:
                        has_normal = True
            if has_normal == False:
                for i in outgoing_edges:
                    if story.causaledges[i].target.rank == min(target_ranks):
                        story.causaledges[i].secondary = False
                        story.causaledges[i].color = "black"
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
        ## Remove unused state nodes from intros.
        #nodes_to_remove = []
        #edges_to_remove = []
        #for i in range(len(story.statenodes)):
        #    statenode = story.statenodes[i]
        #    if statenode not in story.rule_outputs:
        #        outgoing_edges = []
        #        for edge in story.causaledges:
        #            if edge.source == statenode:
        #                outgoing_edges.append(edge)
        #        if len(outgoing_edges) == 0: # Remove state node
        #            nodes_to_remove.insert(0, i)
        #            for j in range(len(story.causaledges)):
        #                edge = story.causaledges[j]
        #                if edge.target == statenode:
        #                    if j not in edges_to_remove:
        #                        edges_to_remove.insert(0, j)               
        #for i in nodes_to_remove:
        #    del(story.statenodes[i])
        #for j in edges_to_remove:
        #    del(story.causaledges[j])
        story.create_hyperedges()
    # Write stories with edited states.
    for i in range(len(stories)):
        stories[i].filename = "edits-{}.dot".format(i+1)
    for story in stories:
        story.occurrence = None
        story.get_maxrank()
        story.build_dot_file(showintro, addedgelabels, showedgelabels,
                             edgeid, edgeocc, edgeprob, statstype, weightedges)
        output_path = "{}/tmp/{}".format(eoi, story.filename)
        outfile = open(output_path, "w")
        outfile.write(story.dot_file)
        outfile.close()

    #time_start = time.perf_counter()
    ## Then write state (context + edit).
    #""" Get the cumulative relevant context for each state node. """
    #for story in stories:
    #    print(story.filename)
    #    story.hyperedges = []
    #    story.rule_outputs = []
    #    for edge in story.causaledges:
    #        if isinstance(edge.source, EventNode):
    #            if edge.source.intro == False:
    #                if isinstance(edge.target, StateNode):
    #                    story.rule_outputs.append(edge.target)
    #    for statenode in story.rule_outputs:
    #        cumul_nodes = []
    #        # Get upstream path of state nodes.
    #        # This part takes a long time on large graphs, may be improved.
    #        paths = story.follow_edges("up", statenode, ignore_conflict=False)
    #        state_paths = []
    #        for path in paths:
    #            state_path = []
    #            for node in path:
    #                if isinstance(node, StateNode):
    #                    state_path.append(node)
    #            state_paths.append(state_path)
    #        # Keep only the first value encountered for each site encountered
    #        # while following path up.
    #        seen_sites = []
    #        for agent in statenode.edit:
    #            if agent["sites"] != None:
    #                for site in agent["sites"]:
    #                    site_str = write_kappa_site(site, "num", True)
    #                    seen_sites.append(site_str)
    #            else:
    #                agent_str = write_kappa_agent(agent, "num", True)
    #                seen_sites.append(agent_str)
    #        for path in state_paths:
    #            for node in path:
    #                # All edited sites of a state node must be unseen before
    #                # we can keep the node.
    #                keep_node = True
    #                for agent in node.edit:
    #                    if agent["sites"] != None:
    #                        for site in agent["sites"]:
    #                            site_str = write_kappa_site(site, "num", True)
    #                            if site_str in seen_sites:
    #                                keep_node = False
    #                            else:
    #                                seen_sites.append(site_str)
    #                    else:
    #                        agent_str = write_kappa_agent(agent, "num", True)
    #                        if agent_str in seen_sites:
    #                            keep_node = False
    #                        else:
    #                            seen_sites.append(agent_str)
    #                if keep_node == True and node not in cumul_nodes:
    #                    cumul_nodes.append(node)
    #        # Also add all the state nodes which have the same event as source.
    #        # !! I won't need this part once I implement parallel context !!
    #        for edge in story.causaledges:
    #            if edge.target == statenode:
    #                source_event = edge.source
    #        for edge in story.causaledges:
    #            if edge.source == source_event:
    #                if edge.target != statenode:
    #                    cumul_nodes.append(edge.target)
    #        # Check which of the cumul nodes are relevant for the future of the
    #        # current state node.
    #        relevant_nodes = []
    #        for cumul_node in cumul_nodes:
    #            # Find all target event nodes.
    #            target_events = []
    #            for edge in story.causaledges:
    #                if edge.source == cumul_node:
    #                    target_events.append(edge.target)
    #            # Check if there is at least one of the target nodes which is
    #            # in the future of current state node (has an upstream path).
    #            downstream_paths = story.follow_edges("down", statenode,
    #                                                  target_events,
    #                                                  ignore_conflict=True,
    #                                                  stop_at_first=True)
    #            if len(downstream_paths) > 0:                
    #                relevant_nodes.append(cumul_node)
    #        # Build current state node context from the state of
    #        # all the relevant_nodes.
    #        full_state = copy.deepcopy(statenode.edit)
    #        for relevant_node in relevant_nodes:
    #            for agent in relevant_node.edit:
    #                context_agent = copy.deepcopy(agent)
    #                if context_agent["type"] == None:
    #                    for context_site in context_agent["sites"]:
    #                        context_site["type"] = "context"
    #                elif context_agent["type"] != None:
    #                    context_agent["type"] = "context"
    #                full_state.append(context_agent)
    #        statenode.state = group_sites_by_agent(full_state)
    #        lbl = write_context_expression(statenode.state)
    #        statenode.label = lbl
    #    for statenode in story.statenodes:
    #        if statenode not in story.rule_outputs:
    #            statenode.state = statenode.edit
    #    story.create_hyperedges()
    #    story.align_vertical()

    #time_stop = time.perf_counter()
    #time_diff = time_stop - time_start
    #print("Time:", time_diff)

    #time_start = time.perf_counter()

    # Add factice agents to do the latter work of agents
    # which where fortuitously involved in several events.
    # !!! Will fail if intro events introduce bonds !!!
    for story in stories:
        print(story.filename)
        story.hyperedges = []
        story.build_adjacency(hyper=False)
        for eventnode in story.eventnodes:
            if eventnode.label == eoi:
                story.eoi_node = eventnode
        story.get_all_reachables()
        # Gather intro_outputs and rule_outputs.
        story.rule_outputs = []
        story.intro_outputs = []
        originals = []
        for event in story.eventnodes:
            if event.intro == True:
                for edge in event.outgoing:
                    story.intro_outputs.append(edge.target)
                    originals.append(edge.target)
            elif event.intro == False:
                for edge in event.outgoing:
                    story.rule_outputs.append(edge.target)
        # Gather state nodes that represent an original state.
        for statenode in story.rule_outputs:
            for introstate in story.intro_outputs:
                are_same = compare_states(statenode.edit, introstate.edit,
                                          ignorevalue=False, ignoreid=True)
                if are_same == True:
                    originals.append(statenode)
                    break
        # Make a dict of all agent ids in the story.
        ids_suffix = {}
        for statenode in story.statenodes:
            for agent in statenode.edit:
                idstr = str(agent["id"])
                if idstr not in ids_suffix.keys():
                    ids_suffix[idstr] = "b"
        # And a dict if eventids for the creation of new intro nodes.
        eventids = {}
        new_intro_nodes = {}
        # Events must be ordered by rank (the order of the nodes within a same
        # rank does not matter). They are already ordered at this point because
        # read_dot reads the node in the order in which they appear in the
        # story file.
        moved_edits = []
        for eventnode in story.eventnodes:
            if eventnode.intro == False:
                # Separate immediate upstream nodes in two groups.
                # upstream_originals are state nodes part of list originals and
                # reach the current event through a transitive edge.
                # upstream_others are all other immediate upstream nodes.
                upstream_originals = []
                upstream_others = []
                for edge in eventnode.incoming:
                    if edge.secondary == True and edge.source in originals:
                         upstream_originals.append(edge.source)
                    else:
                         upstream_others.append(edge.source)
                # Gather agent ids from upstream_originals.
                original_ids = []
                for upstream_original in upstream_originals:
                    for agent in upstream_original.edit:
                        if str(agent["id"]) not in original_ids:
                            original_ids.append(str(agent["id"]))
                            # Get eventid of upstream intro node.
                            intronode = upstream_original.incoming[0].source
                            eventids[str(agent["id"])] = [intronode.label,
                                                          intronode.nodeid]
                # Gather agent ids from upstream_others.
                other_ids = []
                for upstream_other in upstream_others:
                    for agent in upstream_other.edit:
                        if str(agent["id"]) not in other_ids:
                            other_ids.append(str(agent["id"]))
                # original_id suffixes are changed to match corresponding
                # other_ids. If no other_id matches an original_id, the new
                # original_id suffix is taken from ids_suffix.
                id_changes = {}
                for original_id in original_ids:
                    original_number = remove_suffix(original_id)
                    original_found = False
                    for other_id in other_ids:
                        other_number = remove_suffix(other_id)
                        if other_number == original_number:
                            original_found = True
                            if other_id != original_id:
                                id_changes[original_id] = other_id
                    if original_found == False:
                        suffix = ids_suffix[original_number]
                        new_id = "{}{}".format(original_number, suffix)
                        id_changes[original_id] = new_id
                        # Create new intro node.
                        eventids[original_id]
                        node_id = "{}{}".format(eventids[original_id][1], suffix)
                        label = eventids[original_id][0]
                        rank = eventnode.rank - 1
                        new_node = EventNode(node_id, label, rank, intro=True)
                        story.eventnodes.append(new_node)
                        new_intro_nodes[new_id] = new_node
                        # Increment suffix.
                        suffix_index = base_alphabet_to_10(suffix)
                        new_suffix = base_10_to_alphabet(suffix_index+1)
                        ids_suffix[original_number] = new_suffix
                # Apply id changes to upstream_originals and downstream state
                # nodes that are reachable from current event.
                reachable_statenodes = []
                for node in eventnode.reachable:
                    if isinstance(node, StateNode):
                        reachable_statenodes.append(node)
                target_nodes = upstream_originals + reachable_statenodes
                changed_nodes = []
                for target_node in target_nodes:
                    for agent in target_node.edit:
                        is_changed = change_agent_id(agent, id_changes)
                        if is_changed == True and target_node in upstream_originals:
                            if target_node not in changed_nodes:
                                changed_nodes.append(target_node)
                    lbl = write_context_expression(target_node.edit)
                    target_node.label = lbl
                # Detach upstream_originals from their initial intro node and attach
                # them to the new one.
                for changed_node in changed_nodes:
                    agent_id = str(changed_node.edit[0]["id"])
                    expected_intro = new_intro_nodes[agent_id]
                    if expected_intro != changed_node.incoming[0].source:
                        for edge in story.causaledges:
                            if edge.target == changed_node:
                                break
                        story.causaledges.remove(edge)
                        new_edge = CausalEdge(expected_intro, changed_node)
                        story.causaledges.append(new_edge) 
                        changed_node.incoming = [new_edge]
                        changed_node.rank = eventnode.rank - 0.5
                        # Set edges as non-transitive.
                        if changed_node not in moved_edits:
                            for edge in changed_node.outgoing:
                                if edge.target == eventnode:
                                    edge.secondary = False
                                    edge.color = "black"
                            moved_edits.append(changed_node)               

        story.create_hyperedges()
    # Write stories with factice agents.
    for i in range(len(stories)):
        stories[i].filename = "facticeagent-{}.dot".format(i+1)
    for story in stories:
        story.build_dot_file(showintro, addedgelabels, showedgelabels,
                             edgeid, edgeocc, edgeprob, statstype, weightedges)
        output_path = "{}/tmp/{}".format(eoi, story.filename)
        outfile = open(output_path, "w")
        outfile.write(story.dot_file)
        outfile.close()


    # Build states with upstream context.
    for story in stories:
        print(story.filename)
        story.hyperedges = []
        story.rule_outputs = []
        story.build_adjacency(hyper=False)
        for edge in story.causaledges:
            if isinstance(edge.source, EventNode):
                if edge.source.intro == False:
                    if isinstance(edge.target, StateNode):
                        story.rule_outputs.append(edge.target)

        for eventnode in story.eventnodes:
            if eventnode.label == eoi:
                story.eoi_node = eventnode

        story.get_all_reachables()
        for statenode in story.statenodes:
            statenode.cumulnodes = []
        for cr in range(1, story.maxrank):
            current_rank = cr + 0.5
            for statenode in story.rule_outputs:
                if statenode.rank == current_rank:
                    #fullcumul, src_rule = get_fullcumul(statenode, story)
                    fullcumul = get_fullcumul(statenode, story)
                    # Check if cumul node is relevant for the future of
                    # statenode. Relevant cumul nodes have at least one
                    # path to the EOI that does not pass through the
                    # statenode.
                    relevantcumul = []
                    remainingcumul = []
                    for cumulnode in fullcumul:
                        relevant = story.reachability_with_block(cumulnode,
                            statenode.reachable, statenode)
                        if relevant == True:
                            relevantcumul.append(cumulnode)
                        #downstream_paths = story.follow_edges("down",
                        #    cumulnode, [story.eoi_node], block=src_rule,
                        #    ignore_conflict=True, stop_at_first=False)
                        ## If there is no path, then downstream_paths = [].
                        #if len(downstream_paths) > 0:
                        #    relevantcumul.append(cumulnode)
                    statenode.cumulnodes = relevantcumul

                    ## Find immediate upstream state nodes.
                    ## Also get neighbors of state node.
                    #for edge in story.causaledges:
                    #    if edge.target == statenode:
                    #        src_rule = edge.source
                    #        break
                    #upstream_nodes = []
                    #neighbors = []
                    #for edge in story.causaledges:
                    #    if edge.secondary == False:
                    #        if edge.target == src_rule:
                    #            upstream_nodes.append(edge.source)
                    #        if edge.source == src_rule:
                    #            if edge.target != statenode:
                    #                neighbors.append(edge.target)
                    ## Add neighbors of upstream state nodes.
                    #up_rules = []
                    #for edge in story.causaledges:
                    #    if edge.target in upstream_nodes:
                    #        if edge.source not in up_rules:
                    #            up_rules.append(edge.source)
                    #upstream_neighbors = []
                    #for edge in story.causaledges:
                    #    if edge.source in up_rules:
                    #        if edge.target not in upstream_nodes:
                    #            upstream_neighbors.append(edge.target)
                    #upstream_nodes += upstream_neighbors
                    ## Add the immediate upstream nodes and their cumulnodes
                    ## to the cumulnodes of the current statenode.
                    #for upstream_node in upstream_nodes:
                    #    statenode.cumulnodes.append(upstream_node)
                    #    for up_cumulnode in upstream_node.cumulnodes:
                    #        if up_cumulnode not in statenode.cumulnodes:
                    #            statenode.cumulnodes.append(up_cumulnode)


                    ## Not even needed if I remove irrelevant cumul correctly.
                    ##
                    ## ! Wrong ! 
                    ## ! Remove from cumulnodes any node that has at least one !
                    ## ! edit site in common with the current node edit sites !
                    ##
                    ## Remove from cumulnodes any node that has at least one
                    ## edit site in common with the current node edit sites or
                    ## a neighbor node edit sites.
                    #current_sites = []
                    #for node in [statenode] + neighbors:
                    #    for agent in node.edit:
                    #        if agent["sites"] != None:
                    #            for site in agent["sites"]:
                    #                site_str = write_kappa_site(site, "num", True)
                    #                current_sites.append(site_str)
                    #        else:
                    #            agent_str = write_kappa_agent(agent, "num", True)
                    #            current_sites.append(agent_str)
                    #if statenode.nodeid == "state53":
                    #    for s in current_sites:
                    #        print(s)
                    #    print("---")
                    #cumul_to_remove = []
                    #for i in range(len(statenode.cumulnodes)):
                    #    prevnode = statenode.cumulnodes[i]
                    #    #if statenode.nodeid == "state53":
                    #    #    print(prevnode)
                    #    for agent in prevnode.edit:
                    #        if agent["sites"] != None:
                    #            for site in agent["sites"]:
                    #                site_str = write_kappa_site(site, "num", True)
                    #                #if statenode.nodeid == "state30":
                    #                #    print(">>", site_str)
                    #                if site_str in current_sites:
                    #                    if i not in cumul_to_remove:
                    #                        cumul_to_remove.insert(0, i)
                    #        else:
                    #            agent_str = write_kappa_agent(agent, "num", True)
                    #            if agent_str in current_sites:
                    #                if i not in cumul_to_remove:
                    #                    cumul_to_remove.insert(0, i)
                    ##if statenode.nodeid == "state30":
                    ##    print(cumul_to_remove)
                    #for i in cumul_to_remove:
                    #    del(statenode.cumulnodes[i])


                    ## Remove cumul nodes that are no more useful for future.
                    ## AKA keep only relevant context.
                    #cumul_to_remove = []
                    #for i in range(len(statenode.cumulnodes)):
                    #    cumulnode = statenode.cumulnodes[i]


                    #    ## Find all target event nodes.
                    #    #target_events = []
                    #    #for edge in story.causaledges:
                    #    #    if edge.source == cumulnode:
                    #    #        target_events.append(edge.target)
                    #    ## Check if there is at least one of the target nodes
                    #    ## which is in the future of current state node
                    #    ## (has an upstream path).
                    #    #downstream_paths = []
                    #    #if len(target_events) > 0:
                    #    #    downstream_paths = story.follow_edges("down",
                    #    #        statenode, target_events, ignore_conflict=True,
                    #    #        stop_at_first=True)
                    #    #if len(downstream_paths) == 0:
                    #    #    cumul_to_remove.insert(0, i)


                    #    # Check if cumul node is relevant for the future of
                    #    # statenode. Relevant cumul nodes have at least one
                    #    # path to the EOI that does not pass through the
                    #    # statenode.
                    #    downstream_paths = story.follow_edges("down",
                    #        cumulnode, [eoi_node], block=src_rule,
                    #        ignore_conflict=True, stop_at_first=False)
                    #    one_path_down = False
                    #    for path in downstream_paths:
                    #        if len(path) > 0:
                    #            one_path_down = True
                    #            break
                    #    if one_path_down == False:
                    #        cumul_to_remove.insert(0, i)
                    #for i in cumul_to_remove:
                    #    del(statenode.cumulnodes[i])

                    ## Add neighbors nodes. This part can probably be removed
                    ## once parallel context is implemented.
                    #neighbors = []
                    #for edge in story.causaledges:
                    #    if edge.source == src_rule:
                    #        if edge.target != statenode:
                    #            neighbors.append(edge.target)
                    # Build current state node context from the state of
                    # all the relevant_nodes.
                    full_state = copy.deepcopy(statenode.edit)
                    for cumulnode in statenode.cumulnodes: # + neighbors:
                        for agent in cumulnode.edit:
                            context_agent = copy.deepcopy(agent)
                            if context_agent["type"] == None:
                                for context_site in context_agent["sites"]:
                                    context_site["type"] = "context"
                            elif context_agent["type"] != None:
                                context_agent["type"] = "context"
                            full_state.append(context_agent)
                    statenode.state = group_sites_by_agent(full_state)
                    lbl = write_context_expression(statenode.state)
                    statenode.label = lbl
        for statenode in story.statenodes:
            if statenode in story.rule_outputs:
                statenode.introstate = False
            #if statenode not in story.rule_outputs:
            else:
                statenode.state = statenode.edit
                statenode.introstate = True
        story.create_hyperedges()
        story.align_vertical()
    #time_stop = time.perf_counter()
    #time_diff = time_stop - time_start
    #print("Time:", time_diff)

    # Write stories with upstream context on state nodes.
    for i in range(len(stories)):
        stories[i].filename = "state-{}.dot".format(i+1)
    for story in stories:
        story.build_dot_file(showintro, addedgelabels, showedgelabels,
                             edgeid, edgeocc, edgeprob, statstype, weightedges)
        output_path = "{}/tmp/{}".format(eoi, story.filename)
        outfile = open(output_path, "w")
        outfile.write(story.dot_file)
        outfile.close()

    # Globally relevant (harmonize).
    # If a given context is relevant in one story, it must be set as
    # relevant in all stories. Must take into account that a same type of
    # edit may be used several times in the context of a node.
    all_cumul_edits = {}
    for story in stories:
        for statenode in story.rule_outputs:
            # Get the edit of each cumul node from the current state node.
            # Since agent ids will not be considered, there may be many
            # equivalent edits.
            cumul_edits = []
            for cn in statenode.cumulnodes:
                cumul_edits.append(cn.edit)
            # Get edit label for the current state node.
            edit_lbl = write_context_expression(statenode.edit, hideid=True)
            # If this edit_lbl was not seen before, assign all cumul_edits.
            if edit_lbl not in all_cumul_edits.keys():
                all_cumul_edits[edit_lbl] = cumul_edits
            # Else, if edit_lbl was already seen, add only new cumul_edits,
            # including dublicates that were not observed before.
            else:
                allcumulcopy = copy.deepcopy(all_cumul_edits[edit_lbl])
                additional_cumul_edits = []
                for cumul_edit in cumul_edits:
                    cumul_found = False
                    for i in range(len(allcumulcopy)):
                        are_same = compare_states(allcumulcopy[i],
                                                  cumul_edit,
                                                  ignorevalue=False,
                                                  ignoreid=True)
                        if are_same == True:
                            cumul_found = True
                            del(allcumulcopy[i])
                            break
                    if cumul_found == False:
                        additional_cumul_edits.append(cumul_edit)
                all_cumul_edits[edit_lbl] += additional_cumul_edits
    #print(all_cumul_edits["<b>FES(act{t})</b>"])
    # Rebuild states with relevant upstream context from all stories.
    for story in stories:
        story.hyperedges = []
        # Reset cumulnodes for each state node.
        for statenode in story.statenodes:
            statenode.cumulnodes = []
        for eventnode in story.eventnodes:
            if eventnode.label == eoi:
                eoi_node = eventnode
        for cr in range(1, story.maxrank):
            current_rank = cr + 0.5
            for statenode in story.rule_outputs:
                if statenode.rank == current_rank:
                    #fullcumul, src_rule = get_fullcumul(statenode, story)
                    fullcumul = get_fullcumul(statenode, story)
                    edit_lbl = write_context_expression(statenode.edit, hideid=True)
                    allcumulcopy = copy.deepcopy(all_cumul_edits[edit_lbl])
                    # Recheck relevance. I have to redo it instead of taking
                    # back the results from the previous round because adding
                    # nodes from allcumulcopy will change the cumul of the
                    # next rank.
                    relevantcumul = []
                    remainingcumul = []
                    #if edit_lbl == "<b>FES(act{t})</b>":
                    #if edit_lbl == "<b>BCR(Y177[1])</b>, <b>FYN(tyr_kin[1])</b>":
                    #    print("----")
                    #    for f in fullcumul:
                    #        print(f)
                    for cumulnode in fullcumul:
                        relevant = story.reachability_with_block(cumulnode,
                            statenode.reachable, statenode)
                        if relevant == True:
                            relevantcumul.append(cumulnode)
                            # Remove one cumulnode type from allcumulcopy.
                            # This is because I have to keep track of cumul
                            # nodes that are relevant in the current stories
                            # and remove them from allcumulcopy. Otherwise, any
                            # duplicate in fullcumul that is of a type found in
                            # allcumulcopy will be both kept as relevant.
                            for i in range(len(allcumulcopy)):
                                are_same = compare_states(allcumulcopy[i],
                                                          cumulnode.edit,
                                                          ignorevalue=False,
                                                          ignoreid=True)
                                if are_same == True:
                                    del(allcumulcopy[i])
                                    break
                        else:
                            remainingcumul.append(cumulnode)
                    #if edit_lbl == "<b>BCR(Y177[1])</b>, <b>FYN(tyr_kin[1])</b>":
                    #    print("----", story.filename)
                    #    for r in relevantcumul:
                    #        print(r)
                    # Put remaining nodes in relevant nodes if they are found
                    # in allcumulcopy.
                    for remainingnode in remainingcumul:
                        for i in range(len(allcumulcopy)):
                            are_same = compare_states(allcumulcopy[i],
                                                      remainingnode.edit,
                                                      ignorevalue=False,
                                                      ignoreid=True)
                            if are_same == True:
                                relevantcumul.append(remainingnode)
                                del(allcumulcopy[i])
                                break
                    statenode.cumulnodes = relevantcumul
                    # Build current state node context from the state of
                    # all the relevant_nodes.
                    full_state = copy.deepcopy(statenode.edit)
                    for cumulnode in statenode.cumulnodes: # + neighbors:
                        for agent in cumulnode.edit:
                            context_agent = copy.deepcopy(agent)
                            if context_agent["type"] == None:
                                for context_site in context_agent["sites"]:
                                    context_site["type"] = "context"
                            elif context_agent["type"] != None:
                                context_agent["type"] = "context"
                            full_state.append(context_agent)
                    statenode.state = group_sites_by_agent(full_state)
                    lbl = write_context_expression(statenode.state)
                    statenode.label = lbl
        for statenode in story.statenodes:
            if statenode in story.rule_outputs:
                statenode.introstate = False
            else:
                statenode.state = statenode.edit
                statenode.introstate = True
        story.create_hyperedges()
        story.align_vertical()
    # Find ubiquitous context sites.
    ubiquitous_edits = {}
    for story in stories:
        for statenode in story.rule_outputs:
            # Get the edit of each cumul node from the current state node.
            # At that point, statenode.cumulnodes is the relevantcumul computed
            # just before. Since agent ids will not be considered, there may be
            # many equivalent edits.
            cumul_edits = []
            for cn in statenode.cumulnodes:
                cumul_edits.append(cn.edit)
            # Get edit label for the current state node.
            edit_lbl = write_context_expression(statenode.edit, hideid=True)
            # If this edit_lbl was not seen before, assign all cumul_edits.
            if edit_lbl not in ubiquitous_edits.keys():
                ubiquitous_edits[edit_lbl] = cumul_edits
            # Else, if edit_lbl was already seen, remove cumul_edits that are
            # not present in the current statenode.
            else:
                ubi = ubiquitous_edits[edit_lbl]
                if len(ubi) > 0:
                    #cumulcopy = copy.deepcopy(cumul_edits)
                    ubiquitous_to_remove = []
                    for i in range(len(ubi)):
                        ubiquitous_found = False
                        for j in range(len(cumul_edits)):
                            are_same = compare_states(ubi[i],
                                                      cumul_edits[j],
                                                      ignorevalue=False,
                                                      ignoreid=True)
                            if are_same == True:
                                ubiquitous_found = True
                                del(cumul_edits[j])
                                break
                        if ubiquitous_found == False:
                            ubiquitous_to_remove.insert(0, i)
                    for i in ubiquitous_to_remove:
                        del(ubiquitous_edits[edit_lbl][i])
    # For each statenode, remove cumulnodes that contain ubiquitous context.
    for story in stories:
        for statenode in story.rule_outputs:
            edit_lbl = write_context_expression(statenode.edit, hideid=True)
            ubicopy = copy.deepcopy(ubiquitous_edits[edit_lbl])
            cumulnodes_to_remove = []
            for i in range(len(statenode.cumulnodes)):
                for j in range(len(ubicopy)):
                    are_same = compare_states(statenode.cumulnodes[i].edit,
                                              ubicopy[j],
                                              ignorevalue=False,
                                              ignoreid=True)
                    if are_same == True:
                        cumulnodes_to_remove.insert(0, i)
                        del(ubicopy[j])
                        break
            for i in cumulnodes_to_remove:
                del(statenode.cumulnodes[i])
            # Build current state node context without ubiquitous context. 
            full_state = copy.deepcopy(statenode.edit)
            for cumulnode in statenode.cumulnodes: # + neighbors:
                for agent in cumulnode.edit:
                    context_agent = copy.deepcopy(agent)
                    if context_agent["type"] == None:
                        for context_site in context_agent["sites"]:
                            context_site["type"] = "context"
                    elif context_agent["type"] != None:
                        context_agent["type"] = "context"
                    full_state.append(context_agent)
            statenode.state = group_sites_by_agent(full_state)
            lbl = write_context_expression(statenode.state)
            statenode.label = lbl

    # Write stories with harmonized context on state nodes.
    for i in range(len(stories)):
        stories[i].filename = "harmonized-{}.dot".format(i+1)
    for story in stories:
        story.build_dot_file(showintro, addedgelabels, showedgelabels,
                             edgeid, edgeocc, edgeprob, statstype, weightedges)
        output_path = "{}/tmp/{}".format(eoi, story.filename)
        outfile = open(output_path, "w")
        outfile.write(story.dot_file)
        outfile.close()

    # ----------

    ## Add parallel context to states.
    #for story in stories:
    #    print(story.filename)
    #    # 0) For each event nodes, build a list of incoming states
    #    #    (excluding intro states).
    #    for eventnode in story.eventnodes:
    #        incoming_states = []
    #        for edge in story.causaledges:
    #            if edge.secondary == False and edge.target == eventnode:
    #                if edge.source.introstate == False:
    #                    incoming_states.append(edge.source)
    #        eventnode.incoming_states = incoming_states
    #    # For each state node ...
    #    for statenode in story.statenodes:
    #        # 1) Get all downstream event nodes.
    #        downstream_paths = story.follow_edges("down",
    #                                              statenode,
    #                                              ignore_conflict=False)
    #        downstream_nodes = [statenode]
    #        for path in downstream_paths:
    #            for node in path:
    #                #if isinstance(node, EventNode):
    #                if node not in downstream_nodes:
    #                    downstream_nodes.append(node)
    #        # 2) Keep only event nodes that have agents from the current state
    #        #    node in their test.
    #        relevant_events = []
    #        for downstream_node in downstream_nodes:
    #            if isinstance(downstream_node, EventNode):
    #                keep_event = False
    #                for test in downstream_node.tests:
    #                    for test_agent in test:
    #                        for state_agent in statenode.state:
    #                            if test_agent["name"] == state_agent["name"]:
    #                                if test_agent["id"] == state_agent["id"]:
    #                                    keep_event = True
    #                                    break
    #                if keep_event == True:
    #                    relevant_events.append(downstream_node)
    #        # 3) Check which relevant event nodes have many incoming states.
    #        relevant_statenodes = []
    #        for relevant_event in relevant_events:
    #            if len(relevant_event.incoming_states) > 1:
    #                # 4) Check which incoming states are not in
    #                #    downstream nodes.
    #                for incoming_state in relevant_event.incoming_states:
    #                    if incoming_state not in downstream_nodes:
    #                        relevant_statenodes.append(incoming_state)
    #        # 5) Add the state of relevant state nodes as parallel context of
    #        #    the current state node.
    #        statenode.parallel_state = []
    #        for relevant_statenode in relevant_statenodes:
    #            relevant_state = relevant_statenode.state
    #            for agent in relevant_state:
    #                para_agent = copy.deepcopy(agent)
    #                if para_agent["type"] == None:
    #                    for para_site in para_agent["sites"]:
    #                        para_site["type"] = "parallel"
    #                elif para_agent["type"] != None:
    #                    para_agent["type"] = "parallel"
    #                statenode.parallel_state.append(para_agent)
    #        # 5.1) Remove any parallel site that is already
    #        #      present from upstream context.
    #        for para_agent in statenode.parallel_state:
    #            sites_to_remove = []
    #            for j in range(len(para_agent["sites"])):
    #                para_site = para_agent["sites"][j]
    #                site_found = False
    #                # Check if this site is already in statenode.state.
    #                for agent in statenode.state:
    #                    for site in agent["sites"]:
    #                        are_same, ig = compare_sites(site, para_site,
    #                                                     ignoretype=True)
    #                        if are_same == True:
    #                            site_found = True
    #                            break
    #                if site_found == True:
    #                    sites_to_remove.insert(0, j)
    #            for j in sites_to_remove:
    #                del(para_agent["sites"][j])
    #    # 6) Put parallel context information into the state of each
    #    #    state node.
    #    for statenode in story.statenodes:
    #        newstate = statenode.state + statenode.parallel_state
    #        statenode.state = group_sites_by_agent(newstate)
    #        lbl = write_context_expression(statenode.state)
    #        statenode.label = lbl

    #        #if statenode.nodeid == "state105":    
    #        #    for agent in statenode.state:
    #        #        print(agent)
    #        #    print("===")
    #        #    for relevant_event in relevant_events:
    #        #        print(relevant_event.label)
    #            
    #

    ## Write stories with parallel context.
    #for i in range(len(stories)):
    #    stories[i].filename = "parallel-{}.dot".format(i+1)
    #for story in stories:
    #    story.build_dot_file(showintro, addedgelabels, showedgelabels,
    #                         edgeid, edgeocc, edgeprob, statstype, weightedges)
    #    output_path = "{}/tmp/{}".format(eoi, story.filename)
    #    outfile = open(output_path, "w")
    #    outfile.write(story.dot_file)
    #    outfile.close()

    # ----------

    # Distinguish events that are applications of a same rule but
    # in a different context.
    rule_names = []
    res_states = []
    for story in stories:
        for eventnode in story.eventnodes:
            # eventnode.output is a list of states, one per target state node.
            eventnode.output = get_output_of_node(story, eventnode)
            # Assign output sites to rule name.
            if eventnode.label not in rule_names:
                rule_names.append(eventnode.label)
                res_states.append([{"output": eventnode.output,
                                    "occurrence": 1}])
            else:
                rule_index = rule_names.index(eventnode.label)
                output_found = False
                for i in range(len(res_states[rule_index])):
                    output = res_states[rule_index][i]
                    # check if output["output"] is the same as the
                    # current output_sites.
                    are_same = compare_outputs(output["output"],
                                               eventnode.output)
                    if are_same == True:
                        output_found = True
                        res_states[rule_index][i]["occurrence"] += 1
                if output_found == False:
                   res_states[rule_index].append({"output": eventnode.output,
                                                  "occurrence": 1})
    # Sort outputs by occurrence for each rule.
    sorted_res_states = []
    for rule_states in res_states:
        sorted_states = sorted(rule_states, key=lambda x: x["occurrence"],
                               reverse=True)
        sorted_res_states.append(sorted_states)
    # Assign new labels to rules that occur in different context.
    for story in stories:
        for eventnode in story.eventnodes:
            if eventnode.intro == False:
                rule_index = rule_names.index(eventnode.label)
                if len(sorted_res_states[rule_index]) > 1:
                    output_states = get_output_of_node(story, eventnode)
                    output_set = None
                    for i in range(len(sorted_res_states[rule_index])):
                        output = sorted_res_states[rule_index][i]
                        are_same = compare_outputs(output["output"],
                                                   output_states)
                        if are_same:
                            output_set = i
                    if output_set != None:
                        prime = " "
                        for i in range(output_set):
                            prime += "'"
                        eventnode.label = "{}{}".format(eventnode.label, prime)
                    else:
                        raise ValueError("Output nodes not found for node {}, "
                                         "label '{}'.".format(eventnode.nodeid,
                                         eventnode.label))
    ## Check which state nodes should be distinguished using context
    ## (since they appear in different contexts)
    ## !!! Remove this part once ubiquitous context is removed !!!
    #editset = []
    #statesets = []
    #for story in stories:
    #    story.rule_outputs = []
    #    for edge in story.causaledges:
    #        if isinstance(edge.source, EventNode):
    #            if edge.source.intro == False:
    #                if isinstance(edge.target, StateNode):
    #                    story.rule_outputs.append(edge.target)
    #    #for statenode in story.statenodes:
    #    for statenode in story.rule_outputs:
    #        edit_found = False
    #        for i in range(len(editset)):
    #            are_same = compare_states(statenode.edit, editset[i],
    #                                      ignorevalue=False, ignoreid=True)
    #            if are_same == True:
    #               edit_found = True
    #               # Check if the full state is in this stateset.
    #               state_found = False
    #               for j in range(len(statesets[i])):
    #                   same = compare_states(statenode.state, statesets[i][j],
    #                                         ignorevalue=False, ignoreid=True)
    #                   if same == True:
    #                       state_found = True
    #                       break
    #               if state_found == False:
    #                   statesets[i].append(statenode.state)
    #        if edit_found == False:
    #            editset.append(statenode.edit)
    #            statesets.append([statenode.state])
    #for story in stories:
    #    #story.rule_outputs = []
    #    #for edge in story.causaledges:
    #    #    if isinstance(edge.source, EventNode):
    #    #        if edge.source.intro == False:
    #    #            if isinstance(edge.target, StateNode):
    #    #                story.rule_outputs.append(edge.target)
    #    for statenode in story.rule_outputs:
    #        for i in range(len(editset)):
    #            are_same = compare_states(statenode.edit, editset[i],
    #                                      ignorevalue=False, ignoreid=True)
    #            if are_same == True:
    #                if len(statesets[i]) > 1:
    #                    statenode.differentiate = True
    #                else:
    #                    statenode.differentiate = False
    #                break

    # Ensure that all state nodes which have the same state have the same
    # label across all stories (since the order in which the agents are
    # written on the label is arbitrary and may depend on the order with
    # which the context elements were gathered along the story).
    possible_edits = []
    standard_edits = []
    possible_states = []
    standard_labels = []
    for story in stories:
        story.intro_outputs = []
        for edge in story.causaledges:
            if isinstance(edge.source, EventNode):
                if edge.source.intro == True:
                    if isinstance(edge.target, StateNode):
                        story.intro_outputs.append(edge.target)
        for statenode in story.statenodes:
            # Build standard edit label.
            edit_index = None
            for i in range(len(possible_edits)):
                are_same = compare_states(statenode.edit, possible_edits[i],
                                          ignorevalue=False, ignoreid=True)
                if are_same == True:
                    edit_index = i
                    break
            # Define standard edit label the first time that an edit is found.
            if edit_index == None:
                possible_edits.append(statenode.edit)
                standard_edit = write_context_expression(statenode.edit,
                                                          hidevalue=False,
                                                          hideid=True)
                standard_edits.append(standard_edit)
                statenode.stdedit = standard_edit
                if statenode in story.intro_outputs:
                    statenode.label = standard_edit
            # Otherwise assign the already chosen standard edit label.
            else:
                statenode.stdedit = standard_edits[edit_index]
                if statenode in story.intro_outputs:
                    statenode.label = standard_edits[edit_index]
            # Build standard state label.
            if statenode in story.rule_outputs:
                state_index = None
                for j in range(len(possible_states)):
                    are_same = compare_states(statenode.state,
                                              possible_states[j],
                                              ignorevalue=False, ignoreid=True)
                    if are_same == True:
                        state_index = j
                        break
                # Define standard label the first time that a state is found.
                if state_index == None:
                    possible_states.append(statenode.state)
                    standard_label = write_context_expression(statenode.state,
                                                              hidevalue=False,
                                                              hideid=True)
                    standard_labels.append(standard_label)
                    statenode.label = standard_label
                # Otherwise assign the already chosen standard label.
                else:
                    statenode.label = standard_labels[state_index]

    # Writes stories with distinguished events.
    for i in range(len(stories)):
        stories[i].filename = "context-{}.dot".format(i+1)
    for story in stories:
        story.build_dot_file(showintro, addedgelabels, showedgelabels,
                             edgeid, edgeocc, edgeprob, statstype, weightedges)
        output_path = "{}/tmp/{}".format(eoi, story.filename)
        outfile = open(output_path, "w")
        outfile.write(story.dot_file)
        outfile.close()

#    # Simplify state node labels by removing the part that is common to all
#    # state labels with a same edit.
#
#    # 1) Gather lists of all nodes with the same edit.
#    edit_set = []
#    same_edit_nodes = []
#    for story in stories:
#        rule_outputs = []
#        for edge in story.causaledges:
#            if isinstance(edge.source, EventNode):
#                if edge.source.intro == False:
#                    if isinstance(edge.target, StateNode):
#                        rule_outputs.append(edge.target)
#        for statenode in rule_outputs:
#            edit_found = False
#            for i in range(len(edit_set)):
#                are_same = compare_states(statenode.edit, edit_set[i],
#                                          ignorevalue=False, ignoreid=True)
#                if are_same == True:
#                    edit_found = True
#                    same_edit_nodes[i].append(statenode)
#                    break
#            if edit_found == False:
#                edit_set.append(statenode.edit)
#                same_edit_nodes.append([statenode])
#    # 2) For each group of nodes, take the state of the first node and find
#    #    equivalent agents among the state of all other nodes
#    for group_of_nodes in same_edit_nodes:
#        print("----")
#        print(group_of_nodes[0].label)
#        similitudes = []
#        if len(group_of_nodes) > 1 :
#            first_node = group_of_nodes[0]
#            for statenode in group_of_nodes[1:]:
#                simil = []
#                for agent_first in first_node.state:
#                    compare = []
#                    for agent_other in statenode.state:
#                        value = 0.0
#                        if agent_first["name"] == agent_other["name"]:
#                            value += 1
#                            # Check sites
#                            for site_first in agent_first["sites"]:
#                                for site_other in agent_other["sites"]:
#                                    # ?? Check the partner first ??
#                                    if site_first["name"] == site_other["name"]:
#                                        value += 1
#                                        if site_first["value"] != None:
#                                            if site_first["value"] == site_other["value"]:
#                                                value += 1
#                                        if site_first["bond"] != None:
#                                            if site_first["bond"] == "." and site_other["bond"] == ".":
#                                                value += 1
#                                            else:
#                                                if site_first["bond"]["partner"]["agentname"] == site_other["bond"]["partner"]["agentname"]:
#                                                    value += 0.5
#                                                    if site_first["bond"]["partner"]["sitename"] == site_other["bond"]["partner"]["sitename"]:
#                                                        value += 0.5
#                        compare.append(value)
#                    simil.append(compare)
#                similitudes.append(simil)
#            for simil in similitudes:
#                print(simil)
#
#
#
#    # 3) For each node in a group, remove all sites from agents for which
#    #    equivalent agents with the forst node have been found.
#
#    # Writes stories with distinguished events.
#    for i in range(len(stories)):
#        stories[i].filename = "simplified-{}.dot".format(i+1)
#    for story in stories:
#        story.build_dot_file(showintro, addedgelabels, showedgelabels,
#                             edgeid, edgeocc, edgeprob, statstype, weightedges)
#        output_path = "{}/tmp/{}".format(eoi, story.filename)
#        outfile = open(output_path, "w")
#        outfile.write(story.dot_file)
#        outfile.close()

    # Build dual stories.
    for story in stories:
        print(story.filename)
        # Remove unused state nodes from intros.
        nodes_to_remove = []
        edges_to_remove = []
        for i in range(len(story.statenodes)):
            statenode = story.statenodes[i]
            #if statenode not in story.rule_outputs:
            outgoing_edges = []
            for edge in story.causaledges:
                if edge.source == statenode:
                    outgoing_edges.append(edge)
            if len(outgoing_edges) == 0: # Remove state node
                nodes_to_remove.insert(0, i)
                for j in range(len(story.causaledges)):
                    edge = story.causaledges[j]
                    if edge.target == statenode:
                        if j not in edges_to_remove:
                            edges_to_remove.insert(0, j)               
        for i in nodes_to_remove:
            del(story.statenodes[i])
        for j in edges_to_remove:
            del(story.causaledges[j])
        # Remove secondary edges.
        secondary_to_remove = []
        for i in range(len(story.causaledges)):
            edge = story.causaledges[i]
            if edge.secondary == True:
                secondary_to_remove.insert(0,i)
        for i in secondary_to_remove:
            del(story.causaledges[i])
        # Remove intro edits.
        nodes_to_remove = []
        edges_to_remove = []
        for j in range(len(story.causaledges)):
            edge1 = story.causaledges[j]
            if isinstance(edge1.source, EventNode):
                if edge1.source.intro == True:
                    edges_to_remove.append(j)
                    introedit = edge1.target
                    for k in range(len(story.statenodes)):
                       if story.statenodes[k] == introedit:
                           nodes_to_remove.insert(0, k) 
                    for l in range(len(story.causaledges)):
                        edge2 = story.causaledges[l]
                        if edge2.source == introedit:
                            edges_to_remove.append(l)
                            # Check if edge already exists.
                            already_exists = False
                            for e in story.causaledges:
                                if e.source == edge1.source:
                                    if e.target == edge2.target:
                                        already_exists = True
                                        break
                            if already_exists == False:
                                new_edge = CausalEdge(edge1.source,
                                                      edge2.target)
                                story.causaledges.append(new_edge)
        sorted_nodes_to_remove = sorted(nodes_to_remove, reverse=True)
        sorted_edges_to_remove = sorted(edges_to_remove, reverse=True)
        for k in sorted_nodes_to_remove:
            del(story.statenodes[k])
        for j in sorted_edges_to_remove:
            del(story.causaledges[j])
        ## Lower the rank of al non-intro nodes by 0.5.
        #for eventnode in story.eventnodes:
        #    if eventnode.intro == False:
        #        eventnode.rank = eventnode.rank - 0.5
        #for statenode in story.statenodes:
        #    statenode.rank = statenode.rank - 0.5

        # Increase the rank of intro nodes by 0.5
        for eventnode in story.eventnodes:
            if eventnode.intro == True:
                eventnode.rank = eventnode.rank + 0.5
        story.get_maxrank()
        story.create_hyperedges()
        story.align_vertical()

    # Writes stories with event and state nodes.
    for i in range(len(stories)):
        stories[i].filename = "dualstory-{}.dot".format(i+1)
        #stories[i].processed = True
    for story in stories:
        story.build_dot_file(showintro, addedgelabels, showedgelabels,
                             edgeid, edgeocc, edgeprob, statstype, weightedges)
        output_path = "{}/tmp/{}".format(eoi, story.filename)
        outfile = open(output_path, "w")
        outfile.write(story.dot_file)
        outfile.close()

    # Write stories with event nodes only.
    # Use a copy of stories to leave the original stories for the subsequent
    # states version.
    for story in stories:
        story.clear_adjacency()
    storiescopy = copy.deepcopy(stories)
    for story in storiescopy:
        fuse_multiple_outputs_with_causaledges(story) 
        keep_events_only(story)
    for i in range(len(storiescopy)):
        storiescopy[i].filename = "modstory-{}.dot".format(i+1)
    for story in storiescopy:
        story.create_hyperedges()
        story.align_vertical()
        story.get_maxrank()
        story.build_dot_file(showintro, addedgelabels, showedgelabels,
                             edgeid, edgeocc, edgeprob, statstype, weightedges)
        output_path = "{}/tmp/{}".format(eoi, story.filename)
        outfile = open(output_path, "w")
        outfile.write(story.dot_file)
        outfile.close()

    # Write stories with state nodes only.
    for story in stories:
        story.build_adjacency(hyper=True)
        keep_states_only(story)
    for i in range(len(stories)):
        stories[i].filename = "statestory-{}.dot".format(i+1)
    for story in stories:
        story.build_dot_file(showintro, addedgelabels, showedgelabels,
                             edgeid, edgeocc, edgeprob, statstype, weightedges)
        output_path = "{}/tmp/{}".format(eoi, story.filename)
        outfile = open(output_path, "w")
        outfile.write(story.dot_file)
        outfile.close()

    # Write state dict in json format.
    for i in range(len(stories)):
        statedict = {"states": {}, "edits": {}}
        for statenode in stories[i].statenodes:
            statedict["states"][statenode.nodeid] = statenode.state
            statedict["edits"][statenode.nodeid] = statenode.edit
        json_path = "{}/tmp/statefile-{}.json".format(eoi, i+1)
        statefile = open(json_path, "w")
        json.dump(statedict, statefile)
        statefile.close()


def keep_events_only(graph):
    """
    Remove state nodes in a causal graph that contains both events and states
    nodes.
    This must be done for each story before folding because folding the
    dualstories loses some information about edge statistics.
    """

    # Rebuild adjacency lists using causal edges.
    graph.build_adjacency()
    # Connect event nodes together.
    edges_to_remove = []
    for edge in graph.causaledges:
        if isinstance(edge.target, EventNode) and edge.source.intro == False:
            up_edge = edge.source.incoming[0]
            edge.source = up_edge.source
            edge_index = graph.causaledges.index(up_edge)
            if edge_index not in edges_to_remove:
                edges_to_remove.append(edge_index)
    sorted_edges = sorted(edges_to_remove, reverse=True)
    for i in sorted_edges:
        del(graph.causaledges[i])
    # Remove all state nodes.
    graph.statenodes = []
    # Move intro nodes up half a rank.
    for eventnode in graph.eventnodes:
        if eventnode.intro == True:
            eventnode.rank = eventnode.rank - 0.5           


def keep_states_only(graph):
    """
    Remove event nodes in a causal graph that contains both events and states
    nodes.
    This must be done on the pathway after the folding of dualstories to keep
    separate edges when they involve different events (rules).
    """

    # The first nodes become the output of the event
    # nodes which had first=True.
    for hyperedge in graph.hyperedges:
        first_in_sources = False
        for source in hyperedge.sources:
            if source.first == True:
                first_in_sources = True
                break
        if first_in_sources == True:
            hyperedge.target.first = True
    # Connect state nodes together.
    hyperedges_to_remove = []
    for hyperedge in graph.hyperedges:
        if hyperedge.target.label != graph.eoi:
            if isinstance(hyperedge.target, EventNode):
                #outgoing_edges = []
                #for hyperedge2 in graph.hyperedges:
                #    if hyperedge.target in hyperedge2.sources:
                #        outgoing_edges.append(hyperedge2)
                outgoing_edges = hyperedge.target.outgoing
                if len(outgoing_edges) == 1:
                    # Reconnect the hyperedge targeting the event
                    # to now target the output state.
                    for subedge in hyperedge.edgelist:
                        subedge.target = outgoing_edges[0].target
                    hyperedge.target = outgoing_edges[0].target
                    he_index = graph.hyperedges.index(outgoing_edges[0])
                    if he_index not in hyperedges_to_remove:
                        hyperedges_to_remove.append(he_index)
                elif len(outgoing_edges) > 1:
                    # Shrink the event node. Nodes that are shrank hide
                    # their label and are shown as a small white circle.
                    hyperedge.target.shrink = True
                    # Also change outgoing edges type to conflict of all
                    # incoming edge are of type conflict.
                    in_all_conflict = True
                    for subedge in hyperedge.edgelist:
                        if subedge.relationtype != "conflict":
                            in_all_conflict = False
                            break
                    if in_all_conflict == True:
                        for out_hyperedge in outgoing_edges:
                            for subedge in out_hyperedge.edgelist:
                                subedge.relationtype = "conflict"
    sorted_hyperedges = sorted(hyperedges_to_remove, reverse=True)
    for i in sorted_hyperedges:
        del(graph.hyperedges[i])
    # Remove events.
    events_to_remove = []
    for j in range(len(graph.eventnodes)):
        if graph.eventnodes[j].intro == False:
            lbl = graph.eventnodes[j].label
            shrk = graph.eventnodes[j].shrink
            if lbl != graph.eoi and shrk == False:
                events_to_remove.insert(0, j)
    for j in events_to_remove:
        del(graph.eventnodes[j])
    # Move ranks up by one half.
    for node in graph.eventnodes + graph.statenodes:
        if node.label != graph.eoi:
            node.rank = node.rank - 0.5
    graph.get_maxrank()
    # Find essential edges, edges that have an intro as source that
    # is not mentioned in the edit of the target.
    for hyperedge in graph.hyperedges:
        for subedge in hyperedge.edgelist:
            subedge.essential = False
            if subedge.source.intro == True:
                if isinstance(hyperedge.target, EventNode):
                    target_edits = hyperedge.target.edits
                else:
                    target_edits = [hyperedge.target.edit]
                # All agents from the intro edits must be seen in
                # the edit of the target.
                intro_agent_in_target = True
                for intro_edit in subedge.source.edits:
                    for intro_agent in intro_edit:
                        agent_found = False
                        for target_edit in target_edits:
                            for target_agent in target_edit:
                                if target_agent["name"] == intro_agent["name"]:
                                    if target_agent["id"] == intro_agent["id"]:
                                        agent_found = True
                                        break
                        if agent_found == False:
                            subedge.essential = True
                            break
    # Remove shrank nodes if needed.
    remove_shrank_nodes(graph)


def remove_suffix(idstr):
    """ Remove the trailing letters on a id string. """

    idnum = ""
    for char in str(idstr):
        if char.isdigit():
            idnum += char

    return idnum


def change_agent_id(agent, id_changes):
    """ Change an agent id along with its appearances in sites and partner. """

    is_changed = False
    if str(agent["id"]) in id_changes.keys():
        agent["id"] = id_changes[str(agent["id"])]
        is_changed = True
    for site in agent["sites"]:
        if str(site["agentid"]) in id_changes.keys():
            site["agentid"] = id_changes[str(site["agentid"])]
            is_changed = True
        if isinstance(site["bond"], dict):
            partnerid = str(site["bond"]["partner"]["agentid"])
            if partnerid in id_changes.keys():
                site["bond"]["partner"]["agentid"] = id_changes[partnerid]
                is_changed = True

    return is_changed


def get_fullcumul(statenode, story):
    """ Assign the cumulative context from the past to given state node. """

    # Find immediate upstream state nodes.
    # Also get neighbors of state node.
    src_rule = statenode.incoming[0].source
    upstream_nodes = []
    neighbors = []
    for inedge in src_rule.incoming:
        if inedge.secondary == False:
            upstream_nodes.append(inedge.source)
    for outedge in src_rule.outgoing:
        if outedge.secondary == False:
            neighbors.append(outedge.source)
    # Add neighbors of upstream state nodes.
    up_rules = []
    for up_node in upstream_nodes:
        for up_edge in up_node.incoming:
            if up_edge.source not in up_rules:
                up_rules.append(up_edge.source)
    upstream_neighbors = []
    for up_rule in up_rules:
        for down_edge in up_rule.outgoing:
            if down_edge.target not in upstream_nodes:
                if down_edge.target not in upstream_neighbors:
                    upstream_neighbors.append(down_edge.target)
    upstream_nodes += upstream_neighbors
    # Add the immediate upstream nodes and their cumulnodes
    # to the cumulnodes of the current statenode.
    fullcumul = []
    for upstream_node in upstream_nodes:
        fullcumul.append(upstream_node)
        for up_cumulnode in upstream_node.cumulnodes:
            if up_cumulnode not in fullcumul:
                fullcumul.append(up_cumulnode)

    return fullcumul # , src_rule
   

def remove_shrank_nodes(story):
    """
    For each hyperedge that targets a shrank node and has only one subedge,
    create new edges from the source of that hyperedge to the targets of the
    shrank node. Remove the hyperedge in question. If this results in the
    shrank node having no incoming edge, remove the shrank node.
    """

    # First, reset incoming and outgoing edges information.
    for node in story.eventnodes + story.statenodes:
        node.incoming = []
        node.outgoing = []
        node.reachable = []
    for hyperedge in story.hyperedges:
        if len(hyperedge.edgelist) == 1:
            if isinstance(hyperedge.target, EventNode):
                if hyperedge.target.shrink == True:
                    shrink_targets = []
                    for hyperedge2 in story.hyperedges:
                        if hyperedge.target in hyperedge2.sources:
                            shrink_targets.append(hyperedge2.target)
                    first = True
                    for shrink_target in shrink_targets:
                        if first == False:
                            #print("edgelist", hyperedge.edgelist)
                            #print("weight", hyperedge.weight)
                            #print("layout_weight", hyperedge.layout_weight)
                            #print("rel_wei", hyperedge.rel_wei)
                            #print("occurrence", hyperedge.occurrence)
                            #print("rel_occ", hyperedge.rel_occ)
                            #print("number", hyperedge.number)
                            #print("rel_num", hyperedge.rel_num)
                            #print("relationtype", hyperedge.relationtype)
                            #print("color", hyperedge.color)
                            #print("midcolor", hyperedge.midcolor)
                            #print("secondary", hyperedge.secondary)
                            #print("underlying", hyperedge.underlying)
                            #print("reverse", hyperedge.reverse)
                            #print("labelcarrier", hyperedge.labelcarrier)
                            #print("indicator", hyperedge.indicator)
                            #print("hyperid", hyperedge.hyperid)
                            #print("pos", hyperedge.pos)
                            #print("labelpos", hyperedge.labelpos)
                            #print("overridewidth", hyperedge.overridewidth)
                            #print("overridelabel", hyperedge.overridelabel)
                            #print("essential", hyperedge.essential)
                            #print("target", hyperedge.target)
                            #print("sources", hyperedge.sources)
                            #print(hyperedge.sources[0].incoming,
                            #       hyperedge.sources[0].outgoing)
                            #print(hyperedge.edgelist[0].source.__dict__.keys())
                            #print(hyperedge.edgelist[0].source.incoming,
                            #      hyperedge.edgelist[0].source.outgoing)
                            new_hedge = copy.deepcopy(hyperedge)
                            new_hedge.target = shrink_target
                            new_hedge.edgelist[0].target = shrink_target
                            story.hyperedges.append(new_hedge)
                        if first == True:
                            hyperedge.target = shrink_target
                            hyperedge.edgelist[0].target = shrink_target
                            first = False
    nodes_to_remove = []
    for i in range(len(story.eventnodes)):
        eventnode = story.eventnodes[i]
        if eventnode.shrink == True:
            incoming_edges = []
            outgoing_edges = []
            for j in range(len(story.hyperedges)):
                hyperedge = story.hyperedges[j]
                if hyperedge.target == eventnode:
                    incoming_edges.append(j)
                if eventnode in hyperedge.sources:
                    outgoing_edges.insert(0, j)
            if len(incoming_edges) == 0:
                # Remove that shrank node.
                #print(eventnode.label, outgoing_edges)
                for j in outgoing_edges:
                    del(story.hyperedges[j])
                nodes_to_remove.insert(0, i)
    for i in nodes_to_remove:
        del(story.eventnodes[i])


def state_from_action(signatures, action, bnd_num):
    """ Find the resulting state of an action from the trace file. """

    # state = [agents]
    # agent = {"name": ,"id": ,"sites": [sites], "type": (None if sites exist),
    #          "action": (None if sites exist)}
    # site = {"name": ,"bond": ,"value", "agentname": ,"agentid": ,
    #         "type": (edit or context), ("action": or "test:")}
    # bond = {"num": ,"partner": }
    # partner = {"agentname": ,"agentid": ,"sitename":}

    state = []
    if action[0] == 0: # Create (I only look at the agent, not the sites).
        ag_n = action[1][1]
        agid_n = action[1][0]
        entry = signatures[ag_n]
        agentname = entry["name"]
        agent = {"name": agentname, "id": agid_n, "sites": None, "action":0}
        state.append(agent)
    if action[0] == 1: # Mod_internal
        ag_n = action[1][0][1]
        agid_n = action[1][0][0]
        site_n = action[1][1]
        val_n = action[2]
        entry = signatures[ag_n]
        agentname = entry["name"]
        sitename = entry["decl"][site_n]["name"]
        value = entry["decl"][site_n]["decl"][0][val_n]["name"]
        site = {"name": sitename, "bond": None, "value": value,
                "agentname": agentname, "agentid": agid_n, "action": 1}
        agent = {"name": agentname, "id": agid_n, "sites": [site],
                 "action": None}
        state.append(agent)
    if action[0] == 2 or action[0] == 3: # Bind or Bind_to
        ag1_n = action[1][0][1]
        agid1_n = action[1][0][0]
        site1_n = action[1][1]
        entry1 = signatures[ag1_n]
        agentname1 = entry1["name"]
        sitename1 = entry1["decl"][site1_n]["name"]
        ag2_n = action[2][0][1]
        agid2_n = action[2][0][0]
        site2_n = action[2][1]
        entry2 = signatures[ag2_n]
        agentname2 = entry2["name"]
        sitename2 = entry2["decl"][site2_n]["name"]
        partner1 = {"agentname":agentname1, "agentid": agid1_n,
                    "sitename":sitename1}
        partner2 = {"agentname":agentname2, "agentid": agid2_n,
                    "sitename":sitename2}
        site1 = {"name": sitename1,
                 "bond": {"num": bnd_num, "partner": partner2},
                 "value": None, "agentname": agentname1, "agentid": agid1_n,
                 "action": action[0]}
        site2 = {"name": sitename2,
                 "bond": {"num": bnd_num, "partner": partner1},
                 "value": None, "agentname": agentname2, "agentid": agid2_n,
                 "action": action[0]}
        agent1 = {"name": agentname1, "id": agid1_n, "sites": [site1],
                  "action": None}
        agent2 = {"name": agentname2, "id": agid2_n, "sites": [site2],
                  "action": None}
        state.append(agent1)
        state.append(agent2)
        bnd_num += 1
    if action[0] == 4: # Free
        ag_n = action[1][0][1]
        agid_n = action[1][0][0]
        site_n = action[1][1]
        entry = signatures[ag_n]
        agentname = entry["name"]
        sitename = entry["decl"][site_n]["name"]
        site = {"name": sitename, "bond": ".", "value": None,
                "agentname": agentname, "agentid": agid_n, "action": 4}
        agent = {"name": agentname, "id": agid_n, "sites": [site],
                 "action": None}
        state.append(agent)
    #if action[0] == 5: # Remove (I still do not have any example).
    for site in state:
        site["type"] = "edit"
    for agent in state:
        if agent["sites"] != None:
            agent["type"] = None
            for site in agent["sites"]:
                site["type"] = "edit"
        elif agent["sites"] == None:
            agent["type"] = "edit"

    return state, bnd_num


def state_from_test(signatures, test, bnd_num):
    """ Find the required state of a test from the trace file. """

    state = []
    if test[0] == 0: # Is_Here
        ag_n = test[1][1]
        agid_n = test[1][0]
        entry = signatures[ag_n]
        agentname = entry["name"]
        agent = {"name": agentname, "id": agid_n, "sites": None, "test":0}
        state.append(agent)
    if test[0] == 1: # Has_Internal
        ag_n = test[1][0][1]
        agid_n = test[1][0][0]
        site_n = test[1][1]
        val_n = test[2]
        entry = signatures[ag_n]
        agentname = entry["name"]
        sitename = entry["decl"][site_n]["name"]
        value = entry["decl"][site_n]["decl"][0][val_n]["name"]
        site = {"name": sitename, "bond": None, "value": value,
                "agentname": agentname, "agentid": agid_n, "test": 1}
        agent = {"name": agentname, "id": agid_n, "sites": [site],
                 "test": None}
        state.append(agent)
    if test[0] == 2 or test[0] == 3: # Is_Free or Is_Bound.
        ag_n = test[1][0][1]
        agid_n = test[1][0][0]
        site_n = test[1][1]
        entry = signatures[ag_n]
        agentname = entry["name"]
        sitename = entry["decl"][site_n]["name"]
        if test[0] == 2:
            bnd = "."
        elif test[0] == 3:
            bnd = "_"
        site = {"name": sitename, "bond": bnd, "value": None,
                "agentname": agentname, "agentid": agid_n, "test": test[0]}
        agent = {"name": agentname, "id": agid_n, "sites": [site], 
                 "test": None}
        state.append(agent)
    #if test[0] == 4: # Has_Binding_type (No example yet).
    if test[0] == 5: # Is_Bound_to
        ag1_n = test[1][0][1]
        agid1_n = test[1][0][0]
        site1_n = test[1][1]
        entry1 = signatures[ag1_n]
        agentname1 = entry1["name"]
        sitename1 = entry1["decl"][site1_n]["name"]
        ag2_n = test[2][0][1]
        agid2_n = test[2][0][0]
        site2_n = test[2][1]
        entry2 = signatures[ag2_n]
        agentname2 = entry2["name"]
        sitename2 = entry2["decl"][site2_n]["name"]
        partner1 = {"agentname":agentname1, "agentid": agid1_n,
                    "sitename":sitename1}
        partner2 = {"agentname":agentname2, "agentid": agid2_n,
                    "sitename":sitename2}
        site1 = {"name": sitename1,
                 "bond": {"num": bnd_num, "partner": partner2},
                 "value": None, "agentname": agentname1, "agentid": agid1_n,
                 "test": test[0]}
        site2 = {"name": sitename2,
                 "bond": {"num": bnd_num, "partner": partner1},
                 "value": None, "agentname": agentname2, "agentid": agid2_n,
                 "test": test[0]}
        agent1 = {"name": agentname1, "id": agid1_n, "sites": [site1],
                  "test": None}
        agent2 = {"name": agentname2, "id": agid2_n, "sites": [site2],
                  "test": None}
        state.append(agent1)
        state.append(agent2)
        bnd_num += 1

    return state, bnd_num


def group_sites_by_agent(state):
    """ Group all the sites of a given state by agents. """

    new_state = []
    for agent in state:
        agent_found = False
        for seen_agent in new_state:
            if agent["name"] == seen_agent["name"]:
                if agent["id"] == seen_agent["id"]:
                    agent_found = True
                    for site in agent["sites"]:
                        seen_agent["sites"].append(site)
        if agent_found == False:
            new_state.append(agent)
    ## Also group bond and value of same sites together
    ## (i.e. A(x[.] x{p}) becomes A(x[.]{p}))
    ## ... This is a bad idea to do this at this point,
    ##leave it for when we write the label.
    #for agent in new_state:
    #    seen_sites = []
    #    for i in range(len(agent["sites"])):
    #        site = agent["sites"][i]
    #        if site["name"] not in seen_sites:
    #            seen_sites.append(site["name"])

    return new_state


def write_kappa_agent(agent, bond="num", hidevalue=False, hideid=False):
    """
    Write an agent as a string using Kappa language.
    The value of bond can be either 'num' or 'partner'.
    """

    agent_str = ""
    # Check if agent contains at least one edited site.
    edited = False
    if agent["type"] == "edit":
        edited = True
    for site in agent["sites"]:
        if site["type"] == "edit":
            edited = True
            break
    # Write agent.
    was_edit = False
    if edited == True:
        agent_str += "<b>"
        was_edit = True
    agent_str += "{}".format(agent["name"])
    if hideid == False:
        agent_str += ":{}".format(agent["id"])
    agent_str += "("
    # Write sites. Bond and value of a same site are written together
    # (i.e. A(x[.] x{p}) is written as A(x[.]{p}))
    first_site = True
    seen_sites = []
    for site in agent["sites"]:
        if site["type"] == "context" and was_edit == True:
            agent_str += "</b>"
            was_edit = False
        if first_site == False and site["name"] not in seen_sites:
            agent_str += "&nbsp;"
        else:
            first_site = False
        if site["name"] not in seen_sites:
            agent_str += "{}".format(site["name"])
            seen_sites.append(site["name"])
        # If site already seen and was edit, close the edit
        if hidevalue == False:
            if site["bond"] != None:
                agent_str += "["
                if isinstance(site["bond"], dict):
                    if bond == "num":
                        agent_str += "{}".format(site["bond"]["num"])
                    elif bond == "partner":
                        partner = site["bond"]["partner"]
                        agent_str += "{}.{}".format(partner["sitename"],
                                                     partner["agentname"])
                        if hideid == False:
                            agent_str += "{}".format(partner["agentid"])
                else:
                    agent_str += site["bond"]
                agent_str += "]"
            if site["value"] != None:
                agent_str += "{{{}}}".format(site["value"])
        elif hidevalue == True:
            if site["bond"] != None:
                agent_str += "[]"
            if site["value"] != None:
                agent_str += "{}"
    # Close agent parenthesis.
    agent_str += ")"
    if edited == True and was_edit == True:
        agent_str += "</b>"

    return agent_str


def write_kappa_site(site, bond="num", hidevalue=False, hideid=False):
    """
    Write an site as a string using Kappa language.
    The value of bond can be either 'num' or 'partner'.
    A site is written as an agent containing a single site.
    """

    site_str = ""
    # Write agentname.
    site_str = "{}".format(site["agentname"])
    if hideid == False:
        site_str += ":{}".format(site["agentid"])
    site_str += "("
    # Write sites.
    site_str += "{}".format(site["name"])
    if hidevalue == False:
        if site["bond"] != None:
            site_str += "["
            if isinstance(site["bond"], dict):
                if bond == "num":
                    site_str += "{}".format(site["bond"]["num"])
                elif bond == "partner":
                    partner = site["bond"]["partner"]
                    site_str += "{}.{}".format(partner["sitename"],
                                                 partner["agentname"])
                    if hideid == False:
                        site_str += "{}".format(partner["agentid"])
            else:
                site_str += site["bond"]
            site_str += "]"
            #if site["bond"] == ".":
            #    site_str += "[.]"
            #else:
            #    if bond == "num":
            #        site_str += "{}".format(site["bond"]["num"])
            #    elif bond == "partner":
            #        partner = site["bond"]["partner"]
            #        site_str += "{}.{}".format(partner["sitename"],
            #                                     partner["agentname"])
            #        if hideid == False:
            #            site_str += "{}".format(partner["agentid"])
        if site["value"] != None:
            site_str += "{{{}}}".format(site["value"])
    elif hidevalue == True:
        if site["bond"] != None:
            site_str += "[]"
        if site["value"] != None:
            site_str += "{}"
    # Close agent parenthesis.
    site_str += ")"

    return site_str


def edit_vs_test(edit, test, return_correspondances=False):
    """
    Determine if what is produced by an edit can be used by a given test.
    The test may contain wildcards _ or #.
    """


    # Gather edit and test sites. Keep agents without sites in a separate list.
    edit_sites = []
    test_sites = []
    edit_empty = []
    test_empty = []
    for edit_agent in edit:
        if edit_agent["sites"] != None:
            for edit_site in edit_agent["sites"]:
                edit_sites.append(edit_site)
        else:
            edit_empty.append(edit_agent)
    for test_agent in test:
        if test_agent["sites"] != None:
            for test_site in test_agent["sites"]:
                test_sites.append(test_site)
        else:
            test_empty.append(test_agent)
    # Compare empty agents.
    ag_related = True
    ag_correspondances = []
    test_indexes = list(range(len(test_empty)))
    for edit_agent in edit_empty:
        for i in test_indexes:
            test_agent = test_empty[i]
            same_name = edit_agent["name"] == test_agent["name"]
            same_id = edit_agent["id"] == test_agent["id"]
            if same_name and same_id:
                are_same = True
                ag_correspondances.append(i)
                test_indexes.remove(i)
                break
            else:
                are_same = False
        if are_same == False:
            ag_related = False
            break
    if len(test_indexes) > 0:
        ag_related = False        
    # Compare sites.
    # Make a first pass to check which edit sites should
    # be ignored because of wildcards in tests.
    edits_to_ignore = []
    test_indexes = list(range(len(test_sites)))
    for edit_site in edit_sites:
        for i in test_indexes:
            test_site = test_sites[i]
            are_rel, ignore_partner = edit_vs_test_sites(edit_site, test_site)
            if ignore_partner == True:
                # Find the site which is bound to current edit_site.
                partner = edit_site["bond"]["partner"]
                partner_found = False
                for j in range(len(edit_sites)):
                    edit_site2 = edit_sites[j]
                    same_site = edit_site2["name"] == partner["sitename"]
                    same_ag = edit_site2["agentname"] == partner["agentname"]
                    same_agid = edit_site2["agentid"] == partner["agentid"]
                    if same_site and same_ag and same_agid:
                        partner_found = True
                        edits_to_ignore.append(j)
                        break
                if partner_found == False:
                    raise ValueError("Could not find partner of site {}:{}({})"
                        .format(edit_site["agentname"], edit_site["agentid"],
                                edit_site["name"]))
    # Now compare sites, ignoring edit sites determined above.
    site_related = True
    site_correspondances = []
    for j in range(len(edit_sites)):
        if j not in edits_to_ignore:
            edit_site = edit_sites[j]
            for i in test_indexes:
                test_site = test_sites[i]
                are_rel, ignore_partner = edit_vs_test_sites(edit_site, test_site)
                if are_rel == True:
                    site_related = True
                    site_correspondances.append(i)
                    test_indexes.remove(i)
                    break
                else:
                    site_related = False
            if are_rel == False:
                site_related = False
                break
    if len(test_indexes) > 0:
        site_related = False
    if ag_related == True and site_related == True:
        are_related = True
    else:
        are_related = False

    if return_correspondances == False:
        return are_related
    elif return_correspondances == True:
        return are_related, ag_correspondances, site_correspondances


def edit_vs_test_sites(edit, test):
    """
    Determine if a site produce by an edit can be used by a given test site.
    The test may contain wildcards _ or #.
    """

    ignore_partner = False
    are_rel = True
    if edit["name"] != test["name"]:
        are_rel = False
    if edit["agentname"] != test["agentname"]:
        are_rel = False
    if edit["agentid"] != test["agentid"]:
        are_rel = False
    if test["value"] != "#":
        if edit["value"] != test["value"]:
            are_rel = False
    # Check bonds.
    if are_rel == True:
        # If test bond is a wildcard #, anything goes.
        if test["bond"] == "#":
            ignore_partner = True
            return are_rel, ignore_partner
        # Check if both sites are free or without binding.
        elif edit["bond"] == "." or test["bond"] == ".":
            if edit["bond"] != "." or test["bond"] != ".":
                are_rel = False
        elif edit["bond"] == None or test["bond"] == None:
            if edit["bond"] != None or test["bond"] != None:
                are_rel = False
        # Else check if they have equivalent bindings.
        else:
            if test["bond"] == "_":
                ignore_partner = True
            explicit_edit = isinstance(edit["bond"], dict)
            explicit_test = isinstance(test["bond"], dict)
            if explicit_edit == True and explicit_test == True:
                edit_partner = edit["bond"]["partner"]
                test_partner = test["bond"]["partner"]
                if edit_partner["agentname"] != test_partner["agentname"]:
                    are_rel = False
                if edit_partner["agentid"] != test_partner["agentid"]:
                    are_rem = False
                if edit_partner["sitename"] != test_partner["sitename"]:
                    are_rel = False

    return are_rel, ignore_partner


def compare_states(state1, state2, ignorevalue=False, ignoreid=False,
                   ignoretype=False):
    """ Determine if two states (two lists of agents) are the same. """

    list1 = state1.copy()
    list2 = state2.copy()
    found1 = []
    found2 = []
    for i in range(len(list1)):
        agent1 = list1[i]
        for agent2 in list2:
            same_agents, ignrs = compare_agents(agent1, agent2, ignorevalue,
                                                ignoreid, ignoretype)
            # Search state for partner of ignrs[0][site k]
            if same_agents == True:
                found1.insert(0, i)
                break
    for j in range(len(list2)):
        agent2 = list2[j]
        for agent1 in list1:
            same_agents, ignrs = compare_agents(agent2, agent1, ignorevalue,
                                                ignoreid, ignoretype)
            if same_agents == True:
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


def compare_agents(agent1, agent2, ignorevalue=False, ignoreid=False,
                   ignoretype=False):
    """ Determine if two agents are the same (with same sites). """

    ignore_sites = [[], []]
    are_same = True
    if agent1["name"] != agent2["name"]:
        are_same = False
    if ignoretype == False:
        if agent1["type"] != agent2["type"]:
           are_same = False
    if ignoreid == False:
        if agent1["id"] != agent2["id"]:
           are_same = False
    if are_same == True:
        # Compare the sites.
        list1 = agent1["sites"].copy()
        list2 = agent2["sites"].copy()
        found1 = []
        found2 = []
        ignore1 = []
        ignore2 = []
        for i in range(len(list1)):
            site1 = list1[i]
            for site2 in list2:
                same_sites, ignr = compare_sites(site1, site2, ignorevalue,
                                                 ignoreid, ignoretype)
                if ignr != 0: 
                    if ignr == 1:
                        # Add ignr site to ignore_sites list to later pass
                        # it to compare_states.
                        ignore_sites[0].append(i)
                        # Search agent1 for the index of partner of site1.
                        index = search_partner_in_agent(agent1, site1)
                    if ignr == 2:
                        index = search_partner_in_agent(agent2, site2)
                    if index != None:
                        ignore1.append(index)
                if same_sites == True:
                    found1.insert(0, i)
                    break
        for j in range(len(list2)):
            site2 = list2[j]
            for site1 in list1:
                same_sites, ignr = compare_sites(site2, site1, ignorevalue,
                                                 ignoreid, ignoretype)
                if ignr != 0:
                    if ignr == 1:
                        # Add ignr site to ignore_sites list to later pass
                        # it to compare_states.
                        ignore_sites[1].append(j)
                        # Search agent2 for the index of partner of site2.
                        index = search_partner_in_agent(agent2, site2)
                    if ignr == 2:
                        index = search_partner_in_agent(agent1, site1)
                    if index != None:
                        ignore2.append(index)
                if same_sites == True:
                    found2.insert(0, j)
                    break
        to_remove1 = sorted(found1+ignore1, reverse=True)
        for i in to_remove1:
            del(list1[i])
        to_remove2 = sorted(found2+ignore2, reverse=True)
        for j in to_remove2:
            del(list2[j])
        if len(list1) == 0 and len(list2) == 0:
            are_same = True
        else:
            are_same = False 

    return are_same, ignore_sites


def search_partner_in_agent(agent, site):
    """ Search agent for the index of partner of given site. """

    partner_index = None
    list1 = agent["sites"]
    bnd_num = site["bond"]["num"]
    for k in range(len(list1)):
        if list1[k] != site:
            if isinstance(list1[k]["bond"], dict):
                if list1[k]["bond"]["num"] == bnd_num:
                    partner_index = k
                    break

    return partner_index


def compare_sites(site1, site2, ignorevalue=False, ignoreid=False,
                  ignoretype=False):
    """ Determine if two sites are the same. """

    ignore_partner = 0
    are_same = True
    if site1["name"] != site2["name"]:
        are_same = False
    if ignoretype == False:
        if site1["type"] != site2["type"]:
            are_same = False
    if ignorevalue == False:
        if site1["value"] != "#" and site1["value"] != "#":
            if site1["value"] != site2["value"]:
                are_same = False
    # Check bonds.
    if are_same == True:
        # If any bond is a wildcard, anything goes.
        # Remember to tell compare_agents to then ignore
        # the partner of the non-wildcard site.
        if site1["bond"] == "#":
            ignore_partner = 2
        if site2["bond"] == "#":
            ignore_partner = 1
        if site1["bond"] == "#" and site2["bond"] == "#":
           raise ValueError("Comparing two # wildcards is not"
                            "expected to happen.") 
        if ignore_partner != 0:
            return True, ignore_partner
        # Check if both sites are free or with no binding.
        elif site1["bond"] == None or site2["bond"] == None:
            if site1["bond"] != None or site2["bond"] != None:
                are_same = False
        elif site1["bond"] == "." or site2["bond"] == ".":
            if site1["bond"] != "." or site2["bond"] != ".":
                are_same = False
        # Else check if they have equivalent bindings.
        else:
            if site1["bond"] == "_":
                ignore_partner = 2
            if site2["bond"] == "_":
                ignore_partner = 1
            if site1["bond"] == "_" and site2["bond"] == "_":
               raise ValueError("Comparing two _ wildcards is not"
                                "expected to happen.")
            has_partner1 = isinstance(site1["bond"], dict)
            has_partner2 = isinstance(site2["bond"], dict)
            if has_partner1 == True and has_partner2 == True:
                partner1 = site1["bond"]["partner"]
                partner2 = site2["bond"]["partner"]
                if partner1["agentname"] != partner2["agentname"]:
                    are_same = False
                if ignoreid == False:
                    if partner1["agentid"] != partner2["agentid"]:
                        are_same = False
                if partner1["sitename"] != partner2["sitename"]:
                    are_same = False
    #has_partner1 = isinstance(site1["bond"], dict)
    #has_partner2 = isinstance(site2["bond"], dict)
    #if has_partner1 != has_partner2:
    #    if site1["bond"] != "_" and site2["bond"] != "_":
    #        are_same = False
    #if are_same == True:
    #    if has_partner1 == False:
    #        if site1["bond"] != site2["bond"]:
    #            are_same = False
    #    elif has_partner1 == True:
    #        if site1["bond"] != "_" and site2["bond"] != "_":
    #            partner1 = site1["bond"]["partner"]
    #            partner2 = site2["bond"]["partner"]
    #            if partner1["agentname"] != partner2["agentname"]:
    #                are_same = False
    #            if ignoreid == False:
    #                if partner1["agentid"] != partner2["agentid"]:
    #                    are_same = False
    #            if partner1["sitename"] != partner2["sitename"]:
    #                are_same = False

    return are_same, ignore_partner


#def compare_states(state1, state2, ignorevalue=False, ignoreid=False):
#    """ Determine if two states (two lists of sites) are the same. """
#
#    list1 = state1.copy()
#    list2 = state2.copy()
#    found1 = []
#    found2 = []
#    for i in range(len(list1)):
#        site1_str = write_kappa_agent(list1[i], bond="partner",
#                                           hidevalue=ignorevalue,
#                                           hideid=ignoreid)
#        for site2 in list2:
#            site2_str = write_kappa_agent(site2, bond="partner",
#                                               hidevalue=ignorevalue,
#                                               hideid=ignoreid)
#            if site1_str == site2_str:
#                found1.insert(0, i)
#                break
#    for j in range(len(list2)):
#        site2_str = write_kappa_agent(list2[j], bond="partner",
#                                           hidevalue=ignorevalue,
#                                           hideid=ignoreid)
#        for site1 in list1:
#            site1_str = write_kappa_agent(site1, bond="partner",
#                                               hidevalue=ignorevalue,
#                                               hideid=ignoreid)
#            if site2_str == site1_str:
#                found2.insert(0, j)
#                break
#    for i in found1:
#        del(list1[i])
#    for j in found2:
#        del(list2[j])
#    if len(list1) == 0 and len(list2) == 0:
#        are_same = True
#    else:
#        are_same = False
#
#    return are_same


def get_output_of_node(story, node):
    """
    Get the output of a given event node. The output is a list of
    states (with context), which are themselves a list of agents.
    There is one state per target state node of a given event node.
    """

    # Get output state nodes.
    output_nodes = []
    for edge in story.causaledges:
        if edge.source == node:
            if isinstance(edge.target, StateNode):
                output_nodes.append(edge.target)
    # Get states from each output node.
    output_states = []
    for statenode in output_nodes:
        if node.intro == True:
            output_states.append(statenode.edit)
        else:
            output_states.append(statenode.state)

    return output_states


#def write_context_expression2(state, hidevalue=False, hideid=False):
#    """
#    Write a Kappa language string with the edit in bold font and context in
#    normal font. The string is made to be read as html to allow special fonts.
#    """
#
#    print(state)
#    # Sort agents by the number of bonds that they have.
#    for agent in state:
#        nbonds = 0
#        for site in agent["sites"]:
#            if site["bond"] != None and site["bond"] != ".":
#                nbonds += 1
#        agent["nbonds"] = nbonds
#    sorted_agents = sorted(state, key=lambda x: x["nbonds"])
#    # Put agents containing the edited sites first.
#    full_context = []
#    used_agents = []
#    for i in range(len(sorted_agents)):
#        agent = sorted_agents[i]
#        if agent["type"] == "edit":
#            full_context.append(agent)
#            used_agents.insert(0, i)
#            break
#        for site in agent["sites"]:
#            if site["type"] == "edit":
#                full_context.append(agent)
#                used_agents.insert(0, i)
#                break
#    for i in used_agents:
#        del(sorted_agents[i])
#    print(full_context)
#    # Iteratively add context agents to the left or right of the full_context
#    # list depending on whether they are bound to an agent that is on the
#    # first or second half of the current full_context list.
#    bond_found = True
#    while bond_found == True:
#        bond_found = False
#        half_index = int(len(full_context)/2)
#        for i in range(len(sorted_agents)):
#            remaining_agent = sorted_agents[i]
#            # Check if remaining agent has one site that is bound to an agent
#            # that is already in full_context.
#            for rem_site in remaining_agent["sites"]:
#                if rem_site["bond"] != None and rem_site["bond"] != ".":
#                    partner = rem_site["bond"]["partner"]
#                    partnerid = "{}:{}".format(partner["agentname"],
#                                               partner["agentid"])
#                    for j in range(len(full_context)):
#                        existing_agent = full_context[j]
#                        fullid = "{}:{}".format(existing_agent["name"],
#                                                existing_agent["id"])
#                        if partnerid == fullid:
#                            # Add this remaining agent to full_context list.
#                            if j < half_index:
#                                full_context.insert(0, remaining_agent)
#                            elif j >= half_index:
#                                full_context.append(remaining_agent)
#                            bond_found = True
#                            break
#                if bond_found == True:
#                    break
#            if bond_found == True:
#                del(sorted_agents[i])
#                break
#    # Add remaining unbound agents.
#    for agent in sorted_agents:
#        full_context.append(agent)
#    # Reassign bond numbers.
#    for agent in full_context:
#        for site in agent["sites"]:
#            if site["bond"] != None and site["bond"] != ".":
#                site["bond"]["num"] = 0
#    bondid = 1
#    for agent1 in full_context:
#        for site1 in agent1["sites"]:
#            if site1["bond"] != None and site1["bond"] != ".":
#                if site1["bond"]["num"] == 0:
#                    # Find partner.
#                    partner = site1["bond"]["partner"]
#                    partnerid = "{}:{}".format(partner["agentname"],
#                                               partner["agentid"])
#                    for agent2 in full_context:
#                        fullid = "{}:{}".format(agent2["name"], agent2["id"])
#                        if fullid == partnerid:
#                            for site2 in agent2["sites"]:
#                                if site2["name"] == partner["sitename"]:
#                                    site1["bond"]["num"] = bondid
#                                    site2["bond"]["num"] = bondid
#                                    bondid += 1
#                                    break
#    print("---")
#    print(full_context)
#    # Sort sites within agents. Put sites used for bindings on the left if
#    # they bind to an agent that appears before in the full_context list.
#    # Put sites that have only a value and no binding at the end.
#    for i in range(len(full_context)):
#        agent1 = full_context[i]
#        ordered_edits = []
#        ordered_context = []
#        for site1 in agent1["sites"]:
#            if site1["bond"] != None and site1["bond"] != ".":
#                partner = site1["bond"]["partner"]
#                partnerid = "{}:{}".format(partner["agentname"],
#                                           partner["agentid"])
#                for j in range(len(full_context)):
#                    agent2 = full_context[j]
#                    fullid = "{}:{}".format(agent2["name"], agent2["id"])
#                    if fullid == partnerid:
#                        if j < i:
#                            if site1["type"] == "edit":
#                                ordered_edits.insert(0, site1)
#                            elif site1["type"] == "context":
#                                ordered_context.insert(0, site1)
#                        elif j >= i:
#                            if site1["type"] == "edit":
#                                ordered_edits.append(site1)
#                            elif site1["type"] == "context":
#                                ordered_context.append(site1)
#                        break
#        for site1 in agent1["sites"]:
#            if site1["bond"] == None or site1["bond"] == ".":
#                if site1["type"] == "edit":
#                    ordered_edits.append(site1)
#                elif site1["type"] == "context":
#                    ordered_context.append(site1)
#        agent1["sites"] = ordered_edits + ordered_context
#    print("---")
#    print(full_context)
#    # Write string with sites in the order determined by the previous sorting.
#    context_str = ""
#    for i in range(len(full_context)):
#        agent_str = write_kappa_agent(full_context[i], "num",
#                                      hidevalue, hideid)
#        context_str += agent_str
#        if i < len(full_context)-1:
#            context_str += ", "
#
#    return context_str


def write_context_expression(state, hidevalue=False, hideid=False):
    """
    Write a Kappa language string with the edit in bold font and context in
    normal font. The string is made to be read as html to allow special fonts.
    """

    # Sort agents by the number of bonds that they have.
    for agent in state:
        nbonds = 0
        for site in agent["sites"]:
            if site["bond"] != None:
                nbonds += 1
        agent["nbonds"] = nbonds
    sorted_agents = sorted(state, key=lambda x: x["nbonds"])
    # Put agents containing the edited sites first.
    full_context = []
    used_agents = []
    for i in range(len(sorted_agents)):
        agent = sorted_agents[i]
        if agent["type"] == "edit":
            full_context.append(agent)
            used_agents.insert(0, i)
            break
        for site in agent["sites"]:
            if site["type"] == "edit":
                full_context.append(agent)
                used_agents.insert(0, i)
                break
    for i in used_agents:
        del(sorted_agents[i])
    # Iteratively add context agents to the left or right of the full_context
    # list depending on whether they are bound to an agent that is on the
    # first or second half of the current full_context list.
    bond_found = True
    while bond_found == True:
        bond_found = False
        half_index = int(len(full_context)/2)
        for i in range(len(sorted_agents)):
            remaining_agent = sorted_agents[i]
            # Check if remaining agent has one site that is bound to an agent
            # that is already in full_context.
            for rem_site in remaining_agent["sites"]:
                if rem_site["bond"] != None and rem_site["bond"] != ".":
                    partner = rem_site["bond"]["partner"]
                    partnerid = "{}:{}".format(partner["agentname"],
                                               partner["agentid"])
                    for j in range(len(full_context)):
                        existing_agent = full_context[j]
                        fullid = "{}:{}".format(existing_agent["name"],
                                                existing_agent["id"])
                        if partnerid == fullid:
                            # Add this remaining agent to full_context list.
                            if j < half_index:
                                full_context.insert(0, remaining_agent)
                            elif j >= half_index:
                                full_context.append(remaining_agent)
                            bond_found = True
                            break
                if bond_found == True:
                    break
            if bond_found == True:
                del(sorted_agents[i])
                break
    # Add remaining unbound agents.
    for agent in sorted_agents:
        full_context.append(agent)
    # Reassign bond numbers.
    for agent in full_context:
        for site in agent["sites"]:
            #if site["bond"] != None and site["bond"] != ".":
            if isinstance(site["bond"], dict):
                site["bond"]["num"] = 0
    bondid = 1
    for agent1 in full_context:
        for site1 in agent1["sites"]:
            #if site1["bond"] != None and site1["bond"] != ".":
            if isinstance(site1["bond"], dict):
                if site1["bond"]["num"] == 0:
                    # Find partner.
                    partner = site1["bond"]["partner"]
                    partnerid = "{}:{}".format(partner["agentname"],
                                               partner["agentid"])
                    for agent2 in full_context:
                        fullid = "{}:{}".format(agent2["name"], agent2["id"])
                        if fullid == partnerid:

                            for site2 in agent2["sites"]:
                                if site2["name"] == partner["sitename"]:
                                    #sb2 = site2["bond"]
                                    #if sb2 != None and sb2 != ".":
                                    if isinstance(site2["bond"], dict):
                                        site1["bond"]["num"] = bondid
                                        site2["bond"]["num"] = bondid
                                        bondid += 1
                                        break
    # Sort sites within agents. Put sites used for bindings on the left if
    # they bind to an agent that appears before in the full_context list.
    # Put sites that have only a value and no binding at the end.
    for i in range(len(full_context)):
        agent1 = full_context[i]
        ordered_edits = []
        ordered_context = []
        ordered_para = []
        for site1 in agent1["sites"]:
            #if site1["bond"] != None and site1["bond"] != ".":
            if isinstance(site1["bond"], dict):
                partner = site1["bond"]["partner"]
                partnerid = "{}:{}".format(partner["agentname"],
                                           partner["agentid"])
                for j in range(len(full_context)):
                    agent2 = full_context[j]
                    fullid = "{}:{}".format(agent2["name"], agent2["id"])
                    if fullid == partnerid:
                        if j < i:
                            if site1["type"] == "edit":
                                ordered_edits.insert(0, site1)
                            elif site1["type"] == "context":
                                ordered_context.insert(0, site1)
                            elif site1["type"] == "parallel":
                                ordered_para.insert(0, site1)
                        elif j >= i:
                            if site1["type"] == "edit":
                                ordered_edits.append(site1)
                            elif site1["type"] == "context":
                                ordered_context.append(site1)
                            elif site1["type"] == "parallel":
                                ordered_para.insert(0, site1)
                        break
        for site1 in agent1["sites"]:
            #if site1["bond"] == None or site1["bond"] == ".":
            if not isinstance(site1["bond"], dict):
                if site1["type"] == "edit":
                    ordered_edits.append(site1)
                elif site1["type"] == "context":
                    ordered_context.append(site1)
                elif site1["type"] == "parallel":
                    ordered_para.append(site1)
        #agent1["sites"] = ordered_edits + ordered_context + ordered_para
        # Put ordered_context sites after ordered_edits of a same site.
        agent1["sites"] = []
        for eds in ordered_edits:
            agent1["sites"].append(eds)
            conts_to_remove = []
            for i in range(len(ordered_context)):
                conts = ordered_context[i]
                if conts["name"] == eds["name"]:
                    agent1["sites"].append(conts)
                    conts_to_remove.insert(0, i)
            for i in conts_to_remove:
                del(ordered_context[i])
        for conts in ordered_context:
            agent1["sites"].append(conts)
  

    # Write string with sites in the order determined by the previous sorting.
    context_str = ""
    for i in range(len(full_context)):
        agent_str = write_kappa_agent(full_context[i], "num",
                                      hidevalue, hideid)
        context_str += agent_str
        if i < len(full_context)-1:
            context_str += ", "

    return context_str


#def read_context_expression(label):
#    """ Obtain the state of a node based on its label. """
#
#
#    return state


def compare_outputs(output1, output2):
    """ Determine if two outputs (two lists of states) are the same. """

    list1 = output1.copy()
    list2 = output2.copy()
    found1 = []
    found2 = []
    for i in range(len(list1)):
        state1 = list1[i]
        for state2 in list2:
            same_states = compare_states(state1, state2, ignorevalue=False,
                                         ignoreid=True)
            if same_states:
                found1.insert(0, i)
                break
    for j in range(len(list2)):
        state2 = list2[j]
        for state1 in list1:
            same_states = compare_states(state2, state1, ignorevalue=False,
                                         ignoreid=True)
            if same_states:
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


#def compare_state_test(state, test):
#    """ Determine if a state and test represent the same species. """
#
#    list1 = state.copy()
#    list2 = test.copy()
#    found1 = []
#    found2 = []
#    for i in range(len(list1)):
#        site1_str = write_kappa_expression(list1[i], bond="partner")
#        for site2 in list2:
#            site2_str = write_kappa_expression(site2, bond="partner")
#            if site1_str == site2_str:
#                found1.insert(0, i)
#                break
#    for j in range(len(list2)):
#        site2_str = write_kappa_expression(list2[j], bond="partner")
#        for site1 in list1:
#            site1_str = write_kappa_expression(site1, bond="partner")
#            if site2_str == site1_str:
#                found2.insert(0, j)
#                break
#    for i in found1:
#        del(list1[i])
#    for j in found2:
#        del(list2[j])
#    if len(list1) == 0 and len(list2) == 0:
#        are_same = True
#    else:
#        are_same = False
#
#    return are_same


#def accumulatecontext(eoi, showintro=True, addedgelabels=False,
#                      showedgelabels=False, edgeid=True, edgeocc=False,
#                      edgeuse=False, statstype="abs", writedot=True,
#                      weightedges=True):
#    """ Get the cumulative relevant context for each state node. """
#
#    # Reading section.
#    story_files = get_dot_files("{}".format(eoi), "edits")
#    stories = []
#    for story_file in story_files:
#        story_path = "{}/{}".format(eoi, story_file)
#        stories.append(CausalGraph(story_path, eoi))
#    # Write stories with context on state nodes.
#    for i in range(len(stories)):
#        stories[i].filename = "context-{}.dot".format(i+1)
#    for story in stories:
#        story.build_dot_file(showintro, addedgelabels, showedgelabels,
#                             edgeid, edgeocc, edgeuse, statstype, weightedges)
#        output_path = "{}/{}".format(eoi, story.filename)
#        outfile = open(output_path, "w")
#        outfile.write(story.dot_file)
#        outfile.close()

# ;;;;;;;;;;;;;;;;;;;;;;;; End of Context Section ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

# ^^^^^^^^^^^^^^^^^^^^^^ Dual Story Merging Section ^^^^^^^^^^^^^^^^^^^^^^^^^^^

def getuniquestories(eoi, causalgraphs=None, siphon=False, showintro=True,
                     addedgelabels=False, showedgelabels=False, edgeid=True,
                     edgeocc=False, edgeprob=False, statstype="abs",
                     weightedges=False, color=True, writedot=True,
                     rmprev=False, msg=True):
    """ Get unique stories, dual stories and state stories. """

    print("modstory")
    mergedualstories(eoi, "modstory", showintro=showintro,
                     addedgelabels=addedgelabels,
                     showedgelabels=showedgelabels, edgeid=edgeid,
                     edgeocc=edgeocc, edgeprob=edgeprob,
                     weightedges=weightedges, color=color, writedot=writedot,
                     rmprev=rmprev)
    print("dualstory")
    mergedualstories(eoi, "dualstory", showintro=showintro,
                     addedgelabels=addedgelabels,
                     showedgelabels=showedgelabels, edgeid=edgeid,
                     edgeocc=edgeocc, edgeprob=edgeprob,
                     weightedges=weightedges, color=color, writedot=writedot,
                     rmprev=rmprev)
    print("statestory")
    mergedualstories(eoi, "statestory", showintro=showintro,
                     addedgelabels=addedgelabels,
                     showedgelabels=showedgelabels, edgeid=edgeid,
                     edgeocc=edgeocc, edgeprob=edgeprob,
                     weightedges=weightedges, color=color, writedot=writedot,
                     rmprev=rmprev)


def mergedualstories(eoi, prefix, causalgraphs=None, siphon=False, showintro=True,
                     addedgelabels=False, showedgelabels=False, edgeid=True,
                     edgeocc=False, edgeprob=False, statstype="abs",
                     weightedges=False, color=True, writedot=True,
                     rmprev=False, msg=True):
    """
    Merge equivalent dual stories into a unique dual story while counting
    occurrence.
    """

    # Reading section.
    if causalgraphs == None:
        story_files = get_dot_files("{}/tmp".format(eoi), prefix)
        stories = []
        for story_file in story_files:
            story_path = "{}/tmp/{}".format(eoi, story_file)
            stories.append(CausalGraph(story_path, eoi))
    else:
        stories = causalgraphs
        story_files = None
    # Doing the work.
    merged_stories = []
    while len(stories) > 0:
        current_story = stories[0]
        equivalent_list = [0]
        for i in range(1, len(stories)):
            same_story, ev, st, equi_edges = equivalent_graphs(current_story,
                                                               stories[i],
                                                               True, True)
            if same_story == True:
                equivalent_list.insert(0, i)
                current_story.occurrence += stories[i].occurrence
                for j in range(len(current_story.hyperedges)):
                    equi_index = equi_edges[j]
                    weight = stories[i].hyperedges[equi_index].weight
                    current_story.hyperedges[j].weight += weight
        # Find the original dual stories from which each unique
        # story comes from.
        original_stories = []
        for index in equivalent_list:
            file_name = stories[index].filename
            dash = file_name.rfind("-")
            period = file_name.rfind(".")
            #if "_node" in file_name:
            #    underscore = file_name.index("_node")
            #    previd = file_name[:period]
            #else:
            previd = file_name[dash+1:period]
            original_stories.append(previd)
        current_story.prevcores = original_stories
        merged_stories.append(current_story)
        for i in equivalent_list:
            del(stories[i])
        #for i in range(len(analogous_list)-1, -1, -1):
        #    index = analogous_list[i]
        #    del(stories[index])
    sorted_stories = sorted(merged_stories, key=lambda x: x.occurrence,
                            reverse=True)
    # Propagate new hyperedge weights to their edge lists.
    for story in sorted_stories:
        for hyperedge in story.hyperedges:
            for subedge in hyperedge.edgelist:
                subedge.weight = hyperedge.weight
        story.align_vertical()
    # Write merged dual stories.
    for i in range(len(sorted_stories)):
        sorted_stories[i].filename = "{}-{}.dot".format(prefix, i+1)
    if writedot == True:
        for story in sorted_stories:
            #story.compute_relstats()
            #story.compute_visuals(showintro, color)
            story.get_maxrank()
            story.build_nointro()
            story.build_dot_file(showintro, addedgelabels, showedgelabels,
                                 edgeid, edgeocc, edgeprob, statstype, weightedges)
            output_path = "{}/unique/{}".format(eoi, story.filename)
            outfile = open(output_path, "w")
            outfile.write(story.dot_file)
            outfile.close()
        if prefix == "dualstory":
            # Write state dict in json format.
            for i in range(len(sorted_stories)):
                statedict = {"states": {}, "edits": {}}
                for statenode in sorted_stories[i].statenodes:
                    statedict["states"][statenode.nodeid] = statenode.state
                    statedict["edits"][statenode.nodeid] = statenode.edit
                json_path = "{}/unique/statefile-{}.json".format(eoi, i+1)
                statefile = open(json_path, "w")
                json.dump(statedict, statefile)
                statefile.close()

    return sorted_stories


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


def equivalent_graphs(graph1, graph2, enforcerank=True,
                      return_correspondances=False):
    """
    Equivalent graphs have equivalent nodes (events and states) and equivalent
    hyperedges.
    """

    corr_st, corr_ed = [], []
    equi_events, corr_ev = equivalent_node_lists(graph1.eventnodes,
                                                 graph2.eventnodes,
                                                 enforcerank, True)
    if equi_events == True:
        equi_states, corr_st = equivalent_node_lists(graph1.statenodes,
                                                     graph2.statenodes,
                                                     enforcerank, True)
    if equi_events == True and equi_states == True:
        equi_edges, corr_ed = equivalent_hyperedge_lists(graph1.hyperedges,
                                                         graph2.hyperedges,
                                                         enforcerank, True)
    if equi_events == True and equi_states == True and equi_edges == True:
        are_equivalent = True
    else:
        are_equivalent = False

    if return_correspondances == False:
        return are_equivalent
    elif return_correspondances == True:
        return are_equivalent, corr_ev, corr_st, corr_ed


def equivalent_hyperedge_lists(edgelist1, edgelist2, enforcerank=True,
                               return_correspondances=False):
    """
    Find whether two lists contain equivalent hyperedges (connected to nodes
    with same labels).
    """

    are_equivalent = True
    correspondances = []
    list2_indexes = list(range(len(edgelist2)))
    for hyperedge1 in edgelist1:
        equi_edges = False
        for i in list2_indexes:
            hyperedge2 = edgelist2[i]
            equi_edges = equivalent_hyperedges(hyperedge1, hyperedge2,
                                               enforcerank)
            if equi_edges == True:
                correspondances.append(i)
                list2_indexes.remove(i)
                break
        # If the previous loop finished with the last hyperedge from list2
        # still not being equivalent.
        if equi_edges == False:
            are_equivalent = False
            break
    if len(list2_indexes) > 0:
        are_equivalent = False

    if return_correspondances == False:
        return are_equivalent
    elif return_correspondances == True:
        return are_equivalent, correspondances


def equivalent_hyperedges(hyperedge1, hyperedge2, enforcerank=True,
                          return_correspondances=False,
                          disregard_duplicates=False):

    """
    Find whether two hyperedges connect to nodes with the same labels and
    optionally at same ranks.
    """

    equi_targets = equivalent_nodes(hyperedge1.target,
                                    hyperedge2.target,
                                    enforcerank)

    if disregard_duplicates == False:
        hyperedge2_sources = hyperedge2.sources
    elif disregard_duplicates == True:
        hyperedge2_sources = []
        seen_labels = []
        for source in hyperedge2.sources:
            if source.label not in seen_labels:
                hyperedge2_sources.append(source)
                seen_labels.append(source.label)

    equi_srcs, corr_srcs = equivalent_node_lists(hyperedge1.sources,
                                                 hyperedge2_sources,
                                                 enforcerank, True)
    if equi_targets == True and equi_srcs == True:
        are_equivalent = True
    else:
        are_equivalent = False

    if return_correspondances == False:
        return are_equivalent
    elif return_correspondances == True:
        return are_equivalent, corr_srcs
    


def equivalent_node_lists(nodelist1, nodelist2, enforcerank=True,
                          return_correspondances=False):
    """
    Find whether two lists contain nodes with same labels and are optionally
    at the same ranks.
    There may be many nodes with the same label in a given list. If so, the
    number of dublicates must match between the two lists.
    Can also return a correspondance list, which for each element of list1
    gives the index of the corresponding element in list2.
    """

    are_equivalent = True
    correspondances = []
    list2_indexes = list(range(len(nodelist2)))
    for node1 in nodelist1:
        equi_nodes = False
        for i in list2_indexes:
            node2 = nodelist2[i]
            equi_nodes = equivalent_nodes(node1, node2, enforcerank)
            if equi_nodes == True:
                correspondances.append(i)
                list2_indexes.remove(i)
                break
        # If the previous loop finished with the last node from list2
        # still not being equivalent.
        if equi_nodes == False:
            are_equivalent = False
            break
    if len(list2_indexes) > 0:
        are_equivalent = False

    if return_correspondances == False:
        return are_equivalent
    elif return_correspondances == True:
        return are_equivalent, correspondances


def equivalent_nodes(node1, node2, enforcerank=True):
    """
    Find whether two nodes have the same label and optionally are at the same
    rank.
    """

    are_equivalent = False
    if node1.label == node2.label:
        are_equivalent = True
    if enforcerank == True:
        if node1.rank != node2.rank:
            are_equivalent = False

    return are_equivalent

# ^^^^^^^^^^^^^^^^^^^ End of Dual Story Merging Section ^^^^^^^^^^^^^^^^^^^^^^^

# ==================== Causal Cores Merging Section ===========================

#def mergecores(eoi, causalgraphs=None, siphon=False, showintro=True,
#               addedgelabels=False, showedgelabels=False, edgeid=True,
#               edgeocc=False, edgeuse=False, statstype="abs",
#               weightedges=False, color=True, writedot=True, rmprev=False,
#               msg=True):
#    """
#    Merge analogous causal cores and count occurrence.
#    Write the final cores as meshed graphs.
#    """
#
#    # Reading section.
#    if causalgraphs == None:
#        if siphon == False:
#            causal_core_files = get_dot_files(eoi, "causalcore")
#        elif siphon == True:
#            causal_core_files = get_dot_files(eoi, "siphon")
#        causal_cores = []
#        for core_file in causal_core_files:
#            core_path = "{}/{}".format(eoi, core_file)
#            causal_cores.append(CausalGraph(core_path, eoi))
#    else:
#       causal_cores = causalgraphs
#       causal_core_files = None
#    # Doing the work.
#    merged_cores = []
#    while len(causal_cores) > 0:
#        current_core = causal_cores[0]
#        analogous_list = [0]
#        for i in range(1, len(causal_cores)):
#            same_core, equi_meshes = analogous_graphs(current_core,
#                                                      causal_cores[i])
#            if same_core == True:
#                analogous_list.append(i)
#                current_core.occurrence += causal_cores[i].occurrence
#                for j in range(len(current_core.meshes)):
#                    equi_index = equi_meshes[j]
#                    uses = causal_cores[i].meshes[equi_index].uses
#                    current_core.meshes[j].uses += uses
#        prevcores = []
#        for index in analogous_list:
#            file_name = causal_cores[index].filename
#            dash = file_name.rfind("-")
#            period = file_name.rfind(".")
#            if "_node" in file_name:
#                underscore = file_name.index("_node")
#                previd = file_name[:period]
#            else:
#                previd = file_name[dash+1:period]
#            prevcores.append(previd)
#        current_core.prevcores = prevcores
#        merged_cores.append(current_core)
#        for i in range(len(analogous_list)-1, -1, -1):
#            index = analogous_list[i]
#            del(causal_cores[index])
#    sorted_cores = sorted(merged_cores, key=lambda x: x.occurrence,
#                          reverse=True)
#    for i in range(len(sorted_cores)):
#        sorted_cores[i].filename = "meshedcore-{}.dot".format(i+1)
#    for graph in sorted_cores:
#        graph.compute_relstats()
#        graph.compute_visuals(showintro, color)
#        graph.build_dot_file(showintro, addedgelabels, showedgelabels,
#                             edgeid, edgeocc, edgeuse, statstype, weightedges)
#    # Writing section.
#    if writedot == True:
#        for graph in sorted_cores:
#            output_path = "{}/{}".format(eoi, graph.filename)
#            outfile = open(output_path, "w")
#            outfile.write(graph.dot_file)
#            outfile.close()
#    if rmprev == True:
#        if causal_core_files == None:
#            causal_core_files = get_dot_files(eoi, "causalcore")
#        for core_file in causal_core_files:
#            file_path = "{}/{}".format(eoi, core_file)
#            os.remove(file_path)
#    if msg == True:
#        print("Merging equivalent causal cores, {} unique cores obtained."
#              .format(len(sorted_cores)))
#
#    return sorted_cores


#def get_dot_files(eoi, prefix=None):
#    """ Get the number of the first and last stories. """
#
#    tmp_file_list = os.listdir("{}".format(eoi))
#    file_list = []
#    for file_name in tmp_file_list:
#        if "dot" in file_name:
#            if prefix == None:
#                file_list.append(file_name)
#            else:
#                dash = file_name.rfind("-")
#                if file_name[:dash] == prefix:
#                    file_list.append(file_name)
#    file_dicts = []
#    for file_name in file_list:
#        dash = file_name.rfind("-")
#        period = file_name.rfind(".")
#        number = int(file_name[dash+1:period])
#        file_dicts.append({"file": file_name, "num": number})
#    sorted_dicts = sorted(file_dicts, key=lambda x: x["num"])
#    sorted_list = []
#    for d in sorted_dicts:
#        sorted_list.append(d["file"])
#
#    return sorted_list
#
#
#def analogous_graphs(graph1, graph2):
#    """
#    Analogous causal graphs have analogous meshes. That is, all their meshes
#    are between nodes with same labels at same ranks.
#    """
#
#    equi_meshes = []
#    if graph1.maxrank == graph2.maxrank:
#        graph2_indexes = list(range(len(graph2.meshes)))
#        all_edges_found = True
#        for mesh1 in graph1.meshes:
#            for i in graph2_indexes:
#                mesh2 = graph2.meshes[i]
#                are_equi = analogous_meshes(mesh1, mesh2, enforcerank=True)
#                if are_equi == True:
#                    equi_meshes.append(i)
#                    graph2_indexes.remove(i)
#                    break
#            if are_equi == False:
#                all_edges_found = False
#                break
#        # All the edges from graph2 should have been used
#        # at this point for both graphs to be equivalent.
#        if all_edges_found == True:
#            if len(graph2_indexes) > 0:
#                equi_graphs = False
#            else:
#                equi_graphs = True
#        else:
#            equi_graphs = False
#    else:
#        equi_graphs = False
#
#    return equi_graphs, equi_meshes


#def analogous_meshes(mesh1, mesh2, enforcerank=True):
#    """
#    Find whether two meshes connect to event nodes with same labels
#    with analogous midedges.
#    Optionally, nodes may be also required to be at same ranks.
#    """
#
#    nn1 = len(mesh1.midnodes)
#    nn2 = len(mesh2.midnodes)
#    ne1 = len(mesh1.midedges)
#    ne2 = len(mesh2.midedges)
#    if nn1 == nn2 and ne1 == ne2:
#        are_equi = True
#    else:
#        are_equi = False
#    if are_equi == True:
#        sources1, targets1 = mesh1.get_events()
#        sources2, targets2 = mesh2.get_events()
#        equi_sources = analogous_nodes(sources1, sources2, enforcerank)
#        equi_targets = analogous_nodes(targets1, targets2, enforcerank)
#        if equi_sources == True and equi_targets == True:
#            are_equi = True
#        else:
#            are_equi = False
#    if are_equi == True:
#        neighbors1 = mesh1.extend_midedges()
#        neighbors2 = mesh2.extend_midedges()
#        equi_midedges = analogous_midedges(neighbors1, neighbors2,
#                                           enforcerank)
#        if equi_midedges == True:
#            are_equi = True
#        else:
#            are_equi = False
#
#    return are_equi
#
#
#def analogous_nodes(nodelist1, nodelist2, enforcerank=True):
#    """
#    Find whether two lists of nodes contain nodes with
#    same labels.
#    Optionally, nodes may be also required to be at same ranks.
#    (This is comparable to the function "same_objects" used in
#    method "equivalent_meshes".)
#    """
#
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
#
#    return are_equi
#
#
#def analogous_midedges(neighbors1, neighbors2, enforcerank=True):
#    """
#    Find whether two lists of midedges (described as their respective
#    neighbors) connect to event nodes with same labels.
#    Optionally, nodes may be also required to be at same ranks.
#    (This is comparable to the function "equivalent_midedges" used in
#    method "equivalent_meshes".)
#    """
#
#    list1 = neighbors1.copy()
#    list2 = neighbors2.copy()
#    found1 = []
#    found2 = []
#    for i in range(len(list1)):
#        s1 = list1[i]["srcs"]
#        t1 = list1[i]["trgs"]
#        for j in range(len(list2)):
#            s2 = list2[j]["srcs"]
#            t2 = list2[j]["trgs"]
#            if list1[i]["reltype"] == list2[j]["reltype"]:
#                if analogous_nodes(s1, s2, enforcerank):
#                    if analogous_nodes(t1, t2, enforcerank):
#                        found1.insert(0, i)
#                        break
#    for j in range(len(list2)):
#        s2 = list2[j]["srcs"]
#        t2 = list2[j]["trgs"]
#        for i in range(len(list1)):
#            s1 = list1[i]["srcs"]
#            t1 = list1[i]["trgs"]
#            if list2[j]["reltype"] == list1[i]["reltype"]:
#                if analogous_nodes(s2, s1, enforcerank):
#                    if analogous_nodes(t2, t1, enforcerank):
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
#
#    return are_equi

# ================ End of Causal Cores Merging Section ========================

# ......................... Folding Section ...................................

def buildpathways(eoi, causalgraphs=None, siphon=False, ignorelist=[],
                  showintro=False, addedgelabels=True, showedgelabels=True,
                  edgeid=True, edgeocc=False, edgeprob=True, statstype="rel",
                  weightedges=True, color=True, writedot=True, rmprev=False,
                  intropos="top", rulepos="top"):
    """
    Build eventpathway, dualpathway and statepathway.
    What is actually done is:
    - Fold stories into an eventpathway that display only events.
    - Fold dualstories into a dualpathway that displays both events and states.
    - Remove events in dualpathway for statepathway that displays only states.
    The eventpathway could not be obtained by removing states in dualpathways
    because some edge statistics would be lost.
    """

    print("modstory")
    foldpathway(eoi, "modstory", showintro=showintro,
                addedgelabels=addedgelabels, showedgelabels=showedgelabels,
                edgeid=edgeid, edgeocc=edgeocc, edgeprob=edgeprob,
                weightedges=weightedges, color=color, writedot=writedot,
                rmprev=rmprev, intropos=intropos, rulepos=rulepos)
    print("dualstory")
    foldpathway(eoi, "dualstory", showintro=showintro,
                addedgelabels=addedgelabels, showedgelabels=showedgelabels,
                edgeid=edgeid, edgeocc=edgeocc, edgeprob=edgeprob,
                weightedges=weightedges, color=color, writedot=writedot,
                rmprev=rmprev, intropos=intropos, rulepos=rulepos)


def foldpathway(eoi, prefix, causalgraphs=None, siphon=False, ignorelist=[],
                showintro=False, addedgelabels=True, showedgelabels=True,
                edgeid=True, edgeocc=False, edgeprob=True, statstype="rel",
                weightedges=True, color=True, writedot=True, rmprev=False,
                intropos="top", rulepos="top", computenum=True):
    """ Build dual pathway by folding (quotienting) all the stories. """

    # Reading section.
    if causalgraphs == None:
        story_files = get_dot_files("{}/unique".format(eoi), prefix)
        stories = []
        for story_file in story_files:
            story_path = "{}/unique/{}".format(eoi, story_file)
            stories.append(CausalGraph(story_path, eoi))
    else:
        stories = causalgraphs
        story_files = None
    pathway = stories[0]
    foldstory(pathway)
    #for i in range(1, 10):
    for i in range(1, len(stories)):
        print(i+1)
        pathway.occurrence += stories[i].occurrence
        pathway.eventnodes += stories[i].eventnodes
        pathway.statenodes += stories[i].statenodes
        pathway.hyperedges += stories[i].hyperedges
        foldstory(pathway)
    #pathway.align_vertical()
    pathway.hypergraph = True
    # Ranking with intropos="top", rulepos="top" is the fastest as it does
    # not require follow_hyperedges, which takes long on large graphs.
    #pathway.rank_sequentially(intropos="top", rulepos="top")
    #pathway.rank_sequentially(intropos="top", rulepos="bot")
    #pathway.rank_sequentially(intropos="bot", rulepos="top") # <--
    #pathway.rank_sequentially(intropos="bot", rulepos="bot")
    #pathway.get_maxrank()
    pathway.build_adjacency(hyper=True)
    pathway.rank_sequentially(intropos=intropos, rulepos=rulepos)
    # Customize ranks of pathway.
    for node in pathway.eventnodes + pathway.statenodes:
        node.rank = node.rank /2 +0.5
    pathway.get_maxrank()

    # Compute the number of unique instances of precedence relationships.
    if computenum == True:
        initial_files = get_dot_files("{}/tmp".format(eoi), prefix)
        initial_stories = []
        for initial_file in initial_files:
            initial_path = "{}/tmp/{}".format(eoi, initial_file)
            initial_stories.append(CausalGraph(initial_path, eoi))
        hyperedge_instances = []
        for i in range(len(pathway.hyperedges)):
            hyperedge_instances.append([])
        for initial_story in initial_stories:
            for story_hedge in initial_story.hyperedges:
                corresponding_hedge = None
                for i in range(len(pathway.hyperedges)):
                    pathway_hedge = pathway.hyperedges[i]
                    are_equi = equivalent_hyperedges(story_hedge,
                                                     pathway_hedge,
                                                     False, False, True)
                    if are_equi == True:
                        corresponding_hedge = i
                        break
                # Check if the precedence relationship denoted by the story
                # hyperedge is an instance that is not already present in
                # hyperedge_instances[i]
                story_hedge_is_unique = True
                for seen_instance in hyperedge_instances[i]:
                    # Check if story_hedge is the same instance as
                    # seen_instance.
                    is_same = same_instance(story_hedge, seen_instance)
                    if is_same == True:
                        story_hedge_is_unique = False
                        break
                if story_hedge_is_unique == True:
                    hyperedge_instances[i].append(story_hedge)
        for i in range(len(pathway.hyperedges)):
            pathway_hedge = pathway.hyperedges[i]
            n = len(hyperedge_instances[i])
            pathway_hedge.number = n
            for subedge in pathway_hedge.edgelist:
                subedge.number = n

    # Resequentialize node ids. This is required because otherwise the method
    # build_dot_file will draw nodes with same id as if they were merged.
    event_number = 1
    state_number = 1
    for rank_int in range(int(2*pathway.maxrank)+1):
        current_rank = rank_int/2
        for statenode in pathway.statenodes:
            if statenode.rank == current_rank:
                statenode.nodeid = "state{}".format(state_number)
                state_number += 1
        for eventnode in pathway.eventnodes:
            if eventnode.rank == current_rank:
                eventnode.nodeid = "event{}".format(event_number)
                event_number += 1

    # Find path dependent hubs (PDH). Path dependent hubs have different
    # contexts (rule node with 's), with each context having outgoing
    # hyperedges with different rules as sources or targets.
    primed_rules = []
    for eventnode in pathway.eventnodes:
        if eventnode.label[-1] == "'":
            rule = eventnode.label.replace("'", "").strip()
            if rule not in primed_rules:
                primed_rules.append(rule)
    # For each primed_rule, find the group of outgoing edges of each node
    # that correspond to that primed_rule.
    outedge_dict = {}
    for eventnode in pathway.eventnodes:
        primed_rule = None
        for rule in primed_rules:
            if rule in eventnode.label:
                primed_rule = rule
                break
        eventnode.rule = primed_rule
        if primed_rule != None:
            hgroup = []
            for hyperedge in eventnode.outgoing:
                if isinstance(hyperedge.target, EventNode):
                    slbls = []
                    for s in hyperedge.sources:
                        slbls.append(s.label.replace("'", "").strip())
                    tlbl = hyperedge.target.label.replace("'", "").strip()
                    h = {"sources": slbls, "target": tlbl}
                    hgroup.append(h)
                elif isinstance(hyperedge.target, StateNode):
                    # Need to reach one level more up and down to get to
                    # the event nodes.
                    for hyperedge2 in hyperedge.target.outgoing:
                        slbls = []
                        for s in hyperedge2.sources:
                            if s != hyperedge.target:
                                for hyperedge3 in s.incoming:
                                    for s2 in hyperedge3.sources:
                                        slbls.append(s2.label.replace("'", "")
                                                     .strip())
                        tlbl = hyperedge2.target.label.replace("'", "").strip()
                        h = {"sources": slbls, "target": tlbl}
                        hgroup.append(h)
            if primed_rule not in outedge_dict.keys():
                outedge_dict[primed_rule] = [hgroup]
            else:
                outedge_dict[primed_rule].append(hgroup)
    # For a given primed_rule, if all the groups of outgoing edges are the
    # same (ignoring 's), then the primed_rule is NOT a path dependent hub.
    # It is a PDH otherwise.
    pdhs = {}
    pdh_num = 1
    for primed_rule in outedge_dict.keys():
        firstgroup = outedge_dict[primed_rule][0]
        all_same = True
        for i in range(1, len(outedge_dict[primed_rule])):
            # Compare firstgroup with the edge group i.
            secondgroup = outedge_dict[primed_rule][i]
            same_hedge_group = True
            list2_indexes = list(range(len(secondgroup)))
            for h1 in firstgroup:
                same_h = False
                src_cnt1 = collections.Counter(h1["sources"])
                for j in list2_indexes:
                    h2 = secondgroup[j]
                    src_cnt2 = collections.Counter(h2["sources"])
                    if h1["target"] == h2["target"] and src_cnt1 == src_cnt2:
                        list2_indexes.remove(j)
                        same_h = True
                        break
                if same_h == False:
                    same_hedge_group = False
                    break
            if len(list2_indexes) > 0:
                same_hedge_group = False
            if same_hedge_group == False:
                all_same = False
                break
        if all_same == True:
            pdhs[primed_rule] = False
        else:
            pdhs[primed_rule] = pdh_num
            pdh_num += 1
    for eventnode in pathway.eventnodes:
        if eventnode.rule != None:
            eventnode.pdh = pdhs[eventnode.rule]
    # Set state nodes that are outputs of PDH events as PDH themselves.
    # Also assign a number to each state PDH so they can be selected later.
    for eventnode in pathway.eventnodes:
        if eventnode.pdh != False:
            output_state_nodes = []
            for hyperedge in eventnode.outgoing:
                if isinstance(hyperedge.target, StateNode):
                    hyperedge.target.pdh = eventnode.pdh

    # Save a copy of the pathway up to that point. Will be used to create
    # the final pathway where context is ignored.
    pathwaycopy = copy.deepcopy(pathway)

    ##compute_mesh_occurrence(eoi, pathway)
    ##pathway.compute_visuals(showintro, color)
    #pathway.build_nointro()
    #pathway.reverse_subedges()
    #pathway.assign_label_carriers(showintro)
    #pathway.align_vertical()
    ##pathway.compute_relstats()   

    # Write pathway with splited nodes when they have different context. 
    if prefix == "modstory":
        pathway.filename = "eventpathway-split.dot"
    if prefix == "dualstory":
        pathway.filename = "dualpathway-split.dot"
        fuse_multiple_outputs_with_hyperedges(pathway)
    pathway.build_nointro()
    pathway.reverse_subedges()
    pathway.assign_label_carriers(showintro)
    #pathway.align_vertical()
    # Remove 's
    for eventnode in pathway.eventnodes:
        eventnode.label = eventnode.label.replace("'", "").strip()
    pathway.build_dot_file(showintro, addedgelabels, showedgelabels,
                           edgeid, edgeocc, edgeprob, statstype,
                           weightedges, edgewidthscale=1.5)
    output_path = "{}/{}".format(eoi, pathway.filename)
    outfile = open(output_path, "w")
    outfile.write(pathway.dot_file)
    outfile.close()

    if prefix == "dualstory":
        # Write state nodes only on split graph.
        pathway.filename = "statepathway-split.dot"
        pathway.build_adjacency(hyper=True)
        keep_states_only(pathway)
        pathway.build_nointro()
        pathway.build_dot_file(showintro, addedgelabels, showedgelabels,
                               edgeid, edgeocc, edgeprob, statstype,
                               weightedges, edgewidthscale=1.5)
        output_path = "{}/{}".format(eoi, pathway.filename)
        outfile = open(output_path, "w")
        outfile.write(pathway.dot_file)
        outfile.close()

    # Merge nodes corresponding to a same rule (now ignoring 's).
    if prefix == "modstory":
        pathwaycopy.filename = "eventpathway.dot"
    if prefix == "dualstory":
        pathwaycopy.filename = "dualpathway.dot"
        for statenode in pathwaycopy.statenodes:
            statenode.label = statenode.stdedit
            statenode.state = statenode.edit
    for eventnode in pathwaycopy.eventnodes:
        eventnode.label = eventnode.label.replace("'", "").strip()
    foldstory(pathwaycopy)
    pathwaycopy.build_adjacency(hyper=True)
    pathwaycopy.rank_sequentially(intropos=intropos, rulepos=rulepos)
    if prefix == "dualstory":
        # Customize ranks of pathway.
        for node in pathwaycopy.eventnodes + pathwaycopy.statenodes:
            node.rank = node.rank / 2 + 0.5
        pathwaycopy.get_maxrank()
        fuse_multiple_outputs_with_hyperedges(pathwaycopy)
    pathwaycopy.build_nointro()
    pathwaycopy.reverse_subedges()
    pathwaycopy.assign_label_carriers(showintro)
    #pathwaycopy.align_vertical()
    pathwaycopy.build_dot_file(showintro, addedgelabels, showedgelabels,
                               edgeid, edgeocc, edgeprob, statstype,
                               weightedges, edgewidthscale=1.5)
    output_path = "{}/{}".format(eoi, pathwaycopy.filename)
    outfile = open(output_path, "w")
    outfile.write(pathwaycopy.dot_file)
    outfile.close()

    if prefix == "dualstory":
        # Write state nodes only.
        pathwaycopy.filename = "statepathway.dot"
        pathwaycopy.build_adjacency(hyper=True)
        keep_states_only(pathwaycopy)
        pathwaycopy.build_nointro()
        pathwaycopy.build_dot_file(showintro, addedgelabels, showedgelabels,
                                   edgeid, edgeocc, edgeprob, statstype,
                                   weightedges, edgewidthscale=1.5)
        output_path = "{}/{}".format(eoi, pathwaycopy.filename)
        outfile = open(output_path, "w")
        outfile.write(pathwaycopy.dot_file)
        outfile.close()


def fuse_multiple_outputs_with_causaledges(graph):
    """
    Merge state nodes that come from a same rule if they have all the same
    target events. To use on a graph with causaledges.
    """

    # Rebuild adjacency lists using causal edges.
    graph.build_adjacency()
    # Get all the groups of nodes coming from a same event.
    node_groups = []
    for eventnode in graph.eventnodes:
        if eventnode.intro == False:
            new_group = []
            if len(eventnode.outgoing) > 1:
                for edge in eventnode.outgoing:
                    new_group.append(edge.target)
            if len(new_group) > 0:
                # Check if this new node_group is already present.
                already_present = False
                for node_group in node_groups:
                    same_nodes = same_objects(node_group, new_group)
                    if same_nodes == True:
                        already_present = True
                        break
                if already_present == False:
                    node_groups.append(new_group)
    # Check if each node within a group have all the upstream events.
    groups_to_remove = []
    for i in range(len(node_groups)):
        node_group = node_groups[i]
        # Get all upstream events from each node within the group.
        all_up_events = []
        for node in node_group:
            for edge in node.incoming:
                if edge.source not in all_up_events:
                    all_up_events.append(edge.source)
        # Check if each node has an edge from each upstream event.
        for node in node_group:
            all_sources = []
            for edge in node.incoming:
                if edge.source not in all_sources:
                    all_sources.append(edge.source)
            are_same = same_objects(all_sources, all_up_events)
            if are_same == False:
                groups_to_remove.insert(0, i)
                break
    for i in groups_to_remove:
        del(node_groups[i])
    # Check if each node within a group have all the downstream events.
    groups_to_remove2 = []
    for j in range(len(node_groups)):
        node_group = node_groups[j]
        # Get all downstream events from each node within the group.
        all_down_events = []
        for node in node_group:
            for edge in node.outgoing:
                if edge.target not in all_down_events:
                    all_down_events.append(edge.target)
        # Check if each node has an edge to each downstream event.
        for node in node_group:
            all_targets = []
            for edge in node.outgoing:
                if edge.target not in all_targets:
                    all_targets.append(edge.target)
            are_same = same_objects(all_targets, all_down_events)
            if are_same == False:
                groups_to_remove2.insert(0, j)
                break
    for j in groups_to_remove2:
        del(node_groups[j])
    # Merge all the nodes within a group.
    # Effectively, put all states in first node and delete all other
    # nodes along with edges that touch them.
    nodes_to_remove = []
    for node_group in node_groups:
        #first_edit = copy.deepcopy(node_group[0].edit)
        first_state = copy.deepcopy(node_group[0].state)
        for i in range(1, len(node_group)):
            nodes_to_remove.append(node_group[i])
            #for agent in node_group[i].edit:
            #    agentcopy = copy.deepcopy(agent)
            #    first_edit.append(agentcopy)
            for agent in node_group[i].state:
                agentcopy = copy.deepcopy(agent)
                first_state.append(agentcopy)
        #node_group[0].edit = group_sites_by_agent(first_edit)
        node_group[0].state = group_sites_by_agent(first_state)
        lbl = write_context_expression(node_group[0].state, hideid=True)
        node_group[0].label = lbl
    # Delete edges that touch a node to remove.
    edge_indexes = []
    for i in range(len(graph.causaledges)):
        edge = graph.causaledges[i]
        if edge.target in nodes_to_remove or edge.source in nodes_to_remove:
            edge_indexes.insert(0, i)
    for i in edge_indexes:
        del(graph.causaledges[i])
    # Remove nodes.
    node_indexes = []
    for k in range(len(graph.statenodes)):
        if graph.statenodes[k] in nodes_to_remove:
            node_indexes.insert(0, k)
    for k in node_indexes:
        del(graph.statenodes[k])


def fuse_multiple_outputs_with_hyperedges(graph):
    """
    Merge state nodes that come from a same rule if they have all the same
    target events. To use on a graph with hyperedges.
    """

    # Rebuild adjacency lists using hyperedges.
    graph.build_adjacency(hyper=True)
    # Get all the groups of nodes coming from a same event.
    node_groups = []
    for eventnode in graph.eventnodes:
        if eventnode.intro == False:
            new_group = []
            if len(eventnode.outgoing) > 1:
                for hyperedge in eventnode.outgoing:
                    new_group.append(hyperedge.target)
            if len(new_group) > 0:
                # Check if this new node_group is already present.
                already_present = False
                for node_group in node_groups:
                    same_nodes = same_objects(node_group, new_group)
                    if same_nodes == True:
                        already_present = True
                        break
                if already_present == False:
                    node_groups.append(new_group)
    # Check if each node within a group have all the upstream events.
    groups_to_remove = []
    for i in range(len(node_groups)):
        node_group = node_groups[i]
        # Get all upstream events from each node within the group.
        all_up_events = []
        for node in node_group:
            for hyperedge in node.incoming:
                for source in hyperedge.sources:
                    if source not in all_up_events:
                        all_up_events.append(source)
        # Check if each node has an edge from each upstream event.
        for node in node_group:
            all_sources = []
            for hyperedge in node.incoming:
                for source in hyperedge.sources:
                    if source not in all_sources:
                        all_sources.append(source)
            are_same = same_objects(all_sources, all_up_events)
            if are_same == False:
                groups_to_remove.insert(0, i)
                break
    for i in groups_to_remove:
        del(node_groups[i])
    # Check if each node within a group have all the downstream events.
    groups_to_remove2 = []
    for j in range(len(node_groups)):
        node_group = node_groups[j]
        # Get all downstream events from each node within the group.
        all_down_events = []
        for node in node_group:
            for hyperedge in node.outgoing:
                if hyperedge.target not in all_down_events:
                    all_down_events.append(hyperedge.target)
        # Check if each node has an edge to each downstream event.
        for node in node_group:
            all_targets = []
            for hyperedge in node.outgoing:
                if hyperedge.target not in all_targets:
                    all_targets.append(hyperedge.target)
            are_same = same_objects(all_targets, all_down_events)
            if are_same == False:
                groups_to_remove2.insert(0, j)
                break
    for j in groups_to_remove2:
        del(node_groups[j])

    ## Check if each node within a group are together in all their
    ## outgoing hyperedges.
    #groups_to_remove2 = []
    #for j in range(len(node_groups)):
    #    node_group = node_groups[j]
    #    # Get all outgoing hyperedges from each node within the group.
    #    all_outgoing = []
    #    for node in node_group:
    #        for hyperedge in node.outgoing:
    #            if hyperedge not in all_outgoing:
    #                all_outgoing.append(hyperedge)
    #    # Check if each node has all outgoing hyperedges.
    #    for node in node_group:
    #        are_same = same_objects(node.outgoing, all_outgoing)
    #        if are_same == False:
    #            groups_to_remove2.insert(0, j)
    #            break
    #for j in groups_to_remove2:
    #    del(node_groups[j])

    # Merge all the nodes within a group.
    # Effectively, put all states in first node and delete all other
    # nodes along with edges that touch them.
    nodes_to_remove = []
    for node_group in node_groups:
        #first_edit = copy.deepcopy(node_group[0].edit)
        first_state = copy.deepcopy(node_group[0].state)
        for i in range(1, len(node_group)):
            nodes_to_remove.append(node_group[i])
            #for agent in node_group[i].edit:
            #    agentcopy = copy.deepcopy(agent)
            #    first_edit.append(agentcopy)
            for agent in node_group[i].state:
                agentcopy = copy.deepcopy(agent)
                first_state.append(agentcopy)
        #node_group[0].edit = group_sites_by_agent(first_edit)
        node_group[0].state = group_sites_by_agent(first_state)
        lbl = write_context_expression(node_group[0].state, hideid=True)
        node_group[0].label = lbl
    # Delete hyperedge branches that touch a node to remove.
    hyperedge_indexes = []
    for i in range(len(graph.hyperedges)):
        hyperedge = graph.hyperedges[i]
        if hyperedge.target in nodes_to_remove:
            hyperedge_indexes.insert(0, i)
        subedge_indexes = []
        for j in range(len(hyperedge.edgelist)):
            if hyperedge.edgelist[j].source in nodes_to_remove:
                subedge_indexes.insert(0, j)
        for j in subedge_indexes:
            del(graph.hyperedges[i].edgelist[j])
    for i in hyperedge_indexes:
        del(graph.hyperedges[i])
    # Remove nodes.
    node_indexes = []
    for k in range(len(graph.statenodes)):
        if graph.statenodes[k] in nodes_to_remove:
            node_indexes.insert(0, k)
    for k in node_indexes:
        del(graph.statenodes[k])


def colorpaths(eoi, causalgraph=None, siphon=False, ignorelist=[],
               showintro=False, addedgelabels=True, showedgelabels=True,
               edgeid=True, edgeocc=False, edgeprob=True, statstype="rel",
               weightedges=True, color=True, writedot=True, rmprev=False,
               intropos="top", rulepos="top", sel=None):
    """
    Color paths that pass through the chosen PDH (Path Dependent Hub).
    If sel is not a PDH, give a warning message and do not write colored graph.
    """

    print("eventpathway")
    colorgraph(eoi, "eventpathway", showintro=showintro,
               addedgelabels=addedgelabels, showedgelabels=showedgelabels,
               edgeid=edgeid, edgeocc=edgeocc, edgeprob=edgeprob,
               weightedges=weightedges, color=color, writedot=writedot,
               rmprev=rmprev, intropos=intropos, rulepos=rulepos, sel=sel)
    print("dualpathway")
    colorgraph(eoi, "dualpathway", showintro=showintro,
               addedgelabels=addedgelabels, showedgelabels=showedgelabels,
               edgeid=edgeid, edgeocc=edgeocc, edgeprob=edgeprob,
               weightedges=weightedges, color=color, writedot=writedot,
               rmprev=rmprev, intropos=intropos, rulepos=rulepos, sel=sel)


def colorgraph(eoi, prefix, causalgraph=None, siphon=False, ignorelist=[],
               showintro=False, addedgelabels=True, showedgelabels=True,
               edgeid=True, edgeocc=False, edgeprob=True, statstype="rel",
               weightedges=True, color=True, writedot=True, rmprev=False,
               intropos="top", rulepos="top", computenum=True, sel=None):
    """
    Color paths that pass through the chosen PDH (Path Dependent Hub).
    If sel is not a PDH, give a warning message and do not write colored graph.
    """

    # Reading section.
    if causalgraph == None:
        pathway_path = "{}/{}-split.dot".format(eoi, prefix)
        pathway = CausalGraph(pathway_path, eoi)
    else:
        pathway = causalgraph
    # Read PDHs.
    pdhs = []
    for node in pathway.eventnodes + pathway.statenodes:
        if ":" in node.label:
            colon = node.label.index(":")
            pdh_num = int(node.label[colon+1:])
            node.pdh = pdh_num
            if pdh_num not in pdhs:
                pdhs.append(pdh_num)
    if len(pdhs) == 0:
        raise ValueError("No edge coloration required, nothing done.")
    elif len(pdhs) == 1 and sel == None:
        sel = pdhs[0]
    elif sel not in pdhs:
        raise ValueError("Please select a valid PDH number (integer "
                         "after colon : on node labels)")
    # Name output graph according to selected PDH.
    if prefix == "eventpathway":
        pathway.filename = "eventpathway-color{}.dot".format(sel)
    if prefix == "dualpathway":
        pathway.filename = "dualpathway-color{}.dot".format(sel)
    # Select PDH nodes. If state nodes are present, select them instead
    # of the events.
    selectednodes = []
    if len(pathway.statenodes) == 0:
        for eventnode in pathway.eventnodes:
            if eventnode.pdh == sel:
                selectednodes.append(eventnode)
    else:
        for statenode in pathway.statenodes:
            if statenode.pdh == sel:
                selectednodes.append(statenode)
    # For each selected node, get the upstream concurrent path of every
    # target node. If, while building the upstream path, an other selected
    # node is reached, remove it from list selectednodes or remove its
    # corresponding upstream path from paths_nodes if it was already
    # processed.
    # Get the targets of selected nodes.
    pathway.build_adjacency(hyper=True)
    selectedtargets = []
    for node in selectednodes:
        for hyperedge in node.outgoing:
            if hyperedge.target not in selectedtargets:
                selectedtargets.append(hyperedge.target)
    # Get upstream paths (concurrent) from selected targets.
    paths_nodes = []
    paths_hedges = []
    for selectedtarget in selectedtargets:
        up_nodes, up_hedges = concurr_paths_up(pathway, selectedtarget)
        for pn in up_nodes:
            paths_nodes.append(pn)
        for ph in up_hedges:
            paths_hedges.append(ph)
    # Make a list of all path dependent hubs present in each path.
    paths_pdhs = []
    for path in paths_nodes:
        path_pdhs = []
        for node in path:
            if node.pdh != False:
                path_pdhs.append(node)
        paths_pdhs.append(path_pdhs)
    # Manually build a set of all combinations of pdhs found.
    sets_pdhs = []
    for path_pdhs in paths_pdhs:
        pdhs_already_present = False
        for set_pdhs in sets_pdhs:
            same_pdhs = False
            if len(set_pdhs) == len(path_pdhs):
                same_pdhs = True
                for i in range(len(path_pdhs)):
                    #if path_pdhs[i].label != set_pdhs[i].label:
                    if path_pdhs[i] != set_pdhs[i]:
                        same_pdhs = False
                        break
            if same_pdhs == True:
                pdhs_already_present = True
                break
        if pdhs_already_present == False:
            sets_pdhs.append(path_pdhs)
    # Find to which element of sets_pdhs corresponds each element
    # of paths_pdhs.
    set_indexes = []
    for path_pdhs in paths_pdhs:
        for i in range(len(sets_pdhs)):
            set_pdhs = sets_pdhs[i]
            same_pdhs = False
            if len(set_pdhs) == len(path_pdhs):
                same_pdhs = True
                for j in range(len(path_pdhs)):
                    #if path_pdhs[j].label != set_pdhs[j].label:
                    if path_pdhs[j] != set_pdhs[j]:
                        same_pdhs = False
                        break
            if same_pdhs == True:
                set_indexes.append(i)
                break
    # Assign color ids to each hyperedge from each path.
    for hyperedge in pathway.hyperedges:
        hyperedge.color_ids = []
    for i in range(len(paths_hedges)):
        color_id = set_indexes[i] + 1
        path_hedges = paths_hedges[i]
        for hedge in path_hedges:
            if color_id not in hedge.color_ids:
                hedge.color_ids.append(color_id)
    # Redo quotient, but ignoring context.
    #for eventnode in pathway.eventnodes:       
    #    eventnode.label = eventnode.label.replace("'", "").strip()
    for statenode in pathway.statenodes:
        statenode.label = statenode.stdedit
    foldstory(pathway, fusedges=False)
    # Fuse identical hyperedges take colors into account.
    pair_found = True
    while pair_found == True:
        pair_found = False
        for i in range(len(pathway.hyperedges)):
            h1 = pathway.hyperedges[i]
            for j in range(i+1, len(pathway.hyperedges)):
                h2 = pathway.hyperedges[j]
                are_equi, corr = equivalent_hyperedges(h1, h2, False, True)
                if are_equi == True:
                    for k in range(len(h1.edgelist)):
                        main_edge = h1.edgelist[k]
                        other_edge = h2.edgelist[corr[k]]
                        main_edge.weight += other_edge.weight
                        main_edge.number += other_edge.number
                    h1.color_ids += h2.color_ids
                    del(pathway.hyperedges[j])
                    pair_found = True
                    break
            if pair_found == True:
                break
    # Remove color from any hyperedge that has all the colors.
    a = list(range(1, len(sets_pdhs)+1))
    for hyperedge in pathway.hyperedges:
        c = hyperedge.color_ids
        if collections.Counter(c) == collections.Counter(a):
            hyperedge.color_ids = []
    #pathway.rank_sequentially(intropos=intropos, rulepos=rulepos)
    pathway.build_nointro()
    # Assign colors.
    assign_colors(pathway, len(sets_pdhs))
    # Write colored pathway.
    pathway.get_maxrank()
    pathway.build_dot_file(showintro, addedgelabels, showedgelabels,
                           edgeid, edgeocc, edgeprob, statstype,
                           weightedges, edgewidthscale=1.5)
    output_path = "{}/{}".format(eoi, pathway.filename)
    outfile = open(output_path, "w")
    outfile.write(pathway.dot_file)
    outfile.close()

    if prefix == "dualpathway":
        # Write state nodes only on colored graph.
        pathway.filename = "statepathway-color{}.dot".format(sel)
        pathway.build_adjacency(hyper=True)
        keep_states_only(pathway)
        pathway.build_nointro()
        pathway.build_dot_file(showintro, addedgelabels, showedgelabels,
                               edgeid, edgeocc, edgeprob, statstype,
                               weightedges, edgewidthscale=1.5)
        output_path = "{}/{}".format(eoi, pathway.filename)
        outfile = open(output_path, "w")
        outfile.write(pathway.dot_file)
        outfile.close()



def concurr_paths_up(graph, from_node):
    """
    Return a list of all acyclic upstream concurrent paths from a given node.
    """

    fringes_nodes = [[from_node]]
    paths_nodes = [[]]
    paths_hedges = [[]]
    #seen_nodes = [] # To detect loops.
    seen_rules = [] # To stop when a second pdh is reached.
    ends_reached = False
    while ends_reached == False:
        ends_reached = True
        # Get new fringe hyperedges and keep them nested per path.
        nested_fringes_hedges = []
        for path_fringe_nodes in fringes_nodes:
            hedges = []
            for node in path_fringe_nodes:
                if len(node.incoming) > 0:
                    hedges.append(node.incoming)
                    if ends_reached == True:
                        ends_reached = False
            combinations = list(itertools.product(*hedges))
            nested_fringes_hedges.append(combinations)
        # Copy paths up to now if there is more than one fringe_hedges
        # combination for a given path.
        offset = 0
        prev_fringes_nodes = copy.deepcopy(fringes_nodes)
        l = len(fringes_nodes)
        for p in range(l):
            local_hedges = nested_fringes_hedges[p]
            for i in range(len(local_hedges)-1):
                fringe_nodes_copy = fringes_nodes[offset].copy()
                path_nodes_copy = paths_nodes[offset].copy()
                path_hedges_copy = paths_hedges[offset].copy()
                fringes_nodes.insert(offset+i, fringe_nodes_copy)
                paths_nodes.insert(offset+i, path_nodes_copy)
                paths_hedges.insert(offset+i, path_hedges_copy)
            offset += len(local_hedges)
            #if len(local_hedges) > 1:
            #    fringe_nodes_copy = fringes_nodes[p].copy()
            #    path_nodes_copy = paths_nodes[p].copy()
            #    path_hedges_copy = paths_hedges[p].copy()
            #    for i in range(1, len(local_hedges)):
            #        fringes_nodes.insert(offset, fringe_nodes_copy)
            #        paths_nodes.insert(offset, path_nodes_copy)
            #        paths_hedges.insert(offset, path_hedges_copy)
        # Unnest fringes_hedges.
        fringes_hedges = []
        for path_fringe_hedges in nested_fringes_hedges:
            for hedges in path_fringe_hedges:
                fringes_hedges.append(hedges)
        # Add fringe nodes to path nodes.
        for p in range(len(fringes_nodes)):
            path_fringe_nodes = fringes_nodes[p]
            for node in path_fringe_nodes:
                paths_nodes[p].append(node)
        # Find new fringe nodes from fringe hyperedges.
        # Ignore nodes that are already in path to avoid loops.
        fringes_nodes = []
        for p in range(len(fringes_hedges)):
            path_fringe_hedges = fringes_hedges[p]
            path_fringe_node = []
            for hedge in path_fringe_hedges:
                for source in hedge.sources:
                    if source not in paths_nodes[p]:
                        path_fringe_node.append(source)
            fringes_nodes.append(path_fringe_node)
        # Add fringe hyperedges to path edges.
        for p in range(len(fringes_hedges)):
            path_fringe_hedges = fringes_hedges[p]
            for hedge in path_fringe_hedges:
                paths_hedges[p].append(hedge)
    #    if from_node.label == "A phos G":
    #        print("Hedges")
    #        for p in fringes_hedges:
    #            print(p)
    #        print("New nodes")
    #        for p in fringes_nodes:
    #            print(p)
    #        print("Path nodes")
    #        for p in paths_nodes:
    #            print(p)
    #        print("Path hedges")
    #        for p in paths_hedges:
    #            print(p)
    #if from_node.label == "A phos G":
    #    for path in paths_nodes:
    #        print("====")
    #        for n in path:
    #            print(n)

    # Backtracking.
    # For each concurrent path, backtrack down to next pdh if a second
    # pdh of the same rule as the selected node was reached while going up.

    return paths_nodes, paths_hedges


#def get_path(graph, direction, from_node, stop_one_rule=True, block=None,
#             ignore_conflict=False, stop_at_first=False):
#    """ Return a list of all acyclic paths from a given node. """
#
#    all_paths = [[from_node]]
#    seen_rules = [[]]
#    ends_reached = False
#    while ends_reached == False:
#        ends_reached = True
#        for i in range(len(all_paths)):
#            path = all_paths[i]
#            next_nodes = []
#            if direction == "up":
#                for hyperedge in path[-1].incoming:
#                    
#            elif direction == "down":
#                for hyperedge in path[-1].outgoing:
#                    next_nodes.append(hyperedge.target)
#
#
#
#    # //////////////////////////
#    def follow_edges(self, direction, from_node, to_nodes=[], block=None,
#                 ignore_conflict=False, stop_at_first=False):
#    """
#    Return a list of all acyclic paths from a given node to the top of the
#    graph (using direction="up") or to the bottom (using direction="down").
#    If to_nodes are provided, return only the paths that go from from_node
#    to any of the to_nodes.
#    """
#    # This method takes a lot of time on large graphs and can most probably
#    # be improved to speed up calculation.
#    
#    all_paths = [[from_node]]
#    ends_reached = False
#    while ends_reached == False:
#        ends_reached = True
#        for i in range(len(all_paths)):
#            path = all_paths[i]
#            next_nodes = []
#            for edge in self.causaledges:
#                skip = False
#                if ignore_conflict == True:
#                    if edge.relationtype == "conflict":
#                        skip = True
#                if skip == False:
#                    if direction == "up":
#                        if edge.target == path[-1]:
#                            next_nodes.append(edge.source)
#                    elif direction == "down":
#                        if edge.source == path[-1]:
#                            next_nodes.append(edge.target)
#            if len(next_nodes) > 0 and path[-1] not in to_nodes:
#                ends_reached = False
#                if len(next_nodes) > 0:
#                    path_copy = path.copy()
#                path.append(next_nodes[0])
#                for i in range(1, len(next_nodes)):
#                    new_path = path_copy.copy()
#                    new_path.append(next_nodes[i])
#                    all_paths.append(new_path)
#        # Remove looping paths.
#        for i in range(len(all_paths)-1, -1, -1):
#            if len(all_paths[i]) != len(set(all_paths[i])):
#                del(all_paths[i])
#        # Remove paths that end with blocking node if defined.
#        if block != None:
#            for i in range(len(all_paths)-1, -1, -1):
#                if all_paths[i][-1] == block:
#                    del(all_paths[i])
#        # Exit prematurely if only one path is sufficient.
#        if stop_at_first == True:
#            for path in all_paths:
#                if path[-1] in to_nodes:
#                    ends_reached = True
#    # Remove paths that do not end with one of the to_nodes if to_nodes
#    # was defined.
#    if len(to_nodes) > 0:
#        for i in range(len(all_paths)-1, -1, -1):
#            if all_paths[i][-1] not in to_nodes:
#                del(all_paths[i])
#    # Remove the from_node in each path (the first node).
#    for i in range(len(all_paths)):
#        del(all_paths[i][0])
#    
#    return all_paths
#    # ////////////////////////////



#        # !!! THIS IS THE BEGINNING OF THE PART TO CHANGE !!!
#
#        # Initalize color id lists.
#        next_col_id = 1
#        selectednodes = []
#        for hyperedge in pathway.hyperedges:
#            if selectedlabel in hyperedge.target.label:
#                hyperedge.color_ids = [next_col_id]
#                next_col_id += 1
#                if hyperedge.target not in selectednodes:
#                    selectednodes.append(hyperedge.target)
#            else:
#                hyperedge.color_ids = []
#        # Propagate color ids downward 1 step (or two steps if dual story).
#        next_nodes = []
#        for selectednode in selectednodes:
#            incoming_edges = []
#            outgoing_edges = []
#            for hyperedge in pathway.hyperedges:
#                if hyperedge.target == selectednode:
#                    incoming_edges.append(hyperedge)
#                if selectednode in hyperedge.sources:
#                    outgoing_edges.append(hyperedge)
#                    if hyperedge.target not in next_nodes:
#                        next_nodes.append(hyperedge.target)
#            # Find color ids from incoming edges.
#            new_ids = []
#            for inedge in incoming_edges:
#                new_ids.append(inedge.color_ids[0])
#            # Assign new color ids to outgoing edges.
#            # Every outgoing edge takes up all the colors from incoming edges.
#            for outedge in outgoing_edges:
#                outedge.color_ids += new_ids
#        if prefix == "dualstory":
#            current_nodes = next_nodes
#            next_nodes = []
#            for selectednode in current_nodes:
#                incoming_edges = []
#                outgoing_edges = []
#                for hyperedge in pathway.hyperedges:
#                    if hyperedge.target == selectednode:
#                        incoming_edges.append(hyperedge)
#                    if selectednode in hyperedge.sources:
#                        outgoing_edges.append(hyperedge)
#                        if hyperedge.target not in next_nodes:
#                            next_nodes.append(hyperedge.target)
#                # Find color ids from incoming edges.
#                new_ids = []
#                for inedge in incoming_edges:
#                    new_ids.append(inedge.color_ids[0])
#                # Assign new color ids to outgoing edges.
#                # Every outgoing edge takes up all the colors from incoming edges.
#                for outedge in outgoing_edges:
#                    outedge.color_ids += new_ids            
#        startingnodes = next_nodes
#        # Propagate color ids upward. Breadth-first.
#
#
#        current_nodes = startingnodes
#        seen_nodes = startingnodes
#        while len(current_nodes) > 0:
#            next_nodes = []
#            for node in current_nodes:
#                incoming_edges = []
#                outgoing_edges = []
#                for hyperedge in pathway.hyperedges:
#                    if hyperedge.target == node:
#                        incoming_edges.append(hyperedge)
#                        for source in hyperedge.sources:
#                            if source not in seen_nodes:
#                                if source not in next_nodes:
#                                    next_nodes.append(source)
#                                    seen_nodes.append(source)
#                    if node in hyperedge.sources:
#                        outgoing_edges.append(hyperedge)
#                # Find new color ids from outgoing edges.
#                existing_ids = []
#                for inedge in incoming_edges:
#                    existing_ids += inedge.color_ids
#                new_ids = []
#                for outedge in outgoing_edges:
#                    for color_id in outedge.color_ids:
#                        if color_id not in existing_ids:
#                            if color_id not in new_ids:
#                                new_ids.append(color_id)
#                # Assign new color ids to incoming edges.
#                if len(incoming_edges) > 0:
#                    incoming_edges[0].color_ids += new_ids
#                    for i in range(1, len(incoming_edges)):
#                        for j in range(len(new_ids)):
#                            incoming_edges[i].color_ids.append(next_col_id)
#                            next_col_id += 1
#            current_nodes = next_nodes
#
#
#        # Make a final passage from top to bottom.
#        current_nodes = []
#        for node in pathway.eventnodes:
#            if node.first == True:
#                current_nodes.append(node)
#        while len(current_nodes) > 0:
#            next_nodes = []
#            for node in current_nodes:
#                incoming_edges = []
#                outgoing_edges = []
#                for hyperedge in pathway.hyperedges:
#                    if hyperedge.target == node:
#                        incoming_edges.append(hyperedge)
#                    if node in hyperedge.sources:
#                        outgoing_edges.append(hyperedge)
#                        if hyperedge.target not in next_nodes:
#                            if hyperedge.target not in startingnodes:
#                                next_nodes.append(hyperedge.target)
#                # Find new color ids from incoming edges.
#                existing_ids = []
#                for outedge in outgoing_edges:
#                    existing_ids += outedge.color_ids
#                new_ids = []
#                for inedge in incoming_edges:
#                    for color_id in inedge.color_ids:
#                        if color_id not in existing_ids:
#                            if color_id not in new_ids:
#                                new_ids.append(color_id)
#                # Add suplementary colors to outgoing edges.
#                if len(outgoing_edges) > 0:
#                    for outedge in outgoing_edges:
#                        outedge.color_ids += new_ids
#            current_nodes = next_nodes
#
#        # Remove color from any edge that has all the colors.
#        a = list(range(1, next_col_id))
#        for hyperedge in pathway.hyperedges:
#            c = hyperedge.color_ids
#            if collections.Counter(c) == collections.Counter(a):
#                hyperedge.color_ids = []
#
#        # Color edges according their color ids and some color palette.
#        assign_colors(pathway, next_col_id)
#        #color_palette = ["black", "blue1", "firebrick2", "chartreuse3",
#        #                 "gold2", "darkviolet", "darkorange1", "deepskyblue1",
#        #                 "springgreen", "brown2", "magenta", "orange"]
#        #for hyperedge in pathway.hyperedges:
#        #    l = len(hyperedge.color_ids)
#        #    if l > 0:
#        #        if l == 1:
#        #            color_str = color_palette[hyperedge.color_ids[0]]
#        #        else:
#        #            col = color_palette[hyperedge.color_ids[0]]
#        #            color_str = '"{}'.format(col)
#        #            for i in range(1, l):
#        #                col = color_palette[hyperedge.color_ids[i]]
#        #                color_str += ':{}'.format(col)
#        #            color_str += ';{:.2}"'.format(1/float(l))
#        #        hyperedge.color = color_str
#        #        hyperedge.midcolor = color_palette[hyperedge.color_ids[0]]
#        #        for subedge in hyperedge.edgelist:
#        #            subedge.color = color_str
#
#        pathway.build_dot_file(showintro, addedgelabels, showedgelabels,
#                               edgeid, edgeocc, edgeprob, statstype,
#                               weightedges, edgewidthscale=1.5)
#        output_path = "{}/{}".format(eoi, pathway.filename)
#        outfile = open(output_path, "w")
#        outfile.write(pathway.dot_file)
#        outfile.close()
#
#        # !!! THIS IS THE END OF THE PART TO CHANGE !!!
#
#    if prefix == "modstory":
#        for eventnode in pathway.eventnodes:       
#            eventnode.label = eventnode.label.replace("'", "").strip()
#        foldstory(pathway, fusedges=False)
#        # Fuse identical hyperedges take colors into account.
#        pair_found = True
#        while pair_found == True:
#            pair_found = False
#            for i in range(len(pathway.hyperedges)):
#                h1 = pathway.hyperedges[i]
#                for j in range(i+1, len(pathway.hyperedges)):
#                    h2 = pathway.hyperedges[j]
#                    are_equi, corr = equivalent_hyperedges(h1, h2, False, True)
#                    if are_equi == True:
#                        for k in range(len(h1.edgelist)):
#                            main_edge = h1.edgelist[k]
#                            other_edge = h2.edgelist[corr[k]]
#                            main_edge.weight += other_edge.weight
#                            main_edge.number += other_edge.number
#                        h1.color_ids += h2.color_ids
#                        del(pathway.hyperedges[j])
#                        pair_found = True
#                        break
#                if pair_found == True:
#                    break
#        pathway.build_nointro()
#        assign_colors(pathway, next_col_id)
#        pathway.filename = "eventpathway-merged.dot"
#        pathway.build_dot_file(showintro, addedgelabels, showedgelabels,
#                               edgeid, edgeocc, edgeprob, statstype,
#                               weightedges, edgewidthscale=1.5)
#        output_path = "{}/{}".format(eoi, pathway.filename)
#        outfile = open(output_path, "w")
#        outfile.write(pathway.dot_file)
#        outfile.close()
#        # Colorless version of the merged event pathway.
#        for hyperedge in pathway.hyperedges:
#            hyperedge.color = "black"
#            hyperedge.midcolor = "black"
#            for subedge in hyperedge.edgelist:
#                subedge.color = "black"
#        pathway.filename = "eventpathway-merged-black.dot"
#        pathway.build_dot_file(showintro, addedgelabels, showedgelabels,
#                               edgeid, edgeocc, edgeprob, statstype,
#                               weightedges, edgewidthscale=1.5)
#        output_path = "{}/{}".format(eoi, pathway.filename)
#        outfile = open(output_path, "w")
#        outfile.write(pathway.dot_file)
#        outfile.close()

    #return pathway


def assign_colors(pathway, num_colors):
    """
    Assign edge colors according to their color ids and some color palette.
    """

    # Color edges according their color ids and some color palette.
    color_palette = ["black", "blue1", "firebrick2", "chartreuse3",
                     "gold2", "darkviolet", "darkorange1", "deepskyblue1",
                     "springgreen", "brown2", "magenta", "orange"]
    maxcol = num_colors - len(color_palette) + 1
    for hyperedge in pathway.hyperedges:
        l = len(hyperedge.color_ids)
        if l > 0:
            if l == 1:
                if hyperedge.color_ids[0] < len(color_palette):
                    color_str = color_palette[hyperedge.color_ids[0]]
                else:
                    colsec = hyperedge.color_ids[0] - len(color_palette) + 1
                    color_str = ':{:.3f} 1 1'.format(float(col_sec/maxcol))
            else:
                col = color_palette[hyperedge.color_ids[0]]
                color_str = '"{}'.format(col)
                for i in range(1, l):
                    if hyperedge.color_ids[i] < len(color_palette):
                        col = color_palette[hyperedge.color_ids[i]]
                        color_str += ':{}'.format(col)
                    else:
                        colsec = hyperedge.color_ids[0]-len(color_palette)+1
                        color_str += ':{:.3f} 1 1'.format(float(colsec/maxcol))
                color_str += ';{:.2}"'.format(1/float(l))
            hyperedge.color = color_str
            hyperedge.midcolor = color_palette[hyperedge.color_ids[0]]
            for subedge in hyperedge.edgelist:
                subedge.color = color_str


def foldstory(story, fusedges=True):
    """
    Fold (quotient) a dual story based on the label of its nodes. Also
    accumulate the weight of the hyperedges that overlap. 
    """

    # Merge event nodes.
    pair_found = True
    while pair_found == True:
        pair_found = False
        for i in range(len(story.eventnodes)):
            for j in range(i+1, len(story.eventnodes)):
                if story.eventnodes[i].label == story.eventnodes[j].label:
                    # Merging two nodes, effectively deletes a node.
                    merge_nodes([i, j], story.eventnodes, story.hyperedges)
                    pair_found = True
                    break
            if pair_found == True:
                break
    # Merge state nodes.
    pair_found = True
    while pair_found == True:
        pair_found = False
        for i in range(len(story.statenodes)):
            for j in range(i+1, len(story.statenodes)):
                if story.statenodes[i].label == story.statenodes[j].label:
                    # Merging two nodes, effectively deletes a node.
                    merge_nodes([i, j], story.statenodes, story.hyperedges)
                    pair_found = True
                    break
            if pair_found == True:
                break
    # Fuse subedges that are identical within hyperedges.
    for hyperedge in story.hyperedges:
        pair_found = True
        while pair_found == True:
            pair_found = False
            for i in range(len(hyperedge.edgelist)):
                edge1 = hyperedge.edgelist[i]
                for j in range(i+1, len(hyperedge.edgelist)):
                    edge2 = hyperedge.edgelist[j]
                    if edge1.source.label == edge2.source.label:
                        # It seems summing the weights here is incorrect.
                        # We want the edge weights to denote how much a
                        # relationship was useful.
                        #edge1.weight += edge2.weight
                        del(hyperedge.edgelist[j])
                        pair_found = True
                        break
                if pair_found == True:
                    break
        hyperedge.update()
    # Fuse hyperedges that are identical.
    if fusedges == True:
        pair_found = True
        while pair_found == True:
            pair_found = False
            for i in range(len(story.hyperedges)):
                for j in range(i+1, len(story.hyperedges)):
                    are_equi, corr = equivalent_hyperedges(story.hyperedges[i],
                                                           story.hyperedges[j],
                                                           False, True)
                    if are_equi == True:
                        for k in range(len(story.hyperedges[i].edgelist)):
                            main_edge = story.hyperedges[i].edgelist[k]
                            other_edge = story.hyperedges[j].edgelist[corr[k]]
                            main_edge.weight += other_edge.weight
                            main_edge.number += other_edge.number
                        del(story.hyperedges[j])
                        pair_found = True
                        break
                if pair_found == True:
                    break
    for hyperedge in story.hyperedges:
        hyperedge.update()


#def merge_nodes(index1, index2, nodelist, hyperedgelist):
#    """
#    Merge two nodes in a given story. In practice, remove the second node
#    and redirect its edges to the first node.
#    """
#
#    node1 = nodelist[index1]
#    node2 = nodelist[index2]
#    for hyperedge in hyperedgelist:
#        for subedge in hyperedge.edgelist:
#            if subedge.source == node2:
#                subedge.source = node1
#            if subedge.target == node2:
#                subedge.target = node1
#        if hyperedge.target == node2:
#            hyperedge.target = node1
#    del(nodelist[index2])


def merge_nodes(index_list, nodelist, hyperedgelist, deledges=False):
    """
    Merge two nodes in a given story. In practice, remove the second node
    and redirect its edges to the first node.
    """

    node1 = nodelist[index_list[0]]
    for i in range(1, len(index_list)):
        index2 = index_list[i]
        node2 = nodelist[index2]
        for hyperedge in hyperedgelist:
            for subedge in hyperedge.edgelist:
                if subedge.source == node2:
                    subedge.source = node1
                if subedge.target == node2:
                    subedge.target = node1
            if hyperedge.target == node2:
                hyperedge.target = node1
        del(nodelist[index2])


def same_instance(hedge1, hedge2):
    """
    Check of two hyperedges from two different stories are actually the same
    instance in the trace.
    """

    are_same_instance = True
    trg_id1 = "{}{}".format(hedge1.target.label, hedge1.target.eventid)
    trg_id2 = "{}{}".format(hedge2.target.label, hedge2.target.eventid)
    if trg_id1 != trg_id2:
        are_same_instance = False
    if are_same_instance == True:
        src_ids1 = []
        for source1 in hedge1.sources:
            src_ids1.append("{}{}".format(source1.label, source1.eventid))
        sorted_ids1 = sorted(src_ids1)
        src_ids2 = []
        for source2 in hedge2.sources:
            src_ids2.append("{}{}".format(source2.label, source2.eventid))
        sorted_ids2 = sorted(src_ids2)
        if sorted_ids1 != sorted_ids2:
            are_same_instance = False

    return are_same_instance


#def folddualstories(eoi, causalgraphs=None, siphon=False, ignorelist=[],
#                    showintro=False, addedgelabels=True, showedgelabels=True,
#                    edgeid=True, edgeocc=False, edgeprob=True, statstype="rel",
#                    weightedges=True, color=True, writedot=True, rmprev=False):
#    """ Fold (quotient) dual stories into one dual pathway. """
#
#    # Reading section.
#    if causalgraphs == None:
#        story_files = get_dot_files("{}".format(eoi), "unique")
#        stories = []
#        for story_file in story_files:
#            story_path = "{}/{}".format(eoi, story_file)
#            stories.append(CausalGraph(story_path, eoi))
#    else:
#       stories = causalgraphs
#       story_files = None
#    # Doing the work.
#    flush_ignored(stories, story_files, ignorelist)
#    dualpathway = CausalGraph(eoi=eoi, processed=True)
#    dualpathway.occurrence = 0
#    event_number = 1
#    state_number = 1
#    event_labels = []
#    state_labels = []
#    for story in stories:
#        dualpathway.occurrence += story.occurrence
#        # Add event nodes.
#        for eventnode in story.eventnodes:
#            if eventnode.label not in event_labels:
#                event_labels.append(eventnode.label)
#                n_id = "node{}".format(event_number)
#                new_event = EventNode(n_id, eventnode.label, eventnode.rank,
#                                      intro=eventnode.intro,
#                                      first=eventnode.first)
#                dualpathway.eventnodes.append(new_event)
#                event_number += 1
#        # Add states nodes.
#        for statenode in story.statenodes:
#            if statenode.label not in state_labels:
#                state_labels.append(statenode.label)
#                n_id = "state{}".format(state_number)
#                new_state = StateNode(n_id, statenode.label, statenode.rank,
#                                      intro=statenode.intro,
#                                      first=statenode.first)
#                dualpathway.statenodes.append(new_state)
#                state_number += 1
#        # Add hyperedges.
#        for hyperedge in story.hyperedges:
#            hyperedge_found = False
#            for pathwayhedge in dualpathway.hyperedges:
#                equi_edges = equivalent_hyperedges(hyperedge, pathwayhedge,
#                                                   enforcerank=False)
#                if equi_edges == True:
#                    hyperedge_found = True
#                    pathwayhedge.weight += hyperedge.weight
#                    break
#            if hyperedge_found == False:
#                new_hyperedge = HyperEdge(hyperedge.edgelist)
#                dualpathway.hyperedges.append(new_hyperedge)
#    # Propagate new hyperedge weights to their edge lists.
#    for hyperedge in dualpathway.hyperedges:
#        hyperedge.update()
#    # Rerank graph.
#    dualpathway.rank_sequentially()
#    #compute_mesh_occurrence(eoi, pathway)
#    #pathway.compute_visuals(showintro, color)
#    #pathway.compute_relstats()
#    # Write dual pathway.
#    dualpathway.filename = "dualpathway.dot"
#    dualpathway.build_dot_file(showintro, addedgelabels, showedgelabels,
#                               edgeid, edgeocc, edgeprob, statstype,
#                               weightedges)
#    output_path = "{}/{}".format(eoi, dualpathway.filename)
#    outfile = open(output_path, "w")
#    outfile.write(dualpathway.dot_file)
#    outfile.close()


#def foldcores(eoi, causalgraphs=None, siphon=False, ignorelist=[],
#              showintro=False, addedgelabels=True, showedgelabels=True,
#              edgeid=True, edgeocc=False, edgeuse=True, statstype="rel",
#              weightedges=True, color=True, writedot=True, rmprev=False):
#    """ Fold meshed cores into a single event pathway. """
#
#    # Reading section.
#    if causalgraphs == None:
#        # Using cores or eventpaths both work. But it can suffle the nodes
#        # horizontally, yielding a different graph, but with same ranks for
#        # all nodes.
#        core_files = get_dot_files(eoi, "unique")
#        meshedcores = []
#        for core_file in core_files:
#            core_path = "{}/{}".format(eoi, core_file)
#            meshedcores.append(CausalGraph(core_path, eoi))
#    else:
#        meshedcores = causalgraphs
#        core_files = None
#    # Doing the work.
#    flush_ignored(meshedcores, core_files, ignorelist)
#    pathway = CausalGraph(eoi=eoi, meshedgraph=True)
#    pathway.occurrence = 0
#    node_number = 1
#    seen_labels = []
#    midid = 1
#    for meshedcore in meshedcores:
#        # Add nodes.
#        pathway.occurrence += meshedcore.occurrence
#        for node in meshedcore.eventnodes:
#            if node.label not in seen_labels:
#                seen_labels.append(node.label)
#                n_id = "node{}".format(node_number)
#                pathway.eventnodes.append(EventNode(n_id, node.label,
#                                                    node.rank,
#                                                    intro=node.intro,
#                                                    first=node.first))
#                node_number += 1
#        # Add meshes (edges).
#        for mesh in meshedcore.meshes:
#            mesh_found = False
#            for pathwaymesh in pathway.meshes:
#                if analogous_meshes(mesh, pathwaymesh, enforcerank=False):
#                    mesh_found = True
#                    pathwaymesh.uses += mesh.uses
#                    pathwaymesh.weight += mesh.uses
#                    break
#            if mesh_found == False:
#                add_mesh(pathway, mesh, midid)
#                midid += len(mesh.midnodes)
#
#    # Uncomment the next 3 lines and comment pathway.rank_sequentially()
#    # to build unranked version of graph
#    #pathway.rank_intermediary()
#    #pathway.get_maxrank()
#    #pathway.sequentialize_ids()
#    pathway.rank_sequentially()
#    pathway.filename = "eventpathway.dot"
#    compute_mesh_occurrence(eoi, pathway)
#    pathway.compute_visuals(showintro, color)
#    pathway.compute_relstats()
#    pathway.build_dot_file(showintro, addedgelabels, showedgelabels, edgeid,
#                           edgeocc, edgeuse, statstype, weightedges)
#    # Writing section.
#    if writedot == True:
#        output_path1 = "{}/{}".format(eoi, pathway.filename)
#        outfile1 = open(output_path1, "w")
#        outfile1.write(pathway.dot_file)
#        outfile1.close()
#    if rmprev == True:
#        if path_files == None:
#            path_files = get_dot_files(eoi, "meshedcore")
#        for path_file in path_files:
#            file_path = "{}/{}".format(eoi, path_file)
#            os.remove(file_path)
#    print("Merging all event paths into one event pathway.")
#
#    return pathway        


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
           # Find analogous mesh in mapped core template.
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


#def complete_req(graph, mod_nodes):
#    """
#    Extend req_species of each mod node with the req of upstream nodes.
#    """
#
#    mod_no_intro = []
#    for mod_node in mod_nodes:
#        if mod_node.intro == False:
#            mod_no_intro.append(mod_node)
#    for mod_node in mod_no_intro:
#        top_nodes = []
#        for top_node in mod_nodes:
#            if top_node != mod_node:
#                top_nodes.append(top_node)
#        paths = graph.follow_edges("up", mod_node, top_nodes)
#        all_reqs = []
#        for path in paths:
#            path_reqs = []
#            for i in range(len(path)-1):
#                current_node = path[i]
#                #for current_req in current_node.req_species:
#                #    if not species_in(current_req, path_reqs):
#                #        path_reqs.append(current_req.copy())
#                for current_req in current_node.req_species:
#                    cur_ag = current_req["agent"]
#                    cur_site = current_req["site"]
#                    cur_bnd = current_req["bound_agent"]
#                    cur_state = current_req["state"]
#                    add_req = True
#                    for path_req in path_reqs:
#                        path_ag = path_req["agent"]
#                        path_site = path_req["site"]
#                        path_bnd = path_req["bound_agent"]
#                        path_state = path_req["state"]
#                        if path_ag == cur_ag and path_site == cur_site:
#                            if path_bnd == None and cur_bnd == None:
#                                add_req = False
#                            if path_state == None and cur_state == None:
#                                add_req = False
#                    if add_req == True:
#                        path_reqs.append(current_req.copy())
#                up_node = path[i+1]
#                for up_res in up_node.res_species:
#                    res_ag = up_res["agent"]
#                    res_site = up_res["site"]
#                    for path_req in path_reqs:
#                        req_ag = path_req["agent"]
#                        req_site = path_req["site"]
#                        if req_ag == res_ag and req_site == res_site:
#                            if up_res["bound_agent"] != None:
#                                if path_req["bound_agent"] == "_":
#                                    bnd_ag = up_res["bound_agent"]
#                                    bnd_site = up_res["bound_site"]
#                                    path_req["bound_agent"] = bnd_ag
#                                    path_req["bound_site"] = bnd_site
#                                    partner = {"agent": bnd_ag,
#                                               "site": bnd_site,
#                                               "bound_agent": req_ag,
#                                               "bound_site": req_site,
#                                               "state": None}
#                                    path_reqs.append(partner)
#            all_reqs.append(path_reqs)
#
#        comp_req_set = []
#        for path_reqs in all_reqs:
#            for path_req in path_reqs:
#                if not species_in(path_req, comp_req_set):
#                    comp_req_set.append(path_req)
#        mod_node.full_req = comp_req_set


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
