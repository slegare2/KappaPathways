#! /usr/bin/python3

"""
Produce the pathway to EOI from given Kappa simulation.

# Example usage:

from kappa_pathways import KappaPathway

# Set EOI and kappa model.
event_of_interest = "EGFR(Y1092{p})"
kappa_model = "ptyr-model2-act.ka"

# Get pathway using KappaPathway class.
pathway = KappaPathway(event_of_interest, kappa_model,
                       sim_time=3600, seed=235866)
output_path = "path-{}.dot".format(event_of_interest)
output_file = open(output_path, "w")
output_file.write(pathway.dot_file)

"""

import os
import shutil
import json
import random

## Custom story range.
#custom_first = None
#custom_last = None
#
## Get the list of stories in current directory.
#file_list = os.listdir(".")
#num_list = []
#for f in file_list:
#    if "story" in f and "loop" not in f:
#        dash = f.index("-")
#        period = f.index(".")
#        num = int(f[dash+1:period])
#        num_list.append(num)
#first_story = min(num_list)
#last_story = max(num_list)
#if custom_first != None:
#    first_story = custom_first
#if custom_last != None:
#    last_story = custom_last



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
        self.run_kaflow()
        self.get_stories()
        #self.extract_linear_paths()
        #self.push_path()
        self.path_nodes()
        self.path_edges()
        self.count_occur()
        #self.sort_by_rank()
        #self.build_new_dot()
        self.new_dot_simple()


    def add_eoi(self):
        """ Create a new Kappa model where the EOI is added. """

        last_dot = self.kappa_path.rfind(".")
        self.prefix = self.kappa_path[:last_dot]
        self.kappa_eoi_path = "{}-eoi.ka".format(self.prefix)
        self.trace_path = "{}-eoi.json".format(self.prefix)
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
            #simple_dict = {}
            #simple_dict["nodes"] = story_obj.cleared_nodes
            #simple_dict["edges"] = story_obj.cleared_edges
            #end_nodes = self.get_end_nodes(simple_dict)
            #simple_dict["end_nodes"] = end_nodes
            #simple_dict["rule_names"] = story_obj.rule_names
            #simple_name = "simple{}".format(story_name)
            #self.simple_stories[simple_name] = simple_dict
            full_dict = {}
            full_dict["nodes"] = story_obj.story_nodes
            full_dict["edges"] = story_obj.story_edges
            full_dict["rule_names"] = story_obj.rule_names
            self.full_stories[story_name] = full_dict
            ## Optionally, write simplified stories to dot file for debugging.
            #output_path = "{}.dot".format(simple_name)
            #output_file = open(output_path, "w")
            #output_file.write(story_obj.dot_file)
            #output_file.close()


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


#    def get_end_nodes(self, sto_dict):
#        """ Find end nodes as nodes that never appear as sources. """
#
#        end_nodes = []
#        for node in sto_dict["nodes"]:
#            rule_is_source = False
#            for edge in sto_dict["edges"]:
#                if node["node_id"] == edge["source"]:
#                    rule_is_source = True
#                    break
#            if rule_is_source == False:
#                end_nodes.append(node)
#        return end_nodes


#    def extract_linear_paths(self):
#        """
#        Transform every story in possibly many linear paths. That is:
#
#            1 5               1  1  5
#            | |               |  |  |
#            2 6               2  2  6
#           /|/     becomes    |  |  |   
#          7 3                 7  3  3
#          | |                 |  |  |
#          8 4                 8  4  4
#
#        """
#
#        self.linear_paths = []
#        for story_name in self.simple_stories.keys():
#            story = self.simple_stories[story_name]
#            for end_point in story["end_nodes"]:
#                paths = [[end_point]]
#                top_reached = False
#                while top_reached == False:
#                    top_reached = True
#                    for i in range(len(paths)):
#                        if paths[i][-1] != "complete":
#                            top_reached = False
#                            sources = []
#                            for edge in story["edges"]:
#                                if edge["target"] == paths[i][-1]["node_id"]:
#                                    for node in story["nodes"]:
#                                        if node["node_id"] == edge["source"]:
#                                            sources.append(node)
#                            if len(sources) == 0:
#                                paths[i].append("complete")
#                            else:
#                                for j in range(1, len(sources)):
#                                    paths.append(paths[i].copy())
#                                paths[i].append(sources[0])
#                                for k in range(1, len(sources)):
#                                    paths[-k].append(sources[k])
#                for path in paths:
#                    reverse_path = []
#                    for i in range(len(path)-2, -1, -1):
#                        reverse_path.append(path[i])
#                    self.linear_paths.append(reverse_path)
#        #for path in self.linear_paths:
#        #    print("---")
#        #    print(path)
#            

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

    def path_nodes(self):
        """ Get all possible node from every story. """

        self.pathway_nodes = []
        self.pathway_node_ids = {}
        self.pathway_rule_names = {}
        pathway_rules = []
        for story_name in self.full_stories.keys():
            story = self.full_stories[story_name]
            for node in story["nodes"]:
                rule_name = story["rule_names"][node]
                if " ina" not in rule_name and " unbinds" not in rule_name:
                    if rule_name not in pathway_rules:
                        pathway_rules.append(rule_name)
        #print(len(pathway_rules))
        #for r in pathway_rules:
        #    print(r)
        node_id = 1
        for rule in pathway_rules:
            self.pathway_nodes.append({"node_id": node_id,
                                       "rule_name": rule})
            self.pathway_node_ids[rule] = node_id
            self.pathway_rule_names[node_id] = rule
            node_id += 1


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


#    def sort_by_rank(self):
#        """ Group story nodes by rank. """
#
#        self.sorted_nodes = {}
#        self.sorted_edges = []
#        for rank in range(1, self.end_rank+1):
#            node_list = []
#            for node in self.pathway_nodes:
#                if node["rank"] == rank:
#                    node_list.append(node)
#                    for edge in self.pathway_edges:
#                        if edge["source"] == node["node_id"]:
#                            self.sorted_edges.append(edge)
#            rank_str = "{}".format(rank)
#            self.sorted_nodes[rank_str] = node_list


#    def build_new_dot(self):
#        """ Build the dot file that contains the new story. """
#
#        dot_str = "digraph G{\n"
#        dot_str += "  ranksep=.3 ;\n"
#        for rank in range(1, self.end_rank+1):
#            rank_str = "{}".format(rank)
#            dot_str += ('{{ rank = same ; "{}" [shape=plaintext] ;\n'
#                        .format(rank_str))
#            node_list = self.sorted_nodes[rank_str]
#            for node in node_list:
#                if node["rule_name"] == '"{}"'.format(self.eoi):
#                    node_shape = "ellipse"
#                    node_color = "indianred2"
#                else:
#                    node_shape = "invhouse"
#                    # Compute node colors based on counts
#                    count = self.rule_counts_tot[node["rule_name"]]
#                    high = self.max_rule_count_tot
#                    frac = count / high
#                    #print("{:40} {}".format(node["rule_name"], count))
#                    #node_color = '"0.55 {} 1.000"'.format(frac)
#                    node_color = '"0.55 0.6 1.000"'.format(frac)
#                dot_str += '"node{}" [label={}, '.format(node["node_id"],
#                                                   node["rule_name"]) 
#                dot_str += 'shape={}, style=filled, '.format(node_shape)
#                dot_str += 'fillcolor={}] ;\n'.format(node_color)
#            dot_str += "}\n"
#        for rank in range(1, self.end_rank):
#            rank_str = "{}".format(rank)
#            next_rank = "{}".format(rank+1)
#            dot_str += ('"{}" -> "{}" [style="invis"];\n'
#                        .format(rank_str, next_rank))
#        for edge in self.sorted_edges:
#            source_rule = self.pathway_rule_names[edge["source"]]
#            target_rule = self.pathway_rule_names[edge["target"]]
#            edge_rules = "{}-{}".format(source_rule, target_rule)
#            #count = self.edge_counts[edge_rules]
#            #low = self.min_edge_count
#            #high = self.max_edge_count
#            count = self.edge_uses[edge_rules]
#            low = self.min_edge_use
#            high = self.max_edge_use
#            minpen = 1
#            maxpen = 10
#            dpen = maxpen-minpen
#            frac = (count-low) / (high-low)
#            weight = minpen + dpen * frac
#            pensize = "{}".format(weight)
#            #pensize = "1"
#            dot_str += ('"node{}" -> "node{}" [penwidth={}] ;\n'
#                        .format(edge["source"], edge["target"], pensize))
#        dot_str += "}"
#        self.dot_file = dot_str


    def new_dot_simple(self):
        """
        Write dot file without caring for node rank to see if graphviz
        makes the correct layout by himself.
        """

        random.shuffle(self.pathway_nodes)
        random.shuffle(self.pathway_edges)

        dot_str = "digraph G{\n"
        dot_str += "  ranksep=.3 ;\n"
        for node in self.pathway_nodes:
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
        for edge in self.pathway_edges:
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
