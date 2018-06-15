# DTLMedian.py
# Written July 2017 by Andrew Ramirez and Eli Zupke
# Utilizes previous code to find the "median" reconciliation of a pair of gene and species trees
# and other related information


# -1. DATA STRUCTURE QUICK REFERENCE:
#
#
#   DTL or Median Reconciliation graph:
#       { mapping_node: [event1, event2, ... eventn, number], ...}
#   Event:
#       ('event_type', child_mapping_node1, child_mapping_node2)
#
#   Mapping node (case indicates capitalization standards):
#       ('gene_node', 'SPECIES_NODE')
#   or in loss or contemporary event nodes:
#       (None, None)
#
#
#   (edge) trees:
#       {('R','N'): ('R','N', ('N','C1'), ('N','C2')) ...}
#       aka:
#       {root_edge: (root_edge[0], root_edge[1], child1_edge, child2_edge) ...}
#
#   vertex_trees:
#       {'N':('C1','C2') ...}
#

import optparse
from operator import itemgetter
import numpy as np
import DTLReconGraph
import Diameter
import csv
import os



def mapping_node_sort(ordered_gene_node_list, ordered_species_node_list, mapping_node_list):
    """
    :param ordered_gene_node_list: an ordered dictionary of the gene nodes, where each key is a node
    and the corresponding values are children of that node in the gene tree. Note that the order (pre-
    or post-) in which this tree is passed determines the final order - a preorder gene node list will
    return mapping nodes sorted in preorder. The species and gene orderings must match.
    :param ordered_species_node_list: same as for the gene node list above, except for the species tree
    :param mapping_node_list: a list of all mapping nodes within a given reconciliation graph
    :return: the given mapping nodes except sorted in the order corresponding to the order in which
    the species and gene nodes are passed in (see description of ordered_gene_node_list for more on this).
    The returned mapping node list is sorted first by gene node and then by species node
    """

    # In order to sort the mapping nodes, we need a way to convert them into numbers. These two lookup tables allow
    # us to achieve a lexicographical ordering with gene nodes more significant than species nodes.
    gene_level_lookup = {}
    species_level_lookup = {}

    # By multiplying the gene node keys by the number of species nodes, we can ensure that a mapping node with a later
    # gene node always comes before one with a later species node, because a gene node paired with the last species
    # node will be one level less than the next gene node paired with the first species node.
    gene_multiplier = len(ordered_species_node_list)

    for i1, gene_node in enumerate(ordered_gene_node_list):
        gene_level_lookup[gene_node] = i1 * gene_multiplier

    for i2, species_node in enumerate(ordered_species_node_list):
        species_level_lookup[species_node] = i2 

    # The lambda function looks up the level of both the gene node and the species nodes and adds them together to
    # get a number to give to the sorting algorithm for that mapping node. The gene node is weighted far more heavily
    # than the species node to make sure it is always more significant.
    sorted_list = sorted(mapping_node_list, key=lambda node: gene_level_lookup[node[0]] + species_level_lookup[node[1]])

    return sorted_list


def generate_scores(preorder_mapping_node_list, dtl_recon_graph, gene_root):
    """
    Computes frequencies for every event
    :param preorder_mapping_node_list: A list of all mapping nodes in DTLReconGraph in double preorder
    :param dtl_recon_graph: The DTL reconciliation graph that we are scoring
    :param gene_root: The root of the gene tree
    :return: 0. A file structured like the DTLReconGraph, but with the lists of events replaced
                with dicts, where the keys are the events and the values are the scores of those events, and
             1. The number of MPRs in DTLReconGraph.
    """

    # Initialize the dictionary that will store mapping node and event counts (which also acts as a memoization
    # dictionary)
    counts = dict()

    # Initialize the very start count, for the first call of countMPRs
    count = 0

    # Loop over all given minimum cost reconciliation roots
    for mapping_node in preorder_mapping_node_list:
        if mapping_node[0] == gene_root:

            # This will also populate the counts dictionary with the number of MPRs each event and mapping node is in
            count += count_mprs(mapping_node, dtl_recon_graph, counts) # RR - this sums up each root node's frequency in reconciliations. Since each reconciliation starts at a root, we are essentially counting the total number of reconciliations in our MPR

    # Initialize the scores dict. This dict contains the frequency score of each mapping node
    scores = dict()
    for mapping_node in preorder_mapping_node_list:
        scores[mapping_node] = 0.0

    # This entry is going to be thrown away, but it seems neater to just let calculateScoresOfChildren
    # add scores to an unused entry than to check to see if they are (None, None) in the first place.
    scores[(None, None)] = 0.0

    # The scored graph is like the DTLReconGraph, except instead of individual events being in a list, they are the
    # keys of a dictionary where the values are the frequency scores of those events. So, event_scores takes event
    # nodes as keys and (after being filled below) has the frequencies of those events in MPRs as the values
    event_scores = {}

    for mapping_node in preorder_mapping_node_list:

        # If we are at the root of the gene tree, then we need to initialize the score entry
        if mapping_node[0] == gene_root:
            scores[mapping_node] = counts[mapping_node]
        # This fills up the event scores dictionary
        calculate_scores_for_children(mapping_node, dtl_recon_graph, event_scores, scores, counts)

    # Normalize all of the event_scores
    for mapping_node in preorder_mapping_node_list:
        for event in dtl_recon_graph[mapping_node]:
            if event == ('C',(None,None),(None,None)):           #we don't serve your kind around here
                continue
            event_scores[event] = event_scores[event] / float(count)
    
    event_scores[('C',(None,None),(None,None))] = event_scores[('C',(None,None),(None,None))]/float(count)

    return event_scores, count


def count_mprs(mapping_node, dtl_recon_graph, counts):
    """
    :param mapping_node: an individual mapping node that maps a node
    for the parasite tree onto a node of the host tree, in the format
    (p, h), where p is the parasite node and h is the host node
    :param dtl_recon_graph: A DTL reconciliation graph (see data structure quick reference at top of file)
    :param counts: a dictionary representing the running memo that is passed
    down recursive calls of this function. At first it is just an empty
    dictionary (see above function), but as it gets passed down calls, it collects
    keys of mapping nodes or event nodes and values of MPR counts. This memo improves runtime
    of the algorithm
    :return: the number of MPRs spawned below the given mapping node in the graph
    """

    # Search the counts dictionary for previously calculated results (this is the memoization)
    if mapping_node in counts:
        return counts[mapping_node] # RR - this method is memoization, NOT dynamic programming. That's why we're going top-down

    # Base case, occurs if being called on a child produced by a loss or contemporary event 
    if mapping_node == (None, None): # RR - we check for losses or contemporary events first so that when we loop through events later (line 159-160) we know those will only be duplications, speciations or transfers
        return 1 

    # Initialize a variable to keep count of the number of MPRs
    count = 0

    # Loop over all event nodes corresponding to the current mapping node
    for eventNode in dtl_recon_graph[mapping_node]: # RR - do we know that all of these events are only duplications, speciations and transfers? because only for these are we guaranteed two children. Edit: see line 153

        # Save the children produced by the current event
        mapping_child1 = eventNode[1]
        mapping_child2 = eventNode[2]

        # Add the product of the counts of both children (over all children) for this event to get the parent's count
        counts[eventNode] = count_mprs(mapping_child1, dtl_recon_graph, counts) * count_mprs(mapping_child2,
                                                                                             dtl_recon_graph, counts)
        count += counts[eventNode]

    # Save the result in the counts
    counts[mapping_node] = count

    return count


def calculate_scores_for_children(mapping_node, dtl_recon_graph, event_scores, mapping_scores, counts):
    """
    This function calculates the frequency score for every mapping node that is a child of an event node that is a
    child of the given mapping node, and stores them in scoredGraph. RR - it also first calculates the frequency score of the event nodes that are children of the given mapping nodes.
    :param mapping_node: The mapping node that is the parent of the two scores we will compute
    :param dtl_recon_graph: The DTL reconciliation graph (see data structure quick reference at top of file)
    :param event_scores: The scored DTL reconciliation graph (see data structure quick reference at top of file)
    :param mapping_scores: The score for each mapping node (which will ultimately be thrown away) that this function
    helps build up
    :param counts: The counts generated in countMPRs (from the bottom-up). Note that the counts are filled during a
    bottom-up traversal, and the scores are filled in during a top-down traversal after the counts
    :return: Nothing, but scoredGraph is built up.
    """

    assert mapping_scores[mapping_node] != 0, "Sorting error! Ensure that parents are calculated before children" # RR - probably irrelevant, but are we sure that the mapping score always !=0 if it's been calculated? Couldn't it be that that node is just never used?

    # This multiplier results in  counts[event_node] / counts[mapping_node] for each event node, which is the % of
    # this mapping node's scores (scores[mapping_node]) that it gives to each event node.
    multiplier = float(mapping_scores[mapping_node]) / counts[mapping_node]

    # Iterate over every event
    for event_node in dtl_recon_graph[mapping_node]:

        event_scores[event_node] = multiplier * counts[event_node]

        # Save the children produced by the current event
        mapping_child1 = event_node[1]
        mapping_child2 = event_node[2]
        mapping_scores[mapping_child1] += event_scores[event_node]
        mapping_scores[mapping_child2] += event_scores[event_node]


def compute_median(dtl_recon_graph, event_scores, postorder_mapping_nodes, mpr_roots):
    """
    :param dtl_recon_graph: A dictionary representing a DTL Recon Graph.
    :param event_scores: A dictionary with event nodes as keys and values corresponding to the frequency of
    that events in MPR space for the recon graph
    :param postorder_mapping_nodes: A list of the mapping nodes in a possible MPR, except sorted first in
    postorder by species node and postorder by gene node
    :param mpr_roots: A list of mapping nodes that could act as roots to an MPR for the species and
    gene trees in question, output from the findBestRoots function in DTLReconGraph.py
    :return: A new dictionary which is has the same form as a DTL reconciliation graph except every
    mapping node only has one event node, along with the number of median reconciliations for the given DTL
    reconciliation graph, as well as the root(s) of the median MPR for the given graph. Thus, this graph will
    represent a single reconciliation: the median reconciliation. (RR - would a better way to phrase this be - our graph will represent a subgraph of our DTL Recon Graph, such that all reconciliations possible to reconcile from the subgraph are median reconciliations)?
    """

    # Note that for a symmetric median reconciliation, each frequency must have 0.5 subtracted from it

    # Initialize a dict that will store the running total frequency sum incurred up to the given mapping node,
    # and the event node that directly gave it that frequency sum. Keys are mapping nodes, values are tuples
    # consisting of a list of event nodes that maximize the frequency - 0.5 sum score for the lower level,
    # and the corresponding running total frequency - 0.5 sum up to that mapping node
    sum_freqs = dict()

    # Loop over all mapping nodes for the gene tree
    for map_node in postorder_mapping_nodes:

        # Contemporaneous events need to be caught from the get-go
        if dtl_recon_graph[map_node] == [('C', (None, None), (None, None))]:
            sum_freqs[map_node] = ([('C', (None, None), (None, None))], 0.5)  # C events have freq 1, so 1 - 0.5 = 0.5
            continue  # Contemporaneous events should be a lone event in a list, so we move to the next mapping node

        # Get the events for the current mapping node and their running (frequency - 0.5) sums, in a list
        events = list()
        for event in dtl_recon_graph[map_node]:

            # Note that 'event' is of the form: ('event ID', 'Child 1', 'Child 2'), so the 0th element is the event
            # ID and the 1st and 2nd elements are the children produced by the event
            if event[0] == 'L':  # Losses produce only one child, so we only need to look to one lower mapping node
                events.append((event, sum_freqs[event[1]][1] + event_scores[event] - 0.5))
            else:  # Only other options are T, S, and D, which produce two children
                events.append((event, sum_freqs[event[1]][1] + sum_freqs[event[2]][1] + event_scores[event] - 0.5))

        # Find and save the max (frequency - 0.5) sum
        max_sum = max(events, key=itemgetter(1))[1]

        # Initialize list to find all events that gives the current mapping node the best (freq - 0.5) sum
        best_events = list()

        # Check to see which event(s) produce the max (frequency - 0.5) sum
        for event in events:
            if event[1] == max_sum:
                best_events.append(event[0])

        # Help out the garbage collector by discarding the now-useless non-optimal events list
        del events

        # Save the result for this mapping node so it can be used in higher mapping nodes in the graph
        sum_freqs[map_node] = (best_events[:], max_sum) # RR - best_events[:] slices the list and is a short cut to making a copy of the list that we want, which is why we think they did it

    # Get all possible roots of the graph and their running frequency scores, in a list, for later use
    possible_root_combos = [(root, sum_freqs[root][1]) for root in mpr_roots]

    # Find the best frequency - 0.5 sum for all of the potential roots for the median
    best_sum = max(possible_root_combos, key=itemgetter(1))[1]

    # Find all of the root combos for a median by filtering out the roots that don't give the best freq - 0.5 sum
    best_root_combos = list(filter(lambda x: x[1] == best_sum, possible_root_combos))

    # Extract just the roots from the previously filtered out list
    best_roots = [root[0] for root in best_root_combos]

    # Adjust the sum_freqs dictionary so we can use it with the buildDTLReconGraph function from DTLReconGraph.py
    for map_node in sum_freqs:

        # We place the event tuples into lists so they work well with the diameter algorithm
        sum_freqs[map_node] = sum_freqs[map_node][0]  # Only use the events, no longer the associated frequency sum

    # Use the buildDTLReconGraph function from DTLReconGraph.py to find the median recon graph
    # Note that build_dtl... requires a list of the best roots for a reconciliation graph, the events for each
    # mapping node that are viable for an MPR (in our case, the median), and an empty dicitonary to populate
    # as the final return value
    med_recon_graph = DTLReconGraph.build_dtl_recon_graph(best_roots, sum_freqs, {})

    # Check to make sure the median is a subgraph of the DTL reconciliation
    assert check_subgraph(dtl_recon_graph, med_recon_graph), 'Median is not a subgraph of the recon graph!'

    # We can use this function to find the number of medians once we've got the final median recon graph
    n_med_recons = DTLReconGraph.count_mprs_wrapper(best_roots, med_recon_graph)

    return med_recon_graph, n_med_recons, best_roots


def check_subgraph(recon_graph, subrecon):
    """
    :param recon_graph: A reconciliation graph
    :param subrecon: Another reconciliation graph, the one which is supposed to be a subgraph of "recon_graph"
    :return: a boolean value: True if the "subrecon" is really a subgraph of "recon_graph",
    False otherwise
    """

    # Loop over all mapping nodes contained in the median reconciliation graph
    for map_node in subrecon:

        # Loop over mapping nodes
        if map_node not in recon_graph:
            return False
        else:

            # Now events for a given mapping node
            for event in subrecon[map_node]:
                if event not in recon_graph[map_node]:
                    return False
    return True


def choose_random_median_wrapper(median_recon, med_roots, count_dict):
    """
    :param median_recon: the median reconciliation
    :param med_roots: the roots (root mapping nodes) for possible median reconciliations
    :param count_dict: a dictionary detailing how many medians can stem from an individual event
    node
    :return: a randomly, uniformly sampled median reconciliation graph
    """

    # Find the total amount of medians that can stem from the roots
    total_meds = 0.0
    for median_root in med_roots:
        total_meds += count_dict[median_root]

    # Create the choice list for the roots we can choose from, weighted to account for median
    # counts each root can produce
    # Note that numpy is so basic that we need to do a convoluted workaround to choose tuples from a list
    final_root = med_roots[np.random.choice(len(med_roots), p=[count_dict[med_root] / total_meds for med_root in
                                                               med_roots])]

    return choose_random_median(median_recon, final_root, count_dict)


def choose_random_median(median_recon, map_node, count_dict):
    """
    :param median_recon: the full median reconciliation graph, as returned by compute_median
    :param map_node: the current mapping node in the median reconciliation that we're trying
    to find a path from. In the first call, this mapping node will be one of the root mapping
    nodes for the median reconciliation graph, randomly selected
    :param count_dict: a dictionary that tells us how many total medians a given event node can spawn
    :return: a single-path reconciliation graph that is a sub-graph of the median. It is chosen
    randomly but randomly in such a way that event node choice will favor event nodes that lead
    to more MPRs so that the data aren't skewed
    """

    # Initialize the dictionary that will store the final single-path median that we choose
    random_submedian = dict()

    # Find the total number of medians we can get from the current mapping node
    total_meds = float(count_dict[map_node])

    # Use a convoluted numpy workaround to select tuples (events) from a list, taking into account
    # how many medians each event can produce
    next_event = median_recon[map_node][np.random.choice(len(median_recon[map_node]),
                                                         p=[count_dict[event] / total_meds for event in
                                                            median_recon[map_node]])]

    random_submedian.update({map_node: [next_event]})

    # Check for a loss
    if next_event[0] == 'L':
        random_submedian.update(choose_random_median(median_recon, next_event[1], count_dict))

    # Check for events that produce two children
    elif next_event[0] in ['T', 'S', 'D']:
        random_submedian.update(choose_random_median(median_recon, next_event[1], count_dict))
        random_submedian.update(choose_random_median(median_recon, next_event[2], count_dict))

    # Make sure our single path median is indeed a subgraph of the median
    assert check_subgraph(median_recon, random_submedian), 'The randomly chosen single-path median is not a subgraph ' \
                                                           'of the full median!'

    return random_submedian

# Functions that are going to be used for the second median reconciliation (with different measure of distance)

#Function that returns the vertices (nodes_) in preorder of a tree in vertex format given the tree and its root as arguments
# Recall that vertex format of a tree is an ordered dictionary where keys are vertices which relate to tuples, where each tuple consists
# of the vertice's two children

def preorder_vertices(tree, root):
    left_child = tree[root][0]
    right_child = tree[root][1]
        # Base case
    if left_child is None:  # Then right_child == None also
        return [root]
    else:  # Recursive call
        return [root] + \
            preorder_vertices(tree,left_child) + \
            preorder_vertices(tree,right_child)
    
#Function that returns the vertices (nodes_) in post order of a tree in vertex format given the tree and its root as arguments
# Recall that vertex format of a tree is an ordered dictionary where keys are vertices which relate to tuples, where each tuple consists
# of the vertice's two children

def postorder_vertices(tree, root):
    left_child = tree[root][0]
    right_child = tree[root][1]
        # Base case
    if left_child is None:  # Then right_child == None also
        return [root]
    else:  # Recursive call
        return postorder_vertices(tree,left_child) + \
            postorder_vertices(tree,right_child) + \
            [root]

# Function that generates scores for each of the events, telling us the frequency of their appearance in MPRs from the DTL Reconciliation Graph
# Essentially a copy of the generating scores function somewhere above except without normalizing
# Will probably eventually combine the two functions by adding a 'to normalize or not' parameter as an argument
def generate_scores_no_normalizing(preorder_mapping_node_list, dtl_recon_graph, gene_root):
    """
    Computes frequencies for every event
    :param preorder_mapping_node_list: A list of all mapping nodes in DTLReconGraph in double preorder
    :param dtl_recon_graph: The DTL reconciliation graph that we are scoring
    :param gene_root: The root of the gene tree
    :return: 0. A file structured like the DTLReconGraph, but with the lists of events replaced
                with dicts, where the keys are the events and the values are the scores of those events, and
             1. The number of MPRs in DTLReconGraph.
    """

    # Initialize the dictionary that will store mapping node and event counts (which also acts as a memoization
    # dictionary)
    counts = dict()

    # Initialize the very start count, for the first call of countMPRs
    count = 0

    # Loop over all given minimum cost reconciliation roots
    for mapping_node in preorder_mapping_node_list:
        if mapping_node[0] == gene_root:

            # This will also populate the counts dictionary with the number of MPRs each event and mapping node is in
            count += count_mprs(mapping_node, dtl_recon_graph, counts) # RR - this sums up each root node's frequency in reconciliations. Since each reconciliation starts at a root, we are essentially counting the total number of reconciliations in our MPR

    # Initialize the scores dict. This dict contains the frequency score of each mapping node
    scores = dict()
    for mapping_node in preorder_mapping_node_list:
        scores[mapping_node] = 0.0

    # This entry is going to be thrown away, but it seems neater to just let calculateScoresOfChildren
    # add scores to an unused entry than to check to see if they are (None, None) in the first place.
    scores[(None, None)] = 0.0

    # The scored graph is like the DTLReconGraph, except instead of individual events being in a list, they are the
    # keys of a dictionary where the values are the frequency scores of those events. So, event_scores takes event
    # nodes as keys and (after being filled below) has the frequencies of those events in MPRs as the values
    event_scores = {}

    for mapping_node in preorder_mapping_node_list:

        # If we are at the root of the gene tree, then we need to initialize the score entry
        if mapping_node[0] == gene_root:
            scores[mapping_node] = counts[mapping_node]
        # This fills up the event scores dictionary
        calculate_scores_for_children(mapping_node, dtl_recon_graph, event_scores, scores, counts)

    return event_scores, count
                                                                                                            
# A function we use to create a dict of the numPlacing function applied to all mapping nodes (g,s).
#numPlacing() gives us the number of reconciliations that have g placed on s (i.e, a loss event doesn't directly follow the node (g,s)) 
# numPlacing can be calculated by first going through all the events of a given mapping node (g,s). 
# Then, we look at only the events that are not losses, and add together all of their frequencies. 
# This will give us the net number of reconciliations that use a given node (g,s) where g is mapped onto s 
# (because if the node (g,s) exists and g is NOT mapped onto s, it must be because there is a loss event following that mapping node in the reconciliation).
# In our arguments, count is the total number of reconciliations that we get from our reconciliation graph. This value is returned in the generate_scores function as well as in the reconcile function
# NOTE: I don't think it matters whether the ordered species node list and gene list are in preorder or in postorder. I'm not even sure it matters if they are ordererd
# Over here, we will use postorder by convention, simply because the mapping node list is also in postorder.

def construct_numPlacing_table(dtl_recon_graph, ordered_species_node_list, ordered_gene_node_list, postorder_mapping_node_list, event_scores, count):
    # Initialize the numPlacing dict. This dict contains the numPlacing score of each mapping node 
    numPlacing = dict()
    # We want values in this dictionary for every possible combination of gene node x species node, not just the mapping nodes. But the values
    #for the combinations which are not valid mapping nodes will be trivially 0 (which makes our calculations simpler). 
    for gene_node in ordered_gene_node_list:
        for species_node in ordered_species_node_list:
            numPlacing[(gene_node, species_node)] = 0

    # Now, we will just change the values for the ones that are mapping nodes, and may have non-zero values
    for map_node in postorder_mapping_node_list:
    # First checking for leaf nodes, which will have only contemporaneous events following them.
        if dtl_recon_graph[map_node] == [('C', (None, None), (None, None))]:    
            numPlacing[map_node] = count # Because leaf nodes of the form (g,s) are present in all reconciliations and g is always placed on s.
            continue
    # Not sure, but I think we should be able to do this step directly without checking for contemporaneous events first. 
    # TO DO: Check if contemporaneous events are given a frequency of count in generate_scores. If so, we could use just the following step:
        for event in dtl_recon_graph[map_node]:
            if event[0] != 'L': # Making sure the relevant event is not a loss
                numPlacing[map_node] += event_scores[event] # Adding the frequencies of events of (g,s) that are NOT loss events.

    return numPlacing

# A function that calculates the distance between any two nodes on a given tree, where distance is defined as the
# number of edges in the path from one node to the other (since this is a tree, we know there is but one unique path from one node to any other)
# This function will be defined recursively
# Recall that species_tree is in vertex tree format (see DTLMedian) and that it is in the form of a dictionary. 
# The key is nodes and the correspondence are the children of the key node.
def distance(node1, node2, species_tree, ancestral_table):
    # Base case; when both nodes are the same, distance between them is 0.
    if node1 == node2 :
        return 0
    # We want to check if node1 is an ancestor of node2 or vice versa
    elif ancestral_table[node1][node2] == 'an': # check if node1 is an ancestor of node2
        return 1 + distance(node1, parent(node2, species_tree), species_tree, ancestral_table) # we need a function to get the parent of a given node in a tree
    elif ancestral_table[node1][node2] == 'des': #check if node2 is an ancestor of node1
        return 1 + distance(parent(node1, species_tree), node2, species_tree, ancestral_table) # ditto 
    # Last case; checking if the two nodes are incomparable
    elif ancestral_table[node1][node2] == 'in':
        return 2 + distance(parent(node1, species_tree), parent(node2, species_tree), species_tree, ancestral_table)

def parent(node, species_tree):
    for key in species_tree:
        if species_tree[key][0] == node or species_tree[key][1] == node:
            return key
    return node # We will only reach this step if the node doesn't have a parent, which only happens when the node in question is the root of the tree.
                # In this case, we simply return the root itself as its own parent; this does not affect the function in regards to its purpose
                # But even so, we do not expect this to ever happen, because the parent function will only be called on a node that has an ancestor
                # which is a condition the root node will never satisfy.

# A function that creates a dict of the nodeCost function applied to all mapping nodes (g,s). nodeCost() gives us the "cost" of a given node being included in a reconciliation through its difference from other reconciliations.
# nodeCost() is calculated as the summation (for all nodes s' from the species tree) of (distance(s', s) x numPlacing[(g,s')]).
# distance(s',s) is the distance between the nodes s' and s in the species tree.
# The species_tree argument is a tree in vertex dictionary format.
# NOTE: Again, I don't think that it matters if ordered_species_node_list is in pre or post order. We will use post order.
# maybe write a function to calculate distance(s', s) separately?
# in the following function, species_distance_table is a dict that gives us the distance between two species nodes in the species graph. Function not written yet - could probably use recursion for it?

def construct_nodeCost_table(dtl_recon_graph, preorder_mapping_node_list, numPlacing_table, species_tree, ordered_species_node_list, species_distance_table):
        
    # Initialize the nodeCost dict. This dict contains the nodeCost score of each mapping node
    nodeCost = dict()
    for mapping_node in preorder_mapping_node_list:
        nodeCost[mapping_node] = 0
# Q: Can we just use the ordered_species_node_list to find the parent of any given species node? If so, we could rewrite the previous function so we wouldn't require species_tree as an argument (maybe?)
    for mapping_node in preorder_mapping_node_list:
        for species_node in ordered_species_node_list:
            nodeCost[mapping_node] += (distance(mapping_node[1], species_node, species_tree, species_distance_table) * (numPlacing_table[((mapping_node[0]),species_node)]))
    # There might be a problem here; (mapping_node[0], species_node) may not necessarily be a mapping_node in dtl_recon_graph and thus won't have
    # a value in the numPlacing dictionary. However equating numPlacing_table[(mapping_node[0],species_node)] to 0 would still give us the right
    # answer. Perhaps we should just initialize all possible combinations of (gene_node, species_node) to 0 in the beginning (regardless of whether
    # they are mapping nodes or not) and the just change the values that are actually mapping nodes (do this in numPlacing function). 
    return nodeCost

# Function for computing the median reconciliation from the dtl_recon_graph using our new definition of distance between 2 reconciliations.
# Uses dynamic programming algorithm.

def compute_median_2(dtl_recon_graph, postorder_mapping_node_list, nodeCost_table, mpr_roots):

    # Initialize the totalCost dict. This dict contains the totalCost score of each mapping node and each mapping node's events
    totalCost = dict()
    for mapping_node in postorder_mapping_node_list:
        # First checking for the base case of mapping nodes which are leaves
        if dtl_recon_graph[mapping_node] == [('C', (None, None), (None, None))]:
            totalCost[mapping_node] = ([('C',(None,None),(None,None))], 0)
            continue # This tells us if the mapping node in question is a leaf, in which case we don't do anything (since the totalCost of a leaf is 0)

        SDTcostList = list()
        lossCostList = list()
        for event in dtl_recon_graph[mapping_node]:
            # Note that 'event' is of the form: ('event ID', 'Child 1', 'Child 2'), so the 0th element is the event
            # ID and the 1st and 2nd elements are the children produced by the event
            if event[0] == 'L':  # Losses produce only one child, so we only need to look to one lower mapping node
                lossCostList.append((event, totalCost[event[1]][1]))
            else:  # Only other options are T, S, and D, which produce two children
                SDTcostList.append((event, totalCost[event[1]][1] + totalCost[event[2]][1]))
            
        # Now we want to calculate totalCost for the map node itself. To do this, first we need to find min{totalCost(e) + nodeCost(g,s)} for all SDT events e of the node and min{totalCost(l)} for all loss events l of the node    
        if lossCostList == []:
            minLossCost = float('inf')
        else:
            minLossCost = min(lossCostList, key=itemgetter(1))[1]

        if SDTcostList == []:
            minSDTcost = float('inf')
        else:
            minSDTcost = min(SDTcostList, key=itemgetter(1))[1]
            minSDTcost += nodeCost_table[mapping_node] #because if we go through a SDT event (meaning g mapped onto s) then we also need to add the nodeCost
        
        min_cost = min(minLossCost, minSDTcost)

        best_events = list()

        # Check to see which event(s) produce the minimum total cost
        for lossEvent in lossCostList:
            if lossEvent[1] == min_cost:
                best_events.append(lossEvent[0])
        for SDTevent in SDTcostList:
            if SDTevent[1] + nodeCost_table[mapping_node] == min_cost:
                best_events.append(SDTevent[0])
        
        del lossCostList
        del SDTcostList

        # Save the result for this mapping node so it can be used in higher mapping nodes in the graph
        totalCost[mapping_node] = (best_events[:], min_cost) # RR - best_events[:] slices the list and is a short cut to making a copy of the list that we want, which is why we think they did it

    # Get all possible roots of the graph and their running frequency scores, in a list, for later use
    possible_root_combos = [(root, totalCost[root][1]) for root in mpr_roots]

    # Find the best minimum cost sum for all of the potential roots for the median
    best_sum = min(possible_root_combos, key=itemgetter(1))[1]

    # Find all of the root combos for a median by filtering out the roots that don't give the best minimum total cost
    best_root_combos = list(filter(lambda x: x[1] == best_sum, possible_root_combos))

    # Extract just the roots from the previously filtered out list
    best_roots = [root[0] for root in best_root_combos]

    # Adjust the sum_freqs dictionary so we can use it with the buildDTLReconGraph function from DTLReconGraph.py
    for mapping_node in totalCost:

        # We place the event tuples into lists so they work well with the diameter algorithm
        totalCost[mapping_node] = totalCost[mapping_node][0]  # Only use the events, no longer the associated frequency sum

    # Use the buildDTLReconGraph function from DTLReconGraph.py to find the median recon graph
    # Note that build_dtl... requires a list of the best roots for a reconciliation graph, the events for each
    # mapping node that are viable for an MPR (in our case, the median), and an empty dictionary to populate
    # as the final return value
    med_recon_graph = DTLReconGraph.build_dtl_recon_graph(best_roots, totalCost, {})

    # Check to make sure the median is a subgraph of the DTL reconciliation
    assert check_subgraph(dtl_recon_graph, med_recon_graph), 'Median is not a subgraph of the recon graph!'

    # We can use this function to find the number of medians once we've got the final median recon graph
    n_med_recons = DTLReconGraph.count_mprs_wrapper(best_roots, med_recon_graph)

    return med_recon_graph, n_med_recons, best_roots

# Function that returns the number of medians given a species (or host) and gene (or parasite) tree, for both different
# definitions of distance between two different MPRs
# Arguments - file that has the gene tree, species tree, and phi function in correct format, and the costs for duplications, transfers and losses
# Output - A tuple giving the number of medians using the first method, and the number of medians using the second method

def compare_median_numbers(file_name, dup_cost, transfer_cost, loss_cost):
    # First we get reconcile the file to get the host and parasite trees in edge format and the DTL reconciliation graph
    host,paras,graph,num_recon,best_roots = DTLReconGraph.reconcile(file_name,dup_cost,transfer_cost,loss_cost)

    # We want the trees in vertex format to use as arguments in future functions
    vhost = Diameter.reformat_tree(host,'hTop')
    vparas = Diameter.reformat_tree(paras,'pTop') 

    # We also want lists of the vertices of these trees  both in preorder and postorder

    preorder_vhost_nodes = preorder_vertices(vhost[0],vhost[1])
    preorder_vparas_nodes = preorder_vertices(vparas[0],vparas[1])
    postorder_vhost_nodes = postorder_vertices(vhost[0],vhost[1])
    postorder_vparas_nodes = postorder_vertices(vparas[0],vparas[1])

    # Now we get a list of the mapping nodes in the DTL graph, also in both preorder and postorder 
    mapping_node_list = []
    for key in graph:
        mapping_node_list.append(key)
    preorder_mapping_node_list = mapping_node_sort(preorder_vparas_nodes,preorder_vhost_nodes,mapping_node_list)
    postorder_mapping_node_list = mapping_node_sort(postorder_vparas_nodes,postorder_vhost_nodes,mapping_node_list)


    # Calculating the event scores (normalized and not) to use as arguments in median functions
    event_scores,count = generate_scores_no_normalizing(preorder_mapping_node_list,graph,vparas[1])
    event_scores1,count1 = generate_scores(preorder_mapping_node_list,graph,vparas[1])


    # Getting the numPlacing table
    numPlacingTable = construct_numPlacing_table(graph,postorder_vhost_nodes,postorder_vparas_nodes,postorder_mapping_node_list,event_scores,count)

    #Getting the nodeCost table
    species_distance_table = Diameter.calculate_ancestral_table(vhost[0])
    nodeCostTable = construct_nodeCost_table(graph,preorder_mapping_node_list,numPlacingTable,vhost[0],postorder_vhost_nodes,species_distance_table)
    print('Done with getting nodeCostTable')
    #Using the first median function and storing useful values
    med_recon_graph1,n_med_recons1,best_roots1 = compute_median(graph,event_scores1,postorder_mapping_node_list,best_roots)

    #Using the second median function and storing useful values
    med_recon_graph2,n_med_recons2,best_roots2 = compute_median_2(graph,postorder_mapping_node_list,nodeCostTable,best_roots)

    #returning the tuple of number of medians for both types of distance metrics
    print('(' + str(n_med_recons1) + ' , ' + str(n_med_recons2) + ')')
    return (count,n_med_recons1,n_med_recons2)

# Function to find random medians in the MPR space and return a list of the frequencies of events in them

def generate_event_freq_random_median(file_name,dup_cost,transfer_cost,loss_cost):
    #Copying code from previous functions so we can calculate similarly useful variable values
    # First we get reconcile the file to get the host and parasite trees in edge format and the DTL reconciliation graph
    host,paras,graph,num_recon,best_roots = DTLReconGraph.reconcile(file_name,dup_cost,transfer_cost,loss_cost)

    # We want the trees in vertex format to use as arguments in future functions
    vhost = Diameter.reformat_tree(host,'hTop')
    vparas = Diameter.reformat_tree(paras,'pTop') 

    # We also want lists of the vertices of these trees  both in preorder and postorder

    preorder_vhost_nodes = preorder_vertices(vhost[0],vhost[1])
    preorder_vparas_nodes = preorder_vertices(vparas[0],vparas[1])
    postorder_vhost_nodes = postorder_vertices(vhost[0],vhost[1])
    postorder_vparas_nodes = postorder_vertices(vparas[0],vparas[1])

    # Now we get a list of the mapping nodes in the DTL graph, also in both preorder and postorder 
    mapping_node_list = []
    for key in graph:
        mapping_node_list.append(key)
    preorder_mapping_node_list = mapping_node_sort(preorder_vparas_nodes,preorder_vhost_nodes,mapping_node_list)
    postorder_mapping_node_list = mapping_node_sort(postorder_vparas_nodes,postorder_vhost_nodes,mapping_node_list)


    # Calculating the event scores (normalized and not) to use as arguments in median functions
    event_scores,count = generate_scores_no_normalizing(preorder_mapping_node_list,graph,vparas[1])
    event_scores1,count1 = generate_scores(preorder_mapping_node_list,graph,vparas[1])


    # Getting the numPlacing table
    numPlacingTable = construct_numPlacing_table(graph,postorder_vhost_nodes,postorder_vparas_nodes,postorder_mapping_node_list,event_scores,count)

    #Getting the nodeCost table
    species_distance_table = Diameter.calculate_ancestral_table(vhost[0])
    nodeCostTable = construct_nodeCost_table(graph,preorder_mapping_node_list,numPlacingTable,vhost[0],postorder_vhost_nodes,species_distance_table)
    print('Done with getting nodeCostTable')
    #Using the first median function and storing useful values
    med_recon_graph1,n_med_recons1,best_roots1 = compute_median(graph,event_scores1,postorder_mapping_node_list,best_roots)

    #Using the second median function and storing useful values
    med_recon_graph2,n_med_recons2,best_roots2 = compute_median_2(graph,postorder_mapping_node_list,nodeCostTable,best_roots)

    #Remember that eventScores holds the frequency of events in the DTL graph unnormalized, whereas 
    #eventScores1 has the normalized frequency of events (no. of reconciliations in which the event occurs/total reconciliations in graph)
    # Can use count_mprs function to get no. of MPRs spawned below a given event node.

    med_counts1 = dict()
    for root in best_roots1:
        count_mprs(root, med_recon_graph1, med_counts1)

    med_counts2 = dict()
    for root in best_roots2:
        count_mprs(root, med_recon_graph2, med_counts2)

    random_median1 = choose_random_median_wrapper(med_recon_graph1, best_roots1, med_counts1)
    random_median2 = choose_random_median_wrapper(med_recon_graph2, best_roots2, med_counts2)

    freq_scores1 = list()
    freq_scores2 = list()

    for mapping_node in random_median1:
        for event in random_median1[mapping_node]:
            if event == ('C',(None,None),(None,None)):
                continue
            freq_scores1.append(event_scores1[event]) #event_scores1 is a dictionary storing normalized frequencies for all the events in the reconciliation graph, which is the same for both median types
    
    for mapping_node in random_median2:
        for event in random_median2[mapping_node]:
            if event == ('C',(None,None),(None,None)):
                continue
            freq_scores2.append(event_scores1[event]) # ditto
    
    return freq_scores1,freq_scores2

def test_event_freq(filename, dup_cost, transfer_cost, loss_cost):
    host, paras, graph, num_recon, best_roots = DTLReconGraph.reconcile(filename, dup_cost, transfer_cost, loss_cost)

    vhost = Diameter.reformat_tree(host, 'hTop')
    vparas = Diameter.reformat_tree(paras, 'pTop')

    preorder_vhost_nodes = preorder_vertices(vhost[0],vhost[1])
    preorder_vparas_nodes = preorder_vertices(vparas[0],vparas[1])
    postorder_vhost_nodes = postorder_vertices(vhost[0],vhost[1])
    postorder_vparas_nodes = postorder_vertices(vparas[0],vparas[1])

    mapping_node_list = []
    for key in graph:  #graph = DTL reconciliation graph
        mapping_node_list.append(key)

    preorder_mapping_node_list = mapping_node_sort(preorder_vparas_nodes,preorder_vhost_nodes,mapping_node_list)
    postorder_mapping_node_list = mapping_node_sort(postorder_vparas_nodes,postorder_vhost_nodes,mapping_node_list)

    event_scores, count = generate_scores_no_normalizing(preorder_mapping_node_list,graph,vparas[1])
    event_scores1, count1 = generate_scores(preorder_mapping_node_list,graph,vparas[1])

    #numPlacingTable = construct_numPlacing_table(graph,postorder_vhost_nodes,postorder_vparas_nodes,postorder_mapping_node_list,event_scores,count)

    #species_distance_table = Diameter.calculate_ancestral_table(vhost[0])

    #nodeCostTable = construct_nodeCost_table(graph, preorder_mapping_node_list, numPlacingTable, vhost[0], postorder_vhost_nodes, species_distance_table)

    med_recon_graph1,n_med_recons1,best_roots1 = compute_median(graph, event_scores1, postorder_mapping_node_list, best_roots)

    return event_scores, event_scores1, best_roots, graph



def usage():
    """
    :return: the usage statement associated with running this file
    """

    return 'usage: DTLMedian filename dup_cost transfer_cost loss_cost [-r] [-n]'


def main():
    """
    :return: nothing. This function will run the main loop for the command line interface.
    """

    p = optparse.OptionParser(usage=usage())

    p.add_option('-r', '--random', dest='random', help='Add a random median reconciliation from the full median'
                                                       ' reconciliation graph of the given file to the output',
                 action='store_true', default=False)
    p.add_option('-c', '--count', dest='count', help='Add the number of median reconciliations to'
                                                     'the output', action='store_true', default=False)

    options, args = p.parse_args()

    if len(args) == 4:
        try:

            # These will be the outputs we eventually return
            output = []

            # Save arg values
            filename = args[0]
            dup = float(args[1])
            transfer = float(args[2])
            loss = float(args[3])

            # Get basic info just about the dtl recon graph
            species_tree, gene_tree, dtl_recon_graph, mpr_count, best_roots = DTLReconGraph.reconcile(filename, dup,
                                                                                                      transfer, loss)

            # Reformat gene tree and get info on it, as well as for the species tree in the following line
            postorder_gene_tree, gene_tree_root, gene_node_count = Diameter.reformat_tree(gene_tree, "pTop")
            postorder_species_tree, species_tree_root, species_node_count = Diameter.reformat_tree(species_tree,
                                                                                                   "hTop")

            # Get a list of the mapping nodes in preorder
            postorder_mapping_node_list = mapping_node_sort(postorder_gene_tree, postorder_species_tree,
                                                            dtl_recon_graph.keys())

            # Find the dictionary for frequency scores for the given mapping nodes and graph, and the given gene root
            scores_dict = generate_scores(postorder_mapping_node_list[::-1], dtl_recon_graph, gene_tree_root)

            # Now find the median and related info
            median_reconciliation, n_meds, roots_for_median = compute_median(dtl_recon_graph, scores_dict[0],
                                                                             postorder_mapping_node_list, best_roots)

            # We'll always want to output the median
            output.append(median_reconciliation)

            # Check if the user wants the number of medians
            if options.count:
                output.append(n_meds)

            # Check if the user wants a random median
            if options.random:

                # Initialize the dictionary that tells us how many medians can be spawned from a particular event node
                med_counts = dict()

                # Now fill it
                for root in roots_for_median:
                    count_mprs(root, median_reconciliation, med_counts)

                # In case we may want it, calculate a random, uniformly sampled single-path median from the median recon
                random_median = choose_random_median_wrapper(median_reconciliation, roots_for_median, med_counts)
                output.append(random_median)

            # Now print all of the output requested by the user
            for i in range(len(output)):
                if i != (len(output) - 1):
                    print(str(output[i]) + '\n')
                else:
                    print(str(output[i]))

        except ValueError:
            print(usage())
    else:
        print usage()


if __name__ == '__main__':

    main()

def interpretData():
    """
    :param filename: csv file containing the event frequencies from a TreeLife.newick file
    
    :return: means of each of averages, minimums, and maximums
    """
    #First we find all the csv files that we are to interpret
    directories = list(os.walk("."))
    textfilesWithPaths = list()

    #find all text files with their paths prepended
    textfilesWithPaths = list()
    for item in directories:
        if item[0][0:2] == "./":
            path = item
            files = [item[0] + "/" + file for file in item[2]]
            textfilesWithPaths.extend(files)
    
    #only take files ending in "***_*.csv"
    symmetric111TestFiles = list(filter(lambda x: x[-9:] == '111_s.csv', textfilesWithPaths))
    path111TestFiles = list(filter(lambda x: x[-9:] == '111_p.csv', textfilesWithPaths))
    symmetric231TestFiles = list(filter(lambda x: x[-9:] == '231_s.csv', textfilesWithPaths))
    path231TestFiles = list(filter(lambda x: x[-9:] == '231_p.csv', textfilesWithPaths))
    

    #Write csv files for each subdivision of csv inputs
    #Symmetric 111
    with open('Symmetric_111.csv', 'a') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['File Name','Average Mean','Average Minimum','Average Maximum'])

        for file in symmetric111TestFiles:
            #First we read the csv file and initialize rows
            rows = []
            with open(file, 'r') as csvfile:
                csvreader = csv.reader(csvfile)
                #Now we extract each data row one by one
                for row in csvreader:
                    rows.append(row)
                numRows = (csvreader.line_num)      #Get total number of rows

            #Finding the average, minimum, and maximum for each row
            rowAverages = []
            rowSums = [0]*numRows
            numColumns = [0]*numRows
            rowMinimums = [None]*numRows
            rowMaximums = [None]*numRows
            for i in range(len(rows)):
                #parsing each column of a row to sum column values and count how many columns per row
                for j in range(len(rows[i])):
                    rowSums[i] += float(rows[i][j])
                    numColumns[i] += 1
                rowAverages.append(rowSums[i]/numColumns[i])    #Calculate average

            rows = [[float(col) for col in row] for row in rows]        #turn each string value into a float
            
            for i in range(len(rows)):
                rowMinimums[i] = min(rows[i])
                rowMaximums[i] = max(rows[i])

            #Next we find the mean of the averages, the minimums, and the maximums
            meanOfAverages = sum(rowAverages)/numRows
            meanOfMinimums = sum(rowMinimums)/numRows
            meanOfMaximums = sum(rowMaximums)/numRows
            
            #write to the output csv file
            filewriter.writerow([file, meanOfAverages, meanOfMinimums, meanOfMaximums])
 
    #Path 111
    with open('Path_111.csv', 'a') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['File Name','Average Mean','Average Minimum','Average Maximum'])

        for file in path111TestFiles:
            #First we read the csv file and initialize rows
            rows = []
            with open(file, 'r') as csvfile:
                csvreader = csv.reader(csvfile)
                #Now we extract each data row one by one
                for row in csvreader:
                    rows.append(row)
                numRows = (csvreader.line_num)      #Get total number of rows

            #Finding the average, minimum, and maximum for each row
            rowAverages = []
            rowSums = [0]*numRows
            numColumns = [0]*numRows
            rowMinimums = [None]*numRows
            rowMaximums = [None]*numRows
            for i in range(len(rows)):
                #parsing each column of a row to sum column values and count how many columns per row
                for j in range(len(rows[i])):
                    rowSums[i] += float(rows[i][j])
                    numColumns[i] += 1
                rowAverages.append(rowSums[i]/numColumns[i])    #Calculate average

            rows = [[float(col) for col in row] for row in rows]        #turn each string value into a float
            
            for i in range(len(rows)):
                rowMinimums[i] = min(rows[i])
                rowMaximums[i] = max(rows[i])

            #Next we find the mean of the averages, the minimums, and the maximums
            meanOfAverages = sum(rowAverages)/numRows
            meanOfMinimums = sum(rowMinimums)/numRows
            meanOfMaximums = sum(rowMaximums)/numRows
            
            #write to the output csv file
            filewriter.writerow([file, meanOfAverages, meanOfMinimums, meanOfMaximums])

    #Symmetric 231
    with open('Symmetric_231.csv', 'a') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['File Name','Average Mean','Average Minimum','Average Maximum'])

        for file in symmetric231TestFiles:
            #First we read the csv file and initialize rows
            rows = []
            with open(file, 'r') as csvfile:
                csvreader = csv.reader(csvfile)
                #Now we extract each data row one by one
                for row in csvreader:
                    rows.append(row)
                numRows = (csvreader.line_num)      #Get total number of rows

            #Finding the average, minimum, and maximum for each row
            rowAverages = []
            rowSums = [0]*numRows
            numColumns = [0]*numRows
            rowMinimums = [None]*numRows
            rowMaximums = [None]*numRows
            for i in range(len(rows)):
                #parsing each column of a row to sum column values and count how many columns per row
                for j in range(len(rows[i])):
                    rowSums[i] += float(rows[i][j])
                    numColumns[i] += 1
                rowAverages.append(rowSums[i]/numColumns[i])    #Calculate average

            rows = [[float(col) for col in row] for row in rows]        #turn each string value into a float
            
            for i in range(len(rows)):
                rowMinimums[i] = min(rows[i])
                rowMaximums[i] = max(rows[i])

            #Next we find the mean of the averages, the minimums, and the maximums
            meanOfAverages = sum(rowAverages)/numRows
            meanOfMinimums = sum(rowMinimums)/numRows
            meanOfMaximums = sum(rowMaximums)/numRows
            
            #write to the output csv file
            filewriter.writerow([file, meanOfAverages, meanOfMinimums, meanOfMaximums])
 
    #Path 231
    with open('Path_231.csv', 'a') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['File Name','Average Mean','Average Minimum','Average Maximum'])

        for file in path231TestFiles:
            #First we read the csv file and initialize rows
            rows = []
            with open(file, 'r') as csvfile:
                csvreader = csv.reader(csvfile)
                #Now we extract each data row one by one
                for row in csvreader:
                    rows.append(row)
                numRows = (csvreader.line_num)      #Get total number of rows

            #Finding the average, minimum, and maximum for each row
            rowAverages = []
            rowSums = [0]*numRows
            numColumns = [0]*numRows
            rowMinimums = [None]*numRows
            rowMaximums = [None]*numRows
            for i in range(len(rows)):
                #parsing each column of a row to sum column values and count how many columns per row
                for j in range(len(rows[i])):
                    rowSums[i] += float(rows[i][j])
                    numColumns[i] += 1
                rowAverages.append(rowSums[i]/numColumns[i])    #Calculate average

            rows = [[float(col) for col in row] for row in rows]        #turn each string value into a float
            
            for i in range(len(rows)):
                rowMinimums[i] = min(rows[i])
                rowMaximums[i] = max(rows[i])

            #Next we find the mean of the averages, the minimums, and the maximums
            meanOfAverages = sum(rowAverages)/numRows
            meanOfMinimums = sum(rowMinimums)/numRows
            meanOfMaximums = sum(rowMaximums)/numRows
            
            #write to the output csv file
            filewriter.writerow([file, meanOfAverages, meanOfMinimums, meanOfMaximums])
 