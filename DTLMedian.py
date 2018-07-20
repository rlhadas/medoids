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
#   event frequency list (returned as elements of lists from functions like generate_event_freq_random_median_symmetric_set_method() and  generate_event_freq_random_reconciliation()) :
#       ['event1_type',event1_frequency,'event2_type',event2_frequency, ... 'eventn_type', eventn_frequency']

import os
import importlib
import optparse
from operator import itemgetter
import numpy as np
import DTLReconGraph
import Diameter
import math
import csv
import sys

def readcsv(filename):
    filename = 'testEventFreq.csv'
    
    fields = []
    rows = []
    
    with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        fields = csvreader.next()

        for row in csvreader:
            rows.append(row)
        
        print('Total no. of rows: %d'%(csvreader.line_num))
        
    print('Field names are:' + ', '.join(field for field in fields))
    print('\nFirst 5 rows are:\n')
    for row in rows[:5]:
        for col in row:
            print('%10s'%col),
        print('\n')
    return fields,rows


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


def generate_scores(preorder_mapping_node_list, dtl_recon_graph, gene_root, isNormalized):
    """
    Computes frequencies for every event
    :param preorder_mapping_node_list: A list of all mapping nodes in DTLReconGraph in double preorder
    :param dtl_recon_graph: The DTL reconciliation graph that we are scoring
    :param gene_root: The root of the gene tree
    :param isNormalized: a Boolean value, which takes the value True if user would like a dict filled with event frequencies (number of MPRs event is in/total number 
    of MPRs and takes value of False if user would like a dict filled with the actual number of MPRs the event is in (non-normalized values)
    :return: 0. A file structured like the DTLReconGraph, but with the lists of events replaced
                with dicts, where the keys are the events and the values are the scores of those events ('score' is defined based off isNormalized, so either the score 
                is event frequency, or it is the actual number of MPRs the event is in), and
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

    # Normalize all of the event_scores, but only if isNormalized == True
    if isNormalized == True:
        for mapping_node in preorder_mapping_node_list:
            for event in dtl_recon_graph[mapping_node]:
                if event == ('C',(None,None),(None,None)):
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
        return counts[mapping_node] 

    # Base case, occurs if being called on a child produced by a loss or contemporary event 
    if mapping_node == (None, None): # we check for losses or contemporary events first so that when we loop through events later we know those will only be duplications, speciations or transfers
        return 1 

    # Initialize a variable to keep count of the number of MPRs
    count = 0

    # Loop over all event nodes corresponding to the current mapping node
    for eventNode in dtl_recon_graph[mapping_node]: 

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
    Finds the medoid reconciliation of a given DTLReconGraph using the symmetric set difference distance metric.
    :param dtl_recon_graph: A dictionary representing a DTL Recon Graph.
    :param event_scores: A dictionary with event nodes as keys and values corresponding to the frequency of
    that events in MPR space for the recon graph
    :param postorder_mapping_nodes: A list of the mapping nodes in a possible MPR, except sorted first in
    postorder by species node and postorder by gene node
    :param mpr_roots: A list of mapping nodes that could act as roots to an MPR for the species and
    gene trees in question, output from the findBestRoots function in DTLReconGraph.py
    :return: 0. A new dictionary which has the same form as a DTL reconciliation graph except this new graph is a subgraph of the original DTLReconGraph, and
                every possible reconciliation in this new reconciliation graph is a medoid as defined by the symmetric set difference distance metric. Thus 
                this dictionary represents the medoid reconciliation using the symmetric set difference distance metric.
             1. the number of medoid reconciliations for the given DTL reconciliation graph, using the symmetric set difference distance metric.
             2. the root(s) of the obtained medoid reconciliation graph (returned in 0.).
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
        sum_freqs[map_node] = (best_events[:], max_sum) # best_events[:] slices the list and is a short cut to making a copy of the list that we want.

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
    # mapping node that are viable for an MPR (in our case, the median), and an empty dictionary to populate
    # as the final return value
    med_recon_graph = DTLReconGraph.build_dtl_recon_graph(best_roots, sum_freqs, {})

    # Check to make sure the median is a subgraph of the DTL reconciliation
    assert check_subgraph(dtl_recon_graph, med_recon_graph), 'Median is not a subgraph of the recon graph!'

    # We can use this function to find the number of medians once we've got the final median recon graph
    n_med_recons = DTLReconGraph.count_mprs_wrapper(best_roots, med_recon_graph)

    return med_recon_graph, n_med_recons, best_roots

def compute_maximum_average(dtl_recon_graph, event_scores, postorder_mapping_nodes, mpr_roots):
    """
    :param dtl_recon_graph: A dictionary representing a DTL Recon Graph.
    :param event_scores: A dictionary with event nodes as keys and values corresponding to the frequency of
    that events in MPR space for the recon graph.
    :param postorder_mapping_nodes: A list of the mapping nodes in a possible MPR, except sorted first in
    postorder by species node and postorder by gene node
    :param mpr_roots: the root(s) of the given DTL reconciliation graph (param dtl_recon_graph).
    :return: 0. a float value, which is the largest maximum average frequency of all events in any of the MPRs that are in the DTLReconGraph.
             1. a list of integers, each of which is the number of events associated with any MPR that has the maximum average event 
                frequency (note that it's possible that there are different MPRs, with different numbers of events, that all have the same 
                maximum average frequency)   
    NOTE: In further implementation we also want to return the reconciliation that corresponds to this maximum average frequency of its events
    """

    # Initialize a dict that will store the list of running maximum total frequency sums incurred up to the given mapping node,
    # and the event node(s) that directly gave it that frequency sum. Keys are mapping nodes, values are two-element tuples.
    # The first element is a list, for which the ith element is itself a list of event nodes (directly below the corresponding mapping node key)
    # which maximize the sum of event frequencies for an MPR, with exactly (i+1) events, rooted at that mapping node.
    # The second element (which will henceforth be refered to as the frequency vector) is a list of floating point values, in which the ith element is the 
    # maximum sum of all event frequencies for an MPR, with exactly (i+1) events, rooted at that mapping node.
    sum_freqs = dict()
    # This will be the length of the frequency vectors that is stored in sum_freqs
    max_list_count = 0
    numLeaves = 0

    # Here we count the total number of events in the entire reconciliation graph to create an upper-bound for the length of the vectors we create later.
    # Counting all event nodes that are not contemporaneous
    for event in event_scores:
        if event != ('C', (None, None), (None, None)):
            max_list_count += 1
    # Counting all event nodes that are contemporaneous, which is also the number of leaf nodes
    for mapping_node in postorder_mapping_nodes:
        for event in dtl_recon_graph[mapping_node]:
            if event == ('C', (None, None), (None, None)):
                max_list_count += 1
                numLeaves += 1

    # Loop over all mapping nodes for the gene tree in postorder - dynamic programming algorithm
    for map_node in postorder_mapping_nodes:
        # Contemporaneous events need to be caught from the get-go
        if dtl_recon_graph[map_node] == [('C', (None, None), (None, None))]: # Indicates the relevant map_node is a leaf node - our base case
            currList = [0.0]* max_list_count
            currMaxEventList = [None] * max_list_count
            currList[0] = 1.0 # C events have freq 1, and only one event is involved, so maximum average frequency for one event will be 1.
            currMaxEventList[0] = [('C',(None,None),(None,None))] # The only event involved in the maximum average frequency for one event is the contemporaneous event.
            sum_freqs[map_node] = (currMaxEventList[:], currList) # Updating the sum_freqs dict for the leaf node
            continue  # Contemporaneous events should be a lone event in a list, so we move to the next mapping node

        # If the map_node is not a leaf node, then we look at all of its events to calculate their sum_freqs, which will in turn be used to calculate
        # the sum_freqs of the map_node itself.
        events = list()
        for event in dtl_recon_graph[map_node]:
            # Note that 'event' is of the form: ('event ID', 'Child 1', 'Child 2'), so the 0th element is the event
            # ID and the 1st and 2nd elements are the children produced by the event
            if event[0] == 'L':  # Losses produce only one child, so we only need to look to one lower mapping node
                # Since a loss event only has one child, the way to get (i+1) events is to have i events in the loss event's child, and the sum of frequencies is then
                # just the sum of frequencies of its child for i events, plus the loss event's own frequency.
                # Essentially, the frequency vector gets shifted over by 1 index, and the loss event's frequency value is added to all elements of the vector
                currEventList = [0.0] * max_list_count
                for i in range(1, max_list_count):
                    currEventList[i] = sum_freqs[event[1]][1][i-1] + event_scores[event]
                events.append(currEventList)
            else:  # Only other options are T, S, and D, which produce two children
                currEventList = [0.0] * max_list_count
                # Since T, S and D events have two children, the way to get (i+1) events in any MPR below is to have x events from one child and (i+1-x) events from the other child
                # where x must be at least 1 (since it's not possible for an MPR to have 0 events)
                # Also note that (i+1) must at least be 3, since the event itself counts as 1 of the events of the MPR, and the number of events below it must atleast be 2.
                # So we go through all possible values of x and record the maximum sum of event frequencies for each value of x, which is just the sum of its children's sum_freqs and its own frequency
                for i in range(2, max_list_count):
                    tempEventScores = []
                    for j in range(0,i-1):
                        if sum_freqs[event[1]][1][j] == 0.0 or sum_freqs[event[2]][1][i-j-2] == 0.0: # If this is true, then it means there are no MPRs for atleast one of the event's children
                            tempEventScores.append(0.0)
                            continue
                        tempEventScores.append(sum_freqs[event[1]][1][j] + sum_freqs[event[2]][1][i-j-2])
                    currEventList[i] = max(tempEventScores) + event_scores[event] # Adding the event's own frequency to the sum of frequencies score
                events.append(currEventList)
        
        # We have the sum_freqs for all of our mapping node's events, so now we calculate the sum_freqs for the mapping node using those values
        
        # Initialising the frequency vector for the current mapping node
        currList = [0.0] * max_list_count 
        for i in range(max_list_count):
            currList[i] = max(event[i] for event in events) # Find the maximum sum event frequency for each number of events from 0 to max_list_count and add it to the vector

        # Initialising the list of list of events corresponding to the frequency vector
        currMaxEventList = [None] * max_list_count 
        # Going through all possible number of events and finding all those that give us the maximum sum of event frequency
        for i in range(max_list_count):
            if currList[i] == 0.0 : # if this is True, then there are no events which lead to MPRs with (i+1) events, so we leave the value as None
                continue
            tempList = [] # Otherwise, we make a list, appending all events that give us maximum sum of event frequency for MPRs with (i+1) events
            for event in events:
                if event[i] == currList[i]:
                    tempList.append(event)
            currMaxEventList[i] = tempList

        del events

        sum_freqs[map_node] = (currMaxEventList[:], currList) # Then update the sum_freqs dict for the current mapping node

    # Now that the sum_freqs dict is completely filled, we look at all the vectors of all the roots of the DTLReconGraph and calculate the average event frequency
    # for each of them, to find the maximum event frequency out of all of them.
    averages = [] # a list to hold all maximum average value of each root
    currMax = 0   # a changing variable that will hold the maximum event frequency average we have encountered so far
    numEvents = [] # a list that will hold the number of events in any MPR that has the maximum event frequency average

    for root in mpr_roots:
        averageList = []
        # We only want to look at averages for MPRs that have at least one non-contemporaneous event (so number of events is greater than the number of leaves)
        for i in range(numLeaves, max_list_count):
            # the vector holds in its ith position, the maximum sum of event frequencies for an MPR with (i + 1) number of events. Therefore, the average event 
            # frequency will be given by the  value in the ith position of the vector divided by the number of events, which is just (i+1)
            # HOWEVER, we are ignoring contemporaneous events (which all have frequency 1), so we subtract the number of contemporaneous events (i.e numLeaves)
            # from both the numerator and denominator.
            currAverage = (sum_freqs[root][1][i] - numLeaves)/ (float(i+1 - numLeaves)) 
            averageList.append(currAverage)
        maxAverage = max(averageList) # Finding the maximum among all the averages for the current root
        averages.append(maxAverage) # Appending the maximum event frequency average found for the current root to averages

    
    maxMaxAverage = max(averages) # Finding the maximum of the maximum event frequency average of all the roots - which is our final answer

    # There may be multiple MPRs which have maxMaxAverage as their average event frequency - this block of code finds all such MPRs and records the number 
    # of events in each of them, putting these values in numEvents list.
    for root in mpr_roots:
        for i in range(numLeaves, max_list_count):
            currAverage = (sum_freqs[root][1][i] - numLeaves)/ (float(i+1 - numLeaves))
            if currAverage == maxMaxAverage:
                numEvents.append(i+1 - numLeaves) # Number of events that are not contemporaneous         

    return maxMaxAverage, numEvents 
            

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


def preorder_vertices(tree, root):
    """
    :param tree: host or parasite tree (in vertex format) - recall that vertex format of a tree is
    an ordered dictionary where keys are vertices which relate to tuples, where each tuple consists of the vertice's two children
    :param root: the root of the tree which is above all relevant nodes of the tree - for ease of use, the root of a host tree is always 'hTop', and  
    the root of a parasite tree is always 'pTop'
    :return: a list of all the vertices of the tree in preorder ('hTop' and 'pTop' are not included in this list)
    """
    left_child = tree[root][0]
    right_child = tree[root][1]
        # Base case
    if left_child is None:  # Then right_child == None also
        return [root]
    else:  # Recursive call
        return [root] + \
            preorder_vertices(tree,left_child) + \
            preorder_vertices(tree,right_child)
    

def postorder_vertices(tree, root):
    """
    :param tree: host or parasite tree (in vertex format) - recall that vertex format of a tree is
    an ordered dictionary where keys are vertices which relate to tuples, where each tuple consists of the vertice's two children.
    :param root: the root of the tree which is above all relevant nodes of the tree - for ease of use, the root of a host tree is always 'hTop', and 
    the root of a parasite tree is always 'pTop'.
    :return: a list of all the vertices of the tree in postorder ('hTop' and 'pTop' are not included in this list).
    """
    left_child = tree[root][0]
    right_child = tree[root][1]
        # Base case
    if left_child is None:  # Then right_child == None also
        return [root]
    else:  # Recursive call
        return postorder_vertices(tree,left_child) + \
            postorder_vertices(tree,right_child) + \
            [root]

                                                                                                            

def construct_numPlacing_table(dtl_recon_graph, ordered_species_node_list, ordered_gene_node_list, postorder_mapping_node_list, event_scores, count):
    """
    :param dtl_recon_graph: A dictionary representing a DTL Recon Graph.
    :ordered_species_node_list: an ordered list of the nodes (postorder) in the species (or host) tree used to make the corresponding DTL Recon Graph. 
    :ordered_gene_node_list: an ordered list of the nodes (postorder) in the gene (or parasite) tree used to make the corresponding DTL Recon Graph. 
    :postorder_mapping_node_list: A list of the possible mapping nodes in any MPR, except sorted first in
    postorder by gene node and postorder by species node.
    :event_scores: A dictionary with event nodes as keys and values corresponding to the total number of MPRs that event occured in.
    :count: the total number of MPRs in DTLReconGraph (param dtl_recon_graph)
    :return: a dictionary in which the keys are all of the form (g,s), where g is a node from the gene tree, and s is a node from the gene tree, and the corresponding
             values are the number of MPRs in DTLReconGraph in which g is placed on s. Note that g is placed on s in a given MPR iff (g,s) is a mapping node in the 
             MPR and the event following it is NOT a loss. If (g,s) is not a mapping node, then its value in this dictionary is trivially 0.
    """

    # Initialize the numPlacing dict. This dict contains the numPlacing score of each mapping node 
    numPlacing = dict()

    # We want values in this dictionary for every possible combination of (gene node, species node), not just the mapping nodes. But the values
    # for the combinations which are not valid mapping nodes will be trivially 0. 
    
    for gene_node in ordered_gene_node_list:
        for species_node in ordered_species_node_list:
            numPlacing[(gene_node, species_node)] = 0

    # Now, we will just change the numPlacing dict values for the keys that are mapping nodes, and may have non-zero values
    for map_node in postorder_mapping_node_list:
    # First checking for leaf nodes, which will have only contemporaneous events following them.
        if dtl_recon_graph[map_node] == [('C', (None, None), (None, None))]:    
            numPlacing[map_node] = count # Because leaf nodes of the form (g,s) are present in all reconciliations and g is always placed on s.
            continue
    # Otherwise, the mapping node is NOT a leaf node, and we want to add the frequencies of its events that are not losses.
        for event in dtl_recon_graph[map_node]:
            if event[0] != 'L': # Making sure the relevant event is not a loss
                numPlacing[map_node] += event_scores[event] # Adding the frequencies of events of (g,s) that are NOT loss events.

    return numPlacing
 

def distance(node1, node2, species_tree, ancestral_table):
    """
    :param node1: a node in the species tree (param species_tree)
    :param node2: another (possibly same) node in the species tree (param species_tree)
    :param species_tree: a species (or host) tree, in vertex format (see documentation at top of file)
    :param ancestral_table: A nested dictionary. The first dictionary has vertices in the species tree as keys and the values are dictionaries. These dictionaries have
                            as keys vertices of the species tree (again) and values which are strings, representing how the first index relates to the second. 
                            The string 'eq' indicates that the first index and the second are the same (equal), 'an' indicates that the first index is an ancestor 
                            of the second. 'des' indicates the first index is a descendant of the second, and 'in' indicates that the first and second index are
                            incomparable (neither ancestor nor descendant).
    :return: an integer value, which is the distance between node1 and node2. Distance is defined as the number of edges in the unique path from node1 to node2.
    """

    # Base case; when both nodes are the same, distance between them is 0.
    if node1 == node2 :
        return 0
    # We want to check if node1 is an ancestor of node2 or vice versa
    elif ancestral_table[node1][node2] == 'an': # check if node1 is an ancestor of node2
        return 1 + distance(node1, parent(node2, species_tree), species_tree, ancestral_table) 
    elif ancestral_table[node1][node2] == 'des': #check if node2 is an ancestor of node1
        return 1 + distance(parent(node1, species_tree), node2, species_tree, ancestral_table) 
    # Last case; checking if the two nodes are incomparable
    elif ancestral_table[node1][node2] == 'in':
        return 2 + distance(parent(node1, species_tree), parent(node2, species_tree), species_tree, ancestral_table)

def parent(node, species_tree):
    """
    :param node: a node in the species tree (param species_tree)
    :param species_tree: a species tree in vertex format (see documentation at top of file)
    :return: the parent of the given node in the species tree
    """

    for key in species_tree:
        if species_tree[key][0] == node or species_tree[key][1] == node: # Checking to see if node is a child of key (which indicates key is the parent of node)
            return key
    return node # We will only reach this step if the node doesn't have a parent, which only happens when the node in question is the root of the tree.
                # In this case, we simply return the root itself as its own parent; this does not affect the function in regards to its purpose
                # But even so, we do not expect this to ever happen, because the parent function (used only in the distance() function) will only
                # be called on a node that has an ancestor, which is a condition the root node will never satisfy.


def construct_nodeCost_table(dtl_recon_graph, preorder_mapping_node_list, numPlacing_table, species_tree, ordered_species_node_list, species_ancestral_table):
    """
    :param dtl_recon_graph: A dictionary representing a DTL Recon Graph.
    :param preorder_mapping_node_list: A list of the possible mapping nodes in any MPR, except sorted first in
    preorder by gene node and preorder by species node.
    :param numPlacing_table: a dictionary in which the keys are all of the form (g,s), where g is a node from the gene tree, and s is a node from the gene tree, and 
    the corresponding values are the number of MPRs in DTLReconGraph in which g is placed on s. (See construct_numPlacing_table() for more info)
    :param species_tree: the species tree (in vertex format) used to make the corresponding DTL Recon Graph.
    :param ordered_species_node_list: an ordered list of the nodes (postorder) in the species (or host) tree used to make the corresponding DTL Recon Graph. 
    :param species_ancestral_table: an ancestral table for the species tree (see Diameter.calculate_ancestral_table() for more details).
    :return: a dictionary where keys are mapping nodes of the DTLReconGraph and the corresponding values are the nodeCost of that mapping node. nodeCost is the "cost"
             of a given mapping node being included in a reconciliation. nodeCost for a mapping node (g,s) is calculated as the summation of 
             {distance(s,s') x numPlacing[(g,s')]} for all nodes s' in the species tree, where distance(s,s') is the distance between the nodes s and s' in the species
             tree (see distance() function) and numPlacing[(g,s')] is the value corresponding to the key (g,s') in the dictionary created in 
             the construct_numPlacing_table() function.  
    """
        
    # Initialize the nodeCost dict. This dict contains the nodeCost score of each mapping node
    nodeCost = dict()
    for mapping_node in preorder_mapping_node_list:
        nodeCost[mapping_node] = 0

    # This block implements the algorithm that defines 'nodeCost' of any given mapping node. See docstring for more info.
    for mapping_node in preorder_mapping_node_list:
        for species_node in ordered_species_node_list:
            nodeCost[mapping_node] += (distance(mapping_node[1], species_node, species_tree, species_ancestral_table) * (numPlacing_table[((mapping_node[0]),species_node)]))

    return nodeCost


def compute_path_median(dtl_recon_graph, postorder_mapping_node_list, nodeCost_table, mpr_roots):
    """
    Finds the medoid reconciliation of a given DTLReconGraph using the new path difference distance metric.
    :param dtl_recon_graph: A dictionary representing a DTL Recon Graph.
    :param postorder_mapping_node_list: A list of the possible mapping nodes in any MPR of DTLReconGraph (param dtl_recon_graph), except sorted first in
    postorder by gene node and postorder by species node.
    :param nodeCost_table: a dictionary where keys are mapping nodes of the DTLReconGraph and the corresponding values are the nodeCost of those mapping nodes. See
    construct_nodeCost_table() function for more details.
    :param mpr_roots: A list of mapping nodes that could act as roots to an MPR for the species and
    gene trees in question, output from the reconcile function in DTLReconGraph.py
    :return: 0. A new dictionary which has the same form as a DTL reconciliation graph except this new graph is a subgraph of the original DTLReconGraph, and
                every possible reconciliation in this new reconciliation graph is a medoid as defined by the path distance metric. Thus this dictionary
                represents the medoid reconciliation using the path distance metric.
             1. the number of medoid reconciliations for the given DTL reconciliation graph, using the path distance metric.
             2. the root(s) of the obtained medoid reconciliation graph (returned in 0.).
    """

    # Initialize the totalCost dict. This dict contains the totalCost score of each mapping node and each mapping node's events
    totalCost = dict()
    for mapping_node in postorder_mapping_node_list:
        # First checking for the base case of mapping nodes which are leaves - for which totalCost is 0. 
        if dtl_recon_graph[mapping_node] == [('C', (None, None), (None, None))]: # This tells us if the mapping node in question is a leaf, in which case we initialize its totalCost to 0.
            totalCost[mapping_node] = ([('C',(None,None),(None,None))], 0)
            continue 

        SDTcostList = list()
        lossCostList = list()

        # Here we calculate totalCost for all the events of the relevant mapping_node - to do this, we just need to add the totalCosts of the event's children
        for event in dtl_recon_graph[mapping_node]:
            # Note that 'event' is of the form: ('event ID', 'Child 1', 'Child 2'), so the 0th element is the event
            # ID and the 1st and 2nd elements are the children produced by the event
            if event[0] == 'L':  # Losses produce only one child, so we only need to look to one lower mapping node
                lossCostList.append((event, totalCost[event[1]][1]))
            else:  # Only other options are T, S, and D, which produce two children
                SDTcostList.append((event, totalCost[event[1]][1] + totalCost[event[2]][1]))
            
        # Now we want to calculate totalCost for the map node itself. To do this, first we need to find min{totalCost(e) + nodeCost(g,s)} for all SDT events 'e' of the node and min{totalCost(l)} for all loss events 'l' of the node    
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

        # best_events is a list we intend to populate with all events of the mapping_node that minimize totalCost for mapping_node
        # This is done so we can produce the final complete median reconciliation
        best_events = list()

        # Check to see which event(s) produce the minimum total cost
        for lossEvent in lossCostList:
            if lossEvent[1] == min_cost:
                best_events.append(lossEvent[0])
        for SDTevent in SDTcostList:
            if SDTevent[1] + nodeCost_table[mapping_node] == min_cost:
                best_events.append(SDTevent[0])
        
        # deleting the lists we no longer need
        del lossCostList
        del SDTcostList

        # Save the result for this mapping node so it can be used in higher mapping nodes in the graph
        totalCost[mapping_node] = (best_events[:], min_cost) # best_events[:] slices the list and is a short cut to making a copy of the list that we want

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

    # Check to make sure the medoid reconciliaton is a subgraph of the DTL reconciliation
    assert check_subgraph(dtl_recon_graph, med_recon_graph), 'Median is not a subgraph of the recon graph!'

    # We can use this function to find the number of medians once we've got the final median recon graph
    n_med_recons = DTLReconGraph.count_mprs_wrapper(best_roots, med_recon_graph)

    return med_recon_graph, n_med_recons, best_roots

def median_helper(file_name, dup_cost, transfer_cost, loss_cost):
    """
    A helper function that returns a number of useful values that are used as arguments in many other functions.
    :param file_name:the file in which the desired data set it stored, passed as
    a string. The data files that can be used are currently exclusively .newick files. 
    :param dup_cost: the cost associated with a duplication event
    :param transfer_cost: the cost associated with a transfer event
    :param loss_cost: the cost associated with a loss event
    :return: 0. the DTLReconGraph obtained from the dataset stored in the inputted
                file (param file_name), with associated duplication, transfer and loss costs (also arguments).
             1. A list of the possible mapping nodes in any MPR of the DTLReconGraph (returned as explained in 0.) except sorted first in
                preorder by gene node and preorder by species node.
             2. A list of the possible mapping nodes in any MPR of the DTLReconGraph, except sorted first in
                postorder by gene node and postorder by species node.
             3. a dictionary that returns the numPlacing values of keys of the form (g,s), where g and s are gene tree nodes and species tree
                nodes respectively. See return value of construct_numPlacing_table() function for more details.
             4. vertex formatted tree of the host (species) tree obtained from the dataset stored in the inputted file (param file_name).
             5. a list of the nodes of the host (species) tree in postorder.
             6. san ancestral table for the species tree (see Diameter.calculate_ancestral_table() for more details).
             7. A dictionary with event nodes as keys and values corresponding to the frequency of
                that events in MPR space for the recon graph
             8. the root(s) of the DTLReconGraph (returned as explained in 0.).
             9. the total number of MPRs in the DTLReconGraph.
    """
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

    # event_scores is the dict with normalized frequencies
    event_scores,count = generate_scores(preorder_mapping_node_list,graph,vparas[1], True)
    # event_scores_not_normalized is the dict which has the actual number of MPRs each event is in - also note count and count1 should have the same value
    event_scores_not_normalized,count1 = generate_scores(preorder_mapping_node_list,graph,vparas[1],False)

    # Getting the numPlacingTable dict so we can use it as an argument in future functions
    numPlacingTable = construct_numPlacing_table(graph, postorder_vhost_nodes, postorder_vparas_nodes, postorder_mapping_node_list, event_scores_not_normalized, count)

    # Getting the ancestral table dict for the species (host) tree
    species_ancestral_table = Diameter.calculate_ancestral_table(vhost[0])
    
    # Just a check to see how fast this function works - can comment out later if desired
    print('Helper done!')

    # Returning all possibly relevant values that could be used in other functions or as arguments within other functions
    return graph, preorder_mapping_node_list, postorder_mapping_node_list, numPlacingTable, vhost, postorder_vhost_nodes, species_ancestral_table, event_scores, best_roots, count

def get_maximum_average(file_name, dup_cost, transfer_cost, loss_cost):
    """
    an implementation of the compute_maximum_average() function, which just returns the average calculated in the aforementioned function, and the number of events
    closest to the expected number of events in an MPR of the related DTLReconGraph
    :param file_name: the file in which the desired data set it stored, passed as
    a string. The data files that can be used are currently exclusively .newick files. 
    :param dup_cost: the cost associated with a duplication event
    :param transfer_cost: the cost associated with a transfer event
    :param loss_cost: the cost associated with a loss event
    :return: 0. the largest maximum average frequency of all events in any of the MPRs that are in the DTLReconGraph.
             1. the number of events associated with any MPR that has the maximum average event frequency, which is 
                closest to the expected number of events in an MPR of the DTLReconGraph
    """

    # Using the median helper function (median_helper()) to get useful variable values that will be used as function arguments later
    graph, preorder_mapping_node_list, postorder_mapping_node_list, numPlacingTable, vhost, postorder_vhost_nodes, species_ancestral_table, event_scores, best_roots, count = median_helper(file_name,dup_cost,transfer_cost,loss_cost)

    # Running the function to get the expected number of events in a randomly chosen MPR from 'graph' (the DTLReconGraph obtained from the given file and event costs)
    fileName, expectedEventNum = calculate_expected_frequency_MPR(file_name,dup_cost,transfer_cost,loss_cost)

    # Running the function to get the maximum average support value in an MPR in 'graph' and a list of all possible no. of events in an MPR with that maximum average support
    average, numEventsList = compute_maximum_average(graph,event_scores,postorder_mapping_node_list,best_roots)


    # Finding the value in numEventsList which is closest to expectedEventNum to return just that value
    absDif = float('inf')

    for numEvent in numEventsList:
        if abs(expectedEventNum - numEvent) < absDif: # if numEvent is closer to expectedEventNum then the last closest this is 'True'.
            absDif = abs(expectedEventNum - numEvent) # then the new closest is this value
            finalEventNum = numEvent                  # and the new closest number of events is this value of numEvent

    # Returning the maximum average support value and the one value of numEvent in the numEventsList which is closest to the expected number of events in an MPR.
    return average, finalEventNum


def compare_median_numbers(file_name, dup_cost, transfer_cost, loss_cost):
    """
    :param file_name: the file in which the desired data set it stored, passed as
    a string. The data files that can be used are currently exclusively .newick files. 
    :param dup_cost: the cost associated with a duplication event
    :param transfer_cost: the cost associated with a transfer event
    :param loss_cost: the cost associated with a loss event
    :return: 0. the total number of MPRs in the DTLReconGraph obtained from the dataset stored in the inputted
                file (param file_name), with associated duplication, transfer and loss costs (also arguments).
             1. the number of medoid reconciliations in the medoid reconciliation graph (using the symmetric set difference distance metric) of DTLRecongraph.
             2. the number of medoid reconciliations in the medoid reconciliation graph (using the path difference distance metric) of DTLReconGraph.
    """
    # Using the median helper function (median_helper()) to get useful variable values that will be used as function arguments later
    graph, preorder_mapping_node_list, postorder_mapping_node_list, numPlacingTable, vhost, postorder_vhost_nodes, species_ancestral_table, event_scores, best_roots, count = median_helper(file_name,dup_cost,transfer_cost,loss_cost)

    # Getting the nodeCostTable. NOTE: this line is not included in the median helper because construct_nodeCost_table() takes a much longer time to run
    # than any other function in median_helper() and is also not used as much. So we only run this function if necessary.
    nodeCostTable = construct_nodeCost_table(graph,preorder_mapping_node_list,numPlacingTable,vhost[0],postorder_vhost_nodes,species_ancestral_table)
    # print('Done with getting nodeCostTable')
    
    #Using the first median function (symmetric set difference distance metric) and storing useful values
    med_recon_graph1,n_med_recons1,best_roots1 = compute_median(graph,event_scores,postorder_mapping_node_list,best_roots)

    #Using the second median function (path difference distance metric) and storing useful values
    med_recon_graph2,n_med_recons2,best_roots2 = compute_path_median(graph,postorder_mapping_node_list,nodeCostTable,best_roots)

    #returning the tuple of number of medians for both types of distance metrics
    print('(' + str(n_med_recons1) + ' , ' + str(n_med_recons2) + ')')
    return (count,n_med_recons1,n_med_recons2)


def generate_event_freq_random_median_symmetric_set_method(file_name,dup_cost,transfer_cost,loss_cost,sample_number):
    """
    :param file_name: the file in which the desired data set it stored, passed as
    a string. The data files that can be used are currently exclusively .newick files. 
    :param dup_cost: the cost associated with a duplication event
    :param transfer_cost: the cost associated with a transfer event
    :param loss_cost: the cost associated with a loss event
    :param sample_number: an integer value, which is the number of randomly, uniformly sampled MPRs the user would like data for.
    :return: a list of length equal to sample_number, which contains lists as elements. Each of these list elements is an event frequency list; each of these event
             frequency lists are created by first randomly selecting a medoid from the median reconciliation graph (using the symmetric set difference distance metric)
             of DTLReconGraph which is obtained from the inputted data set (param file_name). Then all the events of this randomly selected medoids are recorded in the 
             list as follows - the type of the event is recorded ('Speciation','Duplication','Loss' or 'Transfer') followed by its frequency (see documentation at top
             for more info on event frequency lists). Also note we do NOT include contemporaneous events in this list as we already know they are in all MPRs.
    """
    # Using the median helper function (median_helper()) to get useful variable values that will be used as function arguments later
    graph, preorder_mapping_node_list, postorder_mapping_node_list, numPlacingTable, vhost, postorder_vhost_nodes, species_ancestral_table, event_scores, best_roots, count = median_helper(file_name,dup_cost,transfer_cost,loss_cost)
    
    #Using the first median function (symmetric set difference distance metric) and storing useful values
    med_recon_graph1,n_med_recons1,best_roots1 = compute_median(graph,event_scores,postorder_mapping_node_list,best_roots)

    # Using the count_mprs function to get no. of MPRs spawned below each of the best roots of the median reconciliation, to use later when we randomly select medians 
    # from the medoid reconciliation graph
    med_counts1 = dict()
    for root in best_roots1:
        count_mprs(root, med_recon_graph1, med_counts1)
    
    # A list that will be populated with a number of randomly selected medians (the number is inputted by the user - param sample_number)
    random_median_list = list()

    # Randomly selecting medians and populating random_median_list
    for x in range(sample_number):
        curr_random_median = choose_random_median_wrapper(med_recon_graph1, best_roots1, med_counts1)
        random_median_list.append(curr_random_median)

    # freq_scores1 will be populated with event frequency lists corresponding to each of the medians in random_median_list
    freq_scores1 = list()
    # temporary list that creates each event frequency list and adds it to freq_scores1
    current_freq_score_list = list()
    
    # Looping through all medians in random_median_list to create an event_frequency_list for each of them
    for x in range(sample_number):
        # Nested loop to go through all events in a given median
        for mapping_node in random_median_list[x]:
            for event in random_median_list[x][mapping_node]:
                # Checking what event the event is so that information can be added to the event frequency list
                if event == ('C',(None,None),(None,None)): # If the event is contemporaneous, we ignore it
                    continue
                if event[0] == 'D':
                    current_freq_score_list.append('Duplication')
                elif event[0] == 'T':
                    current_freq_score_list.append('Transfer')
                elif event[0] == 'L':
                    current_freq_score_list.append('Loss')
                elif event[0] == 'S':
                    current_freq_score_list.append('Speciation')
                # After adding what the event is to the list, we add what the frequency of that event is to the list.
                current_freq_score_list.append(event_scores[event]) #event_scores is a dictionary storing normalized frequencies for all the events in the reconciliation graph, which is the same for both median types
        freq_scores1.append(current_freq_score_list) # Adding the current event frequency list to freq_scores1
        current_freq_score_list = [] # Reinitializing current_freq_score_list to an empty list so we can use it again.
    
    return freq_scores1

def generate_event_freq_random_median_path_difference_method(file_name,dup_cost,transfer_cost,loss_cost,sample_number):
    """
    :param file_name: the file in which the desired data set it stored, passed as
    a string. The data files that can be used are currently exclusively .newick files. 
    :param dup_cost: the cost associated with a duplication event
    :param transfer_cost: the cost associated with a transfer event
    :param loss_cost: the cost associated with a loss event
    :param sample_number: an integer value, which is the number of randomly, uniformly sampled MPRs the user would like data for.
    :return: a list of length equal to sample_number, which contains lists as elements. Each of these list elements is an event frequency list; each of these event
             frequency lists are created by first randomly selecting a medoid from the median reconciliation graph (using the path difference distance metric)
             of DTLReconGraph which is obtained from the inputted data set (param file_name). Then all the events of this randomly selected medoids are recorded in the 
             list as follows - the type of the event is recorded ('Speciation','Duplication','Loss' or 'Transfer') followed by its frequency (see documentation at top
             for more info on event frequency lists). Also note we do NOT include contemporaneous events in this list as we already know they are in all MPRs.
    """

    # Using the median helper function (median_helper()) to get useful variable values that will be used as function arguments later
    graph, preorder_mapping_node_list, postorder_mapping_node_list, numPlacingTable, vhost, postorder_vhost_nodes, species_ancestral_table, event_scores, best_roots, count = median_helper(file_name,dup_cost,transfer_cost,loss_cost)
    
    # Getting the nodeCostTable. NOTE: this line is not included in the median helper because construct_nodeCost_table() takes a much longer time to run
    # than any other function in median_helper() and is also not used as much. So we only run this function if necessary.
    nodeCostTable = construct_nodeCost_table(graph,preorder_mapping_node_list,numPlacingTable,vhost[0],postorder_vhost_nodes,species_ancestral_table)

    #Using the second median function (path difference distance metric) and storing useful values
    med_recon_graph2,n_med_recons2,best_roots2 = compute_path_median(graph,postorder_mapping_node_list,nodeCostTable,best_roots)

    # Using the count_mprs function to get no. of MPRs spawned below each of the best roots of the median reconciliation, to use later when we randomly select medians 
    # from the medoid reconciliation graph.
    med_counts2 = dict()
    for root in best_roots2:
        count_mprs(root, med_recon_graph2, med_counts2)

    # A list that will be populated with a number of randomly selected medians (the number is inputted by the user - param sample_number)
    random_median_list = list()

    # Randomly selecting medians and populating random_median_list
    for x in range(sample_number):
        curr_random_median = choose_random_median_wrapper(med_recon_graph2, best_roots2, med_counts2)
        random_median_list.append(curr_random_median)

    # freq_scores2 will be populated with event frequency lists corresponding to each of the medians in random_median_list
    freq_scores2 = list()
    # temporary list that creates each event frequency list and adds it to freq_scores1
    current_freq_score_list = list()

    # Looping through all medians in random_median_list to create an event_frequency_list for each of them
    for x in range(sample_number):
        # Nested loop to go through all events in a given median
        for mapping_node in random_median_list[x]:
            for event in random_median_list[x][mapping_node]:
                # Checking what event the event is so that information can be added to the event frequency list
                if event == ('C',(None,None),(None,None)):
                    continue # If the event is contemporaneous, we ignore it
                if event[0] == 'D':
                    current_freq_score_list.append('Duplication')
                elif event[0] == 'T':
                    current_freq_score_list.append('Transfer')
                elif event[0] == 'L':
                    current_freq_score_list.append('Loss')
                elif event[0] == 'S':
                    current_freq_score_list.append('Speciation')
                # After adding what the event is to the list, we add what the frequency of that event is to the list.
                current_freq_score_list.append(event_scores[event]) #event_scores is a dictionary storing normalized frequencies for all the events in the reconciliation graph, which is the same for both median types
        freq_scores2.append(current_freq_score_list) # Adding the current event frequency list to freq_scores1
        current_freq_score_list = [] # Reinitializing current_freq_score_list to an empty list so we can use it again.
    
    return freq_scores2

def generate_event_freq_random_reconciliation(file_name,dup_cost,transfer_cost,loss_cost,sample_number):
    """
    :param file_name: the file in which the desired data set it stored, passed as
    a string. The data files that can be used are currently exclusively .newick files. 
    :param dup_cost: the cost associated with a duplication event
    :param transfer_cost: the cost associated with a transfer event
    :param loss_cost: the cost associated with a loss event
    :param sample_number: an integer value, which is the number of randomly, uniformly sampled MPRs the user would like data for.
    :return: a list of length equal to sample_number, which contains lists as elements. Each of these list elements is an event frequency list; each of these event
             frequency lists are created by first randomly selecting a reconciliation from the DTLReconGraph obtained from the inputted data set (param file_name).
             Then all the events of this randomly selected medoids are recorded in the list as follows - the type of the event is recorded 
             ('Speciation','Duplication','Loss' or 'Transfer') followed by its frequency (see documentation at top for more info on event frequency lists). Also note 
             we do NOT include contemporaneous events in this list as we already know they are in all MPRs.
    """

    # Using the median helper function (median_helper()) to get useful variable values that will be used as function arguments later
    graph, preorder_mapping_node_list, postorder_mapping_node_list, numPlacingTable, vhost, postorder_vhost_nodes, species_ancestral_table, event_scores, best_roots, count = median_helper(file_name,dup_cost,transfer_cost,loss_cost)
    
    # Using the count_mprs function to get no. of MPRs spawned below each of the best roots of the DTLReconGraph ('graph'), to use later when we randomly select 
    #reconciliations from the DTL reconciliation graph.
    counts = dict()
    for root in best_roots:
        count_mprs(root, graph, counts)

    # A list that will be populated with a number of randomly selected reconciliations (the number is inputted by the user - param sample_number)
    random_reconciliation_list = list()

    # Randomly selecting reconciliations and populating random_reconciliation_list
    for x in range(sample_number):
        curr_random_reconciliation = choose_random_median_wrapper(graph, best_roots, counts)
        random_reconciliation_list.append(curr_random_reconciliation)

    # freq_scores will be populated with event frequency lists corresponding to each of the reconciliations in random_reconciliation_list
    freq_scores = list()
    # temporary list that creates each event frequency list and adds it to freq_scores1
    current_freq_score_list = list()

    # Looping through all reconciliations in random_reconciliation_list to create an event_frequency_list for each of them
    for x in range(sample_number):
        # Nested loop to go through all events in a given reconciliation
        for mapping_node in random_reconciliation_list[x]:
            for event in random_reconciliation_list[x][mapping_node]:
                # Checking what event the event is so that information can be added to the event frequency list
                if event == ('C',(None,None),(None,None)): # If the event is contemporaneous, we ignore it
                    continue
                if event[0] == 'D':
                    current_freq_score_list.append('Duplication')
                elif event[0] == 'T':
                    current_freq_score_list.append('Transfer')
                elif event[0] == 'L':
                    current_freq_score_list.append('Loss')
                elif event[0] == 'S':
                    current_freq_score_list.append('Speciation')
                # After adding what the event is to the list, we add what the frequency of that event is to the list.
                current_freq_score_list.append(event_scores[event]) #event_scores is a dictionary storing normalized frequencies for all the events in the reconciliation graph
        freq_scores.append(current_freq_score_list) # Adding the current event frequency list to freq_scores
        current_freq_score_list = [] # Reinitializing current_freq_score_list to an empty list so we can use it again.
    
    return freq_scores

def calculate_expected_frequency_MPR(file_name,dup_cost,transfer_cost,loss_cost):
    """
    :param file_name: the file in which the desired data set it stored, passed as
    a string. The data files that can be used are currently exclusively .newick files. 
    :param dup_cost: the cost associated with a duplication event
    :param transfer_cost: the cost associated with a transfer event
    :param loss_cost: the cost associated with a loss event
    :return: 0. the file name as inputted (param file_name)
             1. the expected number of events that will be in a randomly selected MPR associated with the DTLReconGraph obtained from the dataset stored in the inputted
                file (param file_name), with associated duplication, transfer and loss costs (also arguments). This expected value is calculated by summing the
                frequencies of all events (except contemporaneous events, which we know are in all MPRs) present in the obtained DTLReconGraph, where frequency 
                of an event is the fraction of MPRs it is present in, out of all MPRs that are obtainable from the DTLReconGraph.
    """
    # Using the median helper function (median_helper()) to get useful variable values that will be used as function arguments later
    graph, preorder_mapping_node_list, postorder_mapping_node_list, numPlacingTable, vhost, postorder_vhost_nodes, species_ancestral_table, event_scores, best_roots, count = median_helper(file_name,dup_cost,transfer_cost,loss_cost)
    
    # totalCount will add together the frequency of all events in the DTLReconGraph ('graph') EXCEPT the contemporaneous events. This will be the expected number of events.
    totalCount = 0
    totalCountcheck = 0 # This is a variable which will be used in calculating expected number of events in another way -  as a check.

    # Going through all events
    for event in event_scores:
        if event == ('C', (None,None), (None,None)): # If the event is contemporaneous, we ignore it
            continue
        # Otherwise, we add its frequency to totalCount
        totalCount +=  event_scores[event]

# Calculating totalCount in a different way to check our answer is correct
    # Nested loop through each mapping node and each of its events
    for mapping_node in graph:
        for event in graph[mapping_node]:
            if event == ('C', (None, None), (None, None)): # If the event is contemporaneous, we ignore it
                continue
            # Otherwise add its frequency to totalCountCheck
            totalCountcheck += event_scores[event]
    
# Note that we did get different values for totalCount and totalCountCheck, but they differed only around the 9th decimal place.
    
    return file_name,totalCount

def calculate_min_and_avg_frequency(file_name,dup_cost,transfer_cost,loss_cost):
    """
    :param file_name: the file in which the desired data set it stored, passed as
    a string. The data files that can be used are currently exclusively .newick files. 
    :param dup_cost: the cost associated with a duplication event
    :param transfer_cost: the cost associated with a transfer event
    :param loss_cost: the cost associated with a loss event
    :return: 0. the file in which the desired data set is stored (same as param file_name)
             1. the average frequency of all events in the DTLReconGraph obtained from the dataset stored in the inputted
                file (param file_name), with associated duplication, transfer and loss costs (also arguments). We disregard contemporaneous events in this calculation.
             2. the minimum frequency out of all event frequencies in the obtained DTLReconGraph.
    """
    # Using the median helper function (median_helper()) to get useful variable values that will be used as function arguments later
    graph, preorder_mapping_node_list, postorder_mapping_node_list, numPlacingTable, vhost, postorder_vhost_nodes, species_ancestral_table, event_scores, best_roots, count = median_helper(file_name,dup_cost,transfer_cost,loss_cost)
    
    # totalCount will ultimately hold the sum of frequency of all events (except contemporaneous) in graph (DTLReconGraph).
    totalCount = 0
    # eventCount counts the number of events in the DTLReconGraph
    eventCount = 0
    # event_frequency_list will hold the frequencies of all events in the DTLReconGraph
    event_frequency_list = list()

    # Looping through all events
    for event in event_scores:
        if event == ('C', (None,None), (None,None)): # If the event is contemporaneous, we ignore it
            continue
        totalCount +=  event_scores[event]
        eventCount += 1
        event_frequency_list.append(event_scores[event])
    
    #Finding the average event frequency by dividing the sum of all event frequencies over the total number of events
    average_frequency = totalCount/float(eventCount)
    # Finding the minimum frequency out of all event frequencies in DTLReconGraph
    min_frequency = min(event_frequency_list)

    return (file_name,average_frequency,min_frequency)


def eventFreqscript(filename,distmetric,dup_cost,transfer_cost,loss_cost,outputFilename,sample_number):
    """
    This function is an implementation of the functions that randomly generate medians/reconciliations and return their event frequency lists. eventFreqscript()
    is used to write these event frequency lists into .csv files.
    :param filename: the file in which the desired data set it stored, passed as
    a string. The data files that can be used are currently exclusively .newick files. 
    :param distmetric: a string, indicating which reconciliation graph the user would like to select random medians from. 'r' means from the DTLReconGraph obtained
    from the dataset stored in the inputted file (param file_name), with associated duplication, transfer and loss costs (also arguments). 's' means the 
    median reconciliation obtained from DTLReconGraph using symmetric set difference distance metric, and 'p' means the median reconcilaition obtained from
    DTLReconGraph using the path difference distance metric.
    :param dup_cost: the cost associated with a duplication event.
    :param transfer_cost: the cost associated with a transfer event.
    :param loss_cost: the cost associated with a loss event.
    :param outputFilename: the name of the .csv file user would like to output into.
    :param sample_number: an integer value, which is the number of randomly, uniformly sampled MPRs the user would like data for.
    :return: Nothing. This function is an implementation of one of generate_event_freq_random_reconciliation (if distmetric == 'r'),
             generate_event_freq_random_median_symmetric_set_method() (if distmetric == 's'), and 
             generate_event_freq_random_median_path_difference_method() (if distmetric == 'p') 
    """

    # csvOutputFilename will be the name of the .csv file data will be written into
    csvOutputFilename = outputFilename + '.csv'
    
    # Checking if the distance metric inputted was valid
    if distmetric!= 's' and distmetric!= 'p'and distmetric!= 'r':
        print('Sorry, that is not a valid distance metric input')
        return

    # Opening the .csv file into which the event frequency lists will be outputted
    with open(csvOutputFilename, 'wb') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        # Running different functions based on what distance metric was input.
        if distmetric == 's':
            freq_scores1 = generate_event_freq_random_median_symmetric_set_method(filename,dup_cost,transfer_cost,loss_cost,sample_number)
            for medianScore in freq_scores1:
                filewriter.writerow(medianScore)
        elif distmetric == 'p':
            freq_scores2 = generate_event_freq_random_median_path_difference_method(filename,dup_cost,transfer_cost,loss_cost,sample_number)
            for medianScore in freq_scores2:
                filewriter.writerow(medianScore)
        elif distmetric == 'r':
            freq_scores3 = generate_event_freq_random_reconciliation(filename,dup_cost,transfer_cost,loss_cost,sample_number)
            for reconScore in freq_scores3:
                filewriter.writerow(reconScore)



def count_leaves(file_name,dup_cost,transfer_cost,loss_cost):
    """
    :param file_name: the file in which the desired data set it stored, passed as
    a string. The data files that can be used are currently exclusively .newick files. 
    :param dup_cost: the cost associated with a duplication event
    :param transfer_cost: the cost associated with a transfer event
    :param loss_cost: the cost associated with a loss event
    :return: the number of leaf nodes in the data set associated with the inputted file (param file_name)
    """

    host,paras,graph,num_recon,best_roots = DTLReconGraph.reconcile(file_name,dup_cost,transfer_cost,loss_cost)

    leafCount = 0

    for mapping_node in graph:
        for event in graph[mapping_node]:
            if event == ('C', (None, None), (None, None)): # Every contemporaneous event corresponds to exactly one leaf-to-leaf mapping
                leafCount += 1
    
    return leafCount



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
            scores_dict = generate_scores(postorder_mapping_node_list[::-1], dtl_recon_graph, gene_tree_root, True)

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
        filewriter.writerow(['File Name','Average Mean','Average Minimum','Average Maximum','Average S. Dev', 'Average # of Columns'])

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
            rowSDevs = [None]*numRows
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
                rowSDevs[i] = np.std(rows[i],ddof=1)        #ddof is 1 because we're calculating a sample sdev

            #Next we find the mean of the averages, the minimums, and the maximums
            meanOfAverages = sum(rowAverages)/numRows
            meanOfMinimums = sum(rowMinimums)/numRows
            meanOfMaximums = sum(rowMaximums)/numRows
            meanOfSDevs = sum(rowSDevs)/numRows
            meanColumns = sum(numColumns)/numRows
            
            #write to the output csv file
            filewriter.writerow([file, meanOfAverages, meanOfMinimums, meanOfMaximums, meanOfSDevs, meanColumns])
 
    #Path 111
    with open('Path_111.csv', 'a') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['File Name','Average Mean','Average Minimum','Average Maximum','Average S. Dev','Average # of Columns'])

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
            rowSDevs = [None]*numRows
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
                rowSDevs[i] = np.std(rows[i],ddof=1)        #ddof is 1 because we're calculating a sample sdev

            #Next we find the mean of the averages, the minimums, and the maximums
            meanOfAverages = sum(rowAverages)/numRows
            meanOfMinimums = sum(rowMinimums)/numRows
            meanOfMaximums = sum(rowMaximums)/numRows
            meanOfSDevs = sum(rowSDevs)/numRows
            meanColumns = sum(numColumns)/numRows
            
            #write to the output csv file
            filewriter.writerow([file, meanOfAverages, meanOfMinimums, meanOfMaximums, meanOfSDevs, meanColumns])

    #Symmetric 231
    with open('Symmetric_231.csv', 'a') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['File Name','Average Mean','Average Minimum','Average Maximum','Average S. Dev', 'Average # of Columns'])

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
            rowSDevs = [None]*numRows
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
                rowSDevs[i] = np.std(rows[i],ddof=1)        #ddof is 1 because we're calculating a sample sdev

            #Next we find the mean of the averages, the minimums, and the maximums
            meanOfAverages = sum(rowAverages)/numRows
            meanOfMinimums = sum(rowMinimums)/numRows
            meanOfMaximums = sum(rowMaximums)/numRows
            meanOfSDevs = sum(rowSDevs)/numRows
            meanColumns = sum(numColumns)/numRows
            
            #write to the output csv file
            filewriter.writerow([file, meanOfAverages, meanOfMinimums, meanOfMaximums, meanOfSDevs, meanColumns])
 
    #Path 231
    with open('Path_231.csv', 'a') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(['File Name','Average Mean','Average Minimum','Average Maximum','Average S. Dev','Average # of Columns'])

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
            rowSDevs = [None]*numRows
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
                rowSDevs[i] = np.std(rows[i],ddof=1)        #ddof is 1 because we're calculating a sample sdev

            #Next we find the mean of the averages, the minimums, and the maximums
            meanOfAverages = sum(rowAverages)/numRows
            meanOfMinimums = sum(rowMinimums)/numRows
            meanOfMaximums = sum(rowMaximums)/numRows
            meanOfSDevs = sum(rowSDevs)/numRows
            meanColumns = sum(numColumns)/numRows
            
            #write to the output csv file
            filewriter.writerow([file, meanOfAverages, meanOfMinimums, meanOfMaximums, meanOfSDevs, meanColumns])
