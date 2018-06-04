# First, we put all the functions that are used to find frequency (number of reconciliations involved) in any given node (event or mapping). To do this, we have three functions in total

def count_mprs(event_node, dtl_recon_graph, counts):
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
    :return: the number of MPRs spawned below the given mapping (event) node in the graph
    """

    # Search the counts dictionary for previously calculated results (this is the memoization)
    if event_node in counts:
        return counts[event_node] # RR - this method is memoization, NOT dynamic programming. That's why we're going top-down

    # Base case, occurs if being called on a child produced by a loss or contemporary event 
    if mapping_node == (None, None): # RR - we check for losses or contemporary events first so that when we loop through events later (line 28-29) we know those will only be duplications, speciations or transfers
        return 1 

    # Initialize a variable to keep count of the number of MPRs
    count = 0

    # Loop over all event nodes corresponding to the current mapping node
    for eventNode in dtl_recon_graph[mapping_node]: # RR - do we know that all of these events are only duplications, speciations and transfers? because only for these are we guaranteed two children. Edit: see line 22

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
    #            event_scores[event] = event_scores[event] / float(count)
    # We remove this above line because at the current time, we don't want the normalized score, we want the absolute score, it is used in a multiplication in construct_numPlacing_table().

    return event_scores, count


# All previously written functions are above. Below are new functions + computeMedian
 

# A function we use to create a dict of the numPlacing function applied to all mapping nodes (g,s). numPlacing() gives us the number of reconciliations that have g placed on s (i.e, a loss event doesn't directly follow the node (g,s)) 
# numPlacing can be calculated by first going through all the events of a given mapping node (g,s). Then, we look at only the events that are not losses, and add together all of their frequencies. 
# This will give us the net number of reconciliations that use a given node (g,s) where g is mapped onto s 
# (because if the node (g,s) exists and g is NOT mapped onto s, it must be because there is a loss event following that mapping node in the reconciliation).
# In our arguments, count is the total number of reconciliations that we get from our reconciliation graph. This value is returned in the generate_scores function

def construct_numPlacing_table(dtl_recon_graph, ordered_species_node_list, ordered_gene_node_list, preorder_mapping_node_list, event_scores, count):


    # Initialize the numPlacing dict. This dict contains the frequency score of each mapping node
    numPlacing = dict()
    # We want values in this dictionary for every possible combination of gene node x species node, not just the mapping nodes. But the values
    #for the combinations which are not valid mapping nodes will be trivially 0 (makes our calculations simpler). 
    for gene_node in ordered_gene_node_list:
        for species_node in ordered_species_node_list:
            numPlacing[(gene_node, species_node)] = 0

    # Now, we will just change the values for the ones that are mapping nodes, and may have non-zero values
    for map_node in preorder_mapping_node_list:
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
#Helper function to convert an ordered node list to a left-oriented(?) list that has an easier structure to work with (maybe)
def convert_node_list(ordered_node_list):
    for key in ordered_node_list:
        if ordered_node_list[key] = (None, None):
            return (ordered_node_list[key], None, None)
        else:
            return (key, convert_node_list(ordered_node_list[key][0]), convert_node_list(ordered_node_list[key][1]))

#Helper function most_recent_ancestor(node1, node2, tree) which gives us the most recent ancestor of two nodes in a tree

def most_recent_ancestor(node1, node2, tree_list):


# A function that calculates the distance between any two nodes on a given tree, where distance is defined as the
# number of edges in the path from one node to the other (since this is a tree, we know there is but one unique path from one node to any other)
# This function will be defined recursively
# Recall that species_tree is in vertex tree format (see DTLMedian) and that it is in the form of a dictionary. 
# The key is nodes and the correspondence are the children of the key node.
def distance(node1, node2, species_tree):
    ancestral_table = calculate_ancestral_table(species_tree)
    # Base case; when both nodes are the same, distance between them is 0.
    if node1 == node2 :
        return 0
    # We want to check if node1 is an ancestor of node2 or vice versa
    elif ancestral_table[node1][node2] == 'an' # check if node1 is an ancestor of node2
        return 1 + distance(node1, parent(node2, species_tree), species_tree) # we need a function to get the parent of a given node in a tree
    elif ancestral_table[node1][node2] == 'des' #check if node2 is an ancestor of node1
        return 1 + distance(parent(node1, species_tree), node2, species_tree) # ditto 
    # Last case; checking if the two nodes are incomparable
    elif ancestral_table[node1][node2] == 'in':
        return 2 + distance(parent(node1, species_tree), parent(node2, species_tree), species_tree)


# to finish this function, we need to create two helper functions; one which tells us if node1 is an ancestor of node2 (or vice versa) and another function that takes a node and a tree, and returns the parent of that node

# Helper function that helps us to find the parent of a given node, given the ordered node list of the tree in which the node resides
# Check if you can do this directly given a (vertex) tree? Seems like it should be possible. 
# Trees most commonly used here are vertex trees. So use that format
# Recall that species_tree in vertex tree format (see DTLMedian) is in the form of a dictionary. 
# The key is nodes and the correspondence are the children of the key node.

def parent(node, species_tree):
    for key in species_tree:
        if species_tree[key][0] == node or species_tree[key][1] == node:
            return key
    return node # We will only reach this step if the node doesn't have a parent, which only happens when the node in question is the root of the tree.
                # In this case, we simply return the root itself as its own parent; this does not affect the function in regards to its purpose
                # But even so, we do not expect this to ever happen, because the parent function will only be called on a node that has an ancestor
                # which is a condition the root node will never satisfy.


# A function that creates a dict of the nodeCost function applied to all mapping nodes (g,s). nodeCost() gives us the "cost" of a given node being included in a reconciliation through its difference from other reconciliations.
# nodeCost() is calculated as the summation (for all nodes s' from the species tree) of (distance(s', s) x numPlacing[(g,s' )]).
# distance(s',s) is the distance between the nodes s' and s in the species tree.
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
            nodeCost[mapping_node] += (distance(mapping_node[1], species_node, species_tree) * (numPlacing_table[((mapping_node[0]),species_node)]))
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
        totalCost[mapping_node] = 0
        for event in dtl_recon_graph[mapping_node]:
            totalCost[event] = 0

    for mapping_node in postorder_mapping_node_list:
        # First checking for the base case of mapping nodes which are leaves
        if dtl_recon_graph[mapping_node] == [('C', (None, None), (None, None))]:
            totalCost[mapping_node] = (('C',(None,None),(None,None)), 0)
            continue # This tells us if the mapping node in question is a leaf, in which case we don't do anything (since the totalCost of a leaf is 0)

        SDTcostList = list()
        lossCostList = list()
        for event in dtl_recon_graph[map_node]:
            # Note that 'event' is of the form: ('event ID', 'Child 1', 'Child 2'), so the 0th element is the event
            # ID and the 1st and 2nd elements are the children produced by the event
            if event[0] == 'L':  # Losses produce only one child, so we only need to look to one lower mapping node
                lossCostList.append((event, totalCost[event[1]][1]))
            else:  # Only other options are T, S, and D, which produce two children
                SDTcostList.append((event, totalCost[event[1]][1] + totalCost[event[2]][1]))
            
        # Now we want to calculate totalCost for the map node itself. To do this, first we need to find min{totalCost(e) + nodeCost(g,s)} for all SDT events e of the node and min{totalCost(l)} for all loss events l of the node    
        if lossCostList == []:
            minLossCost = math.inf
        else:
            minLossCost = min(lossCostList, key=itemgetter(1))[1]

        if SDTcostList == []:
            minSDTcost = math.inf
        else:
            minSDTcost = min(SDTcostList, key=itemgetter(1))[1]
            minSDTcost += nodeCost_Table[mapping_node] #because if we go through a SDT event (meaning g mapped onto s) then we also need to add the nodeCost
        
        min_cost = min(minLossCost, minSDTcost)

        best_events = list()

        # Check to see which event(s) produce the minimum total cost
        for lossEvent in lossCostList:
            if lossEvent[1] == min_cost:
                best_events.append(lossEvent[0])
        for SDTevent in SDTcostList:
            if SDTevent[1] + nodeCost_Table[mapping_node] == min_cost:
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

    # Adjust the sum_freqs dictionary so we can use it with the buildDTLReconGraph function from DTLReconGraph.py
    for map_node in totalCost:

        # We place the event tuples into lists so they work well with the diameter algorithm
        totalCost[map_node] = totalCost[map_node][0]  # Only use the events, no longer the associated frequency sum

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



    
# Can we reuse this code for finding the median reconciliation? If we could convert the frequencies (nodeCosts) to weighted frequencies somehow
# perhaps the same formula would apply. However, note that in the original compute_median, our aim was to maximize the sum of the frequencies 
# of the nodes in the reconciliations, whereas in the new method, we are aiming to minimize the sum of nodeCost of the nodes in the reconciliations   
# so these may be irreconcilable (ha!) or incomparable problems. Reusing the code would be nice though. 

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
    represent a single reconciliation: the median reconciliation. (would a better way to phrase this be - our graph will represent a subgraph of our DTL Recon Graph, such that all reconciliations possible to reconcile from the subgraph are median reconciliations)?
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
