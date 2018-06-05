# dtldiameter and dtlmedian

This code is from Summers 2017 and 2018.  
In Summer 2017, Eli Zupke and Andrew Ramirez developed much of this code to implement Jordan Haack's DTL diameter algorithm (see the resulting APBC 2018 paper for more details)
In Summer 2018, Ramya Ramilingam and Kalani Duran added Melissa Grueter's path-based distance algorithm to compute medoids wrt to path distance.

Note that all of the included code was created using Python 2.7. Users wishing to run these programs on a machine that has a later version of Python may experience several syntax-related errors.

## How to use DTLReconGraph.py

### From the Command Line:

After typing `python DTLReconGraph.py`, one should enter the values for:

  * The name of the file in which the data are stored and
  * The costs of duplication, transfer, and loss, 
  
in that order. If the user doesn't give input in the expected form, a usage statement will be displayed showing the user what inputs to enter and how to enter them. Note that the data file should include: 
  * Two newick trees (formatted, for example, as ((C,D)B,E)A if A has children B and E and B has children C and D) on separate lines, with the species/host tree given first and the gene/parasite tree second, and
  * Mappings between all of the tips of the two trees (formatted, for example, c:C; d:D; e:E if c maps to C, d to D and e to E), with a different one on each line.

This code then prints: 
  * Two Python dictionaries in the format {(parent, child):(parent, child, child 1 of child, child 2 of child):, ...}, where parent is a node in the tree, child is one of parent's children, and children 1 and 2 of child represent children of the child of parent (in other words, grandchildren of the parent) - one for each of the two trees used as input. This is essentially a dictionary representation of the two trees given as input;
  * A Python dictionary with mapping nodes as keys and lists of event nodes that apply to those mapping nodes in a Maximum Parsimony Reconciliation as values, which has the format {(gene node, species node): [(event string, resulting mapping node 1, resulting mapping node 2), minimum cost], ...}. In this representation, the gene and species node values represent the nodes in the given newick trees that are being mapped to each other in the reconciliation; the event string is a character ('D', 'T', 'L', 'C', or 'S', representing duplication, transfer, loss, contemporary event, and speciation, respectively) representing the event that those mapped nodes undergo; lastly, the resulting mapping nodes are the mappings that are induced in the next 'level' of the reconciliation graph as a result of the given event;
  * The number of Maximum Parsimony Reconciliations for the given data set and costs, as an integer
  * A list of mapping nodes we could use to root the gene tree to produce a maximum parsimony reconciliation

Note each of these printed values is separated by a blank line, so the output is easily readable by a user.

### Via Interactive Mode:

One should simply type `python -i DTLReconGraph.py` to access the contents of the code. From there, all major and helper functions will be available, however the two most important are DP and reconcile. 

#### DP

DP is the workhorse of the code - it utilizes several helper functions to actually perform computations and implement the algorithm given in the technical report titled "HMC CS Technical Report CS-2011-1: Faster Dynamic Programming Algorithms for the Cophylogeny Reconstruction Problem". For more details on the algorithm itself, see that report and/or the comments/docstrings included in the file. Given host and parasite trees, tip mappings, and costs for duplication, transfer, and loss, DP implements this algorithm and returns the DTL maximum parsimony reconciliation graph as a dictionary, the total cost associated with the best reconciliation for the given D, T, and L costs, the number of reconciliations for the given trees and tip mappings, and a list of all possible mapping nodes for the root of the gene tree that could be used in an MPR (see the outputs discussed in the section on using this code from the command line for more detail on the output format).

#### Reconcile

Reconcile is the more practically useful function. Since the data are implemented as newick trees and mappings, reconcile utilizes a separate module that both handles getting the data from a separate file and reformats the inputs to work nicely with DP. Reconcile reformat the species and gene trees to match the output format given in the section on running from the command line, and prints out these trees along with the reconciliation graph, the number of Maximum Parsimony Reconciliations as an integer, and a list of mappings (as tuples) of gene nodes onto species nodes that could produce an MPR. Although one may play with any/all functions included in this file and those on which it depends, DP and reconcile are the most important functions.

## How to use RunTests.py

The main program used for calculating diameters and medians is RunTests.py.

### The RunTests.py CLI

Running RunTests from the command line has the following usage pattern:

> RunTests.py [options] file d t l

In this command, `file` is the relative path to the newick file you wish to compute the diameter of, `d`, `t`, and `l` are the event costs. Specifically, `d` is the event cost of duplications, `t` is the event cost of transfers, and `l` is the event cost of losses. This command will calculate the diameter of the file, and print the result to the screen.

#### Optional Flags

There are several flags availible:

* `-l` or `--log` takes one argument, a path to a csv file. It will log the results of the calculation to a row in the provided csv file, creating the file if necessary. Be sure to specify a path without a file extension, as this program will add one for you. 

* `-q` or `--quiet` supresses all text output (except that generated by the `-i` flag).

* `-i` or `--iterate` takes two arguments, and tells the program to iterate over several numbered files and perform the calculation. This first argument is the starting number, and the second argument is the ending number. Note that when this flag is active, the `file` required argument will replace any `#` characters in the file name with whatever number it is currently on, and use that file in its current calculation.

* `-L` or `--loud` prints the bell character after each failed iteration when using the `-i` flag.

* `-d` or `--debug` prints out every (reasonably sized) dynamic programming table generated by this program.

Some flags add extra computations that are output to the screen or logged to a csv file, in addition to the diameter. Those flags are:

* `-z` or `--zero-loss` calculates the diameter when losses do not contribute to the total.

* `-m` or `--median-count` calculates the number of median reconciliations

* `-M` or `--median-diameter` calculates the diameter of the median space

* `-w` or `--worst-median` calculates the distance of a median reconciliation to the furthest MPR, when the median reconciliation is chosen to maximize that distance.

* `-c` or `--cluster` takes one argument, a number. It will select that many median reconciliations at random, and find the distance from each one to the furthest MPR. Usage of the log (`-l`) flag is mandatory when using this mode as it will generate a log file for each file that is calculated, in addton to recording the average, mean, best, and worst distances computed, and the total number of unique medians seen, in the regular logfile.

There is also limited capability to output reconciliation graphs to SIF files, for viewing in external tools, like Cytoscape. These are controlled by two flags:

* `-s` or `--save-graph` saves the entire MPR DTL reconciliation graph to an SIF file in the folder `graphs/`.

* `-S` or `--save-median-graph` saves the median-space DTL reconciliation graph to an SIF file in the folder `graphs/`.


#### Examples

To calculate a single diameter from file `tree1.newick` with costs `D=2`, `T=3`, and `L=1`, use

> python RunTests.py tree1.newick 2 3 1

To log the result to the end of a log file `Results.csv`, use

> python RunTests.py -l Results tree1.newick 2 3 1 

To calculate the diameters for files `tree1.newick`, `tree2.newick`, and `tree4.newick`, use

>python RunTests.py -i 1 5 tree#.newick 2 3 1

The program will skip over the non-existent `tree3.newick`.

To calculate the diameters for every file in the TreeLifeData folder, logging the result to COG_results.csv, while suppressing text and showing tables, use

>python RunTests.py -i 1 6000 -l COG_results -qd TreeLifeData/COG####.newick 2 3 1

To calculate the zero loss diameter of `tree1.newick`, use

> python RunTests.py -z tree1.newick 2 3 1

To calculate the number of medians, the diameter of the median space, and the distance of the worst median from the furthest MPR of `tree1.newick` and log it to `logs/median.csv`, use

> python RunTests.py -mMw -l "logs/median" tree1.newick 2 3 1

To calculate 100 random medians of the entire TreeLifeData dataset, make a warning sound any time the process fails, and log the results to a folder called `results` use

> python RunTests.py -c 100 -L -l results TreeLifeData/COG####.newick 2 3 1

### Via Interactive Mode:

#### Single Files

To calculate the diameter for a single file in interactive mode, call the following function:
> calculate_diameter_from_file(filename, D, T, L, log=None, debug=False, verbose=True, zero_loss=False, median=False,
                                 worst_median=False, median_cluster=0, median_diameter=False)

The arguments in this function are similar to the command line option `calc`. `filename` is the path to the file that will be reconciled by `DTLReconGraph` and the DTL Diameter calculated for. `D`, `T`, and `L` are the event costs, and `log` is the a string containing path to the csv file used for logging (use `None` for no logging).

The rest of the flags and options are identical to the command line.

#### Many Files

To calculate a set of numbered files, use this function:
> repeatedly_calculate_diameter(file_pattern, start, end, d, t, l, log=None, debug=False, verbose=True, loud=False, zero_loss=False, median_count=False, worst_median=False, cluster=0, median_diameter=False)

Where `file_pattern` is the pattern used to find the right files (as described in the command line `-i` flag section), `start` is the starting file number, `end` is the exclusive ending file number, and the rest of the parameters are the same as `calculate_diameter_from_file()`.

## How to Use DTLMedian.py

### From the Command Line

DTLMedian.py has the following usage pattern:

> DTLMedian.py filename dup_cost transfer_cost loss_cost [-r] [-c]

So, the user must enter the name of the file from which to take the data (in .newick format) as well as the costs associated with a duplication event, transfer event, and loss event, respectively. By default, all valid calls will return the full median reconciliation, however the user may add to the output by including the options shown below, with a newline separating each different value returned.

#### Optional Flags

* -r or --random: adds a random single-path median reconciliation, drawn from the full median reconciliation, to the output

* -c or --count: adds the number of single-path medians that could be derived from the full median reconciliation to the output

If both of these options are selected, the order of the printed output is: the full median reconciliation, the number of medians, a randomly selected median.

### Via Interactive Mode

Although the user could access any and all functions within this file via interactive mode, there are no stand-alone functions that easily integrate basic information, such as a data file name, and can produce information about the median from the get-go. In order to make good use of the functions contained in this file, the user would need to have information from DTLReconGraph.py, such as a DTL reconciliation graph and the roots for the DTL reconciliation graph, since the functions in this file take the raw results in those forms.

## How To Use DataAnalysis.py

### From the Command Line

#### Required Argument:

Running DataAnalysis from the command line has the following usage pattern:

> DataAnalysis.py [options] file

In this command, `file` is the relative path to the csv file (output by Diameter.py) you wish to analyze. By default DataAnalysis will find the minimum, maximum, median, and mean of:
* The number of MPRs
* The diameter found
* The size of the gene tree
* The normalized diameter

and the running time of

* Building the DTL reconciliation graph
* Calculating the diameter
* Both combined

and display them to the screen.

#### Optional Flags
These flags can be used to coax some additional functionality out of DataAnalysis.py, most importantly plots!

`-p` will tell the program to show plots. By default, DataAnalysis will plot the normalized diameter against gene tree size and MPR count.

`-n` tells the program to include another plot that uses the non-normalized diameter.

`-t` tells the program to include another plot that plots the timings of the different portions of Diameter.

`-z` tells the program to also make all plots for a zero loss version of the logfile (it will look for one with the same name but ending in `\_zl`)

`-l` tells the program to use LaTeX for plot text rendering

In addition to those plot options, there are a couple more usable flags:

`-c COMPARE_FILE` compares the given csv file with another one (`COMPARE_FILE`) that has the same calculated files in the same order. It will note any mismatches between the two file's diameters.

`-k CHECK_PATH` compares the given csv file with the path that the newick files you originally reconciled reside (`CHECK_PATH`), and notifies you of any duplicate files in your log, or files in the path that are not in the log.
# medoids
