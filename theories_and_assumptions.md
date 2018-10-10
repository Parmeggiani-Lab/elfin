# elfin theories and assumptions

[D] is a definition

[T] is a given truths or unproven theories derived from given truths

[A] is an assumption we use to bound the problem that elfin attempts to solve

[S] is a speculation

## Module chaacteristics
 - [D] A "module" is a protein that has two or more uncapped termini (interfaces)
 - [T] Modules are slightly flexible
 - [A] Each module "pose" can be treated as a rigid snapshot
 - [T] In general modules are of different sizes
 - [T] In general modules have different termini transformations (delta transformation between two termini)
 - [T] Termini are either type N or type C
 - [T] Modules maybe single- or multi-chain
 - [D] A "cap" is an one-terminus protein that can close up a dangling interface
 - [D] A "module sequence" consists of multiple modules that are connected to one another via valid interfaces
 - [D] A "hub" is a multi-chain module that can join different module sequences
 - [D] A "symmetric" hub must have identical module sequences stemming from each of its termini
 - [D] An "asymmetric" hub is a hub without the symmetric hub restriction


## Algorithm characteristics
 - [T] Module interface connectivity can be modelled as an incomplete directed cyclic graph with an uneven distribution of node degrees
 - [D] A "module sequence" is a valid walk in the connectivity graph
 - [D] A "network" is a collection of module sequences joined by hubs
 - [A] Two networks can only be separate networks if they cannot or should not be joined together
 - [A] An elfin design output is a collection of networks
 - [A] Elfin attempts to solve a general shape optimization problem: build a structure that fits a desired shape defined as a 3D guiding graph using smaller modules as components
 - [T] The user does not know beforehand the optimal module choice. It follows that elfin must insert or remove nodes where necessary so modules can actually connect together. Hence a node in the guiding graph does not necessarily correspond to one module in the solution
 - [T] Since elfin may need to insert or remove nodes, the solution size becomes another search parameter
 - [T] Since modules may not be able to satify certain shapes, there may or may not exist a perfect solution
 - [T] In the extremely rare special case where there is an exact optimal solution, it can be checked in O(n) time if the solution is size n (e.g. an artificial case where each node in the guiding path sits exactly on their corresponding optimal module COM) 
 - [T] In the general case, given one input guiding graph and one solution, checking optimality boils down to searching the solution space to see if a better solution exists
 - [T] Searching for the optimally fitting module sequence in the connectivity graph is a special case of the Precedence Constrained [Knapsack Problem](https://en.wikipedia.org/wiki/Knapsack_problem) (PCKP).
 - [T] In addition to PCKP, elfin's problem is unbounded (nodes are repeatable) and immediate-precendence (i.e. last picked node has an arc to the current choice, as opposed to any-one-precedence and all-precendence).
 - [S] Could not confirm complexity lower bound 
 - [T] Solution space grows exponentially with solution size (from my thesis)
 - [T] Many Knapsack problems can be solved with dynamic programming, but that relies on the problem being a 0-1 Knapsack problem. 
 - [S] Elfin is not a 0-1 Knapsack problem because each solution might visit a node multiple times, and the order of visits matter. These seem to be a hint that dynamic programming cannot help elfin solve its problem

## Beyond the general guiding graph
On top of the using a guiding graph (in elfin-ui terms, this would be a pure path guide network) as input, we want to allow the user to manually specify parts of the solution. This comes in the following forms:
 1. Placing specific modules in 3D space
 2. Connecting specific modules to path guides
 3. Specify whether/how much a specific module can translate/rotate in case a more optimal solution can be generated

Numbers 1 and 2 help reduce the search space by reducing the solution size and also by restricting interface choices. However, number 3 <em>adds</em> complexity. If user-specified modules can move/rotate around, I can think of two ways that elfin can effectively use this information:
 1. During candidate evaluation, take into account the spatial tolerance by Kabsch'ing only on nodes that have no spatial tolerance (tolerance needs to be "floored" in preprocessing). Because elfin's scoring function doesn't care about rotation, we need to add rotation tolerance check in this case
 2. Elfin generates a discrete set of possible poses for each spatially tolerant module and tries to solve each case and find the optimal between the poses

 Method 1 seems more straightforward, but might need to violate user-defined tolerances. Method 2 gurantees that user-defined tolerances are observed, but is compute intensive and feels ugly.