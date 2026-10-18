!!! Note
    In this page, we details the implementation of Barnes-Hut algorithm using linear octree,
    which is useful for parallel computing. If you are only interested at Barnes-Hut algorithm,
    you should skip the parts on linear octree and Morton index.

The classic Barnes-Hut algorithm[@barnes_hierarchical_1986]
provides a way to approximate forces in $\mathcal{O}(N \log N)$ without losing accuracy at close range.
Because gravity decays at a quadratic rate, the accuracy of long range interactions
are less important. Therefore, it is reasonable to approximate a far cluster of 
particles as a single particle with mass \(m = m_{\textnormal{cluster}} \)
and coordinate \(x = x_{\textnormal{com, cluster} } \).
One simple choice of criterion is the opening angle \(\theta = l / d\),
where \(l\) is the length of the cubical cell enclosing the cluster
and \(d\) is the distance between the target particle and the
center of mass of the cluster (see figure 1).
This is purely geometric and does not depends on the 
mass or number of particles in the cluster.

<figure style="text-align: center;">
  <img src="../../../examples/media/illustration_barnes_hut.png" alt="Barnes-Hut algorithm" width="800" style="display: block; margin: auto;" />
  <figcaption>Figure 1: Illustration of Barnes-Hut algorithm.</figcaption>
</figure>

## Linear octree construction

??? Note "Source code (Click to expand)"
    ```c linenums="1"
    --8<-- "src/linear_octree.c"
    ```

<figure style="text-align: center;">
  <img src="../../../examples/media/morton_curve.png" alt="Morton curve" width="400" style="display: block; margin: auto;" />
</figure>
Figure 2: Morton curve for Morton index 0 to 7. This is the full curve for level 1,
        or the first \(\frac{1}{8}\) of the full curve at level 2, etc.

Here, we provide a detailed description on constructing a static linear octree
using only linear arrays and Morton indices.
To build a linear octree, we utilize the idea of Morton code (also known as Z-order
or Morton space filling curve). Figure 2 shows a Morton curve at level 1
of the tree. The Morton index is calculated by encoding the spatial coordinate in binary
format. For example, at level 1 we have \(x, y, z \in \{0, 1\}\). For \(x = 0, y = 0, z = 1\),
we have the Morton index \(\underset{z}{1}\underset{y}{0}\underset{x}{0} \textnormal{ (binary)} = 4\). For \(x = 1, y = 1, z = 0\),
we have the Morton index \(\underset{z}{0}\underset{y}{1}\underset{x}{1} \textnormal{ (binary)} = 3\).
As we traverse into deeper level, we stack another three bits after the Morton index for each level. 
For example, a particle has a local Morton index 5 at level 1 and local Morton index 3
at level 2. The full Morton index at level 2 is obtained by

\begin{equation}
    %\label{}
    \overset{5}{\overbrace{\underset{z}{1}\underset{y}{0}\underset{x}{1}}}
    \overset{3}{\overbrace{\underset{z}{0}\underset{y}{1}\underset{x}{1}}}
    \textnormal{ (binary)}
    = 43.
\end{equation}

Unlike a tree data structure, there is a limit for the depth of the octree. For 64-bit
integer, there can only be \(\lfloor64 / 3 \rfloor = 21\) levels (excluding level 0) 
as it takes three bits for each level. This could become an issue for exascale simulations,
but this could be resolved by using integers with more bits. 

The Morton index could be calculated easily using bit-shift operations and loops.
In our project, we fixed the Morton indices to 64-bit integers by default and uses
magic numbers to compute the Morton indices at level 21 directly without a loop.
The magic numbers are generated using the script by [@morton_index_stackoverflow].
Morton index on each level can then be retrieved with bit-shift operations.

<figure style="text-align: center;">
  <img src="../../../examples/media/linear_octree.png" alt="Linear octree" width="800" style="display: block; margin: auto;" />
</figure>
Figure 3: Graphical illustration of linear octree. On the top, there are multiple aligned
        arrays. Each index represent one tree node, and each array represent a piece
        of information stored by the tree node. On below, we have the sorted Morton
        indices and the particle indices sorted with Morton indices.
        Since they are sorted, only the first index and the number of particles are needed
        to obtain the full particle list of each tree node. Tree node 0 is the
        root node with $N$ particles. Tree node 1 and 2 are the proper successor of the root node,
        with 4 and 8 particles respectively.

Now, with the knowledge of Morton index, we can construct a linear octree building algorithm.
The tree is represented with multiple aligned arrays, where each index to the arrays
corresponds to one internal node. An illustration of the linear octree is provided in figure 3.

1. Compute the Morton index for all particles at the deepest
    level (level 21 for 64-bit integers)
2. Sort the Morton index (e.g. radix sort) along with the particle indices, so that
    we have an array of sorted Morton indices, and particle indices
    that corresponds to each Morton index.
3. For each particle, starts from the root node,
    - Check if there are any particles in the corresponding suboctant of the current node.
        This can be done with binary search of the suboctant's Morton index at that level.
        (The binary search also tells us how many particles are in each child node.)
        
        - If not, instantiate a new child node for that suboctant by assigning a tree index.
            Backpropogate the mass and coordinate all the way to the root node.
        - Otherwise, traverse deeper into the child node.

Side note: We do not know beforehand how many internal nodes there will be.
Therefore, the arrays might become full during construction. By using a dynamic array
(one that doubles in size whenever it is full), we can build an octree with as many internal
nodes as needed.

## Tree traversal

??? Note "Source code (Click to expand)"
    ```c linenums="1"
    --8<-- "src/acceleration_barnes_hut.c"
    ```

To compute the acceleration, we use a stack to traverse the tree recursively:

1. Push the root node to the stack.
2. While the stack is not empty:
    - 2.1 Pop the top node from the stack.
    - 2.2 If the node is a leaf node, compute the acceleration directly.
    - 2.3 If the node is not a leaf node, check if it passes the opening angle criterion.
        If yes, compute the acceleration using the center of mass of the node.
        Otherwise, push all child nodes to the stack.

!!! danger
    We need to check whether the current node is included in the branch to avoid
    self-interaction. This can be done by a simple check of the Morton index.

## Benchmark

Here we provide a simple benchmark of the Barnes-Hut algorithm compared to the brute-force 
algorithm. It was done on Macbook Air M1. The brute-force algorithm spent 18.6 seconds 
on \(N = 10^5\), while the Barnes-Hut algorithm
with \(\theta = 1\) and \(\theta = 0.5\) spent \(7.94\) s on \(N = 10^6\) and
\(18.2\) s on \(N = 10^7\) respectively. This shows that Barnes-Hut algorithm could handle 10-100 times more particles
than the brute force algorithm.

(For the benchmark, the particles are uniformly distributed where $\{x, y, z\} \sim U(-1, 1)$)

<figure style="text-align: center;">
  <img src="../../../examples/media/barnes_hut_benchmark.png" alt="Baranes-Hut benchmark" width="600" style="display: block; margin: auto;" />
</figure>
