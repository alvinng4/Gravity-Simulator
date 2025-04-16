The classic Barnes-Hut algorithm \cite{barnes_hierarchical_1986}
provides a way to approximate forces without losing accuracy at close range.
Because gravity decays at a quadratic rate, the accuracy of long range interactions
are less important. Therefore, it is reasonable to approximate a far cluster of 
particles as a single particle with mass \(m = m_{\textnormal{cluster}} \)
and coordinate \(x = x_{\textnormal{com, cluster} } \).
One simple choice of criterion is the opening angle \(\theta = l / d\),
where \(l\) is the length of the cubical cell enclosing the cluster
and \(d\) is the distance between the target particle and the
center of mass of the cluster \mbox{(see figure \ref{fig:barnes_hut}).}
This is purely geometric and does not depends on the 
mass or number of particles in the cluster.

\begin{figure}[!h]
    \begin{subfigure}{0.31\textwidth}
        \centering
        \includegraphics[width=\textwidth]{figures/barnes_hut_a/barnes_hut_a.pdf}
        \caption{Decomposing the spatial domain}
        \label{subfig:barnes_hut_a}
    \end{subfigure}
    \begin{subfigure}{0.31\textwidth}
        \centering
        \includegraphics[width=\textwidth]{figures/barnes_hut_b/barnes_hut_b.pdf}
        \caption{Checking opening angle}
        \label{subfig:barnes_hut_b}
    \end{subfigure}
    \begin{subfigure}{0.31\textwidth}
        \centering
        \includegraphics[width=\textwidth]{figures/barnes_hut_c/barnes_hut_c.pdf}
        \caption{Taking approximation}
        \label{subfig:barnes_hut_c}
    \end{subfigure}
    \caption{Illustration of Barnes-Hut algorithm.}
    \label{fig:barnes_hut}
\end{figure}

At each force evaluation, we perform the following operations:
\begin{itemize}%[leftmargin=*]
    \item[(1)] Construct an octree by decomposing the spatial domain hierarchically.
    \item[(2)] For each particle, traverse the octree from the root into child nodes:
    \begin{itemize}%[leftmargin=*]
        \item[(i)] If the current node has only one particle, evaluate the force directly. Return to parent node.
        \item[(ii)] Otherwise, check the opening criterion. If passed, evaluate the
                    force with the center of mass precomputed in tree construction. Return to parent node.
        \item[(iii)] If not passed, traverse deeper.
    \end{itemize}
\end{itemize}
The technical details for constructing an octree are provided in Appendix \ref{app:tree_construction}.
For a reasonable choice of \(\theta\) \mbox{(e.g. \(\theta \sim 1\))}, the force evaluation of
one particle only takes \(\mathcal{O}(\log N)\) as the expected tree depth to traverse
is \(\mathcal{O}(\log N)\). This gives the time complexity \(\mathcal{O}(N \log N)\)
for \(N\) particles.

Despite the simplicity of the algorithm, it provides an immense speed up over the brute force
\(\mathcal{O}(N^2)\) algorithm. In figure \ref{fig:barnes_hut_benchmark}, we provides a benchmark
of the two algorithms ran on Macbook Air M1 with opening angle \(\theta = 0.5, 1\). As shown from the
result, at the same runtime, Barnes-Hut algorithm could handle \(N \sim 10^6 - 10^7\) particles
while the brute force algorithm could only do \(N \sim 10^5\).
This enables us to run simulations with much larger scales.
As a demonstration, we ran a galaxy collision simulation with 60000 particles using the initial conditions
from Gadget-2 \cite{gadget2}. The result is shown on figure \ref{fig:galaxy_collision}. 
It is done within 30 minutes on Macbook Air M1 with \(\theta = 0.5\).

In modern cosmological simulation code, the tree algorithm
is adopted with one-sided multipole expansion that \linebreak provides higher efficiency and accuracy.
In addition, fast mutlipole algorithm (FMM) \cite{Dehnen_2000} \cite{DEHNEN200227} with dual-sided 
\linebreak multipole expansion provides even higher efficiency
with \(\mathcal{O}(N)\) time complexity.

\begin{figure}[h]
    \centering
    \includegraphics[width=0.5\textwidth]{figures/barnes_hut_benchmark/barnes_hut_benchmark.pdf}
    \caption{
        Benchmarks of direct pairwise computation versus Barnes-Hut algorithm.
        The particles are chosen
        randomly with \((x, y, z) \in [-1, 1)\)
        obtained from a uniform distribution using a random number generator.
        The pairwise algorithm spent 18.6 seconds on \(N = 10^5\), while Barnes-Hut algorithm
        with \(\theta = 1\) and \(\theta = 0.5\) spent \(7.94\) s on \(N = 10^6\) and
        \(18.2\) s on \(N = 10^7\) respectively.     
    }
    \label{fig:barnes_hut_benchmark}
\end{figure}

\begin{figure}[!h]
    \centering
    \includegraphics[width=\textwidth]{figures/galaxy_collision.png}
    \caption{
        Galaxy collision simulation with 20000 disk and 40000 halo particles. 
        The system is evolved for 4 Gyr using the Leapfrog algorithm, 
        with \(\,\Delta t = 2 \textnormal{ Myr}  \) (i.e. 2000 time steps).
        Barnes-Hut algorithm with \(\theta = 0.5\) was used to obtain the 
        gravitational acceleration. Visualizations were done with gadgetviewer.
    }
    \label{fig:galaxy_collision}
\end{figure}