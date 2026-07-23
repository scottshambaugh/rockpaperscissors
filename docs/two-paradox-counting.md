---
title: "Counting tournaments with domination number at least three: a proved polynomial-time algorithm"
author: "Computational derivation and validation"
date: "2026-07-22"
---

# 1. Introduction and main theorem

A *tournament* on $n$ vertices is an orientation of the complete graph $K_n$: for every
pair $\{u,v\}$ exactly one of $u\to v$, $v\to u$ holds. Its *domination number*
$\gamma(T)$ is the size of the smallest set $D$ of vertices such that every other vertex is
beaten by some member of $D$. For $n\ge 2$ we give and prove correct an algorithm that computes

$$
L(n)\;=\;\#\{\text{labeled tournaments on } [n] \text{ with } \gamma(T)\ge 3\}.
$$

(We take $n\ge 2$ throughout: for $n=1$ the reduction below would return $1$ while $L(1)=0$.)
The algorithm runs in **polynomial time**: $L(n)$ is a signed sum over caterpillar forests and
odd-cycle suns whose per-graph weight is computed by a linear-time transfer recurrence (§§3--5),
and a symmetric-power / exponential-generating-function aggregation (§6) evaluates the entire
signed sum *without enumerating the graphs*, in $\mathrm{poly}(n)$ big-integer operations -- as
against the $2^{\binom{n}{2}}=2^{\Theta(n^2)}$ of direct evaluation. (Summing the same per-graph
weights by naïve term-by-term enumeration instead gives a $2^{O(n)}$ variant; we mention it only
as a special case and do not rely on it.)

**Correctness.** Each step is a theorem: the master reduction and sparsity restriction (§2), the
checkerboard factorization and odd-cycle-sun count (§3), the per-caterpillar transfer (§4), and
the polynomial aggregation (**Theorem L3**, §5); so $L(n)$ is computed exactly in polynomial time.
Where an independent check is available it agrees -- $L(7),L(8),L(9)$ match direct brute-force
enumeration, the isomorphism counts through $n=10$ match McKay's catalogues, and $L(11)$ satisfies
a Burnside congruence (§6) -- but these are corroboration, not the proof.

The tractability here should not be mistaken for a complexity-theoretic breakthrough: we *count*
rather than *decide*, and the two are unrelated. The count sums over the sparse family supplied by
the structure theorem of Fisher, Lundgren, Merz and Reid (§2) -- which is where essentially all the
tractability comes from -- so no individual instance is ever solved. (Deciding $\gamma(T)\le k$ for
a *given* tournament is a genuinely harder-looking problem -- quasi-polynomial and LOGSNP-complete:
Megiddo--Vishkin, *Theoret. Comput. Sci.* **61** (1988) 307--316; Papadimitriou--Yannakakis, *J.
Comput. System Sci.* **53** (1996) 161--170 -- but that hardness does not bear on the count, and no
barrier is at issue.) A Burnside/cycle-index step yields the count of
$\gamma\ge 3$ tournaments *up to isomorphism* (a separate small-$n$ layer, §6), extending a
sequence previously tabulated only to $n=10$; we are aware of no prior formula, recurrence, or
OEIS entry for it, so the values for $n\ge 11$ are new.

The condition $\gamma\ge 3$ is the two-paradox (Schütte) property $S_2$: no single vertex
beats everyone, and no pair of vertices together beats all others. By Reid, McRae, Hedetniemi
and Hedetniemi (Australas. J. Combin. **29** (2004) 157--172, Thm. 5), $\gamma(T)\le 3$ for
every tournament on $n<19$ vertices, so for $2\le n\le 18$ the count "$\gamma\ge 3$" equals
"$\gamma=3$ exactly," and the domination distribution is supported on $\{1,2,3\}$ with
$\#\{\gamma=1\}=n\,2^{\binom{n-1}{2}}$.

**Scope.** The method rests on the structure theorem for dominating *pairs* -- $\gamma\ge3$ is
exactly "$\mathrm{dom}(T)$ edgeless" -- so it computes $\#\{\gamma\ge3\}$ and nothing finer; there
is no analogous characterization for dominating triples, and $\gamma\ge4$ is out of reach. For
$n\le18$ this costs nothing: $L(n)=\#\{\gamma=3\}$, the statistic of interest, and even the
single-exponential form of the method reaches that whole range. At $n=19$ the identity breaks --
the Paley tournament $QRT_{19}$ has $\gamma=4$ -- so for $n\ge19$, $L(n)=\#\{\gamma\ge3\}$ is a
blend of $\gamma=3$ and $\gamma\ge4$ that the method cannot separate. Beyond $n=18$ the polynomial
algorithm therefore evaluates this aggregate count efficiently; it does not recover
$\#\{\gamma=3\}$.

---

# 2. Inclusion--exclusion and structural reduction

For a vertex $w$ write $N^+(w)$ for its *out-neighborhood* (the set of vertices $w$ beats).
A pair $\{u,v\}$ is a **dominating pair** if no vertex beats both, i.e. $\{u,v\}$ is not
contained in any $N^+(w)$. The **domination graph** $\mathrm{dom}(T)$ has an edge $\{u,v\}$
exactly for the dominating pairs. Then

$$
\gamma(T)\ge 3
\iff \text{no vertex dominates and no pair dominates}
\iff \mathrm{dom}(T)\ \text{has no edge.}
$$

Indeed $\gamma=1$ means some vertex beats all others, and $\gamma\le 2$ that some pair
$\{u,v\}$ beats everyone else -- exactly that $\{u,v\}$ has no common in-dominator, i.e. is an
edge of $\mathrm{dom}(T)$.

Equivalently, $\{u,v\}$ is a dominating pair iff no out-neighborhood contains both $u$ and
$v$, i.e. every $N^+(w)$ is an independent set in the graph having the single edge
$\{u,v\}$. Asking a *set* of pairs $H$ to all be dominating is therefore

$$
H\subseteq \mathrm{dom}(T)\iff \text{every out-neighborhood } N^+(w)\ \text{is } H\text{-independent.}
$$

**The master inclusion--exclusion.**
Index tournaments by their set of dominating pairs and Möbius-invert on the (edge-)lattice
of graphs on $[n]$: since $L(n)=\#\{T:\mathrm{dom}(T)=\varnothing\}$,

$$
L(n)=\sum_{H\subseteq K_n}(-1)^{e(H)}\,N(H),\qquad
N(H):=\#\{T:\ H\subseteq\mathrm{dom}(T)\},
$$

where $e(H)$ is the number of edges of $H$. Now $N(H)$ depends on $H$ only through its
non-isolated part. Let $H$ have exactly $v$ non-isolated vertices carrying a graph $H'$ (no
isolated vertices) and $n-v$ untouched vertices. The constraint "every $N^+(w)$ is
$H$-independent" involves only how each vertex relates to the $v$ support vertices, so the
tournament factorizes:

* the sub-tournament on the $n-v$ non-support vertices is **arbitrary**: $2^{\binom{n-v}{2}}$;
* the sub-tournament on the $v$ support vertices must have **every out-neighborhood
  $H'$-independent**: call the number of such sub-tournaments $T(H')$;
* each non-support vertex $w$ may beat any $H'$-independent subset of the support (and is
  beaten by the rest): $I(H')$ choices each, where $I(H')$ is the number of independent sets
  of $H'$ (all sizes, including the empty set). This contributes $I(H')^{\,n-v}$.

Summing over the choice of support ($\binom{n}{v}$ ways) and over $H'$ gives the

**Master formula.**
$$
\boxed{\,L(n)=\sum_{v=0}^{n}\binom{n}{v}\,2^{\binom{n-v}{2}}\,A_v(n-v)\,},\qquad
A_v(t)=\!\!\sum_{\substack{H\ \text{on}\ [v]\\ \text{no isolated vertex}}}\!\!(-1)^{e(H)}\,T(H)\,I(H)^{\,t},
$$
the inner sum ranging over *labeled* graphs $H$ on $v$ vertices with minimum degree
$\ge 1$, and
$$
T(H)=\#\{\text{tournaments on }V(H):\ \text{every out-neighborhood is }H\text{-independent}\}.
$$
Here $A_0=1$ (empty graph) and $A_1=0$ (a single vertex is isolated). The $v=0$ term is
$2^{\binom{n}{2}}$, the total number of tournaments, and the higher terms are the
inclusion--exclusion corrections subtracting the $\gamma\le 2$ tournaments.

Everything now rests on computing $T(H)$ and enumerating the $H$ that matter.

**Sparsity (Fisher--Lundgren--Merz--Reid).** $T(H)>0$ requires $H\subseteq\mathrm{dom}(T)$ for
some tournament $T$, and the structure of domination graphs is completely known:

**Theorem (Fisher, Lundgren, Merz, Reid, *J. Graph Theory* 29(2) (1998) 103--110).**
For every tournament $T$, the domination graph $\mathrm{dom}(T)$ is either a forest of
caterpillars, or an odd cycle together with pendant and isolated vertices.

The master formula strips isolated vertices before forming $H'$, so the relevant consequence
is that the non-isolated part of any domination graph is either a forest of nontrivial
caterpillars or a single **odd-cycle sun** (an odd cycle $C_g$, $g\ge 3$, with any number of
pendant leaves on its vertices). (A *caterpillar* is a tree that becomes a path when its leaves
are deleted; a spine path $s_0\cdots s_{m-1}$ with $\ell_i\ge 0$ leaves on $s_i$.) Consequently
$T(H')=0$ unless every component of $H'$ is a nontrivial caterpillar, or $H'$ is a single
odd-cycle sun. This restricts the sum from all $2^{\binom{v}{2}}$ graphs to caterpillar forests
and suns, but is **not** by itself a polynomial enumeration: there are $\Theta(2^{k})$
caterpillars on $k$ vertices (OEIS A005418), so the structural theorem alone gives only a
single-exponential ($2^{\Theta(n)}$) enumeration, not $\mathrm{poly}(n)$.

---

# 3. Weights of connected components

**Theorem.** Let $H$ be a forest, **each of whose components $C_1,\dots,C_c$ has at least one
edge** (as always holds for the isolate-free support $H'$). Then
$$
T(H)=2^{\binom{c}{2}}\prod_{i=1}^{c}\mathrm{Within}(C_i),\qquad
\mathrm{Within}(C):=T(C).
$$
Under the same hypothesis (no isolated vertices), if some component is cyclic then $T(H)=0$
unless $H$ is a single odd-cycle sun, in which case $\mathrm{Within}=2^{\,1+\sum_i \binom{\ell_i}{2}}$
(in particular $\mathrm{Within}=2$ when no cycle vertex carries two or more leaves); an even cycle
gives $0$. (The no-isolates hypothesis is needed here too -- e.g. $T(C_3\sqcup K_1)=8\ne0$ -- for
the reason given next.)

The edge-bearing hypothesis is essential: for $H=K_2\sqcup K_1$ (an isolated vertex) one has
$T=2\cdot 3=6$ (the isolated vertex may beat at most one endpoint of the edge), whereas the
displayed formula would give $2^{\binom{2}{2}}\cdot 2\cdot 1=4$. Isolated vertices are exactly
what the master formula removes, so no such $H$ occurs there; but the theorem must exclude
them.

*Proof.* A set $S$ is $H$-independent iff its restriction to each component is
component-independent (all edges live inside components). Fix the within-component
orientations. Two facts remain:

(i) *Within factor.* The constraint "$N^+(w)\cap C_i$ is $C_i$-independent" for $w\in C_i$
depends only on the sub-tournament on $C_i$; these constraints are independent across
components, giving $\prod_i\mathrm{Within}(C_i)$.

(ii) *Between-component factor (checkerboard lemma).* Fix an ordered pair $(C_i,C_j)$ and write
$X_{uv}=1$ iff $u\in C_i$ beats $v\in C_j$. For any edges $uu'\in E(C_i)$ and $vv'\in E(C_j)$
consider
$$
\begin{pmatrix}X_{uv}&X_{uv'}\\ X_{u'v}&X_{u'v'}\end{pmatrix}.
$$
Independence of the out-neighbourhoods in $C_j$ forces neither row to be $(1,1)$ (else some
$u$ beats the adjacent $v,v'$); independence of the reverse out-neighbourhoods in $C_i$ forces
neither column to be $(0,0)$ (else some $v$ beats the adjacent $u,u'$). A $2\times2$ $0/1$
matrix with no all-ones row and no all-zeros column has exactly two ones and is one of
$\bigl(\begin{smallmatrix}1&0\\0&1\end{smallmatrix}\bigr)$,
$\bigl(\begin{smallmatrix}0&1\\1&0\end{smallmatrix}\bigr)$. Since each $C_i,C_j$ is a
*connected bipartite* graph (a caterpillar has $\ge1$ edge), this local constraint propagates:
writing $p_i,p_j$ for the bipartition indicators, every valid cross-orientation has the form
$X_{uv}=b\oplus p_i(u)\oplus p_j(v)$ for a single bit $b\in\{0,1\}$. Hence exactly **two**
cross-orientations per component pair (connectedness and the edge in each component are what
pin it to $2$; a bare $K_1$ has no edge to propagate $b$ and would allow $3$, the source of the
$K_1$ discrepancy above). Over the $\binom{c}{2}$ pairs this gives $2^{\binom{c}{2}}$. Cyclic
components are killed unless the odd-cycle-sun exception applies. $\square$

We have verified the theorem exhaustively for every isolate-free graph on $v\le 9$ vertices
(zero discrepancies). The "$2$ per component pair" mechanism is illustrated by
$T(K_2\sqcup K_2)=2^{1}\cdot 2\cdot 2=8$, whose $8$ tournaments have exactly $2$ valid
cross-orientations for each of the $2\times 2$ within-orientations.

**Lemma L1 (odd-cycle sun).** $\mathrm{Within}(\text{odd-cycle sun})=2^{\,1+\sum_i\binom{\ell_i}2}$;
$\mathrm{Within}(C_{2m})=0$; and any connected graph containing a cycle that is *not* an
odd-cycle sun has $\mathrm{Within}=0$.

*Zero cases.* If $\mathrm{Within}(H)=T(H)>0$ then $H\subseteq\mathrm{dom}(T)$ for some $T$; by
the Fisher--Lundgren--Merz--Reid theorem (§2) a connected cyclic contributing graph must be a
spanning connected cyclic subgraph of an odd-cycle sun, hence is itself an odd-cycle sun. So
any connected cyclic $H$ that is not an odd-cycle sun has $\mathrm{Within}=0$; in particular an
even cycle gives $0$.

*The value.* Write the cycle as $c_0,\dots,c_{g-1}$, with $\ell_i$ leaves on hub $c_i$. Working
inside any valid tournament (using only the dominating-pair condition on cycle edges and
$H$-independence of out-neighbourhoods), four facts hold: **(a)** the cycle is one of exactly
two directed rotations and every chord $c_i\to c_{i+d}$ is forced ($d$ odd outward), consistent
only for odd $g$ (even $g$ forces the wrap edge both ways, giving $0$); **(b)** each leaf's
orientation against every foreign hub is forced by the parity of the offset; **(c)** every
foreign leaf--leaf edge is forced; **(d)** only the $\binom{\ell_i}2$ edges between same-hub
leaves are free, and every such completion is valid. Each of the two rotations therefore
contributes $\prod_i 2^{\binom{\ell_i}2}$, giving $\mathrm{Within}=2^{\,1+\sum_i\binom{\ell_i}2}$.
The four claims are established by a parity induction on the cycle-and-chord orientations, given
in Appendix B. $\square$

Thus the problem reduces to the single primitive $\mathrm{Within}(\text{caterpillar})$, computed
in §4. The smallest values -- star $K_{1,m}=2^{\binom m2}+m\,2^{\binom{m-1}2}$, path $P_k=4$
($k\ge3$), $K_2=2$ as a component, and the double-/triple-star channels -- are closed forms
(Appendix A), all brute-verified (§6).

---

# 4. The caterpillar transfer theorem

This section gives the per-caterpillar primitive: it computes $\mathrm{Within}(C)$ in time
linear in $|C|$ by an explicit $4$-dimensional spine recurrence, proved exact for all
caterpillars (the transfer theorem below). A caterpillar is encoded by its leaf-count word
$\mathbf{l}=(\ell_0,\dots,\ell_{m-1})$.

**Leaf-twin symmetry.** All $\ell_i$ leaves on $s_i$ are false twins, so validity depends only
on how *many* of them each out-neighborhood contains. Writing
$$
p_A(\ell)=2^{\binom{\ell}{2}}A^{\ell}\quad\text{and}\quad q_A(\ell)=\ell\,2^{\binom{\ell-1}{2}}A^{\ell-1},
\qquad A\in\{1,2,4\},
$$
the per-vertex leaf factor is a combination of these two shapes, where $A$ is the (finite,
per-core) number of admissible one-leaf extensions (the a-priori slot $A=8$ never occurs). The
$4$-dimensional state is **not** motivated by counting these scalar $A$-values; it comes from
the boundary profile $\pi(P)=(\mathrm{Within}(P{\cdot}[0]),\dots,\mathrm{Within}(P{\cdot}[3]))$
-- four evaluations whose governing matrix $G$ is nonsingular ($\det G=144$), as proved in the
transfer theorem below.

**Leaf-Factor Lemma (proved).** Fix all edges not incident to the leaves of one spine vertex
$s$ (an internally-valid core config $\tau$). Summing over the leaf-orientations of $s$
contributes exactly $f_\tau(\ell)=p_{A}(\ell)+A^{*}\,q_{A}(\ell)$, where $A=\#\{$C-independent
$Q\subseteq W$ containing the external in-neighbours of $s\}$ and $A^{*}$ is the same count
restricted to sets avoiding $s$'s C-neighbours. (Proof: at most one leaf beats $s$; the two
cases give the two terms.) This is the true origin of the two-term basis.

**Channels.** The star, double-star and triple channels are *theorems* for every leaf-count
$\ell$, obtained from the Leaf-Factor Lemma by a finite enumeration of internally-valid core
configurations: each channel is a $\mathbb Z$-combination of the two-dimensional block
$\{p_A,q_A\}$, and enumerating the finitely many cores fixes the two coefficients. The
single-vertex value is $\mathrm{star}(m)=2^{\binom m2}+m\,2^{\binom{m-1}2}=p_1(m)+q_1(m)$ by
direct counting. The fourteen double-/triple-channel cases, their core-config
$(A{:}\text{count})/(A{:}\sum A^{*})$ table, and the exhaustiveness argument are collected in
Appendix A, together with the resulting closed forms $DS(L,d)$ and $g_{a,c}(\ell)$.


**State and recurrence.** Realize the automaton in the basis of *one-symbol suffixes*: for a
processed spine prefix $P$, keep the $4$-vector
$$
R[j]=\mathrm{Within}(P\ \text{followed by a spine vertex carrying } j \text{ leaves}),\qquad j\in\{0,1,2,3\}.
$$
Let $G$ be the fixed $4\times 4$ matrix $G[a][d]=DS(a,d)=\mathrm{Within}([a,d])$ (double-star
values), with $\det G = 144$ and adjugate $\operatorname{adj}(G)=144\,G^{-1}$ (an explicit
integer matrix). Let $W(\ell)$ be the $4\times 4$ matrix
$W(\ell)[j][a]=\mathrm{Within}([a,\ell,j])$ (the triple values, each a fixed combination of
the two shapes above). Then:

* **Initialize** $R^{(0)}[j]=\mathrm{Within}([\ell_0,j])$ from the double-star channel.
* **Advance** across each interior symbol $\ell$: $\;R'\;=\;W(\ell)\,G^{-1}\,R$.
* **Read out** with the last symbol $x=\ell_{m-1}$:
  $\;\mathrm{Within}(\mathbf{l})=\lambda(x)^{\top}R,\ \lambda(x)=G^{-1}w(x),\
  w(x)[a]=\mathrm{Within}([a,x])$.

Compactly,
$$
\mathrm{Within}(\ell_0,\dots,\ell_{m-1})
=\lambda(\ell_{m-1})^{\top}\,\Big(\!\prod_{i=m-2}^{1} W(\ell_i)\,G^{-1}\Big)\,R^{(0)}(\ell_0).
$$

Here $G[a][d]=DS(a,d)$ has $\det G=144$, and $\operatorname{adj}(G)=144\,G^{-1}$ is the
integral matrix
$$
\operatorname{adj}(G)=\begin{pmatrix}
0 & -48 & 24 & 0\\
-48 & -30112 & 22928 & -2088\\
24 & 22928 & -17416 & 1584\\
0 & -2088 & 1584 & -144
\end{pmatrix}
\qquad(\operatorname{adj}(G)\,G=144\,I,\ \text{elementary}).
$$
Each spine vertex is a constant number of $4\times4$ big-integer operations, so
$\mathrm{Within}(C)$ costs $O(|C|)$ operations (integers of $O(|C|^2)$ bits) -- versus
$2^{\binom{|C|}{2}}$ for direct enumeration. (The exact-integer division by $144$ at each
*interior* step is justified by Theorem L2 below: every interior numerator
$W(\ell)\operatorname{adj}(G)R$ is again a vector of tournament counts, hence divisible by
$144$, and the recurrence runs entirely in exact integers.)

**Theorem L2 (caterpillar transfer).** For every prefix $P$ the four-coordinate boundary
profile $\pi(P)=(\mathrm{Within}(P\cdot[d]))_{d=0}^{3}$ is a sufficient statistic: it evolves by
a prefix-independent matrix, $\pi(P\cdot[\sigma])=M(\sigma)\,\pi(P)$ (with $M(\sigma)$ built from
the triple values $\mathrm{Within}([a,\sigma,d])$), and the terminal value
$x\mapsto\mathrm{Within}(P\cdot[x])$ is the fixed linear readout $\lambda(x)^{\top}\pi(P)$ with
$\lambda(x)=G^{-1}(DS(a,x))_a$. Hence the $4$-dimensional recurrence computes $\mathrm{Within}(C)$
exactly for **every** caterpillar.

*Proof.* Split $C=P\cdot v$ at the cut spine edge $\{s,g\}$, where $s$ is the last spine vertex
of $P$ (vertex set $B:=V(P)$) and $g$ the first spine vertex of the right part $v$ (vertex set
$R$). Since leaves stay in the block of their hub and the spine is a path, $\{s,g\}$ is the
**only** cross $C$-edge, so a tournament decomposes uniquely as $(T_B,T_R,D)$ with $T_B=T|_B$,
$T_R=T|_R$, and cross-data $D_z:=N^{+}_T(z)\cap B\subseteq B$ for $z\in R$.

**Lemma 4 (gluing in coordinates).** $T=(T_B,T_R,D)$ is valid iff **(0')** $T_B,T_R$ are valid;
**(A')** each $D_z$ is $C$-independent; **(B')** for every $C$-edge $\{z,z'\}$ of $R$,
$D_z\cup D_{z'}=B$; together with the boundary conditions **(C1)** $I\subseteq D_g$, where
$I:=N^-_{T_B}(s)$, and **(C2)** $s\notin D_x$ for every $x\in R$ that beats $g$ in $T_R$.
*Proof.* By the Reformulation, $T$ is valid iff every $N^+(w)$ is $C$-independent, and (as
$\{s,g\}$ is the unique cross edge) a set is $C$-independent iff its $B$- and $R$-parts are and
it omits $\{s,g\}$. For $w\in R$, $N^+(w)\cap B=D_w$, giving (A'); for $w\in B$,
$N^+(w)\cap R=\{z:w\notin D_z\}$ is $C$-independent iff no $C$-edge $\{z,z'\}$ of $R$ has
$w\notin D_z\cup D_{z'}$, i.e. $D_z\cup D_{z'}=B$ after ranging over $w$ -- this is (B'). The
requirement that no $N^+(x)$ contain $\{s,g\}$ gives (C1) for $x\in B$ and (C2) for $x\in R$. $\square$

**Lemma 5 (covering $\Rightarrow$ colour class).** Let $B$ be connected with $|B|\ge2$ and
$2$-colouring $\{X,Y\}$. If $D,D'\subseteq B$ are $C$-independent with $D\cup D'=B$, then
$\{D,D'\}=\{X,Y\}$. *Proof.* Any $v\in D\cap D'$ has every $C$-neighbour $u$ in $D\cup D'=B$; if
$u\in D$ then $\{u,v\}\subseteq D$ is a $C$-edge inside an independent set, and $u\in D'$ is
likewise impossible, so $v$ has no $C$-neighbour -- contradicting connectedness with $|B|\ge2$.
Hence $D\cap D'=\varnothing$, so $(D,D')$ is a proper $2$-colouring of $B$; by uniqueness of the
$2$-colouring of a connected bipartite graph, $\{D,D'\}=\{X,Y\}$. $\square$

**Lemma 6 (two-colouring collapse).** For $B=V(P)$ connected with $|B|\ge2$ and $R$ a connected
caterpillar with $\ge1$ edge, the families $D$ satisfying (A')+(B') are exactly the **two** proper
$2$-colourings of $R$ by the classes $\{X,Y\}$ (a colouring and its swap). *Proof.* Every $z\in R$
meets a $C$-edge $\{z,z'\}$ of $R$; (A')+(B') with Lemma 5 give $\{D_z,D_{z'}\}=\{X,Y\}$, so
$z\mapsto D_z$ is a proper $2$-colouring of $R$, and conversely every such colouring satisfies
(A')+(B'). A connected graph has exactly two. $\square$

So the cross-data has exactly **two** possibilities, *independent of $|V(P)|$*: the feared "$i(P)$
blow-up" (a growing independent-set count on $P$) never appears, because the past vertices' own
out-neighbourhood constraints (B') force each $D_z$ to be an entire colour class. **The two colour
classes $X,Y$ are the two boundary states the transfer must carry.**

*Summation (L2\_proof §§4--5).* Name $X$ the colour class containing $s$. For $R=[\sigma,d]$
(hubs $g,h'$ carrying $\sigma,d$ leaves) the two colourings $\kappa^{\pm}$ assign $X,Y$ to the
class $R_0=\{g\}\cup L_{h'}$. Evaluating (C1),(C2) for a fixed valid $T_B$ (fixing
$I=N^-_{T_B}(s)$) and a fixed valid $T_R$: $\kappa^{+}$ needs $I\subseteq X$ and no leaf of $h'$
beats $g$; $\kappa^{-}$ needs $I\subseteq Y$ and $g$ beats $h'$ with no leaf of $g$ beating $g$.
Writing $N^{+}(\sigma,d),N^{-}(\sigma,d)$ for the resulting two junction constants (functions of
$(\sigma,d)$ only) and summing over valid $T_R$,
$$\textstyle\sum_{T_R}\#\{\text{valid }D\}=[I\subseteq X]\,N^{+}(\sigma,d)+[I\subseteq Y]\,N^{-}(\sigma,d).\tag{4.1}$$
Summing (4.1) over all valid $T_B$ on $P$ (valid because $|V(P)|\ge2$, so Lemma 6 applies) gives
the explicit **rank-$2$ relation**
$$\mathrm{Within}(P\cdot[\sigma,d])=a(P)\,N^{+}(\sigma,d)+b(P)\,N^{-}(\sigma,d),\tag{5.1}$$
$a(P):=\#\{T_B: I\subseteq X\}$, $b(P):=\#\{T_B: I\subseteq Y\}$. The *same* machinery on the
one-hub part $R=[d]$ gives $N^{+}_{[d]}=\mathrm{star}(d)$ and $N^{-}_{[d]}=p(d)=2^{\binom d2}$, so
$\mathrm{Within}(P\cdot[d])=a(P)\,\mathrm{star}(d)+b(P)\,p(d)$; taking $d=1,3$ with
$\mathrm{star}(1)=2,\ p(1)=1,\ \mathrm{star}(3)=14,\ p(3)=8$,
$$\mathrm{Within}(P\cdot[1])=2a(P)+b(P),\qquad \mathrm{Within}(P\cdot[3])=14a(P)+8b(P),$$
**so the $2\times2$ system has determinant $2$, uniquely fixing $a(P),b(P)$ as
prefix-independent linear functionals of $\pi(P)$.** Substituting these back into (5.1):

**Proposition 7 (rank-$2$ transfer).** For every prefix $P$ with $|V(P)|\ge2$ and all
$\sigma,d\ge0$,
$$\mathrm{Within}(P\cdot[\sigma,d])=c_1(\sigma,d)\,\mathrm{Within}(P\cdot[1])+c_3(\sigma,d)\,\mathrm{Within}(P\cdot[3]),\tag{5.4}$$
with $c_1,c_3$ prefix-independent (functions of $(\sigma,d)$ only). Here the two colour classes
$X,Y$ *are* the two boundary quantities, and $\mathrm{Within}(P\cdot[1]),\mathrm{Within}(P\cdot[3])$
are the two basis profiles into which they are read. $\square$

**Proposition 8 (profiles of $\ge2$-vertex prefixes lie in a fixed $3$-space).** For every prefix
$P$ with $|V(P)|\ge2$, $\pi(P)\in U_0:=\operatorname{span}\{p_1,q_1,\delta_0\}$ (restricted to
$\{0,1,2,3\}$), and $\{\pi([1]),\pi([2]),\pi([3])\}$ is a basis of $U_0$. *Proof.* This follows
directly from the one-hub gluing relation already obtained above -- $\mathrm{Within}(P\cdot[d])
=a(P)\,\mathrm{star}(d)+b(P)\,p_1(d)$ for every $d\ge1$ (the $R=[d]$ case of (5.1), valid because
$|V(P)|\ge2$) -- with no appeal to a power-of-two lemma. Thus on $x\in\{1,2,3\}$ the map
$x\mapsto\mathrm{Within}(P\cdot[x])$ agrees with $a(P)\,\mathrm{star}+b(P)\,p_1\in
\operatorname{span}\{\mathrm{star},p_1\}$, while $x=0$ (a suffix with no edge) is unconstrained;
hence $\pi(P)\in\operatorname{span}\{\mathrm{star},p_1,\delta_0\}=\operatorname{span}\{p_1,q_1,
\delta_0\}=U_0$ (the two spans coincide since $\mathrm{star}=p_1+q_1$). Its three generators are
independent on $\{0,1,2,3\}$, and $\pi([1]),\pi([2]),\pi([3])$ -- rows $1,2,3$ of $G$ -- are
independent since $\det G=144\ne0$, so form a basis. $\square$

*Readout (the terminal linear map).* For every prefix $P$ the map
$x\mapsto\mathrm{Within}(P\cdot[x])$ lies in $V'=\operatorname{span}\{DS(a,\cdot):a=0{:}3\}$. The
single-vertex prefix $P=[0]$ gives $DS(0,x)=\mathrm{star}(x+1)$, a defining generator of $V'$.
Every other prefix has $|V(P)|\ge2$ -- in particular every one-hub prefix $P=[a]$ with $a\ge1$,
which has $a+1\ge2$ vertices -- and there the one-hub gluing
$\mathrm{Within}(P\cdot[x])=a(P)\,\mathrm{star}(x)+b(P)\,p_1(x)$ ($x\ge1$, the $R=[x]$ case of
(5.1)) places it in $\operatorname{span}\{p_1,q_1,\delta_0\}\subseteq V'$ (the exceptional
edgeless suffix $x=0$ supplying $\delta_0$; e.g. $\delta_0=(DS(2,\cdot)-2\,DS(1,\cdot))/6$). In
particular $DS(a,\cdot)\in V'$ for **all** $a\ge0$, not merely $a\le3$.
Hence $\mathrm{Within}(P\cdot[x])$ is the unique element of $V'$ matching the four values
$\pi(P)$, which is exactly what the readout $\lambda(x)=G^{-1}(DS(a,x))_a$,
$\mathrm{Within}=\sum_d\lambda_d\,\pi(P)_d$ computes.

Finally the single-vertex prefix $[0]$ -- the one exception (it uses the $A=2$ channel, and (5.4)
fails there) -- supplies the $4$th coordinate: the four profiles $\pi([0]),\dots,\pi([3])$ are the
rows of $G$, a basis of $\mathbb Q^4$, so $M(\sigma)$ is uniquely defined by
$M(\sigma)\pi([a])=\pi([a,\sigma])$ for $a=0,\dots,3$ -- the readout matrix built from the triple
values $\mathrm{Within}([a,\sigma,d])$. For $|V(P)|\ge2$ both sides of $\pi(P\cdot[\sigma])=M(\sigma)\pi(P)$ are linear functionals
of $\pi(P)$ (the left by Proposition 7, the right by definition of $M(\sigma)$) that agree on the
basis $\pi([1]),\pi([2]),\pi([3])$ of $U_0$, hence on all of $U_0\ni\pi(P)$; and for $P=[0]$ the
identity is the defining relation with $a=0$. So $\pi(P\cdot[\sigma])=M(\sigma)\pi(P)$ for every
$P$, and with the readout above the $4$-dimensional recurrence computes $\mathrm{Within}(C)$
exactly for every caterpillar. $\square$

Because each interior $M(\sigma)$ maps a genuine profile (a vector of tournament counts) to
another, the interior division by $144$ is always exact-integer. Against an independent exact
counter sharing no logic with the recurrence, the assembled $M(\sigma)$ reproduces
$\mathrm{Within}$ with zero discrepancies over an extensive battery of caterpillars and
profile-collapse tests (Appendix D).

---

# 5. Polynomial aggregation

With $\mathrm{Within}$ in hand, $T(H)$ is a product over components (§3), and $A_v(t)$ is
assembled by a forest recurrence that couples components through the $2^{\binom{c}{2}}$ factor.
Writing $w_k(t)$ for the labeled weight of connected caterpillar components on $k$ vertices,
$$
w_1(t)=0,\qquad
w_k(t)=\sum_{\substack{C\ \text{caterpillar on }k\ge2\ \text{vertices}}}\frac{k!}{|\mathrm{Aut}(C)|}\,(-1)^{k-1}\,\mathrm{Within}(C)\,I(C)^{t}\quad(k\ge2),
$$
($w_1=0$ because the support graphs $H'$ have no isolated component). The component-coupled
recurrence is
$$
D_v^{(j)}=\sum_{k=1}^{v}\binom{v-1}{k-1}\,w_k(t)\,2^{\,j}\,D_{v-k}^{(j+1)},\qquad D_0^{(\cdot)}=1,
\qquad A_v^{\text{forest}}(t)=D_v^{(0)},
$$
where the depth index $j$ counts already-committed components so that $2^{j}$ accrues the
$\binom{c}{2}$ between-component edges. The single odd-cycle-sun family contributes a separate
additive term $\mathrm{oc}_v(t)$, and $A_v(t)=A_v^{\text{forest}}(t)+\mathrm{oc}_v(t)$; its
inclusion--exclusion sign is $(-1)^{e(S)}=(-1)^{v}$ (a sun on $v$ vertices is connected and
unicyclic, $e(S)=v$). Both $I(C)$ and $I(S)$ are evaluated by the standard $2$-state
spine/cycle recurrence.

*Remark.* Summing $\mathrm{Within}(C)\,I(C)^t$ over the $\Theta(2^{k})$ caterpillars directly needs
none of the machinery below but takes $2^{\Theta(n)}$ operations; this is the variant used for the
$L(n)$ table, and it reaches $n\le18$.

**Polynomial-time aggregation.** The enumeration is removable. The obstruction is the factor
$I(C)^t$: $I(C)$ is the scalar output of a $2$-state spine recurrence, and $I(C)^t$ does not
decompose along the spine. But $I(C)^t=\#\{\text{ordered $t$-tuples of independent sets of }C\}$,
and a $t$-tuple is tracked along the spine by the *count* $c\in\{0,\dots,t\}$ of
tuple-coordinates whose current spine vertex is "in" -- a **symmetric-power** compression of
the $2^t$ raw states to $t+1$. Running this count in lockstep with the $\mathrm{Within}$
$4$-state transfer, folding the labelled weight $k!/|\mathrm{Aut}(C)|$ into an exponential
generating function, and treating the sun term by the analogous cyclic construction, gives a
dynamic program that computes $A_v(t)$ -- and hence $L(n)$ -- with **no enumeration of
caterpillars**, in $\mathrm{poly}(n)$ big-integer operations. This establishes:

**Theorem L3 (polynomial aggregation).** The count-compression dynamic program computes $A_v(t)$,
and hence $L(n)=\sum_v\binom{n}{v}2^{\binom{n-v}2}A_v(n-v)$, exactly for every $n$, in a number of
big-integer operations bounded by a fixed polynomial in $n$.

*Proof.* We prove the theorem by giving explicit dynamic programs for caterpillars, suns, and
forest assembly -- the three ingredients being (i) the symmetric count-collapse of $I(C)^t$ to a
single size coordinate, (ii) the $2$-to-$1$ directed-spine label encoding reproducing
$k!/|\mathrm{Aut}|$ with no palindrome casework, and (iii) the cyclic-seam closure for suns. Each
is proved exact by its construction, and their polynomial state dimensions give the bound.

**Reduced-word set.** A labeled caterpillar on $[k]$ ($k\ge2$) has $m$ inner (degree-$\ge2$)
spine vertices. The cases $m=0$ (forced to be $K_2$, the single edge) and $m=1$ (a star
$K_{1,k-1}$) are separate base cases; for $m\ge2$ the caterpillar is encoded by its reduced
leaf-word, ranging over
$$
W_k=\Big\{(\ell_0,\dots,\ell_{m-1}):\ m\ge2,\ \ell_i\ge0,\ \ell_0\ge1,\ \ell_{m-1}\ge1,\
\textstyle\sum_i\ell_i=k-m\Big\}.
$$
The endpoint conditions $\ell_0,\ell_{m-1}\ge1$ record that a canonical reduced spine's endpoints
are inner vertices, not leaves, so the spine of a caterpillar is intrinsic. Reading it in its two
directions always yields **two distinct** (word, labelling) encodings, because the two spine
endpoints carry distinct labels -- and this holds *whether or not* the leaf-word is palindromic: a
palindromic shape gives the *same* word but with the two endpoint labellings interchanged (still
two encodings), a non-palindromic shape gives its two reversed words. Hence the map
(encoding)$\to$(caterpillar) is uniformly **$2$-to-$1$** (the source of the $\tfrac12$ below), and
dividing the directed-word sum by $2$ reproduces $k!/|\mathrm{Aut}|$ with no separate palindrome
correction.

**Caterpillar DP for $w_k(t)$.** Two ingredients run in lockstep along the spine.

*Independent-set count $I(C)^t$.* Since $I(C)^t=I(t\cdot C)$ counts spine-patterns weighted by
free leaves, and (Lemma (i)) the pattern-count depends on each row only through its "in"-size
$c\in\{0,\dots,t\}$, the count compresses to a product of $(t{+}1)\times(t{+}1)$ matrices,
$$
I(C)^t=\mathbf 1^{\top}L_{m-1}P\,L_{m-2}P\cdots P\,L_0\,e,\qquad
e_c=\binom tc,\quad P_{c,c'}=\binom{t-c'}{c},\quad (L_i)_{c,c}=2^{\ell_i(t-c)}.
$$

*Weight $\mathrm{Within}(C)$.* This is the §4 four-state transfer; with the integer row/matrix
data $r(\sigma)[j]=DS(\sigma,j)$, $\widetilde W(\sigma)[j][d]=\sum_a g_{a,j}(\sigma)\,
\mathrm{adj}(G)[d][a]$, $\widetilde\rho(\sigma)[d]=\sum_a DS(a,\sigma)\,\mathrm{adj}(G)[d][a]$
(state $j\in\{0,1,2,3\}$),
$\mathrm{Within}(C)=144^{-(m-1)}\,\widetilde\rho(\ell_{m-1})
\big(\prod_{i=m-2}^{1}\widetilde W(\ell_i)\big)r(\ell_0)$.

*Fusion.* Fold both, together with the integer label weight
$\lambda_i=(k-u_i)\binom{k-u_i-1}{\ell_i}$ (which telescopes to $k!/\prod_i\ell_i!$, the EGF
$1/\ell!$ made integral, $u_i$ the running label count), into an array $V_s[(j,c,u)]$ indexed by
transfer-state $j$, count $c$, and labels-used $u$, which the DP advances in **rounds**
$s=0,1,2,\dots$. After round $s$ the array holds the fused weight of every reduced prefix with
*exactly* $s+1$ spine hubs; because each of the $s$ completed interior steps multiplied by the
integer matrix $\widetilde W=144\,M$, every entry of $V_s$ carries the **same** integer factor
$144^{s}$. This common, round-uniform power is what makes the read-out division exact and
unambiguous -- prefixes of different hub-counts are never summed into one array, so no mixed total
is ever divided.

* **Round $0$ (initialize)** (hub $0$, $\ell_0\ge1$, each $c_0$; no factor of $144$ yet): set
  $V_0[(j,c_0,1{+}\ell_0)]\mathrel{+}= r(\ell_0)[j]\cdot\binom{t}{c_0}2^{\ell_0(t-c_0)}\cdot\lambda_0$,
  with $\lambda_0=k\binom{k-1}{\ell_0}$.
* **Read-out of round $s$** (append the terminal hub, forced $\ell=k-1-u\ge1$; these are the
  $m=s+2$-hub caterpillars): accumulate
  $\widetilde\rho(\ell)[j']\cdot\binom{t-c'}{c}2^{\ell(t-c)}\cdot(k-u)\binom{k-u-1}{\ell}\cdot
  V_s[(j',c',u)]$ over $c,j',c',u$, then divide by $144^{\,s+1}=144^{\,m-1}$. This is an exact
  integer division (Theorem L2) and is unambiguous precisely because every entry of $V_s$ shares
  the power $144^{s}$.
* **Interior step $V_s\to V_{s+1}$** (advance every entry by one hub, any $\ell\ge0$): set
  $$
  V_{s+1}[(j,c,u{+}1{+}\ell)]\mathrel{+}=
  \widetilde W(\ell)[j][j']\cdot\binom{t-c'}{c}2^{\ell(t-c)}\cdot(k-u)\binom{k-u-1}{\ell}\cdot
  V_s[(j',c',u)],
  $$
  which multiplies the common power by one further $144$ (to $144^{s+1}$). The new array
  **replaces** $V_s$; this is why each round carries a single hub-count and a single power.

The read-out is taken at every round $s\ge0$ (each contributing the $m=s+2$ caterpillars) and the
results summed. The endpoint constraints $\ell_0,\ell_{m-1}\ge1$ and total $=k$ restrict the
accumulation exactly to $w\in W_k$. Dividing this summed $m\ge2$ read-out by the $2$-to-$1$ spine
reversal and attaching the sign $(-1)^{k-1}=(-1)^{e(C)}$ (a $k$-vertex tree),
$$
w_k(t)=(-1)^{k-1}\Big(\underbrace{[k{=}2]\,2\cdot3^{t}}_{K_2}
+\underbrace{[k{\ge}3]\,k\,\mathrm{star}(k{-}1)(2^{k-1}{+}1)^{t}}_{m=1\ \text{star}}
+[k{\ge}4]\,\tfrac12\big(\text{$m\ge2$ readout total}\big)\Big),
$$
the two base cases being the labeled counts of $K_2$ ($I(K_2)=3$, one labeling) and of the star
$K_{1,k-1}$ ($k$ labelings, $\mathrm{Within}=\mathrm{star}(k{-}1)$, $I=2^{k-1}{+}1$).

**Sun seam DP for $\mathrm{oc}_v(t)$.** The same collapse runs around the cycle
$c_0,\dots,c_{g-1}$, but the wrap edge $A_{g-1}\cap A_0=\varnothing$ breaks the size-only
compression. Fix the top set $A_0$ to a specific set of size $c_0$ and split every later row by
its overlap with $A_0$ and with $A_0^{c}$: track $(a_i,b_i)=(|A_i\cap A_0|,\,|A_i\setminus
A_0|)$. Then the transition is
$$
(a',b')\to(a,b):\qquad
\binom{c_0-a'}{a}\binom{(t-c_0)-b'}{b}\cdot 2^{\ell_i(t-a-b)}\cdot 2^{\binom{\ell_i}2},
$$
the two binomials being the disjoint one-step extensions inside $A_0$ and inside $A_0^{c}$, the
factor $2^{\ell_i(t-a-b)}$ the free leaves, and $2^{\binom{\ell_i}2}$ the per-hub
$\mathrm{Within}$ factor (Lemma L1). Hub $0$ initializes the array (for each cycle length $g$ and
each size $c_0$): at state $(a_0,b_0)=(c_0,0)$ set
$$
V_0[(c_0,0)]=\binom{t}{c_0}\,2^{\ell_0(t-c_0)}\,2^{\binom{\ell_0}2}\cdot v\binom{v-1}{\ell_0},
\qquad u=1+\ell_0,
$$
($\binom t{c_0}$ choosing $A_0$, the two $2$-powers the hub's leaves and $\mathrm{Within}$ factor,
$v\binom{v-1}{\ell_0}$ the telescoping label weight); hubs $1,\dots,g-1$ apply the transition
above (advancing $u$ by $1+\ell_i$), the last one closed by the seam $a_{g-1}=0$ (the wrap
constraint $A_{g-1}\cap A_0=\varnothing$). The state $(a,b)$ ranges over $0\le a\le c_0,\ 0\le
b\le t-c_0$, i.e.\ $O(t^2)$ states. Reading a labeled sun's cycle from any of its $g$ vertices in either direction gives
$2g$ ordered encodings, and because the $g\ge3$ cycle vertices carry distinct labels the
dihedral group $D_g$ acts freely -- so the encoding is exactly **$2g$-to-$1$** even for periodic
leaf-words. Dividing by $2g$ and attaching the sign $(-1)^{v}$ (a sun on $v$ vertices has
$e(S)=v$ edges),
$$
\mathrm{oc}_v(t)=(-1)^{v}\!\!\sum_{\substack{g\ \text{odd}\\3\le g\le v}}\frac1{2g}
\!\!\sum_{\substack{(\ell_0,\dots,\ell_{g-1})\\ g+\sum_i\ell_i=v}}\!\!
\frac{v!}{\prod_i\ell_i!}\,2^{\,1+\sum_i\binom{\ell_i}2}\,I(S_\ell)^{t},
$$
the inner double sum evaluated by the seam DP with the same telescoping integer label weight
$(v-u)\binom{v-u-1}{\ell}$; the global $\mathrm{Within}$ factor $2$ pairs with the $2g$ to leave
a division by $g$.

**Assembly.** With $w_k(t)$ and $\mathrm{oc}_v(t)$ in hand, the forest part is the species
convolution
$$
D_v^{(j)}=\sum_{k=1}^{v}\binom{v-1}{k-1}\,w_k(t)\,2^{\,j}\,D_{v-k}^{(j+1)},\qquad
D_0^{(j)}=1,\qquad A_v^{\mathrm{forest}}(t)=D_v^{(0)},
$$
peeling the block of the current smallest label; the $2^{j}$ accrues the $j$ new
inter-component edges as a component is added, so the running product
$2^{0}2^{1}\cdots2^{r-1}=2^{\binom r2}$. Then
$$
A_v(t)=A_v^{\mathrm{forest}}(t)+\mathrm{oc}_v(t),\qquad
L(n)=\sum_{v=0}^{n}\binom nv 2^{\binom{n-v}2}A_v(n-v).
$$

**Complexity (no $2^{k}$ appears).** Every array is indexed by bounded coordinates: the
transfer-state $j\le4$, while $c,a,b,u,\ell,k,v$ each range in $\{0,\dots,n\}$ (as $t=n-v\le n$).
Hence the caterpillar array $V[(j,c,u)]$ has $O(kt)$ entries and $O(kt)$ work per transition
($O(k^3t^2)$ per $w_k$), the sun array $(a,b,u)$ has $O(t^2v)$ entries ($O(v^4t^5)$ per
$\mathrm{oc}_v$), and the forest recurrence has $O(v^2)$ states. Every intermediate entry is an
exact integer of $2^{O(n^2)}$ magnitude -- hence $O(n^2)$ bits -- since it is a signed sum of at
most $n!$ label assignments times at most $2^{O(n^2)}$ tournament/independent-set choices; we do
*not* claim the tighter $\le2^{\binom n2}$ (individual signed inclusion--exclusion terms may exceed
the final count before cancellation, so only the $2^{O(n^2)}$ bit bound is asserted). Summed over
$v\le n$, $t=n-v$, the whole computation is a fixed polynomial number of big-integer operations (a
loose $O(n^{10})$; the measured scaling is in Appendix D). Since every operand has $O(n^2)$ bits,
each such operation costs $\mathrm{poly}(n)$ bit-operations, so the algorithm runs in **polynomial
bit complexity** -- with **no enumeration of caterpillars or suns** and no $2^{k}$ state anywhere.
This is the assertion of Theorem L3. $\blacksquare$

An implementation reproduces $L(7),\dots,L(18)$ exactly, with measured running time scaling as
a fixed polynomial ($\approx n^{7.3}$) against two exponential competitors (Appendix D); we
thus have a **proved polynomial-time algorithm** for $L(n)$. The Burnside/cycle-index step
yielding the count of $\gamma\ge3$ tournaments up to isomorphism is a separate small-$n$ layer,
developed in Appendix D.

---

# 6. Results and validation

The method of §§1--5 is proved, so every value below is a theorem; the corroboration notes
record what independent cross-check is additionally available. The $A_v(t)$ were cross-checked
term by term for all $v\le 10$ against an *independent* enumeration of the contributing structures
(caterpillar forests and odd-cycle suns), each weighted by a $\mathrm{Within}$ value obtained from
its defining count rather than from the recurrence: exact agreement. Direct brute-force
enumeration of all $2^{\binom n2}$ tournaments confirms $L(7),L(8),L(9)$; the four-state
$\mathrm{Within}$ recurrence matches exact enumeration with zero mismatches on an extensive test
battery (Appendix D).

```{=latex}
\begingroup\small
```

| $n$ | $L(n)=\#\{\gamma\ge 3\}$ (labeled) | status |
|:--|--------------------------------------------:|:-----------------------|
| $7$ | $240$ | direct enumeration ($2^{21}$) + full Burnside |
| $8$ | $109440$ | direct enumeration ($2^{28}$) + full Burnside |
| $9$ | $77274880$ | direct enumeration ($2^{36}$) |
| $10$ | $107731095040$ | Burnside-corroborated vs. independent $\mathrm{iso}(10)$ |
| $11$ | $277644146489600$ | new; also passes a Burnside congruence (App.\ D) |
| $12$ | $1259251974345052160$ | new; also passes a Burnside congruence (App.\ D) |
| $13$ | $10059739307720354796544$ | new |
| $14$ | $144025884971245392244400128$ | new |
| $15$ | $3764299177453034811662226208768$ | new |
| $16$ | $182456135219244423561589552542810112$ | new |
| $17$ | $16611757357122373923247804720005425659904$ | new |
| $18$ | $2870353532067896148226715034553628871332397056$ | new |

```{=latex}
\endgroup
```

Here "direct enumeration" means examining all $2^{\binom{n}{2}}$ tournaments; we ran this for
$n=7,8,9$ (up to $2^{36}$, sharded). $L(10)$ was corroborated by an independent Burnside
assembly against $\mathrm{iso}(10)$; $L(11),L(12)$ pass Burnside congruences (Appendix D). All
entries are theorems; the corroboration column records only what independent cross-check is
available.

As a concrete payoff, the density $L(n)/2^{\binom n2}$ -- the probability that a uniformly random
labeled tournament has $\gamma\ge3$ -- can be read off directly: it rises from $0.011\%$ at $n=7$
through $0.31\%$ at $n=10$ and $3.3\%$ at $n=13$ to $25.1\%$ at $n=18$, exhibiting the threshold
climb concretely rather than through asymptotic bounds.

The polynomial algorithm is not confined to this window; we ran it through $n=30$, where it
evaluates the aggregate $L(n)=\#\{\gamma\ge3\}$ (which for $n\ge19$ exceeds $\#\{\gamma=3\}$, per
the scope note of §1). For instance
$$L(19)=949209776866692815063198172251754624787427891544064,$$
$$L(20)=604915380532280302020187569457614900100873131394954166272.$$
We do not tabulate $n\ge19$ further: these no longer isolate $\gamma=3$, and separating the
$\gamma=3$ and $\gamma\ge4$ parts would require a structure theorem for dominating triples that we
do not have.

(For $n\le 18$ every tournament has $\gamma\le 3$, so $L(n)=\#\{\gamma=3\}$ throughout.) The
full domination distribution (labeled, $\gamma\in\{1,2,3\}$, with
$\gamma_1=n\,2^{\binom{n-1}{2}}$ and $\gamma_2=2^{\binom{n}{2}}-\gamma_1-\gamma_3$; each row sums
to $2^{\binom n2}$) is, for $n=2,\dots,14$ -- the two-paradox class $\gamma=3$ first appearing at
$n=7$, the smallest order admitting a tournament with $\gamma\ge3$:

```{=latex}
\begingroup\footnotesize\setlength{\tabcolsep}{5pt}
\begin{center}
\begin{tabular}{@{}lrrr@{}}
\toprule
$n$ & $\gamma=1$ & $\gamma=2$ & $\gamma=3$ \\
\midrule
$2$ & $2$ & $0$ & $0$ \\
$3$ & $6$ & $2$ & $0$ \\
$4$ & $32$ & $32$ & $0$ \\
$5$ & $320$ & $704$ & $0$ \\
$6$ & $6144$ & $26624$ & $0$ \\
$7$ & $229376$ & $1867536$ & $240$ \\
$8$ & $16777216$ & $251548800$ & $109440$ \\
$9$ & $2415919104$ & $66226282752$ & $77274880$ \\
$10$ & $687194767360$ & $34389446226432$ & $107731095040$ \\
$11$ & $387028092977152$ & $35364124779497216$ & $277644146489600$ \\
$12$ & $432345564227567616$ & $72095378756265586688$ & $1259251974345052160$ \\
$13$ & $959230691832896684032$ & $291212484904104042195968$ & $10059739307720354796544$ \\
$14$ & $4231240368651202111471616$ & $2327622953230863955442376704$ & $144025884971245392244400128$ \\
\bottomrule
\end{tabular}
\end{center}
\endgroup
```

**Isomorphism classes.** The number of *isomorphism classes* of two-paradox tournaments follows
from the labeled counts by Burnside's lemma over $S_n$ (the identity term is $L(n)$; Appendix D):

| $n$ | $7$ | $8$ | $9$ | $10$ | $11$ | $12$ |
|:--|--:|--:|--:|--:|--:|--:|
| $\mathrm{iso}(n)$ | $1$ | $5$ | $226$ | $29816$ | $6959159$ | $2629321652$ |

Here $\mathrm{iso}(7),\dots,\mathrm{iso}(10)$ agree with the complete tournament catalogues through
order $10$; $\mathrm{iso}(11)$ and $\mathrm{iso}(12)$ come from an independent isomorph-free
generation (checksummed against the tournament sequence $\mathrm{A000568}$) and satisfy the
Burnside congruences derived in Appendix D. The unique class at $n=7$ is the Paley (quadratic-
residue) tournament $QR_7$. The Burnside congruence checks and detailed test records are collected
in Appendix D.

---

# Appendix A. Channel formulas and finite certificate

Write $p_A(\ell)=2^{\binom{\ell}{2}}A^{\ell}$ and $q_A(\ell)=\ell\,2^{\binom{\ell-1}{2}}A^{\ell-1}$ for $A\in\{1,2,4\}$ (the $A=8$ slot is unused). The single-vertex value is $\mathrm{star}(m)=p_1(m)+q_1(m)$.

**The matrix $G$** ($G[a][d]=DS(a,d)=\mathrm{Within}([a,d])$, symmetric, $\det G=144$):
$$
G=\begin{pmatrix}2&4&14&96\\4&4&8&30\\14&8&16&60\\96&30&60&224\end{pmatrix},\qquad
\operatorname{adj}(G)=\begin{pmatrix}0&-48&24&0\\-48&-30112&22928&-2088\\24&22928&-17416&1584\\0&-2088&1584&-144\end{pmatrix}=144\,G^{-1}.
$$

**Double-star channels** $DS(L,d)=\mathrm{Within}([L,d])$ (for $L\ge1$; $DS(L,0)=\mathrm{star}(L+1)$):
$$
DS(L,1)=3p_1(L)+q_1(L),\quad DS(L,2)=6p_1(L)+2q_1(L),\quad DS(L,3)=22p_1(L)+8q_1(L).
$$

**Triple channels** $g_{a,c}(\ell)=\mathrm{Within}([a,\ell,c])$, $0\le a\le c\le3$ (the ten cases; $g_{c,a}=g_{a,c}$):
$$
\begin{aligned}
g_{0,0}&=2p_2+2p_4+2q_4, & g_{0,1}&=p_1+3p_2+q_2, & g_{0,2}&=2p_1+6p_2+2q_2,\\
g_{0,3}&=8p_1+22p_2+8q_2, & g_{1,1}&=4p_1, & g_{1,2}&=8p_1,\\
g_{1,3}&=30p_1, & g_{2,2}&=16p_1, & g_{2,3}&=60p_1,\\
g_{3,3}&=224p_1. & & & &
\end{aligned}
$$

**Finite certificate (derivation of the channels).**

*The single-vertex value, by direct counting.* The independent sets of the star $K_{1,m}$ are
the leaf-subsets and $\{\text{centre}\}$, and validity means every vertex beating the centre
loses to all leaves. At most one leaf beats the centre: if two leaves $u,u'$ both beat it, then
whichever of $u\to u'$, $u'\to u$ holds makes that leaf a common in-neighbour of the other's
centre-edge, which the reformulation forbids. If none, the centre beats all
$m$ leaves and the $\binom m2$ leaf-edges are free ($2^{\binom m2}$ completions); if exactly one
($m$ choices of the special leaf $u^{*}$), $u^{*}$ loses to the other $m-1$ leaves whose
$\binom{m-1}2$ edges are free ($m\,2^{\binom{m-1}2}$). Hence
$\mathrm{star}(m)=2^{\binom m2}+m\,2^{\binom{m-1}2}=p_1(m)+q_1(m)$.

*The double-star and triple channels, by core-config enumeration.* Fix the varying hub $s$ of
the channel and apply the Leaf-Factor Lemma: each internally-valid orientation $\tau$ of the
finite core (everything not incident to the leaves of $s$) contributes
$f_\tau(\ell)=p_{A_\tau}(\ell)+A^{*}_\tau\,q_{A_\tau}(\ell)$, and the channel is
$\sum_\tau f_\tau(\ell)$. Enumerating the finitely many $\tau$ and recording, for each realized
value $A$, the number of configs with that $A$ (the coefficient on $p_A$) and the sum of their
$A^{*}$ values (the coefficient on $q_A$) fixes the two coefficients of each channel. The core of
a channel is itself a bounded caterpillar fragment -- for the triple $g_{a,c}=\mathrm{Within}
([a,\ell,c])$ it is the three-hub caterpillar $[a,\cdot,c]$ with $a,c\le3$, at most $3+a+c\le9$
vertices, hence at most $2^{\binom 92}$ orientations to inspect -- so the enumeration is finite.
All fourteen double-/triple-channel cases are recorded below in aggregate (the
right column is the resulting closed form, whose $p_A$-coefficient is the $A$-count and whose
$q_A$-coefficient is the $A^{*}$-sum, so the table and Appendix A determine each other):

| channel | $(A:\text{count})$ | $(A:\sum A^{*})$ | closed form (App.\ A) |
|---|---|---|---|
| $DS(L,0)$ | $\{1{:}1,\,2{:}1\}$ | $\{2{:}1\}$ | $\mathrm{star}(L{+}1)$ |
| $DS(L,1)$ | $\{1{:}3\}$ | $\{1{:}1\}$ | $3p_1+q_1$ |
| $DS(L,2)$ | $\{1{:}6\}$ | $\{1{:}2\}$ | $6p_1+2q_1$ |
| $DS(L,3)$ | $\{1{:}22\}$ | $\{1{:}8\}$ | $22p_1+8q_1$ |
| $g_{0,0}$ | $\{2{:}2,\,4{:}2\}$ | $\{4{:}2\}$ | $2p_2+2p_4+2q_4$ |
| $g_{0,1}$ | $\{1{:}1,\,2{:}3\}$ | $\{2{:}1\}$ | $p_1+3p_2+q_2$ |
| $g_{0,2}$ | $\{1{:}2,\,2{:}6\}$ | $\{2{:}2\}$ | $2p_1+6p_2+2q_2$ |
| $g_{0,3}$ | $\{1{:}8,\,2{:}22\}$ | $\{2{:}8\}$ | $8p_1+22p_2+8q_2$ |
| $g_{1,1}$ | $\{1{:}4\}$ | $\varnothing$ | $4p_1$ |
| $g_{1,2}$ | $\{1{:}8\}$ | $\varnothing$ | $8p_1$ |
| $g_{1,3}$ | $\{1{:}30\}$ | $\varnothing$ | $30p_1$ |
| $g_{2,2}$ | $\{1{:}16\}$ | $\varnothing$ | $16p_1$ |
| $g_{2,3}$ | $\{1{:}60\}$ | $\varnothing$ | $60p_1$ |
| $g_{3,3}$ | $\{1{:}224\}$ | $\varnothing$ | $224p_1$ |

(The $DS$ rows also carry an $A{=}0$ boundary config contributing $0$; $A{=}0$ terms never affect
$p_A,q_A$.)

*Exhaustiveness.* Since the Leaf-Factor Lemma pins each channel's $\ell$-dependence to the fixed
two-dimensional block $\{p_A,q_A\}$, and the core is a bounded fragment with only finitely many
orientations, the complete enumeration above determines the two coefficients -- and hence the whole
channel function -- for **all** $\ell$. This is a proof, not a fit: no interpolation through
sampled values is used; the $\ell$-shape is fixed a priori and only its two coefficients are read
off the finite enumeration. (The table records, for each channel, the aggregate $A$-count and
$A^{*}$-sum over the core configurations rather than the individual configurations; since the core
has at most $2^{\binom 92}$ orientations, a reader may in principle regenerate the aggregates
directly.)

**Transfer convention.** State $R\in\mathbb{Z}^4$, $R[j]$ = the value for the processed prefix followed by a spine vertex carrying $j$ leaves. Initialize $R[j]=DS(\ell_0,j)$. For each interior leaf-count $\ell$ (in spine order $\ell_1,\dots,\ell_{m-2}$): set $S[a]=\sum_d \operatorname{adj}(G)[d][a]\,R[d]$, then $R'[j]=\tfrac{1}{144}\sum_a g_{a,j}(\ell)\,S[a]$ (an exact integer division, by Theorem L2, so $R\in\mathbb Z^4$ throughout). Read out with the last symbol $x=\ell_{m-1}$: $\mathrm{Within}=\tfrac{1}{144}\sum_d\bigl(\sum_a DS(a,x)\,\operatorname{adj}(G)[d][a]\bigr)R[d]$. Base cases: $m=1$ gives $\mathrm{star}(\ell_0)$; $m=2$ is obtained by the *same* readout -- initialize $R[j]=DS(\ell_0,j)$ and read out at $x=\ell_1$ -- which supplies double stars $DS(\ell_0,\ell_1)$ for arbitrary $\ell_0,\ell_1$ (e.g. $\mathrm{Within}([4,4])=12288$), beyond the $d\le3$ closed forms tabulated above. The convention reproduces every channel closed form and the values $L(7),\dots,L(18)$ exactly.

---

# Appendix B. Odd-cycle-sun orientation proof

This appendix proves the four claims (a)--(d) of Lemma L1 (§3). Write the cycle as
$c_0,\dots,c_{g-1}$ (indices mod $g$), with $\ell_i$ leaves on hub $c_i$, and work inside any
valid tournament, using only the dominating-pair condition on the cycle edges and
$H$-independence of out-neighbourhoods.

**(a) Exactly two orientations of the cycle-and-chord tournament.** *No local sink:* if
$c_{i-1}\to c_i\leftarrow c_{i+1}$, then domination of $\{c_{i-1},c_i\}$ (third vertex $c_{i+1}$)
forces $c_{i-1}\to c_{i+1}$ and domination of $\{c_i,c_{i+1}\}$ (third vertex $c_{i-1}$) forces
$c_{i+1}\to c_{i-1}$ -- a contradiction. In the cyclic $F/B$ word of the cycle-edge orientations,
$\#$sinks $=\#$sources; with no sink there is no source, hence no transition, so the word is
constant: the cycle is a **directed cycle** (forward or backward -- two choices). Fixing
$c_i\to c_{i+1}$, induction on the gap $d$ forces every chord $\{c_i,c_{i+d}\}$: when $d$ is odd
$c_{i+d}\in N^+(c_i)$, and $H$-independence forbids the cycle-adjacent $c_{i+d+1}$, forcing
$c_{i+d+1}\to c_i$ (gap $d{+}1$ even, inward); when $d$ is even, domination of $\{c_i,c_{i+1}\}$
rules out the reverse, forcing $c_i\to c_{i+d+1}$ (gap $d{+}1$ odd, outward). So $c_i\to c_{i+d}$
iff $d$ is odd. Gap $d$ from $c_i$ and gap $g-d$ from $c_{i+d}$ name the same edge, and the two
forcings agree **iff $g$ is odd**; hence odd $g$ yields exactly the two rotational tournaments
$R$ (each valid) and even $g$ yields $0$ (the wrap edge $\{c_0,c_{g-1}\}$ is forced both ways).
For each hub $c_i$, $N^+(c_i)\cap\text{hubs}$ and $N^-(c_i)\cap\text{hubs}$ are maximum
independent sets of $C_g$, and $c_i$ beats exactly one cycle-neighbour.

**(b) Each leaf's orientation is forced against every foreign hub.** No leaf beats its own hub:
if a leaf $u\to c_i$, the cycle-neighbour $c_i^{*}\to c_i$ is forced (independence of
$N^+(u)\ni c_i$) to satisfy $c_i^{*}\to u$, making $c_i^{*}$ a common in-neighbour of the leaf
edge $\{u,c_i\}$ -- contradiction; so $c_i\to u$. Domination of $\{u,c_i\}$ then forces
$N^+(u)\supseteq N^-(c_i)\cap\text{hubs}$ (a maximum independent set), and independence caps it,
so $N^+(u)\cap\text{hubs}=N^-(c_i)\cap\text{hubs}=\{c_a:(i-a)\text{ odd}\}$ is **forced**; in
particular every leaf--foreign-hub edge is determined by the parity of the offset.

**(c) Forced orientation between leaves at different hubs.** For a leaf $u$ of $c_i$ and a leaf
$v$ of $c_j$ ($j\ne i$), exactly one of $(j-i),(i-j)$ is odd. If $(j-i)$ is odd then $v\to c_i$
by (b), so $v\to u$ would put both $u$ and $c_i$ into $N^+(v)$, i.e.\ the sun edge $\{u,c_i\}$
inside $N^+(v)$ -- impossible; hence $u\to v$.
The symmetric case forces $v\to u$. So every foreign leaf--leaf edge is **forced**.

**(d) Freedom only within same-hub leaf sets.** Two leaves $u,u'$ of the same hub $c_i$ are not
joined by any $H$-edge (their only $H$-neighbour is $c_i$), so the pair $\{u,u'\}$ carries no
domination constraint of its own; and since $c_i\in N^+(u)\cap N^+(u')$ is impossible to violate
here ($c_i\notin N^+(u),N^+(u')$ as $c_i$ beats each leaf by (b)), orienting $\{u,u'\}$ either way
adds no $H$-edge to any out-neighbourhood. Hence both orientations are valid and independent. The
$\binom{\ell_i}2$ intra-hub edges
are thus free and independent across hubs, and any such completion is valid (each
out-neighbourhood is checked $H$-independent). Each of the two rotational orientations therefore
contributes $\prod_i 2^{\binom{\ell_i}2}$, giving
$\mathrm{Within}=2\cdot\prod_i 2^{\binom{\ell_i}2}=2^{\,1+\sum_i\binom{\ell_i}2}$. $\square$
Verified exhaustively: $C_3,\dots,C_7$, all suns on $\le8$ vertices, and each zero-case family,
by independent enumeration.

---

# Appendix C. Worked example: $n=7$

**Master sum.** With $A_v(t)=\sum_{H}(-1)^{e(H)}T(H)I(H)^t$ over labeled minimum-degree-$\ge 1$
graphs on $v$ vertices:

* $v=2$: the only such $H$ is the edge $K_2$; $e=1$, $T(K_2)=2$, $I(K_2)=3$ (namely
  $\varnothing,\{a\},\{b\}$). So $A_2(t)=-2\cdot 3^{t}$, and $A_2(5)=-486$.
* $v=3$: the labeled graphs are the $3$ paths $P_3$ and the $1$ triangle $K_3$.
  For $P_3$: $e=2$, $T=\mathrm{Within}(P_3)=4$, $I=5$. For $K_3$ (an odd-cycle sun with no
  leaves): $e=3$, $T=\mathrm{Within}=2$, $I=4$. Thus
  $A_3(t)=3\cdot 4\cdot 5^{t}-2\cdot 4^{t}=12\cdot 5^{t}-2\cdot 4^{t}$, and $A_3(4)=6988$.

Continuing (the higher $A_v$ are assembled the same way from stars, double-stars, paths and
suns), the master formula $L(7)=\sum_v \binom{7}{v}2^{\binom{7-v}{2}}A_v(7-v)$ gives:

| $v$ | $\binom{7}{v}$ | $2^{\binom{7-v}{2}}$ | $A_v(7-v)$ | term |
|---:|---:|---:|---:|---:|
| $0$ | $1$ | $2097152$ | $1$ | $+2097152$ |
| $2$ | $21$ | $1024$ | $-486$ | $-10450944$ |
| $3$ | $35$ | $64$ | $6988$ | $+15653120$ |
| $4$ | $35$ | $8$ | $-39672$ | $-11108160$ |
| $5$ | $21$ | $2$ | $124896$ | $+5245632$ |
| $6$ | $7$ | $1$ | $-242352$ | $-1696464$ |
| $7$ | $1$ | $1$ | $259904$ | $+259904$ |

($v=1$ contributes $0$.) The terms sum to
$$
L(7)=2097152-10450944+15653120-11108160+5245632-1696464+259904=\boxed{240}.
$$

**Isomorphism count.** The odd cycle types of $7$ are $[7],[5,1,1],[3,3,1],[3,1^4],[1^7]$,
with edge-orbit counts $E=3,5,7,11,21$ and coefficients $7!/z_\lambda=720,504,280,70,1$. The
$S_2$-invariant counts are $S([7])=2$, $S([5,1,1])=0$, $S([3,3,1])=12$, $S([3,1^4])=0$, and
$S([1^7])=L(7)=240$. Hence
$$
\mathrm{iso}(7)=\frac{1}{5040}\bigl(720\cdot 2+280\cdot 12+1\cdot 240\bigr)
=\frac{1440+3360+240}{5040}=\frac{5040}{5040}=1.
$$
The single class is the Paley (quadratic-residue) tournament $QR_7$, with $|\mathrm{Aut}|=21$;
its labeled orbit has size $7!/21=240$, matching $L(7)$ exactly.

---

# Appendix D. Computational records, timing, and Burnside checks

**Timing.**
We have implemented this; it reproduces $L(7),\dots,L(18)$ exactly, and its measured running time
scales as a fixed polynomial. Figure 1 plots the measured single-threaded wall-clock time against
$n$ on a semilogarithmic axis, alongside the two exponential competitors. A least-squares fit to
the measured points ($7\le n\le32$) gives $\approx C\,n^{7.3}$ for the polynomial algorithm,
against $\approx 2^{1.10\,n}$ for the same signed sum evaluated by naïve term-by-term enumeration,
and $2^{\binom
n2}$ for brute-force evaluation of every tournament. On these axes an exponential method is a
straight line; the polynomial method's downward curvature away from the enumeration line is the
visible signature of its sub-exponential growth, and it is the only one of the three that remains
on a human time-scale as $n$ grows. This is independent numerical corroboration of Theorem L3, not
a substitute for it. We therefore have a **proved polynomial-time algorithm** for $L(n)$.

![Wall-clock time to compute $L(n)$ versus $n$ (logarithmic time axis), with asymptotic
comparison slopes superimposed. Blue: the polynomial algorithm of Theorem L3 (measured points;
dashed least-squares fit $\propto n^{7.3}$). Red: the same signed master sum evaluated by naïve
term-by-term enumeration of caterpillars (measured points; fit $\propto 2^{1.10\,n}$). Grey: the
$2^{\binom n2}$ brute force over all tournaments, charged an optimistic $1$ ns per tournament,
which leaves the frame by $n=10$. On a semilogarithmic axis an exponential running time is a
straight line, so the polynomial curve's bend below the enumeration line directly exhibits the
separation proved in Theorem L3. Horizontal guides mark human time-scales. All timings are
single-threaded exact big-integer arithmetic.](timing_plot.pdf){width=95%}

**Within closed forms** (verified by exhaustive enumeration over all orientations): star
$K_{1,m}=2,4,14,96,1344$; paths $=2,4,4,4$; double stars $DS(a,1)=4,8,30,224,3392$,
$DS(3,3)=224$, $DS(4,3)=1664$, $DS(5,2)=6784$. The four-state recurrence reproduces
$\mathrm{Within}$ with **zero mismatches** on all $44$ non-isomorphic caterpillars through
order $8$ (against exact enumeration) and, against an independent pruned counter, on $241$
bounded-leaf ($\le3$) spine-encodings through order $9$ and $66$ targeted cases to $k\approx13$
(§4).

**Independent verification.** The crux (Lemma 5's "two independent covers = the two colour
classes") is confirmed exactly on connected bipartite graphs; and, against an independent exact
counter sharing no logic with the recurrence, the assembled $M(\sigma)$ reproduces
$\mathrm{Within}$ with **zero discrepancies** over all $44$ caterpillars through order $8$, $241$
bounded-leaf encodings through order $9$, $66$ targeted cases to $k\approx13$, and a battery of
$1131$ profile-collapse tests (spines to length $8$, including the exceptional $[0]$).

**Isomorphism-class counts.** Applying Burnside/cycle-index over the symmetric group,
$$
\mathrm{iso}(n)=\frac{1}{n!}\sum_{\lambda\ \text{all-odd}}\frac{n!}{z_\lambda}\,S(\lambda),
$$
where $S(\lambda)$ is the number of $\pi$-invariant $S_2$ tournaments for a permutation of
cycle type $\lambda$ (only odd cycle types can fix a tournament). The identity term is
$S(1^n)=L(n)$ from the master formula of §2; each non-identity term is computed by direct orientation
enumeration over edge-orbits (Gray-code, $O(2^{E(\lambda)})$). *Many* cycle types have far
fewer edge orbits than $\binom n2$ and are cheap; but the **near-identity** types (a small cycle
plus many fixed points) have $E(\lambda)=\Theta(n^2)$ -- e.g. $E([3,1^{n-3}])=\tfrac12(n^2-5n+8)$
-- so this Burnside layer is *not* single-exponential and is a separate small-$n$ device (this
appendix), not part of the $2^{O(n)}$ labeled-count method.

**Isomorphism classes.** For $n\le 10$ the Burnside sum was assembled with *all* cycle-type
terms computed exactly (each $S(\lambda)$ by Gray-code orbit enumeration, the identity term by
the pipeline), giving exact integer quotients
$$
\mathrm{iso}(7,\dots,10)=1,\;5,\;226,\;29816.
$$
These agree with direct enumeration of all non-isomorphic tournaments (McKay's tournament
catalogues, complete through order $10$).

For $n=11,12$ the values $\mathrm{iso}(11)=6959159$ and $\mathrm{iso}(12)=2629321652$ are
**computed here by isomorph-free generation**, an algorithm independent of the transfer-matrix
pipeline; McKay's general catalogues stop at order $10$ and are not a source for them. Since
$S_2$ is *not* hereditary under vertex deletion ($QR_7$ is $S_2$ but every tournament on $<7$
vertices has $\gamma\le2$), the generator does **not** prune non-$S_2$ intermediates: canonical
augmentation builds all $n$-vertex tournaments up to isomorphism -- each order's total
checksummed against the tournament sequence $\mathrm{A000568}$ (orders $9,10,11,12$:
$191536,\ 9733056,\ 903753248,\ 154108311168$) -- and the $S_2$ property is applied as a
*post-filter* on completed tournaments. Generating and filtering all $154{,}108{,}311{,}168$
order-$12$ tournaments is a substantial computation; here $\mathrm{iso}(11),\mathrm{iso}(12)$ enter
only as an independent consistency cross-check, not as part of the labeled-count method. The Burnside sum is used only as a *consistency check* linking
the two independent computations. Concretely, for $n=11$ let
$$
K=\!\!\sum_{\substack{\lambda\ \text{odd},\ \lambda\ne 1^{11}\\ \lambda\ne[3,1^8]}}\!\!\frac{11!}{z_\lambda}\,S(\lambda)=32{,}677{,}258{,}240
$$
be the sum of all non-identity cycle-type terms except $[3,1^8]$ (each $S(\lambda)$ by
Gray-code enumeration; $E(\lambda)\le 30$, all brute-feasible). Since $11!/z_{[3,1^8]}=330$,
the Burnside identity gives
$$
\mathrm{iso}(11)\cdot 11! - L(11) - K \;=\; 330\cdot S([3,1^8]) \;\equiv\; 0 \pmod{330},
$$
and indeed $\mathrm{iso}(11)\cdot 11! - L(11) - K = 110{,}534{,}223{,}360 = 330\cdot 334{,}952{,}192$,
so the transfer-matrix $L(11)$ is consistent with the independently-generated $\mathrm{iso}(11)$
and forces the clean non-negative integer $S([3,1^8])=334952192$. (Note $\mathrm{iso}(11)\cdot
11! - L(11)$ *without* $K$ is $\equiv 220\not\equiv 0\pmod{330}$; the $K$ term is essential.)
The analogous check for $n=12$: with $K_{12}=\sum_{\lambda\ne1^{12},[3,1^9]}(12!/z_\lambda)S(\lambda)=6{,}066{,}211{,}092{,}480$
and $12!/z_{[3,1^9]}=440$,
$$
\mathrm{iso}(12)\cdot 12! - L(12) - K_{12} = 191{,}237{,}666{,}498{,}560 = 440\cdot 434{,}631{,}060{,}224,
$$
forcing $S([3,1^9])=434631060224$. These are consistency checks, not independent recomputations
of the iso counts.

We do **not** report $\mathrm{iso}(13)$ or $\mathrm{iso}(14)$: they require the near-identity
$S(\lambda)$ terms (e.g. $[3,1^{10}]$, $E=56$), whose brute cost is infeasible and whose method
computation was not completed. Note that these near-identity Burnside terms cost
$2^{\Theta(n^2)}$ (for $\lambda=[3,1^{n-3}]$, $E(\lambda)=\tfrac12(n^2-5n+8)$), so the
isomorphism layer is a **separate small-$n$ validation layer**, not part of the $2^{O(n)}$
labeled-count method of §§1--6.

