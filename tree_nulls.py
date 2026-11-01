"""Random-tree null models and comparison metrics for the Sundarbans channel network.

Null models (all on the SAME nodes as the real network, so size-dependent numbers are comparable):

* random_channel_tree  - random spanning forest of the real network rooted at the outlets
                          (Wilson's algorithm = uniform over all such forests). Loops are cut at
                          random, only real channels are used. -> "what do the loops add?"
* random_parent_tree   - every junction / inlet picks a random parent among the nodes already
                          placed, anywhere on the map (outlets are roots, inlets are leaves,
                          junctions take at most 2 children so degree <= 3 like real junctions).
                          -> "what do the loops AND the geography add?"

Every tree is returned as a directed MultiDiGraph (child -> parent, i.e. towards the sea)
plus its undirected MultiGraph, so the same metric code runs on real and null graphs.
"""
import networkx as nx
import numpy as np

SEA = 'SEA'


# ---------------------------------------------------------------- null models
def _tree_from_parent(G, parent):
    D = nx.MultiDiGraph()
    D.add_nodes_from(G.nodes(data=True))
    for child, (p, key) in parent.items():
        d = G.edges[child, p, key] if key is not None and G.has_edge(child, p, key) else {}
        D.add_edge(child, p, key=key, **d)
    return D, nx.MultiGraph(D.to_undirected())


def random_channel_tree(G, outlets, rng):
    """Uniform random spanning forest of G whose trees are rooted at the outlets (Wilson 1996)."""
    nbrs = {n: [(v, k) for _, v, k in G.edges(n, keys=True)] for n in G}
    in_tree, parent = set(outlets), {}
    order = list(G.nodes)
    rng.shuffle(order)
    for start in order:
        u, step = start, {}
        while u not in in_tree:                      # loop-erased random walk until it hits the forest
            step[u] = nbrs[u][rng.integers(len(nbrs[u]))]
            u = step[u][0]
        u = start
        while u not in in_tree:
            in_tree.add(u)
            parent[u] = step[u]
            u = step[u][0]
    return _tree_from_parent(G, parent)


def random_parent_tree(G, inlets, outlets, rng, root_children=1, junction_children=2):
    """Random tree on the same nodes: outlets = roots, inlets = leaves, junctions get random parents."""
    junctions = [n for n in G if n not in set(inlets) | set(outlets)]
    free = {o: root_children for o in outlets}
    parent = {}

    def pick():
        nodes = [n for n, c in free.items() if c > 0]
        w = np.array([free[n] for n in nodes], float)
        p = nodes[rng.choice(len(nodes), p=w / w.sum())]
        free[p] -= 1
        return p

    for j in rng.permutation(junctions):
        parent[j] = (pick(), None)
        free[j] = junction_children
    for i in rng.permutation(inlets):
        parent[i] = (pick(), None)
    D, U = _tree_from_parent(G, parent)
    return D, U


# ---------------------------------------------------------------- metrics
def _with_sea(U, outlets):
    """Simple weighted graph (weight = number of parallel links = conductance) + one SEA node joined to every outlet."""
    S = nx.Graph()
    S.add_nodes_from(U.nodes)
    for u, v in U.edges():
        if u != v:
            S.add_edge(u, v, weight=S[u][v]['weight'] + 1 if S.has_edge(u, v) else 1)
    for o in outlets:
        S.add_edge(SEA, o, weight=1)
    return S


def _n_bridges(U):
    S = nx.Graph()
    seen = set()
    for u, v, k in U.edges(keys=True):
        if (u, v) in seen or (v, u) in seen or u == v:   # parallel link: subdivide so it forms a loop
            m = ('virtual', u, v, k)
            S.add_edge(u, m); S.add_edge(m, v)
        else:
            S.add_edge(u, v); seen.add((u, v))
    br = [(u, v) for u, v in nx.bridges(S)]
    return len(br)


def _n_paths(D, outlets):
    """Number of distinct downstream routes from every node to any outlet (DAG dynamic programming)."""
    npth = {}
    for n in reversed(list(nx.topological_sort(D))):
        npth[n] = 1 if n in outlets else sum(npth[v] for _, v in D.out_edges(n))
    return npth


def gini(x):
    x = np.sort(np.asarray(x, float))
    n = len(x)
    return (2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()) if x.sum() > 0 else 0.0


def robustness_curve(U, inlets, outlets, fracs, reps, rng):
    """Mean share of inlets still connected to the sea after removing a random fraction of links."""
    edges = list(U.edges(keys=True))
    out = []
    for f in fracs:
        vals = []
        for _ in range(reps):
            H = U.copy()
            k = int(round(f * len(edges)))
            if k:
                H.remove_edges_from([edges[i] for i in rng.choice(len(edges), k, replace=False)])
            H.add_edges_from((SEA, o) for o in outlets)
            reach = nx.node_connected_component(H, SEA)
            vals.append(np.mean([i in reach for i in inlets]))
        out.append(np.mean(vals))
    return np.array(out)


def width_function(S, maxhop=40):
    """Number of nodes at each hop distance from the sea (Tejedor's delta width function, topological)."""
    d = nx.single_source_shortest_path_length(S, SEA)
    w = np.bincount([v for n, v in d.items() if n != SEA], minlength=maxhop + 1)[:maxhop + 1]
    return w


def sea_oriented(U, hops):
    """Direction-free flow orientation: every link points from the node farther from the sea (in links)
    to the closer one. Links joining two nodes at the same distance are left out. Always a DAG."""
    D = nx.MultiDiGraph()
    D.add_nodes_from(U.nodes(data=True))
    for u, v, k, d in U.edges(keys=True, data=True):
        if hops[u] > hops[v]:
            D.add_edge(u, v, key=k, **d)
        elif hops[v] > hops[u]:
            D.add_edge(v, u, key=k, **d)
    return D


def compare_metrics(U, D, inlets, outlets, orientation='given'):
    """Scalar parameters used for the real-vs-tree comparison.
    orientation='given' uses D (e.g. RivGraph flow directions); 'sea_distance' ignores D and orients links towards the sea."""
    inlets, outlets = list(inlets), list(outlets)
    V, E, C = U.number_of_nodes(), U.number_of_edges(), nx.number_connected_components(U)
    deg = np.array([d for _, d in U.degree()])
    S = _with_sea(U, outlets)

    hops = nx.single_source_shortest_path_length(S, SEA)
    if orientation == 'sea_distance':
        D = sea_oriented(U, hops)
    inlet_hops = np.array([hops[i] - 1 for i in inlets])       # links from inlet to its outlet

    # effective resistance inlet -> sea (1 ohm per link, parallel links in parallel)
    nodes = list(S.nodes)
    ix = {n: i for i, n in enumerate(nodes)}
    L = nx.laplacian_matrix(S, nodelist=nodes, weight='weight').toarray().astype(float)
    lam = np.linalg.eigvalsh(L)
    Lp = np.linalg.pinv(L, hermitian=True)
    s = ix[SEA]
    R = np.array([Lp[ix[i], ix[i]] + Lp[s, s] - 2 * Lp[ix[i], s] for i in inlets]) - 1   # minus the outlet->sea link

    npth = _n_paths(D, set(outlets))
    ebc = nx.edge_betweenness_centrality(nx.Graph(S))

    return {
        'independent cycles': E - V + C,
        'meshedness': (E - V + C) / (2 * V - 5),
        'bridge share': _n_bridges(U) / E,
        'alt. paths per inlet (log10)': np.log10(np.mean([npth[i] for i in inlets])),
        'outlets reachable per inlet': np.mean([len(nx.descendants(D, i) & set(outlets)) for i in inlets]),
        'fork share (water splits)': np.mean([D.out_degree(n) >= 2 for n in D]),
        'dead-end junctions': np.mean([(d == 1) and (n not in inlets) and (n not in outlets) for n, d in U.degree()]),
        'mean inlet→sea path (links)': inlet_hops.mean(),
        'longest inlet→sea path (links)': inlet_hops.max(),
        'inlet path spread (CV)': inlet_hops.std() / inlet_hops.mean(),
        'resistance inlet→sea (Ω)': R.mean(),
        'resistance / path length': np.mean(R / inlet_hops),
        'global efficiency': nx.global_efficiency(nx.Graph(S)),
        'algebraic connectivity λ2': lam[1],
        'edge betweenness Gini': gini(list(ebc.values())),
        'share of degree-3 nodes': np.mean(deg == 3),
    }
