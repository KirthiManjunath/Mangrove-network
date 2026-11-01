"""Weighted, simplified Sundarbans network: image -> cleaned RivGraph graph -> weights -> simplification.

Keep this file next to sundarbans_rg.py and run it with the `rivgraph` conda environment.

    from sundarbans_simplified_network import build_simplified_network
    net = build_simplified_network('18th Jan 2016.tiff', 'example_2016_simplified')
    net['G']          # simplified nx.MultiGraph (edge attrs: weight, conductance, length, width, link_id, idx)
    net['D']          # simplified nx.MultiDiGraph (same contractions, RivGraph directions)
    net['G_raw']      # graph before simplification, with weights added
    net['inlets'], net['outlets'], net['POS'], net['history_df'], net['spatial_history_df']

    to_unweighted(net['G'])                       # same graph with every width/length/weight removed
    load_graphs('.../all_graphs.pkl', weighted=False)   # saved G1..G10 without weights

From the command line (saves a before/after figure into the output folder):

    python sundarbans_simplified_network.py "18th Jan 2016.tiff" example_2016_simplified

Rules that keep the graph valid:
* Inlets and outlets are never contracted. They define the network's boundary, and the
  random-tree baselines need them to survive.
* A link with a parallel partner is only contracted if every parallel link is itself short.
  Contracting it otherwise would silently delete a real loop around an island.
* Merged nodes keep kind='junction' and also get merged=True and merged_from=[...].
* No node has a self-edge: any in the input are removed first, and contraction never creates one.

tree_nulls.compare_metrics works on the result with orientation='sea_distance'. RivGraph's
directions in net['D'] can contain loops after contraction, so check nx.is_directed_acyclic_graph
before using orientation='given'.
"""
import math
import pickle
import sys

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

from sundarbans_rg import extract_network

KIND_COL = {'inlet': '#2ca02c', 'outlet': '#ff7f0e', 'junction': '#444444'}


def add_geometric_weights(G):
    """Add max-normalised width/length weights to every edge of a MultiGraph, in place.

    raw_weight = width / length, weight = raw_weight / max(raw_weight), so weights lie in (0, 1].
    The raw value is kept as `conductance`. `weight` is relative to the graph it was computed on,
    so compare `conductance`, not `weight`, between years.
    """
    if not G.is_multigraph():
        raise TypeError("Expected a NetworkX MultiGraph, such as net['G'] from extract_network().")

    raw_weights = {}
    for u, v, key, data in G.edges(keys=True, data=True):
        if 'width' not in data or 'length' not in data:
            raise KeyError(f"Edge ({u}, {v}, {key}) must contain 'width' and 'length' attributes.")
        width, length = float(data['width']), float(data['length'])
        if not math.isfinite(width) or width <= 0:
            raise ValueError(f"Edge ({u}, {v}, {key}) has invalid width: {data['width']!r}")
        if not math.isfinite(length) or length <= 0:
            raise ValueError(f"Edge ({u}, {v}, {key}) has invalid length: {data['length']!r}")
        raw_weights[(u, v, key)] = width / length

    if raw_weights:
        max_weight = max(raw_weights.values())
        for (u, v, key), raw_weight in raw_weights.items():
            G[u][v][key]['conductance'] = raw_weight
            G[u][v][key]['weight'] = raw_weight / max_weight
    return G


def _can_contract(H, u, v, terminals, min_length=None):
    """False for inlet/outlet links, and for links whose removal would delete a real loop."""
    if u in terminals or v in terminals:
        return False
    between = H.get_edge_data(u, v)
    if len(between) == 1:
        return True
    return min_length is not None and all(d['length'] < min_length for d in between.values())


def _contract(H, keep, drop, pos=None):
    """Merge `drop` into `keep` in place. `pos` is the merged position (default: their midpoint)."""
    if pos is None:
        pos = tuple((np.array(H.nodes[keep]['pos'], float) + np.array(H.nodes[drop]['pos'], float)) / 2)
    merged_from = H.nodes[keep].get('merged_from', []) + [drop] + H.nodes[drop].get('merged_from', [])
    nx.contracted_nodes(H, keep, drop, self_loops=False, copy=False)
    H.nodes[keep].pop('contraction', None)          # nested dict added by networkx; breaks GraphML export
    H.nodes[keep].update(pos=pos, merged=True, merged_from=merged_from)
    return pos


def remove_self_edges(G):
    """Remove every link that starts and ends at the same node, in place. Returns how many were removed."""
    self_edges = list(nx.selfloop_edges(G, keys=True))
    G.remove_edges_from(self_edges)
    return len(self_edges)


def simplify_network(G, inlets, outlets, factor=2, max_contractions=7, min_nodes=None, min_length=120):
    """Contract redundant links. G is never modified.

    Pass 1: repeatedly contract the highest-weight allowed link above
            median + factor * 1.4826 * MAD, at most max_contractions times.
            These are short, wide links: one junction split in two inside a wide channel.
    Pass 2: contract every allowed link shorter than min_length metres, shortest first.
            120 m is 4 pixels at 30 m; scale it if the pixel size changes.

    Both passes stop at min_nodes (default: 90% of the original node count).
    Returns the simplified graph and one history table per pass.
    H.graph['self_edges_removed'] counts the self-edges taken out before contraction.
    """
    H = G.copy()
    terminals = set(inlets) | set(outlets)
    missing = terminals - set(H)
    if missing:
        raise ValueError(f'Inlets/outlets not in the graph: {sorted(missing)}')
    if min_nodes is None:
        min_nodes = max(1, int(0.9 * H.number_of_nodes()))
    H.graph['self_edges_removed'] = remove_self_edges(H)
    for u, v, k, d in H.edges(keys=True, data=True):
        d.setdefault('link_id', k)                  # contraction renumbers edge keys

    history = []
    while len(history) < max_contractions and H.number_of_nodes() > min_nodes:
        edges = [(u, v, d) for u, v, d in H.edges(data=True) if np.isfinite(d.get('weight', np.nan))]
        if not edges:
            break
        weights = np.array([d['weight'] for *_, d in edges])
        median = np.median(weights)
        threshold = median + factor * 1.4826 * np.median(np.abs(weights - median))
        candidates = [e for e in edges if e[2]['weight'] > threshold and _can_contract(H, e[0], e[1], terminals)]
        if not candidates:
            break
        u, v, d = max(candidates, key=lambda e: e[2]['weight'])
        pos = _contract(H, u, v)
        history.append({'node_1': u, 'node_2': v, 'pos': pos, 'weight': d['weight'], 'threshold': threshold,
                        'length': d['length'], 'width': d['width'], 'n_candidates': len(candidates)})

    spatial_history = []
    while H.number_of_nodes() > min_nodes:
        short = [(u, v, d) for u, v, d in H.edges(data=True)
                 if d['length'] < min_length and _can_contract(H, u, v, terminals, min_length)]
        if not short:
            break
        u, v, d = min(short, key=lambda e: e[2]['length'])
        pos = _contract(H, u, v)
        spatial_history.append({'node_1': u, 'node_2': v, 'pos': pos, 'length': d['length'],
                                'width': d['width'], 'weight': d.get('weight')})

    assert nx.number_of_selfloops(H) == 0
    return H, pd.DataFrame(history), pd.DataFrame(spatial_history)


def simplify_directed(D, history_df, spatial_history_df):
    """Apply the same contractions, in the same order, to the directed graph."""
    H = D.copy()
    remove_self_edges(H)
    for table in (history_df, spatial_history_df):
        for row in table.itertuples():
            _contract(H, row.node_1, row.node_2, row.pos)
    return H


def count_cycles(G):
    return G.number_of_edges() - G.number_of_nodes() + nx.number_connected_components(G)


def to_unweighted(G, simple=False):
    """Same nodes and links with every edge attribute removed.
    simple=True gives an nx.Graph, which merges parallel links and so hides their loops."""
    H = nx.MultiGraph()
    H.add_nodes_from(G.nodes(data=True))
    H.add_edges_from((u, v, k) for u, v, k in G.edges(keys=True))
    return nx.Graph(H) if simple else H


def assign_baseline_weights(tree, G, rng):
    """Give a baseline tree real channel geometry, in place: every tree link gets the width and length
    of a different real link of G, chosen at random. Random-tree links join arbitrary nodes, so they have
    no geometry of their own; this keeps the year's width/length distribution and nothing else."""
    real = [(d['width'], d['length']) for *_, d in G.edges(data=True)]
    picks = rng.choice(len(real), tree.number_of_edges(), replace=False)
    for (u, v, k), i in zip(tree.edges(keys=True), picks):
        tree[u][v][k].update(width=real[i][0], length=real[i][1])
    return add_geometric_weights(tree)


def load_graphs(path, weighted=True, baseline=False):
    """Read the saved {'G1': graph, ...} from all_graphs.pkl. weighted=False strips the weights;
    baseline=True returns the B1..B10 trees instead."""
    with open(path, 'rb') as fh:
        graphs = pickle.load(fh)['B' if baseline else 'G']
    return graphs if weighted else {name: to_unweighted(g) for name, g in graphs.items()}


def build_simplified_network(src, out, grid=None, shoreline_path=None, **simplify_kwargs):
    """Image -> cleaned RivGraph graph -> weights -> simplified graph -> weights renormalised."""
    net = extract_network(src, out, grid=grid, shoreline_path=shoreline_path)
    rn = net['rn']
    inlets, outlets = list(rn.nodes['inlets']), list(rn.nodes['outlets'])

    G_raw = add_geometric_weights(net['G'])
    G, history_df, spatial_history_df = simplify_network(G_raw, inlets, outlets, **simplify_kwargs)
    add_geometric_weights(G)                        # the old maximum may have been contracted away
    D = simplify_directed(net['D'], history_df, spatial_history_df)

    T = net['transform']
    H, W = rn.imshape
    POS = {n: ((d['pos'][0] - T.c) / 1e3, (d['pos'][1] - (T.f + H * T.e)) / 1e3) for n, d in G.nodes(data=True)}
    return {**net, 'G_raw': G_raw, 'D_raw': net['D'], 'G': G, 'D': D, 'inlets': inlets, 'outlets': outlets, 'POS': POS,
            'EXT': [0, W * T.a / 1e3, 0, H * (-T.e) / 1e3],
            'history_df': history_df, 'spatial_history_df': spatial_history_df, 'self_edges_removed': G.graph['self_edges_removed']}


def plot_before_after(net, path):
    """Node-edge drawing of the graph before and after simplification; merged nodes are ringed in red."""
    T = net['transform']
    H = net['rn'].imshape[0]
    fig, axs = plt.subplots(1, 2, figsize=(16, 7))
    for ax, g, title in zip(axs, (net['G_raw'], net['G']), ('before simplification', 'after simplification')):
        pos = {n: ((d['pos'][0] - T.c) / 1e3, (d['pos'][1] - (T.f + H * T.e)) / 1e3) for n, d in g.nodes(data=True)}
        nx.draw_networkx_edges(nx.Graph(g), pos, ax=ax, edge_color='#555', width=0.9)
        nx.draw_networkx_nodes(g, pos, ax=ax, node_size=[26 if g.nodes[n]['kind'] != 'junction' else 8 for n in g],
                               node_color=[KIND_COL[g.nodes[n]['kind']] for n in g],
                               edgecolors=['red' if g.nodes[n].get('merged') else 'k' for n in g],
                               linewidths=[1.5 if g.nodes[n].get('merged') else 0.2 for n in g])
        ax.set_title(f'{title}: {g.number_of_nodes()} nodes, {g.number_of_edges()} links, {count_cycles(g)} cycles')
        ax.set_aspect('equal')
        ax.axis('off')
    plt.tight_layout()
    plt.savefig(path, dpi=130)
    plt.close(fig)


if __name__ == '__main__':
    src = sys.argv[1] if len(sys.argv) > 1 else '18th Jan 2016.tiff'
    out = sys.argv[2] if len(sys.argv) > 2 else 'example_2016_simplified'
    net = build_simplified_network(src, out)
    G_raw, G = net['G_raw'], net['G']
    print(f"real network:       {G_raw.number_of_nodes()} nodes, {G_raw.number_of_edges()} links, {count_cycles(G_raw)} cycles")
    print(f"simplified network: {G.number_of_nodes()} nodes, {G.number_of_edges()} links, {count_cycles(G)} cycles, "
          f"{len(net['inlets'])} inlets, {len(net['outlets'])} outlets")
    print(f"contractions: {len(net['history_df'])} by weight, {len(net['spatial_history_df'])} by length")
    plot_before_after(net, f'{out}/before_after.png')
    print(f'figure saved to {out}/before_after.png')
