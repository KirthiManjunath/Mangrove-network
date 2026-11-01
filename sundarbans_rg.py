"""Image -> cleaned RivGraph network for the Sundarbans.

Shared by the notebooks so every date is processed with exactly the same rules.

    from sundarbans_rg import extract_network
    net = extract_network('18th Jan 2016.tiff', 'example_2016')
    net['rn']        # rivgraph delta object, pruned + flow directions assigned
    net['G']         # nx.MultiGraph  (undirected, link attrs: length, width, idx; node attr: pos)
    net['D']         # nx.MultiDiGraph (upstream -> downstream)
"""
import copy
import os

import geopandas as gpd
import networkx as nx
import numpy as np
import rasterio
from rasterio import features
from rasterio.warp import Resampling, reproject
from shapely.geometry import Point, box, shape
from shapely.affinity import translate
from shapely.ops import linemerge, unary_union
from skimage import morphology, transform

from rivgraph import geo_utils as gu
from rivgraph import ln_utils as lnu
from rivgraph.classes import delta
from rivgraph.deltas import delta_utils as du


def water_mask(rgb, min_px=30):
    """Copernicus NDWI visualisation (blue = water, green = land) -> boolean water mask."""
    R, G, B = rgb.astype(np.int16)
    water = (B - G) > 0
    water = morphology.remove_small_objects(water, min_px)
    return morphology.remove_small_holes(water, min_px)


def closed_land(water, r=80, f=4):
    """Close the land mask with an r-px disk (at 1/f resolution) so river mouths are bridged."""
    small = transform.downscale_local_mean((~water).astype(float), (f, f)) > 0.5
    rs = r // f
    pad = np.pad(small, rs, mode='reflect')            # reflect so the image edge isn't treated as land or sea
    cl = morphology.binary_closing(pad, morphology.disk(rs))[rs:-rs, rs:-rs]
    cl = morphology.remove_small_holes(cl, 3000)       # fill inland water bodies only
    cl = (transform.resize(cl.astype(float), water.shape, order=0) > 0.5) | ~water
    return morphology.remove_small_objects(cl, 20000)


def write_shoreline(land, T, crs, path, frame_m=90):
    """Boundary of the closed land, minus the image frame -> shoreline line file."""
    H, W = land.shape
    west, south, east, north = rasterio.transform.array_bounds(H, W, T)
    landpoly = unary_union([shape(g) for g, v in features.shapes(land.astype('uint8'), transform=T) if v == 1])
    inner = box(west + frame_m, south + frame_m, east - frame_m, north - frame_m)
    shore = landpoly.boundary.intersection(inner)
    if shore.geom_type == 'MultiLineString':
        shore = linemerge(shore)
    # shoreline vertices sit on pixel corners, which diagonal skeleton steps pass through exactly;
    # RivGraph's clip_by_shoreline crashes on the resulting zero-length intersections, so nudge by a fraction of a pixel
    shore = translate(shore, 0.137 * T.a, 0.211 * T.e)
    # RivGraph's sjoin looks for an 'id_left' column, which only appears if the shoreline also has an 'id'
    gpd.GeoDataFrame({'id': [0]}, geometry=[shore], crs=crs).to_file(path)
    return shore, landpoly


def write_inlets(rn, land, crs, path):
    """Inlets = channel ends on the top/left/right image edge, on the land side."""
    H, W = rn.imshape
    G = nx.MultiGraph()
    G.add_nodes_from(rn.nodes['id'])
    for a, b in rn.links['conn']:
        G.add_edge(a, b)
    r, c = np.unravel_index(rn.nodes['idx'], rn.imshape)
    on_edge = (r <= 2) | (c <= 2) | (c >= W - 3)
    pick = [i for i in range(len(r))
            if on_edge[i] and G.degree(rn.nodes['id'][i]) == 1 and land[r[i], c[i]]]
    x, y = gu.idx_to_coords([rn.nodes['idx'][i] for i in pick], rn.imshape, rn.gt)
    gpd.GeoDataFrame(geometry=[Point(a, b) for a, b in zip(x, y)], crs=crs).to_file(path)
    return len(pick)


def prune_sundarbans(rn, shoreline_path, inlets_path):
    """RivGraph's delta.prune_network, made to work when the shoreline splits the scene into
    several estuary systems: the stock version keeps the component of the *first* inlet and then
    crashes on the inlets that were dropped. Here every inlet's component is tried and the largest kept."""
    links, nodes = rn.links, du.find_inlet_nodes(rn.nodes, inlets_path, rn.imshape, rn.gt, rn.wkt)
    links, nodes = lnu.remove_all_spurs(links, nodes, dontremove=list(nodes['inlets']))
    all_inlets = [i for i in dict.fromkeys(nodes['inlets']) if i in nodes['id']]
    best, seen = None, set()
    for first in all_inlets:
        if first in seen:
            continue
        L, N = copy.deepcopy(links), copy.deepcopy(nodes)
        N['inlets'] = [first] + [i for i in all_inlets if i != first]
        L, N = du.clip_by_shoreline(L, N, shoreline_path, rn.imshape, rn.gt, rn.wkt)
        N['inlets'] = [i for i in N['inlets'] if i in N['id']]
        seen |= set(N['inlets']) | {first}
        if best is None or len(L['id']) > len(best[0]['id']):
            best = (L, N)
    links, nodes = best
    links, nodes = lnu.remove_all_spurs(links, nodes, dontremove=list(nodes['inlets']) + list(nodes['outlets']))
    links, nodes = lnu.remove_disconnected_bridge_links(links, nodes)
    links, nodes = lnu.remove_single_pixel_links(links, nodes)
    links, nodes = lnu.find_parallel_links(links, nodes)
    rn.links, rn.nodes = lnu.add_link_conn(links, nodes), nodes
    rn.finalize_ids()


def node_xy(rn):
    """Node id -> (x, y) map coordinates in metres."""
    r, c = np.unravel_index(rn.nodes['idx'], rn.imshape)
    x = rn.gt[0] + (c + 0.5) * rn.gt[1]
    y = rn.gt[3] + (r + 0.5) * rn.gt[5]
    return {nid: (x[i], y[i]) for i, nid in enumerate(rn.nodes['id'])}


def to_graphs(rn):
    """RivGraph links/nodes -> undirected MultiGraph and directed MultiDiGraph (same attributes)."""
    pos = node_xy(rn)
    G, D = nx.MultiGraph(), nx.MultiDiGraph()
    for i, lid in enumerate(rn.links['id']):
        a, b = rn.links['conn'][i]
        attrs = dict(length=rn.links['len'][i], width=rn.links['wid_adj'][i], idx=rn.links['idx'][i])
        G.add_edge(a, b, key=lid, **attrs)
        D.add_edge(a, b, key=lid, **attrs)
    for g in (G, D):
        nx.set_node_attributes(g, {n: pos[n] for n in g.nodes}, 'pos')
        nx.set_node_attributes(g, {n: 'inlet' for n in rn.nodes['inlets']}, 'kind')
        nx.set_node_attributes(g, {n: 'outlet' for n in rn.nodes['outlets']}, 'kind')
        nx.set_node_attributes(g, {n: 'junction' for n in g.nodes if 'kind' not in g.nodes[n]}, 'kind')
    return G, D


def common_grid(srcs):
    """Grid covering the overlap of all images, at the first image's pixel size."""
    with rasterio.open(srcs[0]) as s:
        xres, yres, crs = s.transform.a, -s.transform.e, s.crs
    bs = [rasterio.open(p).bounds for p in srcs]
    west, north = max(b.left for b in bs), min(b.top for b in bs)
    east, south = min(b.right for b in bs), max(b.bottom for b in bs)
    W, H = int((east - west) // xres), int((north - south) // yres)
    return dict(transform=rasterio.Affine(xres, 0, west, 0, -yres, north), width=W, height=H, crs=crs)


def read_rgb(src, grid=None):
    """Read an RGB image, optionally resampled (nearest) onto a common grid."""
    with rasterio.open(src) as s:
        prof = {k: v for k, v in s.profile.items() if k not in ('photometric', 'tiled', 'blockxsize', 'blockysize', 'interleave')}
        if grid is None:
            return s.read(), prof, s.transform, s.crs
        rgb = np.zeros((s.count, grid['height'], grid['width']), s.dtypes[0])
        reproject(s.read(), rgb, src_transform=s.transform, src_crs=s.crs,
                  dst_transform=grid['transform'], dst_crs=grid['crs'], resampling=Resampling.nearest)
    prof.update(transform=grid['transform'], width=grid['width'], height=grid['height'], crs=grid['crs'])
    return rgb, prof, grid['transform'], grid['crs']


def common_shoreline(srcs, grid, path, shore_disk=80):
    """One shoreline for all dates, from the water mask that at least half of the images agree on."""
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    votes = sum(water_mask(read_rgb(p, grid)[0]).astype(np.uint8) for p in srcs)
    water = votes >= len(srcs) / 2
    return write_shoreline(closed_land(water, r=shore_disk), grid['transform'], grid['crs'], path)


def extract_network(src, out, shore_disk=80, grid=None, shoreline_path=None):
    """Full pipeline for one image. Writes mask/shoreline/inlets/rivgraph files into `out`.
    grid: resample onto this common grid first. shoreline_path: use this shoreline instead of one from this image."""
    os.makedirs(out, exist_ok=True)
    rgb, prof, T, crs = read_rgb(src, grid)
    water = water_mask(rgb)
    with rasterio.open(f'{out}/mask.tif', 'w', **{**prof, 'count': 1, 'compress': 'lzw'}) as d:
        d.write(water.astype('uint8'), 1)

    land = closed_land(water, r=shore_disk)
    if shoreline_path is None:
        shore, landpoly = write_shoreline(land, T, crs, f'{out}/shoreline.gpkg')
    else:
        shore = gpd.read_file(shoreline_path).geometry[0]
        landpoly = None
        gpd.read_file(shoreline_path).to_file(f'{out}/shoreline.gpkg')

    rn = delta(os.path.basename(out), f'{out}/mask.tif', results_folder=f'{out}/rivgraph', verbose=False)
    rn.skeletonize()
    rn.compute_network()
    n_raw = len(rn.links['id'])
    write_inlets(rn, land, crs, f'{out}/inlets.gpkg')
    prune_sundarbans(rn, f'{out}/shoreline.gpkg', f'{out}/inlets.gpkg')
    rn.compute_link_width_and_length()
    rn.assign_flow_directions()
    rn.save_network()

    G, D = to_graphs(rn)
    return dict(rn=rn, G=G, D=D, rgb=rgb, water=water, land=land, shore=shore, landpoly=landpoly,
                transform=T, crs=crs, n_raw_links=n_raw)
