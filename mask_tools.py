"""Turn a Sentinel-2 image into the binary, georeferenced GeoTIFF mask RivGraph needs."""
import numpy as np
import rasterio
from PIL import Image
from pyproj import Transformer
from rasterio.transform import from_origin
from skimage import filters, morphology


def water_mask_from_rgb(path, top=45, bottom=45, min_size=20, threshold=None):
    """Water mask (uint8, 1 = water) from a true-colour image.

    In this scene the turbid channels are brighter in red than the dark-green
    mangrove, so the red channel is thresholded (Otsu unless `threshold` is given).
    `top` / `bottom` crop the Copernicus Browser banner and credit strip.
    Returns (mask, valid, threshold); `valid` marks real data pixels.
    """
    rgb = np.array(Image.open(path).convert('RGB')).astype(float)
    valid = rgb.sum(axis=2) > 15  # drops the black no-data wedge
    valid[:top] = False
    valid[rgb.shape[0] - bottom:] = False
    valid = morphology.binary_erosion(valid, morphology.disk(4))

    red = filters.gaussian(rgb[..., 0], 0.7)
    if threshold is None:
        threshold = filters.threshold_otsu(red[valid])
    water = (red > threshold) & valid
    water = morphology.remove_small_objects(water, min_size)
    return water.astype(np.uint8), valid, threshold


def water_mask_from_ndwi_render(path, out_path, min_size=20):
    """Binary water GeoTIFF from a Copernicus Browser NDWI render (3-band uint8).

    The render is blue over water and green over land, so blue minus green is
    thresholded with Otsu (per image, since each date's colour balance differs).
    The output keeps the source CRS and transform. Returns (mask, threshold).
    """
    with rasterio.open(path) as src:
        rgb = src.read().astype(float)
        profile = src.profile
    score = filters.gaussian(rgb[2] - rgb[1], 0.7)
    threshold = filters.threshold_otsu(score)
    water = score > threshold
    water = morphology.remove_small_objects(water, min_size)
    water = morphology.remove_small_holes(water, min_size).astype(np.uint8)

    profile.update(count=1, dtype='uint8', nodata=None, photometric=None)
    with rasterio.open(out_path, 'w', **profile) as dst:
        dst.write(water, 1)
    return water, threshold


def save_geotiff(mask, out_path, center_lonlat, m_per_px, epsg=32645):
    """Write `mask` as a single-band GeoTIFF in UTM 45N, centred on `center_lonlat`.

    A screenshot carries no georeferencing, so the position is approximate
    (centre from the Copernicus Browser link, pixel size from the scale bar).
    """
    cx, cy = Transformer.from_crs(4326, epsg, always_xy=True).transform(*center_lonlat)
    h, w = mask.shape
    transform = from_origin(cx - w * m_per_px / 2, cy + h * m_per_px / 2, m_per_px, m_per_px)
    with rasterio.open(out_path, 'w', driver='GTiff', height=h, width=w, count=1,
                       dtype='uint8', crs=f'EPSG:{epsg}', transform=transform) as dst:
        dst.write(mask, 1)
    return out_path
