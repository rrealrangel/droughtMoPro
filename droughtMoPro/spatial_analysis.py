#!/usr/bin/env python2
# -*- coding: utf-8 -*-
__author__ = 'Roberto A. Real-Rangel (Institute of Engineering UNAM)'
__license__ = 'GNU General Public License version 3'

import numpy as np

import gdal
import ogr
import xarray as xr


def vector2array(layer, xmin, xmax, ymin, ymax, res, nodata):
    cols = int((xmax - xmin) / res) + 1
    rows = int((ymax - ymin) / res) + 1

    # Create the destination data source
    targetDS = (
        gdal.GetDriverByName('MEM').Create(
            '', cols, rows, gdal.GDT_Byte
            )
        )

    targetDS.SetGeoTransform((xmin, res, 0, ymax, 0, - res))
    band = targetDS.GetRasterBand(1)
    band.SetNoDataValue(nodata)

    # Rasterize
    err = gdal.RasterizeLayer(
        targetDS,
        [1],
        layer,
        burn_values=[1],
        options=['ALL_TOUCHED=TRUE']
        )

    if err != 0:
        print("error:", err)

    # Read as array
    array = band.ReadAsArray()
    array = array.astype(float)
    array[np.where(array == 0)] = np.nan
    return(np.flipud(array))


def mask_data(data, layer, mmap_layer):
    mask = vector2array(
        layer=layer,
        xmin=min(data.lon.values),
        xmax=max(data.lon.values),
        ymin=min(data.lat.values),
        ymax=max(data.lat.values),
        res=data.attrs['LatitudeResolution'],
        nodata=-32768
        )

    masked = data * mask
    layer_xmin, layer_xmax, layer_ymin, layer_ymax = mmap_layer.GetExtent()
    trimmed = masked.sel(
        {'lat': (masked.lat >= layer_ymin) & (masked.lat <= layer_ymax),
         'lon': (masked.lon >= layer_xmin) & (masked.lon <= layer_xmax)}
        )

    return(trimmed.load())


def trim_data(data, vmap, res, nodata):
    """
    Parameters:
        data : xarray.Dataset
        vmap : string
        res : float
        nodata : float

    Source:
        https://bit.ly/2HxeOng
    """
    # Open the data source and read in the extent
    source_ds = ogr.Open(vmap)
    source_layer = source_ds.GetLayer(0)
    xmin, xmax, ymin, ymax = source_layer.GetExtent()

    def round_mult(num, mult, to):
        if to == 'up':
            return(mult * round(float(num + mult) / mult))

        elif to == 'down':
            return(mult * round(float(num - mult) / mult))

    xmin = round_mult(num=xmin, mult=res, to='down')
    xmax = round_mult(num=xmax, mult=res, to='up')
    ymin = round_mult(num=ymin, mult=res, to='down')
    ymax = round_mult(num=ymax, mult=res, to='up')
    x_coords = np.arange((xmin + (res / 2)), xmax, res)
    y_coords = np.arange((ymin + (res / 2)), ymax, res)

    # Create the destination data source
    cols = int((xmax - xmin) / res)
    rows = int((ymax - ymin) / res)

    output_source = gdal.GetDriverByName('MEM').Create(
        '', cols, rows, gdal.GDT_Byte
        )

    output_source.SetGeoTransform((xmin, res, 0, ymax, 0, -res))
    output_band = output_source.GetRasterBand(1)
    output_band.SetNoDataValue(nodata)

    # Rasterize
    gdal.RasterizeLayer(output_source, [1], source_layer, burn_values=[1])

    mask = xr.DataArray(
        data=np.flipud(output_band.ReadAsArray()),
        coords={'lat': y_coords, 'lon': x_coords},
        dims=['lat', 'lon']
        )

    trimmed_data = data.loc[dict(lat=mask.lat, lon=mask.lon)]
    masked_data = trimmed_data.where(mask)
    return(masked_data)


def interpolate_idw(distances, values, power=2):
    nominator = 0
    denominator = 0

    for i in range(len(distances)):
        nominator += (values[i] / pow(distances[i], power))
        denominator += (1 / pow(distances[i], power))

    # Return NODATA if the denominator is zero
    if denominator > 0:
        return(nominator/denominator)

    else:
        return(np.nan)
