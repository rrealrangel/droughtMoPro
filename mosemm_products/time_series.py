#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 15:42:52 2019

@author: realrangel

Functions:
    * Open datasets
    * Open spatial filters
    * Apply spatial filters
    * Apply values filter
    * Aggregate by date
    * Export reports
"""
#import sys

from pathlib2 import Path
import gdal
import numpy as np
import ogr
import pandas as pd
import xarray as xr

import mosemm_products.data_manager as dmgr


#def progress_message(current, total, message="- Processing", units=None):
#    """Issue a messages of the progress of the process.
#
#    Generates a progress bar in terminal. It works within a for loop,
#    computing the progress percentage based on the current item
#    number and the total length of the sequence of item to iterate.
#
#    Parameters:
#        current : integer
#            The last item number computed within the for loop. This
#            could be obtained using enumerate() in when calling the for
#            loop.
#        total : integer
#            The total length of the sequence for which the for loop is
#            performing the iterations.
#        message : string (optional; default = "- Processing")
#    """
#    if units is not None:
#        progress = float(current)/total
#        sys.stdout.write(
#            "\r    {} ({:.1f} % of {} processed)".format(
#                message, progress * 100, units
#                )
#            )
#
#    else:
#        progress = float(current)/total
#        sys.stdout.write(
#            "\r    {} ({:.1f} % processed)".format(
#                message, progress * 100
#                )
#            )
#
#    if progress < 1:
#        sys.stdout.flush()
#
#    else:
#        sys.stdout.write('\n')
#
#
#def list_files(parent_dir, pattern):
#    """List all files in a directory with a specified pattern.
#
#    Parameters
#        parent_dir: string
#            Full path of the directory of which the files are to be listed.
#        pattern: string or list of strings
#            Pattern(s) of the files to be listed.
#    """
#    parent_dir = Path(parent_dir)
#    files_list = []
#
#    if isinstance(pattern, str):
#        pattern = [pattern]
#
#    for patt in pattern:
#        files_list.extend(parent_dir.glob(pattern=patt))
#
#    return(files_list)


#def vector2array(layer, data, nodata):
#    res = data.attrs['LatitudeResolution']
#    xmin = min(data.lon.values) - (res / 2)  # Xmin edge (not centroid)
#    ymax = max(data.lat.values) + (res / 2)  # Ymax edge (not centroid)
#    cols = len(data.lon)
#    rows = len(data.lat)
#
#    # Create the destination data source
#    raster_datasource = (gdal.GetDriverByName('MEM').Create(
#        '', cols, rows, gdal.GDT_Byte
#        ))
#
#    raster_datasource.SetGeoTransform((xmin, res, 0, ymax, 0, -res))
#    raster_band = raster_datasource.GetRasterBand(1)
#    raster_band.SetNoDataValue(nodata)
#
#    # Rasterize
#    err = gdal.RasterizeLayer(
#        raster_datasource,
#        [1],
#        layer,
#        burn_values=[1],
#        options=['ALL_TOUCHED=TRUE']
#        )
#
#    if err != 0:
#        print("error:", err)
#
#    # Read as array
#    mask_array = raster_band.ReadAsArray()
#    mask_array = mask_array.astype(float)
#    mask_array[np.where(mask_array == 0)] = np.nan
#    return(np.flipud(mask_array))


#def groupby_intensity(data, nodata):
#    # data[np.where(data == nodata)] = np.nan
#    categories = {
#        'd4': data <= -2,
#        'd3': (-2 < data) & (data <= -1.6),
#        'd2': (-1.6 < data) & (data <= -1.3),
#        'd1': (-1.3 < data) & (data <= -0.8),
#        'd0': (-0.8 < data) & (data <= -0.5),
#        'normal': (-0.5 < data) & (data < 0.5),
#        'w0': (0.5 <= data) & (data < 0.8),
#        'w1': (0.8 <= data) & (data < 1.3),
#        'w2': (1.3 <= data) & (data < 1.6),
#        'w3': (1.6 <= data) & (data < 2),
#        'w4': data >= 2
#        }
#
#    return(categories)
#
#
#def groupby_magnitude(magnitude, nodata):
#    # magnitude[np.where(magnitude == nodata)] = np.nan
#    categories = {
#        'not_drought': magnitude < 1,
#        'm1': (1 <= magnitude) & (magnitude < 3),
#        'm2': (3 <= magnitude) & (magnitude < 6),
#        'm3': (6 <= magnitude) & (magnitude < 9),
#        'm4': (9 <= magnitude) & (magnitude < 12),
#        'm5': magnitude >= 12
#        }
#
#    return(categories)
#
#
#def magnitude_area(categories, mask):
#    not_drought = categories['not_drought'].sum(
#        dim=['lat', 'lon'],
#        skipna=True
#        )
#
#    m1 = categories['m1'].sum(dim=['lat', 'lon'], skipna=True)
#    m2 = categories['m2'].sum(dim=['lat', 'lon'], skipna=True)
#    m3 = categories['m3'].sum(dim=['lat', 'lon'], skipna=True)
#    m4 = categories['m4'].sum(dim=['lat', 'lon'], skipna=True)
#    m5 = categories['m5'].sum(dim=['lat', 'lon'], skipna=True)
#    data_sum = xr.Dataset(
#        data_vars={
#            'not_drought': not_drought,
#            'm1': m1,
#            'm2': m2,
#            'm3': m3,
#            'm4': m4,
#            'm5': m5
#            }
#        )
#
#    return((data_sum / np.nansum(mask)) * 100)
#
#
#
#config = {
#    'datasets_dir': 'C:/Users/rreal/OneDrive/datasets/pysdi',
#    'vmaps_dir': (
#        'C:/Users/rreal/Mega/projects/multiannual/0000-inv-001_pysdi/main'
#        '/pysdi/regions'),
#    'datasets_patt': ['**/*PRESMORUN-01*.nc', '**/*PRESMORUN-01*.nc4'],
#    'vmaps_patt': '**/*.shp',
#    'output_dir': (
#        'C:/Users/rreal/OneDrive/datasets/conagua_mosemm/v2/time_series'
#        )
#    }

#datasets = list_files(
#    parent_dir=config['datasets_dir'],
#    pattern=config['datasets_patt']
#    )
#
#vmaps = list_files(
#    parent_dir=config['vmaps_dir'],
#    pattern=config['vmaps_patt']
#    )
#
#data = xr.open_mfdataset(paths=datasets, concat_dim='time')


def make_mask(array_reference, feature, layer_definition, nodata=-32768):
    # Create a vector region map stored in memory.
    region_datasource = (
        ogr.GetDriverByName('Memory').CreateDataSource(utf8_path='out')
        )
    region_layer = region_datasource.CreateLayer(
        name='',
        geom_type=ogr.wkbPolygon
        )
    region_layer.CreateField(field_def=layer_definition)
    region_layer_definition = region_layer.GetLayerDefn()
    region_feature = ogr.Feature(region_layer_definition)
    region_geometry = feature.GetGeometryRef()
    region_feature.SetGeometry(region_geometry)
    region_layer.CreateFeature(region_feature)

    # Create a raster region map stored in memory
    res = list(set(array_reference.lat.diff(dim='lat').values))[0]
    xmin = min(array_reference.lon.values) - (res / 2)  # Xmin edge
    ymax = max(array_reference.lat.values) + (res / 2)  # Ymax edge
    cols = len(array_reference.lon)
    rows = len(array_reference.lat)
    raster_datasource = (gdal.GetDriverByName('MEM').Create(
        '', cols, rows, gdal.GDT_Byte
        ))
    raster_datasource.SetGeoTransform([
        xmin,  # X coord of the upper left corner of the upper left pixel
        res,  # pixel width
        0,
        ymax,  # Y coord of the upper left corner of the upper left pixel
        0,
        -res  # pixel height
        ])
    raster_band = raster_datasource.GetRasterBand(1)
    raster_band.SetNoDataValue(nodata)

    # Rasterize
    gdal.RasterizeLayer(
        dataset=raster_datasource,
        bands=[1],
        layer=region_layer,
        burn_values=[1],
        options=['ALL_TOUCHED=TRUE']
        )

    # Read as array
    mask_array = raster_band.ReadAsArray().astype(bool)
    array_reference.values = np.flipud(mask_array)
    array_reference.name = 'mask'
    mask = array_reference.loc[{
        'lon': array_reference.sum(dim='lat') > 0,
        'lat': array_reference.sum(dim='lon') > 0
        }]
    region_datasource = None
    raster_datasource = None
    region_feature.Destroy
    region_geometry.Destroy
    return(mask)


def export_dint_area(data, output_file):
    categories = {
        'd4': data <= -2,
        'd3': (-2 < data) & (data <= -1.6),
        'd2': (-1.6 < data) & (data <= -1.3),
        'd1': (-1.3 < data) & (data <= -0.8),
        'd0': (-0.8 < data) & (data <= -0.5),
        'normal': (-0.5 < data) & (data < 0.5),
        'w0': (0.5 <= data) & (data < 0.8),
        'w1': (0.8 <= data) & (data < 1.3),
        'w2': (1.3 <= data) & (data < 1.6),
        'w3': (1.6 <= data) & (data < 2),
        'w4': data >= 2
        }
    d4 = categories['d4'].sum(dim=['lat', 'lon'], skipna=True)
    d3 = categories['d3'].sum(dim=['lat', 'lon'], skipna=True)
    d2 = categories['d2'].sum(dim=['lat', 'lon'], skipna=True)
    d1 = categories['d1'].sum(dim=['lat', 'lon'], skipna=True)
    d0 = categories['d0'].sum(dim=['lat', 'lon'], skipna=True)
    normal = categories['normal'].sum(dim=['lat', 'lon'], skipna=True)
    w0 = categories['w0'].sum(dim=['lat', 'lon'], skipna=True)
    w1 = categories['w1'].sum(dim=['lat', 'lon'], skipna=True)
    w2 = categories['w2'].sum(dim=['lat', 'lon'], skipna=True)
    w3 = categories['w3'].sum(dim=['lat', 'lon'], skipna=True)
    w4 = categories['w4'].sum(dim=['lat', 'lon'], skipna=True)
    data_sum = xr.Dataset(
        data_vars={
            'd4': d4,
            'd3': d3,
            'd2': d2,
            'd1': d1,
            'd0': d0,
            'normal': normal,
            'w0': w0,
            'w1': w1,
            'w2': w2,
            'w3': w3,
            'w4': w4
            }
        )
    cells = float(
        data[{'time': -1}].notnull().sum().values
        )
    dint_area = (data_sum / cells) * 100
    fields = [
        'd4', 'd3', 'd2', 'd1', 'd0', 'normal', 'w0', 'w1', 'w2', 'w3',
        'w4'
        ]
    dint_area = dint_area.to_dataframe().reindex(fields, axis=1)
    dint_area.index = [
        str(i.year) + str(i.month).zfill(2) for i in dint_area.index
        ]
    dint_area.to_csv(output_file)


def export_quant(data, output_file):
    data.load()
    p0 = data.quantile(q=0.00, dim=['lat', 'lon']).drop('quantile')
    p25 = data.quantile(q=0.25, dim=['lat', 'lon']).drop('quantile')
    p50 = data.quantile(q=0.50, dim=['lat', 'lon']).drop('quantile')
    p75 = data.quantile(q=0.75, dim=['lat', 'lon']).drop('quantile')
    p100 = data.quantile(q=1.00, dim=['lat', 'lon']).drop('quantile')
    quant = xr.Dataset(
        data_vars={
            'p0': p0,
            'p25': p25,
            'p50': p50,
            'p75': p75,
            'p100': p100
            }
        )
    fields = ['p0', 'p25', 'p50', 'p75', 'p100']
    quant = quant.to_dataframe().reindex(fields, axis=1)
    quant.index = [
        str(i.year)+str(i.month).zfill(2) for i in quant.index
        ]
    quant.to_csv(output_file)


def export_dmag_area(data, output_file):
    categories = {
        'not_drought': data > -1,
        'm1': (data <= -1) & (-3 < data),
        'm2': (data <= -3) & (-6 < data),
        'm3': (data <= -6) & (-9 < data),
        'm4': (data <= -9) & (-12 < data),
        'm5': data <= -12
        }
    not_drought = categories['not_drought'].sum(
        dim=['lat', 'lon'],
        skipna=True
        )
    m1 = categories['m1'].sum(dim=['lat', 'lon'], skipna=True)
    m2 = categories['m2'].sum(dim=['lat', 'lon'], skipna=True)
    m3 = categories['m3'].sum(dim=['lat', 'lon'], skipna=True)
    m4 = categories['m4'].sum(dim=['lat', 'lon'], skipna=True)
    m5 = categories['m5'].sum(dim=['lat', 'lon'], skipna=True)
    data_sum = xr.Dataset(
        data_vars={
            'not_drought': not_drought,
            'm1': m1,
            'm2': m2,
            'm3': m3,
            'm4': m4,
            'm5': m5
            }
        )
    cells = float(
        data[{'time': -1}].notnull().sum().values
        )
    dmag_area = (data_sum / cells) * 100
    fields = ['not_drought', 'm1', 'm2', 'm3', 'm4', 'm5']
    dmag_area = dmag_area.to_dataframe().reindex(fields, axis=1)
    dmag_area.index = [
        str(i.year) + str(i.month).zfill(2) for i in dmag_area.index
        ]
    dmag_area.to_csv(output_file)


def export_ts(data_files, map_file, nodata, output_dir):
    # Open the input datasets.
    data = xr.open_mfdataset(
        paths=data_files,
        chunks=None,
        concat_dim='time',
        )

    # Open the vector map that contains the masking regions.
    map_datasource = ogr.Open(utf8_path=str(map_file))
    map_layer = map_datasource.GetLayer(0)
    map_auxiliar_feature = map_layer.GetFeature(0)
    map_auxiliar_field = map_auxiliar_feature.GetFieldDefnRef(0)
    map_layer.ResetReading()
    REGIONS_IN_VMAP = map_layer.GetFeatureCount()

    # Create the results subdirectory.
    output_subdir = (
        Path(output_dir) / '/'.join(map_file.parts[-2:]).split('.')[-2]
        )
    output_subdir.mkdir(parents=True, exist_ok=True)

    # Mask and analyze data with each region in the vector map.
    for mr, map_region in enumerate(map_layer):
        id01 = map_region.GetField('ID_01')
        id02 = map_region.GetField('ID_02')

        # Mask the input data
        region_mask = make_mask(
            array_reference=data.Drought_intensity[{'time': -1}].drop('time'),
            feature=map_region,
            layer_definition=map_auxiliar_field,
            nodata=nodata
            )
        data_masked = data.where(region_mask)

        # Export drought intensity: fraction of area.
        output_dint_area_file = output_subdir / 'dint_area_{}_{}.csv'.format(id01, id02)
        export_dint_area(
            data=data_masked.Drought_intensity,
            output_file=output_dint_area_file
            )

        # Export drought intensity: quantiles.
        output_dint_quant_file = output_subdir / 'dint_quant_{}_{}.csv'.format(id01, id02)
        export_quant(
            data=data_masked.Drought_intensity,
            output_file=output_dint_quant_file
            )

        # Export droght magnitude: fraction of area.
        output_dmag_area_file = output_subdir / 'dmag_area_{}_{}.csv'.format(id01, id02)
        export_dmag_area(
            data=data_masked.Drought_magnitude,
            output_file=output_dmag_area_file
            )

        # Export drought magnitude: quantiles.
        output_dmag_quant_file = output_subdir / 'dmag_quant_{}_{}.csv'.format(id01, id02)
        export_quant(
            data=data_masked.Drought_magnitude,
            output_file=output_dmag_quant_file
            )

        dmgr.progress_message(
            current=(mr + 1),
            total=REGIONS_IN_VMAP,
            message="- Processing" '{}'".format(map_file.stem)",
            units='feature'
            )

    # Clear things up
    map_layer.ResetReading()

#def export_time_series(data_files, map_file):
#    map_datasource = ogr.Open(str(map_file))
#    map_layer = map_datasource.GetLayer(0)
#    feature = map_layer.GetFeature(0)
#    map_auxiliar_field = feature.GetFieldDefnRef(0)
#    map_layer.ResetReading()
#    REGIONS_IN_VMAP = map_layer.GetFeatureCount()
#
#    for mr, map_region in enumerate(map_layer):
#        # We define an output in-memory OGR dataset
#        region_datasource = ogr.GetDriverByName('Memory').CreateDataSource(
#            utf8_path='out'
#            )
#        region_layer = region_datasource.CreateLayer(
#            name='',
#            geom_type=ogr.wkbPolygon
#            )
#
#        # Apply the map_auxiliar_field definition from the original to the output
#        region_layer.CreateField(map_auxiliar_field)
#        region_layer_definition = region_layer.GetLayerDefn()
#
#        # For each feature, get the geometry
#        region_geometry = map_region.GetGeometryRef()
#
#        # Create an output feature
#        region_feature = ogr.Feature(region_layer_definition)
#        region_feature.SetGeometry(region_geometry)
#        region_layer.CreateFeature(region_feature)
#
#        # Mask the input data
#        data_masked = mask_data(layer=region_layer)
#
#        # Quantify drought intensity
#        dint_area = intensity_area(
#            data_masked=data_masked,
#            mask=vector2array(
#                layer=region_layer,
#                xmin=min(data.lon.values),
#                xmax=max(data.lon.values),
#                ymin=min(data.lat.values),
#                ymax=max(data.lat.values),
#                res=data.attrs['LatitudeResolution'],
#                nodata=-32768
#                )
#            )
#
#        dint_quant = quantiles(
#            data=data_masked.Drought_intensity,
#            nodata=-32768
#            )
#
#        # Quantify droght magnitude
#        categories = groupby_magnitude(
#            magnitude=data_masked.Drought_magnitude,
#            nodata=-32768
#            )
#
#        dmag_area = magnitude_area(
#            categories=categories,
#            mask=vector2array(
#                layer=region_layer,
#                xmin=min(data.lon.values),
#                xmax=max(data.lon.values),
#                ymin=min(data.lat.values),
#                ymax=max(data.lat.values),
#                res=data.attrs['LatitudeResolution'],
#                nodata=-32768
#                )
#            )
#
#        dmag_quant = quantiles(
#            data=data_masked.Drought_magnitude,
#            nodata=-32768
#            )
#
#        dint_cats_header = [
#            'd4', 'd3', 'd2', 'd1', 'd0', 'normal', 'w0', 'w1', 'w2', 'w3',
#            'w4']
#
#        categories = ['not_drought', 'm1', 'm2', 'm3']
#        fields = ['p0', 'p25', 'p50', 'p75', 'p100']
#        feature_id_01 = id01
#        feature_id_02 = id02
#        output_dir = (
#            Path(config['output_dir']) / '/'.join(
#                map_datasource.parts[-2:]
#                ).split('.')[-2]
#            )
#
#        output_dir.mkdir(parents=True, exist_ok=True)
#
#        # Export drought intensity areas.
#        output_file = output_dir / 'area_{}_{}.csv'.format(
#            feature_id_01, feature_id_02
#            )
#
#        dint_area = dint_area.to_dataframe().reindex(dint_cats_header, axis=1)
#        dint_area.index = [
#            str(i.year)+str(i.month).zfill(2) for i in dint_area.index]
#
#        dint_area.to_csv(output_file)
#
#        # Export drought intensity quantiles.
#        output_dint_quant_file = output_dir / 'dint_{}_{}.csv'.format(
#            feature_id_01, feature_id_02)
#
#        dint_quant = dint_quant.to_dataframe().reindex(fields, axis=1)
#        dint_quant.index = [
#            str(i.year)+str(i.month).zfill(2) for i in dint_quant.index]
#
#        dint_quant.to_csv(output_dint_quant_file)
#
#        # Export drought magnitude quantiles.
#        dmag_quant_fname = output_dir / 'dmag_{}_{}.csv'.format(
#            feature_id_01, feature_id_02)
#
#        dmag_quant = dmag_quant.to_dataframe().reindex(fields, axis=1)
#        dmag_quant.index = [
#            str(i.year)+str(i.month).zfill(2) for i in dmag_quant.index]
#
#        dmag_quant.to_csv(dmag_quant_fname)
#
#        dmgr.progress_message(
#            current=(mr + 1),
#            total=REGIONS_IN_VMAP,
#            message="- Processing '{}'".format(map_datasource.stem),
#            units='feature'
#            )
#
#    # Clear things up
#    region_feature.Destroy
#    region_geometry.Destroy
#
#    # Reset the output layer to the 0th geometry
#    region_layer.ResetReading()
