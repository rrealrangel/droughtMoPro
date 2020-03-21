#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 15:42:52 2019

@author: realrangel

Functions:
    * Open data_files
    * Open spatial filters
    * Apply spatial filters
    * Apply values filter
    * Aggregate by date
    * Export reports
"""
# TODO: Merge intensity_area() and magnitude_area() into one single
# function.
# TODO: (In doubt) Make a function (maybe in data_manager.py) to read
# the map_file and iterate between its features.

from pathlib2 import Path
import ogr
import xarray as xr

from droughtMoPro import data_manager as dmgr


def intensity_area(input_data, header):
    """
    Parameters
    ----------
    input_data
    header

    Returns
    -------
    """
    dint_cats = {
        'not_drought': input_data >= -0.5,
        'd0': (input_data <= -0.5) & (-0.8 < input_data),
        'd1': (input_data <= -0.8) & (-1.3 < input_data),
        'd2': (input_data <= -1.3) & (-1.6 < input_data),
        'd3': (input_data <= -1.6) & (-2 < input_data),
        'd4': (input_data <= -2)
        }
    input_data_sum = xr.Dataset(
        data_vars={
            'not_drought': dint_cats['not_drought'].sum(
                dim=['lat', 'lon'],
                skipna=True
                ),
            'd0': dint_cats['d0'].sum(dim=['lat', 'lon'], skipna=True),
            'd1': dint_cats['d1'].sum(dim=['lat', 'lon'], skipna=True),
            'd2': dint_cats['d2'].sum(dim=['lat', 'lon'], skipna=True),
            'd3': dint_cats['d3'].sum(dim=['lat', 'lon'], skipna=True),
            'd4': dint_cats['d4'].sum(dim=['lat', 'lon'], skipna=True),
            }
        )
    cells = float(input_data[{'time': -1}].notnull().sum().values)
    dint_area = (input_data_sum / cells) * 100
    dint_area = dint_area.to_dataframe().reindex(
        labels=header,
        axis=1
        )
    dint_area.index = [
        str(i.year)+str(i.month).zfill(2) for i in dint_area.index
        ]
    return(dint_area)


def magnitude_area(input_data, header):
    """
    Parameters
    ----------
    input_data
    header

    Returns
    -------
    """
    dmag_cats = {
        'not_drought': input_data > -1,
        'm1': (input_data <= -1) & (-3 < input_data),
        'm2': (input_data <= -3) & (-6 < input_data),
        'm3': (input_data <= -6) & (-9 < input_data),
        'm4': (input_data <= -9) & (-12 < input_data),
        'm5': (input_data <= -12)
        }
    not_drought = dmag_cats['not_drought'].sum(
        dim=['lat', 'lon'],
        skipna=True
        )
    input_data_sum = xr.Dataset(
        data_vars={
            'not_drought': not_drought,
            'm1': dmag_cats['m1'].sum(dim=['lat', 'lon'], skipna=True),
            'm2': dmag_cats['m2'].sum(dim=['lat', 'lon'], skipna=True),
            'm3': dmag_cats['m3'].sum(dim=['lat', 'lon'], skipna=True),
            'm4': dmag_cats['m4'].sum(dim=['lat', 'lon'], skipna=True),
            'm5': dmag_cats['m5'].sum(dim=['lat', 'lon'], skipna=True)
            }
        )
    cells = float(input_data[{'time': -1}].notnull().sum().values)
    dmag_area = (input_data_sum / cells) * 100
    dmag_area = dmag_area.to_dataframe().reindex(
        labels=header,
        axis=1
        )
    dmag_area.index = [
        str(i.year)+str(i.month).zfill(2) for i in dmag_area.index
        ]
    return(dmag_area)


def quantiles(input_data, header):
    """
    Parameters
    ----------
    input_data
    header

    Returns
    -------
    """
    quant = xr.Dataset(
        data_vars={
            'min': input_data.quantile(
                q=0.00,
                dim=['lat', 'lon']
                ).drop('quantile'),
            'p25': input_data.quantile(
                q=0.25,
                dim=['lat', 'lon']
                ).drop('quantile'),
            'p50': input_data.quantile(
                q=0.50,
                dim=['lat', 'lon']
                ).drop('quantile'),
            'p75': input_data.quantile(
                q=0.75,
                dim=['lat', 'lon']
                ).drop('quantile'),
            'max': input_data.quantile(
                q=1.00,
                dim=['lat', 'lon']
                ).drop('quantile')
            }
        )
    quant = quant.to_dataframe().reindex(
        labels=header,
        axis=1
        )
    quant.index = [
        str(i.year)+str(i.month).zfill(2) for i in quant.index
        ]
    return(quant)


def export_ts(data_files, map_files, nodata, output_dir):
    """
    Parameters
    ----------
    data_files : list
    map_files : list
    nodata : int
    output_dir : str

    Returns
    -------
    """
    print("    - Importing data. This may take a while. Please wait.")
    data = xr.open_mfdataset(
        paths=data_files,
        concat_dim='time'
        )
    dint_header = ['not_drought', 'd0', 'd1', 'd2', 'd3', 'd4']
    dmag_header = ['not_drought', 'm1', 'm2', 'm3', 'm4', 'm5']
    quantile_header = ['min', 'p25', 'p50', 'p75', 'max']

    for map_file in map_files:
        map_ = ogr.Open(str(map_file))
        map_layer = map_.GetLayer(0)
        feature = map_layer.GetFeature(0)
        map_field = feature.GetFieldDefnRef(0)
        map_layer.ResetReading()
        output_subdir = (
            Path(output_dir) / '/'.join(map_file.parts[-2:]).split('.')[-2]
            )
        output_subdir.mkdir(parents=True, exist_ok=True)

        for mf, map_feature in enumerate(map_layer):
            # Define an output in-memory OGR dataset
            region = ogr.GetDriverByName('Memory').CreateDataSource('out')
            region_layer = region.CreateLayer('', geom_type=ogr.wkbPolygon)

            # Apply the map_field definition from original to output
            region_layer.CreateField(map_field)
            region_feature_definition = region_layer.GetLayerDefn()

            # For each feature, get the geometry
            map_feature_geometry = map_feature.GetGeometryRef()

            # Create an output feature
            region_feature_geometry = ogr.Feature(region_feature_definition)
            region_feature_geometry.SetGeometry(map_feature_geometry)
            region_layer.CreateFeature(region_feature_geometry)

            # Create a mask
            res_y = data.attrs['LatitudeResolution']
            res_x = data.attrs['LongitudeResolution']
            mask = dmgr.vector2array(
                res_x=res_x,
                res_y=res_y,
                xmin=min(data.lon.values) - (res_x / 2),
                ymax=max(data.lat.values) + (res_y / 2),
                cols=len(data.lon.values),
                rows=len(data.lat.values),
                layer=region_layer,
                nodata=nodata
                )

            # Mask and trim the input data
            data_masked = data * mask
            notnulls_cols = (
                data_masked.Drought_intensity.isel(
                    time=-1
                    ).notnull().sum('lat') > 0
                ).values
            notnulls_rows = (
                data_masked.Drought_intensity.isel(
                    time=-1
                    ).notnull().sum('lon') > 0
                ).values
            data_trimmed = data_masked[{
                'lat': notnulls_rows,
                'lon': notnulls_cols
                }].load()

            # ---- Drought intensity: area fraction ----
            dint_area = intensity_area(
                input_data=data_trimmed.Drought_intensity,
                header=dint_header
                )
            dint_area_fname = output_subdir / (
                '{id01}_{id02}_{sdi}_dint_area.csv'.format(
                    id01=map_feature.GetField('ID_01'),
                    id02=map_feature.GetField('ID_02'),
                    sdi=data.attrs['Title'].lower()
                    )
                )
            dint_area.to_csv(dint_area_fname)

            # ---- Drought intensity: quantiles ----
            dint_quant = quantiles(
                input_data=data_trimmed.Drought_intensity,
                header=quantile_header
                )
            dint_quant_fname = output_subdir / (
                '{id01}_{id02}_{sdi}_dint_quant.csv'.format(
                    id01=map_feature.GetField('ID_01'),
                    id02=map_feature.GetField('ID_02'),
                    sdi=data.attrs['Title'].lower()
                    )
                )
            dint_quant.to_csv(dint_quant_fname)

            # ---- Droght magnitude: area fraction ----
            if 'Drought_magnitude' in data_trimmed.data_vars.keys():
                dmag_area = magnitude_area(
                    input_data=data_trimmed.Drought_magnitude,
                    header=dmag_header
                    )
                dmag_area_fname = output_subdir / (
                    '{id01}_{id02}_{sdi}_dmag_area.csv'.format(
                        id01=map_feature.GetField('ID_01'),
                        id02=map_feature.GetField('ID_02'),
                        sdi=data.attrs['Title'].lower()
                        )
                    )
                dmag_area.to_csv(dmag_area_fname)

                # ---- Drought magnitude: quantiles ----
                dmag_quant = quantiles(
                    input_data=data_trimmed.Drought_magnitude,
                    header=quantile_header
                    )
                dmag_quant_fname = output_subdir / (
                    '{id01}_{id02}_{sdi}_dmag_quant.csv'.format(
                        id01=map_feature.GetField('ID_01'),
                        id02=map_feature.GetField('ID_02'),
                        sdi=data.attrs['Title'].lower()
                        )
                    )
                dmag_quant.to_csv(dmag_quant_fname)

            # Progress message.
            dmgr.progress_message(
                current=mf + 1,
                total=len(map_layer),
                message="- Processing '{}'".format(map_file.stem),
                units='feature'
                )

        # Clear things up
        region_feature_geometry.Destroy
        map_feature_geometry.Destroy

        # Reset the output layer to the 0th geometry
        region_layer.ResetReading()
