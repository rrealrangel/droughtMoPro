# -*- coding: utf-8 -*-
"""
Report generator - (proto)pySDI script

Author
------
    Roberto A. Real-Rangel (Institute of Engineering UNAM; Mexico)

License
-------
    GNU General Public License
"""
from pathlib2 import Path
import ogr
import pandas as pd
import xarray as xr

from droughtMoPro import data_manager as dmgr


def ratio_per_category(input_data, indicator):
    cells = float(input_data[{'time': -1}].notnull().sum().values)
    ratio = {}

    if indicator == 'Drought_intensity':
        count = {
            'D0': ((-0.5 >= input_data) & (input_data > -0.8)).sum(
                dim=['lat', 'lon'],
                skipna=True
                ).values[0],
            'D1': ((-0.8 >= input_data) & (input_data > -1.3)).sum(
                dim=['lat', 'lon'],
                skipna=True
                ).values[0],
            'D2': ((-1.3 >= input_data) & (input_data > -1.6)).sum(
                dim=['lat', 'lon'],
                skipna=True
                ).values[0],
            'D3': ((-1.6 >= input_data) & (input_data > -2.0)).sum(
                dim=['lat', 'lon'],
                skipna=True
                ).values[0],
            'D4': (-2.0 >= input_data).sum(
                dim=['lat', 'lon'],
                skipna=True
                ).values[0]
            }

    elif indicator == 'Drought_magnitude':
        count = {
            'M1': ((-3 < input_data) & (input_data <= -1)).sum(
                dim=['lat', 'lon'],
                skipna=True
                ).values[0],
            'M2': ((-6 < input_data) & (input_data <= -3)).sum(
                dim=['lat', 'lon'],
                skipna=True
                ).values[0],
            'M3': ((-9 < input_data) & (input_data <= -6)).sum(
                dim=['lat', 'lon'],
                skipna=True
                ).values[0],
            'M4': ((-12 < input_data) & (input_data <= -9)).sum(
                dim=['lat', 'lon'],
                skipna=True
                ).values[0],
            'M5': ((input_data <= -12)).sum(
                dim=['lat', 'lon'],
                skipna=True
                ).values[0]
            }

    for category, input_data_sum in count.iteritems():
        ratio[category] = (input_data_sum / cells) * 100

    return(ratio)


def make_report(data_files, map_files, output_dir, nodata=-32768):
    """
    Source:
    https://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html
    """
    data = xr.open_mfdataset(paths=data_files, concat_dim='time')

    for v, variable in data.data_vars.items():
        for map_file in map_files:
            map_ = ogr.Open(str(map_file))
            map_name = map_file.stem
            print("    - The map '{}' was loaded.".format(map_name))
            map_layer = map_.GetLayer(0)
            feature = map_layer.GetFeature(0)
            map_field = feature.GetFieldDefnRef(0)
            map_layer.ResetReading()
            print("    - Rasterizing features and exporting reports.")
            entries = []

            for mf, map_feature in enumerate(map_layer):
                # We define an output in-memory OGR dataset
                region = ogr.GetDriverByName('Memory').CreateDataSource('out')
                region_layer = region.CreateLayer('', geom_type=ogr.wkbPolygon)

                # Apply the map_field definition from original to output
                region_layer.CreateField(map_field)
                region_feature_definition = region_layer.GetLayerDefn()

                # For each feature, get the geometry
                map_feature_geometry = map_feature.GetGeometryRef()

                # Create an output feature
                region_feature_geometry = ogr.Feature(
                    region_feature_definition
                    )
                region_feature_geometry.SetGeometry(map_feature_geometry)
                region_layer.CreateFeature(region_feature_geometry)

                # Create a mask
                res = data.attrs['LatitudeResolution']
                mask = dmgr.vector2array(
                    res=res,
                    xmin=min(variable.lon.values) - (res / 2),
                    ymax=max(variable.lat.values) + (res / 2),
                    cols=len(variable.lon.values),
                    rows=len(variable.lat.values),
                    layer=region_layer,
                    nodata=nodata
                    )

                # Mask and trim the input data
                data_masked = variable * mask
                notnulls_cols = (
                    data_masked.isel(
                        time=-1
                        ).notnull().sum('lat') > 0
                    ).values
                notnulls_rows = (
                    data_masked.isel(
                        time=-1
                        ).notnull().sum('lon') > 0
                    ).values
                data_trimmed = data_masked[{
                    'lat': notnulls_rows,
                    'lon': notnulls_cols
                    }]

                # Compute entry values.
                ratio = ratio_per_category(
                    input_data=data_trimmed,
                    indicator=variable.name
                    )

                ratio['Entidad'] = map_feature.GetField('NOM_ENT')
                ratio['Nombre'] = map_feature.GetField('NOM_MUN')
                ratio['Clave'] = (
                    map_feature.GetField('NUM_ENT') +
                    map_feature.GetField('ID_02')
                    )
                entries.append(ratio)

                dmgr.progress_message(
                    current=(mf + 1),
                    total=len(map_layer),
                    message="- Processing",
                    units=None
                    )

            report = pd.DataFrame(entries)
            report = report[
                ['Clave', 'Nombre', 'Entidad'] +
                [c for c in report if c not in ['Clave', 'Nombre', 'Entidad']]
                ]
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            report.to_excel(
                excel_writer=(
                    output_dir + '/' + v + '.xlsx'
                    ),
                sheet_name=v,
                float_format="%.1f",
                index=False,
                engine='openpyxl'
                )
            region_layer.ResetReading()
