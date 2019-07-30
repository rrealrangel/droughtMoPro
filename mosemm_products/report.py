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
import numpy as np
import ogr
import pandas as pd

import lib.data_manager as dmgr
import lib.spatial_analysis as span


def ratio_per_category(input_array, indicator):
    cells_total = float(np.isfinite(input_array).sum())
    ratio = {}

    if indicator == 'intensity':
        count = {
            'D0': ((-0.8 < input_array) & (input_array <= -0.5)).sum(),
            'D1': ((-1.3 < input_array) & (input_array <= -0.8)).sum(),
            'D2': ((-1.6 < input_array) & (input_array <= -1.3)).sum(),
            'D3': ((-2.0 < input_array) & (input_array <= -1.6)).sum(),
            'D4': (input_array <= -2.0).sum()
            }

    elif indicator == 'magnitude':
        count = {
            'M1': ((-3 < input_array) & (input_array <= -1)).sum(),
            'M2': ((-6 < input_array) & (input_array <= -3)).sum(),
            'M3': ((-9 < input_array) & (input_array <= -6)).sum(),
            'M4': ((-12 < input_array) & (input_array <= -9)).sum(),
            'M5': ((input_array <= -12)).sum()
            }

    else:
        print(
            "Please, enter a valid severity indicator ('intensidad' or "
            "'magnitud' only)."
            )

    for category, cells in count.iteritems():
        ratio[category] = (cells / cells_total) * 100

    return(ratio)


def make_report(data, indicator, regions, fielddef, output_dir, nodata=-32768):
    """
    Source:
    https://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html
    """
    entries = []

    for f, feat in enumerate(regions):
        # We define an output in-memory OGR dataset
        drv = ogr.GetDriverByName('Memory')
        dst_ds = drv.CreateDataSource('out')
        dst_layer = dst_ds.CreateLayer('', geom_type=ogr.wkbPolygon)

        # Apply the field definition from the original to the output
        dst_layer.CreateField(fielddef)
        feature_defn = dst_layer.GetLayerDefn()

        # Get the geometry of each feature
        geom = feat.GetGeometryRef()

        # Create an output feature
        out_geom = ogr.Feature(feature_defn)
        out_geom.SetGeometry(geom)

        # Add the feature with its geometry to the output layer
        dst_layer.CreateFeature(out_geom)

        # Apply mask to input data
        mask = span.vector2array(
            layer=dst_layer,
            xmin=min(data.lon.values),
            xmax=max(data.lon.values),
            ymin=min(data.lat.values),
            ymax=max(data.lat.values),
            res=data.lat.values[1] - data.lat.values[0],
            nodata=nodata
            )

        input_array = (data.values * mask)[-1]

        # Compute entry values.
        ratio = ratio_per_category(
            input_array=input_array,
            indicator=indicator
            )

        ratio['Entidad'] = feat.GetField('NOM_ENT')
        ratio['Nombre'] = feat.GetField('NOM_MUN')
        ratio['Clave'] = feat.GetField('NUM_ENT') + feat.GetField('ID_02')
        entries.append(ratio)

        dmgr.progress_message(
            current=(f + 1),
            total=len(regions),
            message="- Processing",
            units=None
            )

    report = pd.DataFrame(entries)
    report = report[
        ['Clave', 'Nombre', 'Entidad'] +
        [c for c in report if c not in ['Clave', 'Nombre', 'Entidad']]
        ]

    report.to_excel(
        excel_writer=(
            output_dir + '/' + indicator + '.xlsx'
            ),
        sheet_name=indicator,
        float_format="%.1f",
        index=False,
        engine='openpyxl'
        )
    dst_layer.ResetReading()

    data = xr.open_mfdataset(paths=data_files, concat_dim='time')
    drought_intensity = data.Drought_intensity
    drought_magnitude = data.Drought_magnitude

    for theme in [map_files[4]]:  # Using only 'municipios'.
        for data in [drought_intensity, drought_magnitude]:
            vector_ds = ogr.Open(str(theme))
            theme_name = theme.stem
            print("    - The theme '{}' was loaded.".format(theme_name))
            lyr = vector_ds.GetLayer(0)

            # Get a field definition from the original vector file
            feature = lyr.GetFeature(0)
            field = feature.GetFieldDefnRef(0)

            # Reset the original layer so we can read all features
            lyr.ResetReading()
            feature_number = lyr.GetFeatureCount()
            print("    - Rasterizing features and exporting reports.")

            rep.make_report(
                data=data,
                indicator=data.attrs['DroughtFeature'].split('_')[-1],
                regions=lyr,
                fielddef=field,
                output_dir=settings.reports['output_dir'],
                nodata=settings.general['NODATA']
                )
