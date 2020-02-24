#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 22:22:42 2019
@author: realrangel
"""
import os

from grass.script import array as garray
from pathlib2 import Path
import grass.script as grass
import numpy as np
import xarray as xr


def rename_indices(old_name):
    names = {
        'SPI': 'SPI',
        'SRI': 'SRI',
        'SSI': 'SSI',
        'MSDI_PRERUN': 'PreRun',
        'MSDI_PRESMO': 'PreSMo',
        'MSDI_PRESMORUN': 'PreSMoRun'
        }
    return(names[old_name])


def set_grass_env(reference):
    grass.run_command('g.gisenv', set="DEBUG=0")
    grass.run_command('g.gisenv', set="GRASS_VERBOSE=-1")
    res = reference.lat.diff(dim='lat').values[0]
    n = max(reference.lat).values + res
    s = min(reference.lat).values - res
    e = max(reference.lon).values + res
    w = min(reference.lon).values - res
    grass.run_command(
        'g.region',
        n=n,
        s=s,
        e=e,
        w=w,
        rows=len(reference.lat),
        cols=len(reference.lon),
        res=res,
        overwrite=True,
        quiet=True
        )


def export_to_grass(data, map_name, nodata=-32768):
    # Export the array as a raster map to GRASS-GIS
    output_map = garray.array()
    output_map[:, :] = np.flipud(data)
    output_map.write(
        mapname=map_name,
        overwrite=True
        )


def raster2vector_dint(input_map):
    # Reclass rasters following the drought index categories
    grass.mapcalc(
            '$output='
            'if($input <= -2.00, 1,'
            'if($input <= -1.60, 2,'
            'if($input <= -1.30, 3,'
            'if($input <= -0.80, 4,'
            'if($input <= -0.50, 5,'
            'if($input < 0.50, 6,'
            'if($input < 0.80, 7,'
            'if($input < 1.30, 8,'
            'if($input < 1.60, 9,'
            'if($input < 2.00, 10, 11))))))))))',
            input=input_map,
            output=input_map + '_reclass',
            overwrite=True,
            quiet=True
            )
    grass.run_command(
            'r.to.vect',
            flags='sv',
            input=input_map + '_reclass',
            output=input_map,
            type='area',
            overwrite=True,
            quiet=True
            )

    # Update the label (W4, W3, W2, etc.) of the categories used
    labels_cond = {
        'cat=1': 'D4',
        'cat=2': 'D3',
        'cat=3': 'D2',
        'cat=4': 'D1',
        'cat=5': 'D0',
        'cat=6': 'Normal',
        'cat=7': 'W0',
        'cat=8': 'W1',
        'cat=9': 'W2',
        'cat=10': 'W3',
        'cat=11': 'W4'
        }
    main_dir = os.path.dirname(os.path.abspath('__file__'))

    with open(main_dir + "/update_db.sql", "w") as text_file:
        text_file.write(
            'ALTER TABLE {} ADD COLUMN rgb varchar(11);\n'.format(input_map)
            )

        for key, val in labels_cond.iteritems():
            text_file.write(
                "UPDATE {} SET label='{}' WHERE {};\n".format(
                        input_map, val, key
                        )
                )

    # Add the color rules for the drought vector map
    colors_cond = {
        'cat=1': '189:0:38',
        'cat=2': '240:59:32',
        'cat=3': '253:141:60',
        'cat=4': '254:204:92',
        'cat=5': '255:255:178',
        'cat=6': '255:255:255',
        'cat=7': '240:249:232',
        'cat=8': '186:228:188',
        'cat=9': '123:204:196',
        'cat=10': '67:162:202',
        'cat=11': '8:104:172'
        }

    with open(main_dir + "/update_db.sql", "a") as text_file:
        for key, val in colors_cond.iteritems():
            text_file.write("UPDATE {} SET rgb='{}' WHERE {};\n".format(
                input_map, val, key
                ))

    grass.run_command('db.execute', input=main_dir+"/update_db.sql")


def raster2vector_dmag(input_map):
    # Reclass rasters following the drought index categories
    grass.mapcalc(
            '$output='
            'if($input > -1, 0,'
            'if($input > -3, 1,'
            'if($input > -6, 2,'
            'if($input > -9, 3,'
            'if($input > -12, 4, 5)))))',
            input=input_map,
            output=input_map + '_reclass',
            overwrite=True,
            quiet=True
            )
    grass.run_command(
            'r.to.vect',
            flags='sv',
            input=input_map + '_reclass',
            output=input_map,
            type='area',
            overwrite=True,
            quiet=True
            )

    # Update the label (W4, W3, W2, etc.) of the categories used
    labels_cond = {
        'cat=0': 'No drought',
        'cat=1': 'M1',
        'cat=2': 'M2',
        'cat=3': 'M3',
        'cat=4': 'M4',
        'cat=5': 'M5'
        }
    main_dir = os.path.dirname(os.path.abspath('__file__'))

    with open(main_dir + "/update_db.sql", "w") as text_file:
        text_file.write(
            'ALTER TABLE {} ADD COLUMN rgb varchar(11);\n'.format(
                input_map
                ))

        for key, val in labels_cond.iteritems():
            text_file.write(
                "UPDATE {} SET label='{}' WHERE {};\n".format(
                    input_map, val, key
                    ))

    # Add the color rules for the drought vector map
    colors_cond = {
        'cat=0': '255:255:255',
        'cat=1': '255:215:0',
        'cat=2': '253:140:0',
        'cat=3': '254:0:0',
        'cat=4': '128:0:128',
        'cat=5': '64:0:64'
        }

    with open(main_dir + "/update_db.sql", "a") as text_file:
        for key, val in colors_cond.iteritems():
            text_file.write(
                "UPDATE {} SET rgb='{}' WHERE {};\n".format(
                    input_map, val, key
                    )
                )

    grass.run_command('db.execute', input=main_dir+"/update_db.sql")


def kml_style_dint(input_file):
    with open(input_file, 'r') as f:
        kml = f.read()

    kml = kml.splitlines()
    pmark_begin = []

    for ln, line in enumerate(kml):
        if "<Placemark>" in line:
            pmark_begin.append(ln)

    pmark_end = []

    for ln, line in enumerate(kml):
        if "</Placemark>" in line:
            pmark_end.append(ln)

    def style(hexcolor):
        style = (
            "\t<Style><PolyStyle><color>" + hexcolor +
            "</color><outline>0</outline></PolyStyle></Style>"
            )

        return(style)

    for b, begin in enumerate(pmark_begin):
        for line in range(pmark_begin[b], pmark_end[b]):
            if """<SimpleData name="label">""" in kml[line]:
                if "D4" in kml[line]:
                    for line2 in range(pmark_begin[b], pmark_end[b]):
                        if "<Style>" in kml[line2]:
                            kml[line2] = style('992600bd')

                    break

                elif "D3" in kml[line]:
                    for line2 in range(pmark_begin[b], pmark_end[b]):
                        if "<Style>" in kml[line2]:
                            kml[line2] = style('99203bf0')

                elif "D2" in kml[line]:
                    for line2 in range(pmark_begin[b], pmark_end[b]):
                        if "<Style>" in kml[line2]:
                            kml[line2] = style('993c8dfd')

                elif "D1" in kml[line]:
                    for line2 in range(pmark_begin[b], pmark_end[b]):
                        if "<Style>" in kml[line2]:
                            kml[line2] = style('995cccfe')

                elif "D0" in kml[line]:
                    for line2 in range(pmark_begin[b], pmark_end[b]):
                        if "<Style>" in kml[line2]:
                            kml[line2] = style('99b2ffff')

                elif "Normal" in kml[line]:
                    for line2 in range(pmark_begin[b], pmark_end[b]):
                        if "<Style>" in kml[line2]:
                            kml[line2] = style('0')

                elif "W0" in kml[line]:
                    for line2 in range(pmark_begin[b], pmark_end[b]):
                        if "<Style>" in kml[line2]:
                            kml[line2] = style('99e8f9f0')

                elif "W1" in kml[line]:
                    for line2 in range(pmark_begin[b], pmark_end[b]):
                        if "<Style>" in kml[line2]:
                            kml[line2] = style('99bce4ba')

                elif "W2" in kml[line]:
                    for line2 in range(pmark_begin[b], pmark_end[b]):
                        if "<Style>" in kml[line2]:
                            kml[line2] = style('99c4cc7b')

                elif "W3" in kml[line]:
                    for line2 in range(pmark_begin[b], pmark_end[b]):
                        if "<Style>" in kml[line2]:
                            kml[line2] = style('99caa243')

                elif "W4" in kml[line]:
                    for line2 in range(pmark_begin[b], pmark_end[b]):
                        if "<Style>" in kml[line2]:
                            kml[line2] = style('99ac6808')

    with open(input_file, 'w') as kml_out:
        for ln, line in enumerate(kml):
            kml_out.write(line)
            kml_out.write("\n")


def kml_style_dmag(input_file):
    with open(input_file, 'r') as f:
        kml = f.read()

    kml = kml.splitlines()
    pmark_begin = []

    for ln, line in enumerate(kml):
        if "<Placemark>" in line:
            pmark_begin.append(ln)

    pmark_end = []

    for ln, line in enumerate(kml):
        if "</Placemark>" in line:
            pmark_end.append(ln)

    def style(hexcolor):
        style = (
            "\t<Style><PolyStyle><color>" + hexcolor +
            "</color><outline>0</outline></PolyStyle></Style>"
            )

        return(style)

    for b, begin in enumerate(pmark_begin):
        for line in range(pmark_begin[b], pmark_end[b]):
            if """<SimpleData name="label">""" in kml[line]:
                if "No drought" in kml[line]:
                    for line2 in range(pmark_begin[b], pmark_end[b]):
                        if "<Style>" in kml[line2]:
                            kml[line2] = style('0')

                elif "M1" in kml[line]:
                    for line2 in range(pmark_begin[b], pmark_end[b]):
                        if "<Style>" in kml[line2]:
                            kml[line2] = style('9900D7FF')

                elif "M2" in kml[line]:
                    for line2 in range(pmark_begin[b], pmark_end[b]):
                        if "<Style>" in kml[line2]:
                            kml[line2] = style('99008CFF')

                elif "M3" in kml[line]:
                    for line2 in range(pmark_begin[b], pmark_end[b]):
                        if "<Style>" in kml[line2]:
                            kml[line2] = style('990000FF')

                elif "M4" in kml[line]:
                    for line2 in range(pmark_begin[b], pmark_end[b]):
                        if "<Style>" in kml[line2]:
                            kml[line2] = style('99800080')

                elif "M5" in kml[line]:
                    for line2 in range(pmark_begin[b], pmark_end[b]):
                        if "<Style>" in kml[line2]:
                            kml[line2] = style('99400040')

    with open(input_file, 'w') as kml_out:
        for ln, line in enumerate(kml):
            kml_out.write(line)
            kml_out.write("\n")


def export_dint(
        input_file, output_dir, nodata=-32768, output_fname_prefix='',
        overwrite=False
        ):
    """
    """
    # Open the data.
    data = xr.open_dataset(input_file).Drought_intensity

    # Vectorize the input data in GRASS GIS.
    set_grass_env(data)
    export_to_grass(data, map_name='map_to_export', nodata=nodata)
    raster2vector_dint(input_map='map_to_export')

    # Export the data as a KML file.
    # - Build the output file name.
    index = data.attrs['DroughtIndex'].lower()
    time_scale = data.attrs['TemporalScale'][:2].strip().zfill(2)
    year = str(data.time.dt.year.values)
    output_subdir = Path(index + time_scale) / year
    output_fname = (
        output_fname_prefix + rename_indices(data.attrs['DroughtIndex'])
        + time_scale + '_' + year + str(data.time.dt.month.values).zfill(2) +
        '_res.kml'
        )

    # - Make the directory, in case it does not exists.
    (Path(output_dir) / output_subdir).mkdir(
        parents=True,
        exist_ok=True
        )

    # Export the KML file and modify its style.
    output_full_path = Path(output_dir) / output_subdir / output_fname
    grass.run_command(
        'v.out.ogr',
        input='map_to_export',
        output=str(output_full_path),
        format='KML',
        overwrite=True,
        quiet=True
        )
    kml_style_dint(
        input_file=str(output_full_path)
        )


def export_dmag(
        input_file, output_dir, nodata=-32768, output_fname_prefix='',
        overwrite=False
        ):
    """
    """
    # Open the data.
    data = xr.open_dataset(input_file).Drought_magnitude

    # Vectorize the input data in GRASS GIS.
    set_grass_env(data)
    export_to_grass(data, map_name='map_to_export', nodata=nodata)
    raster2vector_dmag(input_map='map_to_export')

    # Export the data as a KML file.
    # - Build the output file name.
    index = data.attrs['DroughtIndex'].lower()
    time_scale = data.attrs['TemporalScale'][:2].strip().zfill(2)
    year = str(data.time.dt.year.values)
    output_subdir = Path(index + time_scale) / year
    output_fname = (
        output_fname_prefix + rename_indices(data.attrs['DroughtIndex'])
        + time_scale + '_dm_' + year +
        str(data.time.dt.month.values).zfill(2) + '_res.kml'
        )

    # - Make the directory, in case it does not exists.
    (Path(output_dir) / output_subdir).mkdir(
        parents=True,
        exist_ok=True
        )

    # Export the KML file and modify its style.
    output_full_path = Path(output_dir) / output_subdir / output_fname
    grass.run_command(
        'v.out.ogr',
        input='map_to_export',
        output=str(output_full_path),
        format='KML',
        overwrite=True,
        quiet=True
        )
    kml_style_dmag(
        input_file=str(output_full_path)
        )
