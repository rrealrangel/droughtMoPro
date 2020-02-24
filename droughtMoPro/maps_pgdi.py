# -*- coding: utf-8 -*-
"""
Created on Sat Nov  3 21:33:05 2018

@author: r.realrangel
"""
import os

from grass.script import array as garray
import grass.script as grass
from scipy.spatial import distance_matrix
import numpy as np
import pandas as pd
import xarray as xr

from droughtMoPro import data_manager as dmgr
from droughtMoPro import spatial_analysis as span

index_cats = {
    # credit: https://www.researchgate.net/profile/Omar_Cenobio-Cruz
    'MSDI_PRESMO01': 1,
    'MSDI_PRESMO03': 2,
    'MSDI_PRESMO06': 3,
    'MSDI_PRESMO09': 4,
    'MSDI_PRESMO12': 5,
    'MSDI_PRERUN01': 6,
    'MSDI_PRERUN03': 7,
    'MSDI_PRERUN06': 8,
    'MSDI_PRERUN09': 9,
    'MSDI_PRERUN12': 10,
    'MSDI_PRESMORUN01': 11,
    'MSDI_PRESMORUN03': 12,
    'MSDI_PRESMORUN06': 13,
    'MSDI_PRESMORUN09': 14,
    'MSDI_PRESMORUN12': 15,
    'SPI01': 16,
    'SPI03': 17,
    'SPI06': 18,
    'SPI09': 19,
    'SPI12': 20,
    'SRI01': 21,
    'SRI03': 22,
    'SRI06': 23,
    'SRI09': 24,
    'SRI12': 25,
    'SSI01': 26,
    'SSI03': 27,
    'SSI06': 28,
    'SSI09': 29,
    'SSI12': 30
    }


def apply_pgdi_filter(data_files, input_filter_fpath):
    dataarrays = {}

    pgdi_filter = xr.open_dataarray(
        filename_or_obj=input_filter_fpath
        )

    for data_file in data_files:
        dataarray = xr.open_dataset(
            filename_or_obj=data_file
            )['Drought_intensity']

        index_name = (
            dataarray.DroughtIndex +
            (dataarray.TemporalScale.split(' ')[0]).zfill(2)
            )

        dataarrays[index_name] = dataarray

    sdi = xr.Dataset(data_vars=dataarrays)
    pgdi = sdi[sdi.data_vars.keys()[0]].copy().rename('PGDI')
    pgdi.values = pgdi.values * np.nan

    # Interpolate sdi filter.
    pgdi_filter = pgdi_filter.interp(
        lat=pgdi.lat,
        lon=pgdi.lon,
        method='nearest'
        )

    for sdi_name, sdi_data in sdi.data_vars.items():
        pgdi = xr.where(
            pgdi_filter == index_cats[sdi_name],
            sdi[sdi_name],
            pgdi
            )

    return(pgdi)


def interp_n_trim_dataset(data, trimmer, output_res, nodata=-32768):
    input_res = data.lat.values[1] - data.lat.values[0]
    ymin = min(data.lat.values) - (input_res / 2)
    ymax = max(data.lat.values) + (input_res / 2)
    xmin = min(data.lon.values) - (input_res / 2)
    xmax = max(data.lon.values) + (input_res / 2)
    output_lat = np.arange(ymin + (output_res / 2), ymax, output_res)
    output_lon = np.arange(xmin + (output_res / 2), xmax, output_res)
    data_stacked = data.stack(x=['lat', 'lon'])
    cell_empty = data_stacked[data_stacked.isnull()].x.values
    cell_data = data_stacked[data_stacked.notnull()].x.values

    distances = pd.DataFrame(
        data=distance_matrix(
            x=list(cell_empty),
            y=list(cell_data)
            ),
        columns=cell_data,
        index=cell_empty
        )

    distances[distances > 1.25] = np.nan
    distances = distances.loc[(distances < 1.25).sum(axis=1) > 0]
    data_filled_edges = data.copy()

    for c, (lat, lon) in enumerate(distances.index.values):
        distances_local = distances.iloc[c][(distances.iloc[c]).notnull()]
        data_local = distances_local.copy()

        for i in data_local.index:
            data_local[i] = data.loc[{'lat': i[0], 'lon': i[1]}].values

        data_filled_edges.loc[{'lat': lat, 'lon': lon}] = span.interpolate_idw(
            distances=distances_local,
            values=data_local,
            power=2
            )

    data_interp = data_filled_edges.interp(lat=output_lat, lon=output_lon)

    data_trimmed = span.trim_data(
        data=data_interp,
        vmap=trimmer,
        res=output_res,
        nodata=nodata
        )

    return(data_trimmed)


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


def export_dint(data, output_fpath, nodata=-32768):
    """
    """
    # Vectorize the input data in GRASS GIS.
    set_grass_env(data)
    export_to_grass(data, map_name='map_to_export', nodata=nodata)
    raster2vector_dint(input_map='map_to_export')

    grass.run_command(
        'v.out.ogr',
        input='map_to_export',
        output=str(output_fpath),
        format='KML',
        overwrite=True,
        quiet=True
        )
    kml_style_dint(
        input_file=str(output_fpath)
        )


def export_pgdi_maps(
        data_files, filter_file, trim_vmap, output_res, output_dir,
        nodata=-32768
        ):
    pgdi_data = apply_pgdi_filter(
        data_files=data_files,
        input_filter_fpath=filter_file
        )

    print("    - Interpolating and trimming PGDI dataset.")
    pgdi_out = interp_n_trim_dataset(
        data=pgdi_data,
        trimmer=trim_vmap,
        output_res=output_res,
        nodata=nodata
        )

    print("    - Exporting the NetCDF file.")
    nc4_file = dmgr.name_pgdi_outfile(
        data=pgdi_out,
        output_dir=output_dir
        )

    nc4_file.parent.mkdir(
        parents=True,
        exist_ok=True
        )

    pgdi_out.to_netcdf(str(nc4_file))

    pgdi_out.attrs['DroughtIndex'] = 'PGDI'
    pgdi_out.attrs['TemporalScale'] = '00'

    print("    - Exporting the KML file.")
    kml_fname = str(nc4_file).replace('.nc4', '.kml')
#
#    export_dint(
#        input_file=str(nc4_file),
#        output_dir=output_dir
#        )
    export_dint(data=pgdi_out, output_fpath=kml_fname, nodata=-32768)
