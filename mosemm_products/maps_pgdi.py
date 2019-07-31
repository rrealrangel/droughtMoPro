# -*- coding: utf-8 -*-
"""
Created on Sat Nov  3 21:33:05 2018

@author: r.realrangel
"""
from scipy.spatial import distance_matrix
import numpy as np
import pandas as pd
import xarray as xr

import lib.spatial_analysis as span

index_cats = {
    1: 'MSDI_PRESMO01',
    2: 'MSDI_PRESMO03',
    3: 'MSDI_PRESMO06',
    4: 'MSDI_PRESMO09',
    5: 'MSDI_PRESMO12',
    6: 'MSDI_PRERUN01',
    7: 'MSDI_PRERUN03',
    8: 'MSDI_PRERUN06',
    9: 'MSDI_PRERUN09',
    10: 'MSDI_PRERUN12',
    11: 'MSDI_PRESMORUN01',
    12: 'MSDI_PRESMORUN03',
    13: 'MSDI_PRESMORUN06',
    14: 'MSDI_PRESMORUN09',
    15: 'MSDI_PRESMORUN12',
    16: 'SPI01',
    17: 'SPI03',
    18: 'SPI06',
    19: 'SPI09',
    20: 'SPI12',
    21: 'SRI01',
    22: 'SRI03',
    23: 'SRI06',
    24: 'SRI09',
    25: 'SRI12',
    26: 'SSI01',
    27: 'SSI03',
    28: 'SSI06',
    29: 'SSI09',
    30: 'SSI12'
    }


def apply_pgdi_filter(data_files, sdi_filter):
    dataarrays = {}

    pgdi_cats = xr.open_dataarray(
        filename_or_obj=sdi_filter
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

    for lat in pgdi.lat.values:
        for lon in pgdi.lon.values:
            try:
                cat = pgdi_cats.sel(
                    indexers={'lat': lat, 'lon': lon}
                    ).values.item()

                pgdi.loc[{'lat': lat, 'lon': lon}] = sdi[index_cats[cat]].sel(
                    indexers={'lat': lat, 'lon': lon}
                    ).values

            except KeyError:
                pass

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


def export_pgdi_maps(
        data_files, filter_file, trim_vmap, output_res, nodata=-32768
        ):
    pgdi_data = apply_pgdi_filter(
        data_files=data_files,
        sdi_filter=filter_file
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
        output_dir=settings.pgdi['output_dir']
        )

    nc4_file.parent.mkdir(
        parents=True,
        exist_ok=True
        )

    pgdi_out.to_netcdf(str(nc4_file))

    if settings.pgdi['export_kml']:
        print("    - Exporting the KML file.")
        kml_file = str(nc4_file).replace('.nc4', '.kml')

        maps.export_kml_dint(
            input_file=str(nc4_file),
            output_file=kml_file,
            array_name=False
            )