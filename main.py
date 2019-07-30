# -*- coding: utf-8 -*-
"""Nonparametric Standardized Drought Indices - pySDI main script
Author
------
    Roberto A. Real-Rangel (Institute of Engineering UNAM; Mexico)

License
-------
    GNU General Public License
"""
import datetime as dt

import mosemm_products.data_manager as dmgr
import mosemm_products.maps_pgdi as pgdi
import mosemm_products.maps_sdi as maps
import mosemm_products.reports as rep
import mosemm_products.time_series as ts

run_start = dt.datetime.now()
print("---- Process started at {start} ----".format(start=run_start))
settings = dmgr.Configurations(config_file='config.toml')

# Export drought intensity maps to KML format files.
if settings.intensity_maps['export']:
    print("- Working on drought intensity maps.")

    # Available source data to generate drought monitoring products.
    data_files = dmgr.list_files(
        parent_dir=settings.general['input_dir'],
        pattern=settings.intensity_maps['datafile_patt'],
        what=settings.intensity_maps['period_to_export']
        )

    for df, data_file in enumerate(data_files):
        maps.export_dint(
            input_file=data_file,
            output_dir=settings.intensity_maps['output_dir'],
            nodata=settings.general['NODATA'],
            output_fname_prefix=settings.general['output_prefix'],
            overwrite=settings.general['overwrite']
            )
        dmgr.progress_message(
            current=(df + 1),
            total=len(data_files),
            message="- Exporting the drought intensity maps",
            units='maps'
            )

# Export drought magnitude maps to KML format files.
if settings.magnitude_maps['export']:
    print("- Working on drought magnitude maps.")

    # Available source data to generate drought monitoring products.
    data_files = dmgr.list_files(
        parent_dir=settings.general['input_dir'],
        pattern=settings.magnitude_maps['datafile_patt'],
        what=settings.magnitude_maps['period_to_export']
        )

    for df, data_file in enumerate(data_files):
        maps.export_dmag(
            input_file=data_file,
            output_dir=settings.magnitude_maps['output_dir'],
            nodata=settings.general['NODATA'],
            output_fname_prefix=settings.general['output_prefix'],
            overwrite=settings.general['overwrite']
            )
        dmgr.progress_message(
            current=(df + 1),
            total=len(data_files),
            message="- Exporting the drought magnitude maps",
            units='maps'
            )

# Export time series.
if settings.time_series['export']:
    print("- Working on time series products.")
    data_files = dmgr.list_files(
        parent_dir=settings.general['input_dir'],
        pattern=settings.time_series['datafile_patt'],
        what=settings.time_series['period_to_export']
        )
    map_files = dmgr.list_files(
        parent_dir=settings.time_series['vmaps_dir'],
        pattern=settings.time_series['mapfile_patt']
        )
    ts.export_ts(
        data_files=data_files,
        map_files=map_files,
        nodata=settings.general['NODATA'],
        output_dir=settings.time_series['output_dir']
        )

# Export reports.
if settings.reports['export']:
    print("- Working on drought reports.")
    print("    - Reading regions maps.")
    data_files = dmgr.list_files(
        parent_dir=settings.general['input_dir'],
        pattern=settings.reports['datafile_patt'],
        what=settings.reports['period_to_export']
        )
    map_files = dmgr.list_files(
        parent_dir=settings.time_series['vmaps_dir'],
        pattern=settings.reports['mapfile_patt']
        )
    rep.make_report(
        data_files=data_files,
        map_files=map_files,
        output_dir=settings.reports['output_dir'],
        nodata=settings.general['NODATA']
        )

# Compute and export Proxy groundwater drought index (PGDI)
# TODO: Apply the PGDI filter to the SDI files with the original resolution.
if settings.pgdi['compute']:
    print("- Computing the PGDI.")
    print("    - Reading SDI files.")
    data_files = dmgr.list_files(
        parent_dir=settings.pgdi['input_data_dir'],
        pattern=settings.pgdi['input_pattern']
        )

    print("    - Applying the SDI filter.")
    pgdi_data = pgdi.apply_pgdi_filter(
        data_files=data_files,
        sdi_filter=settings.pgdi['input_filter_dir'])

    print("    - Interpolating and trimming PGDI dataset.")
    pgdi_out = pgdi.interp_n_trim_dataset(
        data=pgdi_data.copy(),
        trimmer=settings.pgdi['trim_vmap'],
        output_res=settings.pgdi['output_res'],
        nodata=settings.general['NODATA']
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

run_end = dt.datetime.now()
print("---- Process ended at {end}. Time elapsed: {eltime} ----".format(
    end=run_end,
    eltime=str(run_end - run_start)
    ))
