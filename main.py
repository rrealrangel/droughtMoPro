# -*- coding: utf-8 -*-
"""Nonparametric Standardized Drought Indices - pySDI main script
Author
------
    Roberto A. Real-Rangel (Institute of Engineering UNAM; Mexico)

License
-------
    GNU General Public License
"""
from pandas import DateOffset
from pandas import date_range
# import ogr
import datetime as dt
# import xarray as xr

import mosemm_products.data_manager as dmgr
# import mosemm_products.pgdi as pgdi
import mosemm_products.nc4_to_kml as tokml
# import mosemm_products.reports_dev as rep

settings = dmgr.Configurations(config_file='config.toml')

# Export drought intensity maps to KML format files.
if settings.intensity_maps['export']:
    print("- Exporting drought intensity maps.")

    # Available source data to generate drought monitoring products.
    input_files = dmgr.list_files(
        parent_dir=settings.intensity_maps['input_dir'],
        pattern='**/*.nc4'
        )

    if settings.general['period_to_export'] == 'last':
        available_periods = sorted(list(set([
            i.stem.split('_')[-1] for i in input_files
            ])))
        files_to_export = sorted([
            i
            for i in input_files
            if available_periods[-1] in str(i)
            ])

    elif settings.general['period_to_export'] == 'all':
        files_to_export = sorted(input_files[:])

    elif isinstance(settings.general['period_to_export'], list):
        first = settings.general['period_to_export'][0]
        first = dt.date(
            year=int(first[:4]),
            month=int(first[-2:]),
            day=1
            )
        last = settings.general['period_to_export'][-1]
        last = dt.date(
            year=int(last[:4]),
            month=int(last[-2:]),
            day=1
            ) + DateOffset(months=1)
        dates = date_range(start=first, end=last, freq='1M')
        datestamps_to_export = [
            str(i.year) + str(i.month).zfill(2)
            for i in dates
            ]
        files_to_export = sorted([
            i
            for i in input_files
            for j in datestamps_to_export
            if j in str(i)
            ])

    elif isinstance(settings.general['period_to_export'], unicode):
        files_to_export = sorted([
            i
            for i in input_files
            if settings.general['period_to_export'] in str(i)
            ])

    else:
        print("- Nothing to export.")

    for f, file_to_export in enumerate(files_to_export):
        tokml.export_kml_dint(
            input_file=file_to_export,
            output_dir=settings.intensity_maps['output_dir'],
            output_fname_prefix=settings.general['output_prefix'],
            overwrite=settings.general['overwrite']
            )
        dmgr.progress_message(
            current=(f + 1),
            total=len(files_to_export),
            message="- Exporting the drought intensity maps",
            units='maps'
            )

# Export drought magnitude maps to KML format files.
if settings.magnitude_maps['export']:
    print("- Exporting drought magnitude maps.")

    # Available source data to generate drought monitoring products.
    input_files = dmgr.list_files(
        parent_dir=settings.magnitude_maps['input_dir'],
        pattern='**/*.nc4'
        )

    if settings.general['period_to_export'] == 'last':
        available_periods = sorted(list(set([
            i.stem.split('_')[-1] for i in input_files
            ])))
        files_to_export = sorted([
            i
            for i in input_files
            if ((available_periods[-1] in str(i)) and '-01_' in str(i))
            ])

    elif settings.general['period_to_export'] == 'all':
        files_to_export = sorted([i for i in input_files if '-01_' in str(i)])

    elif isinstance(settings.general['period_to_export'], list):
        first = settings.general['period_to_export'][0]
        first = dt.date(
            year=int(first[:4]),
            month=int(first[-2:]),
            day=1
            )
        last = settings.general['period_to_export'][-1]
        last = dt.date(
            year=int(last[:4]),
            month=int(last[-2:]),
            day=1
            ) + DateOffset(months=1)
        dates = date_range(start=first, end=last, freq='1M')
        datestamps_to_export = [
            str(i.year) + str(i.month).zfill(2)
            for i in dates
            ]
        files_to_export = sorted([
            i
            for i in input_files
            for j in datestamps_to_export
            if ((j in str(i)) and ('-01_' in str(i)))
            ])

    elif isinstance(settings.general['period_to_export'], unicode):
        files_to_export = sorted([
            i
            for i in input_files
            if ((settings.general['period_to_export'] in str(i)) and
                ('-01_' in str(i)))
            ])

    else:
        print("- Nothing to export.")

    for f, file_to_export in enumerate(files_to_export):
        tokml.export_kml_dmag(
            input_file=file_to_export,
            output_dir=settings.magnitude_maps['output_dir'],
            output_fname_prefix=settings.general['output_prefix'],
            overwrite=settings.general['overwrite']
            )
        dmgr.progress_message(
            current=(f + 1),
            total=len(files_to_export),
            message="- Exporting the drought magnitude maps",
            units='maps'
            )

# Export time series.



#
## Export reports (csv; or xlsx?).
#if settings.reports['compute']:
#    print("- Computing the drought reports.")
#    print("    - Reading regions maps.")
#
#    regions_files = dmgr.list_files(
#        parent_dir=settings.reports['input_regions_dir'],
#        pattern='*.shp'
#        )
#
#    print("    - Reading SDI files.")
#
#    sdi_files = dmgr.list_files(
#        parent_dir=settings.reports['input_data_dir'],
#        pattern=settings.reports['input_pattern']
#        )
#
#    data = xr.open_mfdataset(paths=sdi_files, concat_dim='time')
#    drought_intensity = data.Drought_intensity
#    drought_magnitude = data.Drought_magnitude
#
#    for theme in [regions_files[4]]:  # Using only 'municipios'.
#        for data in [drought_intensity, drought_magnitude]:
#            vector_ds = ogr.Open(str(theme))
#            theme_name = theme.stem
#            print("    - The theme '{}' was loaded.".format(theme_name))
#            lyr = vector_ds.GetLayer(0)
#
#            # Get a field definition from the original vector file
#            feature = lyr.GetFeature(0)
#            field = feature.GetFieldDefnRef(0)
#
#            # Reset the original layer so we can read all features
#            lyr.ResetReading()
#            feature_number = lyr.GetFeatureCount()
#            print("    - Rasterizing features and exporting reports.")
#
#            rep.make_report(
#                data=data,
#                indicator=data.attrs['DroughtFeature'].split('_')[-1],
#                regions=lyr,
#                fielddef=field,
#                output_dir=settings.reports['output_dir'],
#                nodata=settings.reports['nodata']
#                )
#
## Compute and export Proxy groundwater drought index (PGDI)
## TODO: Apply the PGDI filter to the SDI files with the original resolution.
#if settings.pgdi['compute']:
#    print("- Computing the PGDI.")
#    print("    - Reading SDI files.")
#    sdi_files = dmgr.list_files(
#        parent_dir=settings.pgdi['input_data_dir'],
#        pattern=settings.pgdi['input_pattern']
#        )
#
#    print("    - Applying the SDI filter.")
#    pgdi_data = pgdi.apply_pgdi_filter(
#        sdi_files=sdi_files,
#        sdi_filter=settings.pgdi['input_filter_dir'])
#
#    print("    - Interpolating and trimming PGDI dataset.")
#    pgdi_out = pgdi.interp_n_trim_dataset(
#        data=pgdi_data.copy(),
#        trimmer=settings.pgdi['trim_vmap'],
#        output_res=settings.pgdi['output_res'],
#        nodata=settings.pgdi['nodata']
#        )
#
#    print("    - Exporting the NetCDF file.")
#    nc4_file = dmgr.name_pgdi_outfile(
#        data=pgdi_out,
#        output_dir=settings.pgdi['output_dir']
#        )
#
#    nc4_file.parent.mkdir(
#        parents=True,
#        exist_ok=True
#        )
#
#    pgdi_out.to_netcdf(str(nc4_file))
#
#    if settings.pgdi['export_kml']:
#        print("    - Exporting the KML file.")
#        kml_file = str(nc4_file).replace('.nc4', '.kml')
#
#        tokml.export_kml_dint(
#            input_file=str(nc4_file),
#            output_file=kml_file,
#            array_name=False
#            )
