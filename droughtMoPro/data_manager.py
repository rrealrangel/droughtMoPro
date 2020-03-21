# -*- coding: utf-8 -*-
"""Utilities for manage datasets
Author
------
    Roberto A. Real-Rangel (Institute of Engineering UNAM; Mexico)

License
-------
    GNU General Public License
"""
import numpy as np
import sys

from collections import OrderedDict
from pandas import DateOffset
from pandas import date_range
from pathlib2 import Path
from six import string_types
import datetime as dt
import gdal
import toml
import xarray as xr
import gc


class Configurations():
    """
    """
    def __init__(self, config_file):
        self.config_file = config_file

        config = toml.load(config_file)

        for key, value in config.items():
            setattr(self, key, value)


def vector2array(xmin, ymax, cols, rows, res_x, res_y, layer, nodata):
    # Create the destination data source
    targetDS = (
        gdal.GetDriverByName('MEM').Create(
            '', cols, rows, gdal.GDT_Byte
            )
        )

    targetDS.SetGeoTransform((
        xmin,  # Leftmost pixel position.
        res_x,  # The pixel width.
        0,
        ymax,  # Upper pixel position.
        0,
        -res_y  # The pixel height.
        ))
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


def load_dir(directory):
    """
    Parameters
    ----------
        directory: string
            Full path of the directory to be loaded.
    """
    if not Path(directory).exists():
        create_dir = raw_input(
                "The directory '{}' does not exist.\n"
                "Do you want to create it? [Y] Yes, [N] No. ".
                format(directory))

        if create_dir.lower() == 'y':
            Path(directory).mkdir(parents=True, exist_ok=True)

        else:
            sys.exit("Cannot continue without this directory. Aborting.")

    return(Path(directory))


def list_files(parent_dir, pattern, what='all'):
    """List all files in a directory with a specified pattern.

    Parameters
        parent_dir: string
            Full path of the directory of which the files are to be listed.
        pattern: string or list of strings
            Pattern(s) of the files to be listed.
    """
    parent_dir = Path(parent_dir)
    files_list = []

    if isinstance(pattern, string_types):
        pattern = [pattern]

    for patt in pattern:
        files_list.extend(parent_dir.glob(pattern='**/' + patt))

    if what == 'last':
        available_periods = sorted(list(set([
            i.stem.split('_')[-1] for i in files_list
            ])))
        return(sorted([
            i
            for i in files_list
            if available_periods[-1] in str(i)
            ]))

    elif what == 'all':
        return(sorted(files_list[:]))

    elif isinstance(what, list):
        first = what[0]
        first = dt.date(
            year=int(first[:4]),
            month=int(first[-2:]),
            day=1
            )
        last = what[-1]
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
        return(sorted([
            i
            for i in files_list
            for j in datestamps_to_export
            if j in str(i)
            ]))

    elif isinstance(what, string_types):
        return(sorted([
            i
            for i in files_list
            if what in str(i)
            ]))

    else:
        print("- Nothing to export.")


def name_pgdi_outfile(data, output_dir):
    tlabel_year = str(data.time.dt.year.values.item()).zfill(4)
    tlabel_month = str(data.time.dt.month.values.item()).zfill(2)
    date_label = (tlabel_year + tlabel_month)

    return(Path(
        output_dir +
        '/' +
        tlabel_year +
        '/PGDI_' +
        date_label +
        '.nc4'
        ))


def load_source(path, xyt_vars, chunks_size=-1):
    """Load data.

    Parameters
    ----------
        path:
        xyt_vars: dictionary
        chunks_size: integer, optional (default -1; i. e., no chunks)
    """
    # TODO: Extend the temporal compatility to daily data. Issue #2.
    files_list = []

    for ftype in ['**/*.nc', '**/*.nc4']:
        files_list.extend(load_dir(path).glob(ftype))

    gc.collect()
    data = xr.open_mfdataset(files_list)
    data.attrs.clear()
    data.time.values = (   # This ensures temporal compatibility.
            data.time.values.astype('datetime64[M]'))
    return(data.chunk({
            xyt_vars['t']: -1,
            xyt_vars['x']: chunks_size,
            xyt_vars['y']: chunks_size}))


def drop_array(data, xyt_vars, keeplst=False, drop_patt=False):
    """
    Parameters
    ----------
        data
        xyt_vars: dictionary
        keeplst: bool, optional (default False)
        drop_patt: bool, optional (default False)
    """
    if keeplst is not False:
        joined_list = ['_'.join(i) for i in keeplst]
        vars_to_drop = [i for i in data.var() if i not in joined_list]

    elif drop_patt is not False:
        vars_to_drop = [var
                        for var in data.var()
                        for pattern in drop_patt
                        if pattern in var]

    return(data.drop(vars_to_drop))


def merge_arrays(data, vars_to_merge):
    for merged_variables in vars_to_merge:
        units = list(set([data[i].units for i in merged_variables]))

        if len(units) == 1:
            new_var = '_'.join(merged_variables)
            data[new_var] = sum(
                        [data[i] for i in merged_variables]).assign_attrs(
                                {'units': units[0]})

        data = data.drop(merged_variables)

    return(data)


def convert_units(data):
    """
    Parameters
    ----------
        data: xarray.Dataset
            Dataset of the values which units are to be
            transformed.

    Returns
    -------
        xarray.Dataset
            Dataset of the transformed values.
    """
    time_m = data.time.values.astype('datetime64[M]')
    time_s = data.time.values.astype('datetime64[s]')
    seconds = ((time_m + np.timedelta64(1, 'M') - time_s).
               astype('datetime64[s]')).astype(int)

    for var in data.var():
        try:
            if data[var].units == 'kg m-2 s-1':
                data_aux = data[var].values

                for i, val in enumerate(data_aux):
                    data_aux[i] = val * seconds[i]

                data[var].values = data_aux
                data[var].attrs['units'] = 'mm'

            elif data[var].units == 'K':
                data[var].values = data[var].values - 273.15
                data[var].attrs['units'] = 'C'

            elif data[var].units == 'kg m-2':
                data[var].attrs['units'] = 'mm'

            else:
                pass

        except AttributeError:
            pass

    gc.collect()
    return(data)


def load_all(imp_dataset, xyt_vars, variable=False, drop_patt=False,
             vars_to_merge=False, chunks_size=-1):
    """
    Parameters
    ----------
        imp_dataset: xarray.Dataset
        xyt_vars: dictionary
        variable: list, optional (default False)
        drop_patt: list, optional (default False)
        vars_to_merge: list, optional (default False)
        chunks_size: integer, optional (default -1)
    """
    data = xr.Dataset()

    for path in imp_dataset.itervalues():
        data_aux = load_source(
                path=path,
                xyt_vars=xyt_vars,
                chunks_size=chunks_size)
        data = data.merge(data_aux)

    if vars_to_merge is not False:
        data = merge_arrays(data, vars_to_merge)

    data = drop_array(
            data=data,
            xyt_vars=xyt_vars,
            keeplst=variable,
            drop_patt=drop_patt)

    return(convert_units(data).chunk({
            xyt_vars['x']: chunks_size,
            xyt_vars['y']: chunks_size}))


def delay_time(data, t_lag):
    """Generates a dataset of time-lagged values.

    Parameters
    ----------
        data: xarray.Dataset
        t_lag: list
    """
    # TODO: This function needs to be modified to consider any temporal
    # resolution of the input data. Currently, it only allows monthly data.
    # Issue #2.
    data_lagged = data.copy()
    data_lagged.time.values = (
            data.time.values.astype(
                    'datetime64[M]') + np.timedelta64(t_lag, 'M'))

    for old_name in data_lagged.var():
        new_name = old_name + '_tlag' + str(t_lag)
        data_lagged[new_name] = data_lagged[old_name].rename(new_name)
        data_lagged = data_lagged.drop(old_name)

    gc.collect()
    return(data_lagged)


def accumulate_time(data, t_acc):
    """Generates a dataset of time-accumulated values.

    Parameters
    ----------
        data: xarray.Dataset
        t_acc: list
    """
    # TODO: This function needs to be modified to consider any temporal
    # resolution of the input data. Currently, it only allows monthly data.
    # Issue #2.
    data_accum = data.copy()

    for var in data_accum.var():
        try:
            if data_accum[var].units == 'mm':
                data_accum[var] = data_accum[var].rolling(time=t_acc).sum()

            else:
                data_accum[var] = data_accum[var].rolling(time=t_acc).mean()

        except (AttributeError):
            data_accum[var] = data_accum[var].rolling(time=t_acc).sum()

        except (KeyError):
            data_accum = data_accum.drop(var)

    gc.collect()

    for old_name in data_accum.var():
        new_name = old_name + '_tacc' + str(t_acc)
        data_accum[new_name] = data_accum[old_name].rename(new_name)
        data_accum = data_accum.drop(old_name)

    gc.collect()
    return(data_accum)


def aggregate(data, xyt_vars, variable, t_acc=False, t_lag=False):
    """Delay and aggregate data in time. Future versions will 'delay' and
        aggregate data in space too.

    Parameters:
        data : xarray.Dataset
        xyt_vars : dictionary
        variable : list or False
        t_acc : list (optional; default: False)
        t_lag : list (optional; default: False)
    """

    # Remove (drop) unwanted variables.
    if variable is not False:
        data = drop_array(data, xyt_vars, keeplst=variable, drop_patt=False)

    # Accumulate values in time.
    if t_acc is not False:
        for i, window in enumerate(t_acc):
            if i == 0:
                data_accum = accumulate_time(data, window)

            else:
                data_accum = data_accum.merge(
                        accumulate_time(data, window))

    # Add a time lag in values.
    if t_lag is not False:
        for i, lag in enumerate(t_lag):
            if i == 0:
                data_lagged = delay_time(data_accum, lag)

            else:
                data_lagged = data_lagged.merge(
                        delay_time(data_accum, lag))

    # Remove (drop) nan values.
    start = 0
    stop = len(data_lagged.time)

    if t_acc is not False:
        start += (max(t_acc) - 1)

    if t_lag is not False:
        start += max(t_lag)
        stop -= max(t_lag)

    return(data_lagged[dict(time=slice(start, stop))])


def check_source(raw_dataset, imp_dataset):
    """Import MERRA-2 dataset.

    Parameters
    ----------
    raw_dataset: dictionary
    imp_dataset: dictionary
    """
    # TODO: Check and update the last year dataset if it is out of date.
    # OR always import any missing AND the last/current year.
    for name, path in raw_dataset.iteritems():
        print("- Checking imported {} dataset in '{}'".format(
                name, imp_dataset[name]))
        files_list = []

        for ftype in ['**/*.nc', '**/*.nc4']:
            files_list.extend(load_dir(path).glob(ftype))

        gc.collect()
        years = sorted(set([str(i).split('.')[2][:4] for i in files_list]))

        for year in years:
            imp_file = load_dir(
                    imp_dataset[name])/(
                            name + '_' + year + '.nc4')

            if not imp_file.is_file():
                print("  - Importing {} dataset of {}".format(name, year))
                sources = [
                        i for i in files_list
                        if str(i).split('.')[2][:4] == year]
                xr.open_mfdataset(sources).to_netcdf(imp_file)
                # TODO: Update attributes and include the temporal resolution.

    gc.collect()


def slice_month(dataset, month):
    """
    Parameters
    ----------
        dataset: xarray.Dataset
        month: integer
    """
    return(dataset.sel(time=dataset['time.month'] == month))


def progress_message(current, total, message="- Processing", units=None):
    """Issue a messages of the progress of the process.

    Generates a progress bar in terminal. It works within a for loop,
    computing the progress percentage based on the current item
    number and the total length of the sequence of item to iterate.

    Parameters:
        current : integer
            The last item number computed within the for loop. This
            could be obtained using enumerate() in when calling the for
            loop.
        total : integer
            The total length of the sequence for which the for loop is
            performing the iterations.
        message : string (optional; default = "- Processing")
    """
    if units is not None:
        progress = float(current)/total
        sys.stdout.write(
            "\r    {} ({:.1f} % of {} processed)".format(
                message, progress * 100, units
                )
            )

    else:
        progress = float(current)/total
        sys.stdout.write(
            "\r    {} ({:.1f} % processed)".format(
                message, progress * 100
                )
            )

    if progress < 1:
        sys.stdout.flush()

    else:
        sys.stdout.write('\n')


def monthly_dataset(date, arrays, title):
    """
    """
    data_vars = {
            i.attrs['DroughtFeature']: i.sel({'time': date}) for i in arrays}
    year = str(date.dt.year.values)
    month = str(date.dt.month.values).zfill(2)
    output_dataset = xr.Dataset(data_vars=data_vars)
    attrs = OrderedDict()
    attrs['Title'] = title
    attrs['TemporalRange'] = year + month
    attrs['SouthernmostLatitude'] = min(output_dataset.lat.values)
    attrs['NorthernmostLatitude'] = max(output_dataset.lat.values)
    attrs['WesternmostLongitude'] = min(output_dataset.lon.values)
    attrs['EasternmostLongitude'] = max(output_dataset.lon.values)
    attrs['LatitudeResolution'] = output_dataset.lat.diff(dim='lat').values[0]
    attrs['LongitudeResolution'] = output_dataset.lon.diff(dim='lon').values[0]
    attrs['SpatialCoverage'] = 'Mexico'
    attrs['History'] = (
            'Original file generated:' + dt.datetime.now().isoformat())
    attrs['Contact'] = (
            'Roberto A. Real-Rangel (rrealr@iingen.unam.mx)')
    attrs['Institution'] = (
            'Institute of Engineering of the '
            'National Autonomous University of Mexico (II-UNAM)')
    attrs['Format'] = 'NetCDF-4/HDF-5'
    attrs['VersionID'] = '1.0.0'
    output_dataset.attrs = attrs
    return(output_dataset)


def export_nc4(dataset, output_dir, prefix=None):
    """ Export a given dataset to a NetCDF-4 file (.nc4).

    Parameters:
        dataset : xarray.Dataset
        output_dir : string
        index : string
        temp_scale : integer
        prefix : string (optional; default is None)
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    output_file = (
            prefix + dataset.attrs['Title'].upper() + '_'
            + dataset.attrs['TemporalRange'] + '.nc4')
    dataset.to_netcdf(Path(output_dir) / output_file)
