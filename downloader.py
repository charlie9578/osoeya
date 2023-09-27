"""
This module provides functions for downloading files, including reanalysis data

This module provides functions for downloading data, including long-term historical atmospheric
data from the MERRA2 and ERA5 reanalysis products and returning as pandas DataFrames and saving
data in csv files. Currently by default the module downloads monthly reanalysis data for a time
period of interest using NASA Goddard Earth Sciences Data and Information Services Center
(GES DISC) for MERRA2 and the Copernicus Climate Data Store (CDS) API for ERA5, but this can be
modified to get hourly data, and indeed other data sources available on GES DISC and CDS.

To use this module to download data users must first create user accounts. Instructions can be
found at https://disc.gsfc.nasa.gov/data-access#python-requests and
https://cds.climate.copernicus.eu/api-how-to

In addition you can download data directly from these source:

* Hourly MERRA2 data can be downloaded directly from NASA GES DISC by selecting the
  "Subset / Get Data" link on the following webpage:
  https://disc.gsfc.nasa.gov/datasets/M2T1NXSLV_5.12.4/summary. Specific dates, variables, and
  coordinates can be selected using the OPeNDAP or GES DISC Subsetter download methods.

* Hourly ERA5 data can be downloaded using either the CDS web interface or the CDS API, as
  explained here: https://confluence.ecmwf.int/display/CKB/How+to+download+ERA5. Data for specific
  dates, variables, and coordinates can be downloaded using the CDS web interface via the "Download
  data" tab here:
  https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview.
  Instructions for using the CDS Toolbox API to download ERA5 data programatically is found here:
  https://cds.climate.copernicus.eu/toolbox/doc/how-to/1_how_to_retrieve_data/1_how_to_retrieve_data.html
  (note that the "reanalysis-era5-single-levels" dataset should generally be used).
"""

from __future__ import annotations

import re
import hashlib
import datetime
from pathlib import Path
from zipfile import ZipFile
import json

import cdsapi
import pandas as pd
import xarray as xr
import requests
from tqdm import tqdm
import numpy as np

import logging


logger = logging.getLogger()


def download_file(url: str, outfile: str | Path) -> None:
    """
    Download a file from the web, based on its url, and save to the outfile.

    Args:
        url(:obj:`str`): Url of data to download.
        outfile(:obj:`str` | :obj:`Path`): File path to which the download is saved.

    Raises:
        HTTPError: If unable to access url.
        Exception: If the request failed for another reason.
    """

    outfile = Path(outfile).resolve()
    result = requests.get(url, stream=True)

    try:
        result.raise_for_status()
        try:
            with outfile.open("wb") as f:
                for chunk in tqdm(
                    result.iter_content(chunk_size=1024 * 1024), desc="MB downloaded"
                ):
                    if chunk:
                        f.write(chunk)

            logger.info(f"Contents of {url} written to {outfile}")

        except Exception as e:
            logger.error(e)
            raise

    except requests.exceptions.HTTPError as eh:
        logger.error(eh)
        raise

    except Exception as e:
        logger.error(e)
        raise


def download_zenodo_data(record_id: int, outfile_path: str | Path) -> None:
    """
    Download data from Zenodo based on the Zenodo record_id.

    The following files will be saved to the asset data folder:

        1. record_details.json, which details the Zenodo api details.
        2. all files available for the record_id.

    Args:
        record_id(:obj:`int`): The Zenodo record id.
        outfile_path(:obj:`str` | :obj:`Path`): Path to save files to.

    """

    url_zenodo = r"https://zenodo.org/api/records/"
    r = requests.get(f"{url_zenodo}{record_id}")

    r_json = r.json()

    logger.info("======")
    logger.info("Title: " + r_json["metadata"]["title"])
    logger.info("Version: " + r_json["metadata"]["version"])
    logger.info("URL: " + r_json["links"]["latest_html"])
    logger.info("Record DOI: " + r_json["doi"])
    logger.info("License: " + r_json["metadata"]["license"]["id"])
    logger.info("======\n")

    # create outfile_path if it does not exist
    outfile_path = Path(outfile_path).resolve()
    if not outfile_path.exists():
        outfile_path.mkdir()

    # save record details to json file
    outfile = outfile_path / "record_details.json"

    with outfile.open("wb") as f:
        f.write(r.content)

    # download all files
    files = r_json["files"]
    for f in files:
        url_file = f["links"]["self"]

        outfile = outfile_path / (file_name := f["key"])

        # check if file exists
        if outfile.exists():
            # if it does check the checksum is correct
            with outfile.open("rb") as f_check:
                file_hash = hashlib.md5()
                while chunk := f_check.read(8192):
                    file_hash.update(chunk)

            if f["checksum"][4:] == file_hash.hexdigest():
                logger.info(f"File already exists: {file_name}")

            # download and unzip if the checksum isn't correct
            else:
                logger.info(f"Downloading: {file_name}")
                logger.info(f"File size: {f['size']/(1024*1024):,.2f} MB")

                download_file(url_file, outfile)

                logger.info(f"Saved to: {outfile}\n")

                if outfile.suffix == ".zip":
                    with ZipFile(outfile) as zipfile:
                        zipfile.extractall(outfile_path)

        # download and unzip if the file doesn't exist
        else:
            logger.info(f"\nDownloading: {file_name}")
            logger.info(f"File size: {f['size']/(1024*1024):,.2f} MB")

            download_file(url_file, outfile)

            logger.info(f"Saved to: {outfile}\n")

            if outfile.suffix == ".zip":
                with ZipFile(outfile) as zipfile:
                    zipfile.extractall(outfile_path)


def download_era5_monthly(
    save_pathname: str | Path = "data",
    save_filename: str = "era5_monthly",
    start_date: str = "2000-01",
    end_date: str = None,
) -> pd.DataFrame:
    """
    Get ERA5 data directly from the CDS service. This requires registration on the CDS service.
    See registration details at: https://cds.climate.copernicus.eu/api-how-to

    This function returns monthly ERA5 data from the "ERA5 monthly averaged data on single levels
    from 1959 to present" dataset. See further details regarding the dataset at:
    https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means.
    Only the 10m wind speed is downloaded here.

    As well as returning the data as a dataframe, the data is also saved as monthly NetCDF files and
    a csv file with the concatenated data. These are located in the "save_pathname" directory, with
    "save_filename" prefix. This allows future loading without download from the CDS service.

    Args:
        save_pathname(:obj:`str` | :obj:`Path`): The path where the downloaded reanalysis data will
            be saved.
        save_filename(:obj:`str`): The file name used to save the downloaded reanalysis data.
        start_date(:obj:`str`): The starting year and month that data is downloaded for. This
            should be provided as a string in the format "YYYY-MM". Defaults to "2000-01".
        end_date(:obj:`str`): The final year and month that data is downloaded for. This should be
            provided as a string in the format "YYYY-MM". Defaults to current year and most recent
            month with full data, accounting for the fact that the ERA5 monthly dataset is released
            around the the 6th of the month.

    Returns:
        df(:obj:`dataframe`): A dataframe containing time series of the requested reanalysis
            variables:
            1. windspeed_ms: the wind speed in m/s at 10m height.

    Raises:
        ValueError: If the start_date is greater than the end_date.
        Exception: If unable to connect to the cdsapi client.
    """

    logger.info("Please note access to ERA5 data requires registration")
    logger.info("Please see: https://cds.climate.copernicus.eu/api-how-to")

    # set up cds-api client
    try:
        c = cdsapi.Client()
    except Exception as e:
        logger.error("Failed to make connection to cds")
        logger.error("Please see https://cds.climate.copernicus.eu/api-how-to for help")
        logger.error(e)
        raise

    # create save_pathname if it does not exist
    save_pathname = Path(save_pathname).resolve()
    if not save_pathname.exists():
        save_pathname.mkdir()

    # get the current date minus 37 days to find the most recent full month of data
    now = datetime.datetime.now() - datetime.timedelta(days=37)

    # assign end_year to current year if not provided by the user
    if end_date is None:
        end_date = f"{now.year}-{now.month:02}"

    # convert dates to datetime objects
    start_date = datetime.datetime.strptime(start_date, "%Y-%m")
    end_date = datetime.datetime.strptime(end_date, "%Y-%m")

    # check that the start and end dates are the right way around
    if start_date > end_date:
        logger.error("The start_date should be less than or equal to the end_date")
        logger.error(f"start_date = {start_date.date()}, end_date = {end_date.date()}")
        raise ValueError("The start_date should be less than or equal to the end_date")

    # list all dates that will be downloaded
    dates = pd.date_range(start=start_date, end=end_date, freq="MS", inclusive="both")

    # See: https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=form
    # for formulating other requests from cds
    cds_dataset = "reanalysis-era5-single-levels-monthly-means"
    cds_request = {
        "product_type": "monthly_averaged_reanalysis",
        "format": "netcdf",
        "variable": [
            "10m_wind_speed",
        ],
        "year": None,
        "month": None,
        "time": ["00:00"],
    }

    # download the data
    for date in dates:
        outfile = save_pathname / f"{save_filename}_{date.year}{date.month:02}.nc"

        if not outfile.is_file():
            logger.info(f"Downloading ERA5: {outfile}")

            try:
                cds_request.update({"year": date.year, "month": date.month})
                c.retrieve(cds_dataset, cds_request, outfile)

            except Exception as e:
                logger.error(f"Failed to download ERA5: {outfile}")
                logger.error(e)

def get_nearest_slice(sel,ds_nc,num_points_to_include=1):

    # Get the longitude and latitude values of the nearest point
    nearest_lon = sel['longitude'].values
    nearest_lat = sel['latitude'].values

    # Find the indices of the nearest coordinates in the original dataset
    nearest_indices  = np.where((ds_nc['longitude'] == nearest_lon) & (ds_nc['latitude'] == nearest_lat))

    # nearest_indices is a tuple, so you can extract the indices
    lon_index, lat_index = nearest_indices[0][0], nearest_indices[1][0]

    # Calculate the start and end indices for the slice
    lon_start = max(0, lon_index - num_points_to_include)
    lon_end = min(ds_nc.sizes['longitude'], lon_index + num_points_to_include + 1)  # Add 1 to include the endpoint

    lat_start = max(0, lat_index - num_points_to_include)
    lat_end = min(ds_nc.sizes['latitude'], lat_index + num_points_to_include + 1)  # Add 1 to include the endpoint

    # Extract the slice from the original dataset
    nearest_slice = ds_nc.isel(longitude=slice(lon_start, lon_end), latitude=slice(lat_start, lat_end))

    return nearest_slice

def get_era5_monthly(
    lat: float,
    lon: float,
    save_filename: str,
    save_pathname: str | Path = "data/assets",
    data_pathname: str | Path = "data",
    data_filename: str = "era5_monthly",
    start_date: str = "2000-01",
    end_date: str = None,
) -> pd.DataFrame:
    """ 
        Args:
        lat(:obj:`float`): Latitude in WGS 84 spatial reference system (decimal degrees).
        lon(:obj:`float`): Longitude in WGS 84 spatial reference system (decimal degrees).
        save_pathname(:obj:`str` | :obj:`Path`): The path where the downloaded reanalysis data will
            be saved.
        save_filename(:obj:`str`): The file name used to save the downloaded reanalysis data.
        start_date(:obj:`str`): The starting year and month that data is downloaded for. This
            should be provided as a string in the format "YYYY-MM". Defaults to "2000-01".
        end_date(:obj:`str`): The final year and month that data is downloaded for. This should be
            provided as a string in the format "YYYY-MM". Defaults to current year and most recent
            month with full data, accounting for the fact that the ERA5 monthly dataset is released
            around the the 6th of the month.

    Returns:
        df(:obj:`dataframe`): A dataframe containing time series of the requested reanalysis
            variables:
            1. windspeed_ms: the wind speed in m/s at 10m height.

    Raises:
        ValueError: If the start_date is greater than the end_date.
        Exception: If unable to connect to the cdsapi client.
    """
    
    # download the monthly era5 data
    download_era5_monthly(save_pathname=data_pathname,
        save_filename=data_filename,
        start_date=start_date,
        end_date=end_date)

    # get the saved data
    ds_nc = xr.open_mfdataset(f"{data_pathname}/{data_filename}*.nc")

    # rename variables to conform with OpenOA
    ds_nc = ds_nc.rename_vars(
        {"si10": "windspeed_ms"}
    )

    # select the central node only for now
    if "expver" in ds_nc.dims:
        sel = ds_nc.sel(expver=1, latitude=lat, longitude=lon, method="nearest")
    else:
        sel = ds_nc.sel(latitude=lat, longitude=lon, method="nearest")

    # now take the surrounding nearest nodes as well
    nearest_slice = get_nearest_slice(sel,ds_nc,num_points_to_include=1)

    # convert to a pandas dataframe
    df = nearest_slice.to_dataframe().unstack(["latitude","longitude"])["windspeed_ms"]

    # rename columns based on their coordinates
    df.columns = [str(item) for item in df.columns.values]

    # rename the index to match other datasets
    df.index.name = "datetime"

    # drop any empty rows
    df = df.dropna()

    # crop time series to only the selected time period
    df = df.loc[start_date:end_date]

    # save to csv for easy loading as required
    # create save_pathname if it does not exist
    save_pathname = Path(save_pathname).resolve()
    if not save_pathname.exists():
        save_pathname.mkdir()

    df.to_csv(save_pathname / f"{save_filename}.csv", index=True)

    return df

def download_merra2_monthly(
    save_pathname: str | Path = "data",
    save_filename: str = "era5_monthly",
    start_date: str = "2000-01",
    end_date: str = None,
) -> pd.DataFrame:
    """
    Get MERRA2 data directly from the NASA GES DISC service, which requires registration on the
    GES DISC service. See: https://disc.gsfc.nasa.gov/data-access#python-requests.

    This function returns monthly MERRA2 data from the "M2IMNXLFO" dataset. See further details
    regarding the dataset at: https://disc.gsfc.nasa.gov/datasets/M2IMNXLFO_5.12.4/summary.
    Only surface wind speed, temperature and surface pressure are downloaded here.

    As well as returning the data as a dataframe, the data is also saved as monthly NetCDF files
    and a csv file with the concatenated data. These are located in the "save_pathname" directory,
    with "save_filename" prefix. This allows future loading without download from the CDS service.

    Args:
        save_pathname(:obj:`str` | :obj:`Path`): The path where the downloaded reanalysis data will
            be saved.
        save_filename(:obj:`str`): The file name used to save the downloaded reanalysis data.
        start_date(:obj:`str`): The starting year and month that data is downloaded for. This
            should be provided as a string in the format "YYYY-MM". Defaults to "2000-01".
        end_date(:obj:`str`): The final year and month that data is downloaded for. This should be
            provided as a string in the format "YYYY-MM". Defaults to current year and most recent
            month.

    Returns:
        df(:obj:`dataframe`): A dataframe containing time series of the requested reanalysis
            variables:
            1. windspeed_ms: the surface wind speed in m/s.

    Raises:
        ValueError: If the start_year is greater than the end_year.
    """

    # base url containing the monthly data set M2IMNXLFO
    base_url = r"https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2_MONTHLY/M2IMNXLFO.5.12.4/"

    # create save_pathname if it does not exist
    save_pathname = Path(save_pathname).resolve()
    if not save_pathname.exists():
        save_pathname.mkdir()

    # get the current date minus 37 days to find the most recent full month of data
    now = datetime.datetime.now() - datetime.timedelta(days=37)

    # assign end_year to current year if not provided by the user
    if end_date is None:
        end_date = f"{now.year}-{now.month:02}"

    # convert dates to datetime objects
    start_date = datetime.datetime.strptime(start_date, "%Y-%m")
    end_date = datetime.datetime.strptime(end_date, "%Y-%m")

    # check that the start and end dates are the right way around
    if start_date > end_date:
        logger.error("The start_date should be less than or equal to the end_date")
        logger.error(f"start_date = {start_date.date()}, end_date = {end_date.date()}")
        raise ValueError("The start_date should be less than or equal to the end_date")

    # list all dates that will be downloaded
    dates = pd.date_range(start=start_date, end=end_date, freq="MS", inclusive="both")

    # check what years need downloading
    years = []
    for date in dates:
        outfile = save_pathname / f"{save_filename}_{date.year}{date.month:02}.nc"
        if not outfile.is_file():
            years.append(date.year)

    # list all years that will be downloaded
    years = list(set(years))

    # download the data
    for year in years:
        # get the file names from the GES DISC site for the year
        result = requests.get(f"{base_url}{year}")

        files = re.findall(r"(>MERRA2_\S+.nc4)", result.text)
        files = list(dict.fromkeys(files))
        files = [x[1:] for x in files]

        # download each of the files and save them
        for f in files:
            outfile = save_pathname / f"{save_filename}_{f.split('.')[-2]}.nc"

            if not outfile.is_file():
                # download each file
                url = f"{base_url}{year}/{f}" + r".nc4?SPEEDLML,time,lat,lon"
                download_file(url, outfile)
                    
def get_merra2_monthly(
    lat: float,
    lon: float,
    save_filename: str,
    save_pathname: str | Path = "data/assets",
    data_pathname: str | Path = "data",
    data_filename: str = "merra2_monthly",
    start_date: str = "2000-01",
    end_date: str = None,
) -> pd.DataFrame:
    """
    Get MERRA2 data directly from the NASA GES DISC service, which requires registration on the
    GES DISC service. See: https://disc.gsfc.nasa.gov/data-access#python-requests.

    This function returns monthly MERRA2 data from the "M2IMNXLFO" dataset. See further details
    regarding the dataset at: https://disc.gsfc.nasa.gov/datasets/M2IMNXLFO_5.12.4/summary.
    Only surface wind speed, temperature and surface pressure are downloaded here.

    As well as returning the data as a dataframe, the data is also saved as monthly NetCDF files
    and a csv file with the concatenated data. These are located in the "save_pathname" directory,
    with "save_filename" prefix. This allows future loading without download from the CDS service.

    Args:
        lat(:obj:`float`): Latitude in WGS 84 spatial reference system (decimal degrees).
        lon(:obj:`float`): Longitude in WGS 84 spatial reference system (decimal degrees).
        save_pathname(:obj:`str` | :obj:`Path`): The path where the downloaded reanalysis data will
            be saved.
        save_filename(:obj:`str`): The file name used to save the downloaded reanalysis data.
        start_date(:obj:`str`): The starting year and month that data is downloaded for. This
            should be provided as a string in the format "YYYY-MM". Defaults to "2000-01".
        end_date(:obj:`str`): The final year and month that data is downloaded for. This should be
            provided as a string in the format "YYYY-MM". Defaults to current year and most recent
            month.

    Returns:
        df(:obj:`dataframe`): A dataframe containing time series of the requested reanalysis
            variables:
            1. windspeed_ms: the surface wind speed in m/s.

    Raises:
        ValueError: If the start_year is greater than the end_year.
    """

    logger.info("Please note access to MERRA2 data requires registration")
    logger.info("Please see: https://disc.gsfc.nasa.gov/data-access#python-requests")

    download_merra2_monthly(save_pathname=data_pathname,
        save_filename=data_filename,
        start_date=start_date,
        end_date=end_date)

    # get the saved data
    ds_nc = xr.open_mfdataset(f"{data_pathname}/{data_filename}*.nc")

    # rename variables to conform with OpenOA
    ds_nc = ds_nc.rename_vars(
        {"SPEEDLML": "windspeed_ms"}
    )

    # rename coords to match across code
    ds_nc = ds_nc.rename(
        {"lon":"longitude",
         "lat":"latitude"}
    )

    # wrap -180..179 to 0..359    
    ds_nc.coords["longitude"] = np.mod(ds_nc["longitude"], 360)

    # sort the data
    ds_nc = ds_nc.reindex({ "longitude" : np.sort(ds_nc["longitude"])})
                   
    # select the central node only for now
    sel = ds_nc.sel(latitude=lat, longitude=lon, method="nearest")

    # now take the surrounding nearest nodes as well
    nearest_slice = get_nearest_slice(sel,ds_nc,num_points_to_include=1)

    # convert to a pandas dataframe
    df = nearest_slice.to_dataframe().unstack(["latitude","longitude"])["windspeed_ms"]

    # rename columns based on their coordinates
    df.columns = [str(item) for item in df.columns.values]
    
    # rename the index to match other datasets
    df.index.name = "datetime"

    # drop any empty rows
    df = df.dropna()

    # crop time series to only the selected time period
    df = df.loc[start_date:end_date]

    # save to csv for easy loading as required
    # create save_pathname if it does not exist
    save_pathname = Path(save_pathname).resolve()
    if not save_pathname.exists():
        save_pathname.mkdir()

    df.to_csv(save_pathname / f"{save_filename}.csv", index=True)

    return df


def get_best_correlated_resource(df_resource,df_generation):
        
    # combine the dataframe with the gross production to determine the best node for use
    correlation_data  = pd.concat([df_generation,df_resource],axis=1)

    # determine R2, note correlation coefficient is R, not R-squared, so hence need to apply power of 2
    all_correlations = correlation_data.corr()[df_generation.name]**2

    # find the best correlation of wind speed to production
    best_correlation_coef = all_correlations.nlargest(2).iloc[1]

    # get the wind speed data for the best correlation
    best_correlation_column = list(all_correlations.index[all_correlations==best_correlation_coef])[0]

    df_resource = df_resource[best_correlation_column]

    df_resource = df_resource.to_frame()

    return df_resource


def request_EIA(base_url,params):
    """
    Make a url request to the EIA API_v2
    See: https://www.eia.gov/opendata/documentation.php
    """

    try:
        r = requests.get(base_url, params=params)
        r.raise_for_status()
    except requests.exceptions.HTTPError as e:
        print(e.response.text)
        
    json_output = json.loads(r.text)

    data = json_output["response"]["data"]

    df = pd.DataFrame.from_records(data)

    return(df)



def get_EIA_plant_information(EIA_API_KEY):
    """
    Build and make a request for EIA plant information
    See: https://www.eia.gov/opendata/documentation.php
    """

    base_url = "https://api.eia.gov/v2/electricity/operating-generator-capacity/data/?"

    params = {"frequency":"monthly",
            "data[0]":"nameplate-capacity-mw",
            "data[1]":"longitude",
            "data[2]":"latitude",
            "data[3]":"operating-year-month",
            "data[4]":"planned-retirement-year-month",
            "facets[energy_source_code][]":"WND",
            "start":"2023-05",
            "end":"2023-05",
            "sort[0][column]":"period",
            "sort[0][direction]":"desc",
            "offset":"0",
            "length":"5000",
            "api_key":EIA_API_KEY}

    df = request_EIA(base_url,params)

    return df



# request EIA time-series data
def get_EIA_plant_generation(EIA_API_KEY,plantCode):
    """
    Build and make a request for EIA plant generation
    See: https://www.eia.gov/opendata/documentation.php
    """
    base_url = "https://api.eia.gov/v2/electricity/facility-fuel/data/?"

    params = {"frequency":"monthly",
            "data[0]":"generation",
            "data[1]":"gross-generation",
            "facets[fuel2002][]":"WND",
            "facets[plantCode][]":str(plantCode),
            "facets[primeMover][]":"ALL",
            "sort[0][column]":"period",
            "sort[0][direction]":"desc",
            "offset":"0",
            "length":"5000",
            "api_key":EIA_API_KEY}
          
    df = request_EIA(base_url,params)

    return df