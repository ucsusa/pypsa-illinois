import numpy as np
import sys
import os
from pathlib import Path
import pandas as pd
from shapely import Point
from dotenv import load_dotenv

load_dotenv("../.env")

personal_data = {'api_key':os.environ.get('NREL_API_KEY'),
                 'name':os.environ.get('NAME').replace(' ', '+'),
                 'reason':os.environ.get('REASON').replace(' ', '+'),
                 'affiliation':os.environ.get('AFFIL').replace(' ', '+'),
                 'email':os.environ.get('EMAIL'),
                 'mailing_list':os.environ.get('MAILING_LIST')
                 }

parameters = {'lon':40.09,
              'lat':-88.26,
              'year':2019,
              'leap_day':'false',
              'selector':'POINT',
              'utc':'true',
              'interval':'60',
              'attr_list':['ghi']}


def make_wkt(selector, lat, lon):
    """
    This function generates a well known text (wkt)
    string for use with the NREL API.

    Parameters
    ----------
    selector : String
        Indicates how you want to access the data.
        Accepts: 'POINT', 'MULTIPOINT', 'POLYGON'
    lat : Float or List
        The latitude, or set of latitudes, of interest.
    lon : Float or List
        The longitude, or set of longitudes, of interest.

    Returns
    -------
    wkt : String
        The well known text string.
    """

    method = selector.upper()

    if method == 'POINT':
        wkt = '{method}({lat}%20{lon})'.format(method=method,
                                               lon = lon,
                                               lat = lat)
    else:
        try:
            combinations = [f'{i}%20{j}' for i, j in zip(lon, lat)]
        except ValueError:
            "Longitude and Latitudes are different sizes."
        coord_list = ('%2C').join(combinations)
        wkt = "{method}({coordinates})".format(method=method,
                                                             coordinates=coord_list)

    return wkt


def make_csv_url(parameters, personal_data=personal_data, kind='solar'):
    """
    This function generates a url to access renewable energy
    data through the NREL API. This function requires your
    personal API key. If you want to sign up with NREL and
    get your own API key, visit this website:
    https://developer.nrel.gov/signup/

    Parameters
    ----------
    parameters : dictionary
        This dictionary contains all of the information about
        the dataset you wish to download. Required values are:
        * 'lon' : The longitude of the location of interest, float or list
        * 'lat' : The latitude of the location of interest, float or list
        * 'year' : The year of interest, integer
        * 'leap_day' : Boolean. If 'true', includes leap day.
        * 'utc' : Boolean. If 'true', uses the time at the location.
                  Else, uses your local time.
        * 'selector' : String. Indicates how you want to access the
                       data. Accepts: 'POINT', 'MULTIPOINT', 'POLYGON'
        * 'interval' : String or integer. The desired resolution in minutes.
                       For solar, resolutions are: 30 or 60
                       For wind, resolutions are: 5, 10, 15, 30, or 60
        * 'attr_list' : List. The list of desired columns in your dataset.

    personal_data : dictionary
        This dictionary contains all of the information about
        the user seeking data. Required values are:
        * 'api_key' : The API key generated for you by NREL
        * 'name' : The name of the user
        * 'reason' : How you will use the data
        * 'affiliation' : The institution you are affiliated with
        * 'mailing_list' : 'true' or 'false.' If true, will add you
                            to the mailing list.
        * 'email' : The email of the user

    kind : string
        Indicates what kind of data should be downloaded. Currently
        accepts 'solar' or 'wind.' Default is 'solar'.

    Returns
    -------
    url : string
        A well formed URL to access NREL data.
    """
    databases = {'wind': 'api/wind-toolkit/v2/wind/wtk-download',
                 'solar': 'api/solar/nsrdb_psm3_download'}

    # Which NREL database should be accessed?
    db_to_access = databases[kind.lower()]

    # Generate the WKT string
    wkt = make_wkt(parameters['selector'],
                   parameters['lon'],
                   parameters['lat'])

    # Make sure the dictionary values are properly formatte
    name = personal_data['name'].replace(' ', '+')
    affiliation = personal_data['affiliation'].replace(' ', '+')
    reason=personal_data['reason'].replace(' ', '+')
    attributes = (',').join(parameters['attr_list'])
    leap_day = str(parameters['leap_day']).lower()
    year = int(parameters['year'])
    interval = int(parameters['interval'])
    utc = str(parameters['utc']).lower()
    mailing_list = str(personal_data['mailing_list']).lower()
    email = personal_data['email']
    key = personal_data['api_key']

    # print(wkt)
    # print(db_to_access)

    # Generate string
    url = ("https://developer.nrel.gov/{db}.csv?wkt={wkt}&names={year}"
           "&leap_day={leap}&interval={interval}&utc={utc}&full_name={name}"
           "&email={email}&affiliation={affiliation}&mailing_list={mailing_list}"
           "&reason={reason}&api_key={api}&attributes={attr}").format(db=db_to_access,
                                                                      wkt=wkt,
                                                                      year=year,
                                                                      leap=leap_day,
                                                                      interval=interval,
                                                                      utc=utc,
                                                                      name=name,
                                                                      email=email,
                                                                      affiliation=affiliation,
                                                                      mailing_list=mailing_list,
                                                                      reason=reason,
                                                                      api=key,
                                                                      attr=attributes)

    return url