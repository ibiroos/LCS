# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 11:18:57 2023

@author: josmo
"""

# Data processing
import numpy as np
# For additional mathematical functions.
import math as m
# For handling data in various formats
import pandas as pd
# For graphical representations
import matplotlib.pyplot as plt
# For dates and times
from datetime import datetime, timedelta
# For directories and system files
import os
# For secure file downloading
import urllib
# To show the progress of the download
from tqdm import tqdm
# For searching for files in a specific directory pattern
import glob
# For handling labeled multidimensional data
import xarray as xr
# For handling regular expressions
import re
# For handling images.
from PIL import Image, ImageTk
# To create graphical user interfaces (GUI).
import tkinter as tk
# To color printed text in the console
from termcolor import colored
# We import functions from other programs
from common.readers.reader_factory import read_factory
from common import read_input
from common.boundarybox import BoundaryBox
from common.drawcurrents import drawcurrents
from common.map import Map
from common.drawmap import*


def distancia_esferica(lat1, lon1, lat2, lon2):
    """
    Calculates the distance between two points on a sphere given their latitude and longitude coordinates.
    
    Parameters:
    -----------
    lat1 : float
        Latitude of the first point in degrees.
    lon1 : float
        Longitude of the first point in degrees.
    lat2 : float
        Latitude of the second point in degrees.
    lon2 : float
        Longitude of the second point in degrees.
    
    Returns:
    --------
    d : float
        Distance between the two points in kilometers.
    """
    
    r_tierra = 6371 # Radius of the Earth in kilometers
    d_lat = np.radians(lat2 - lat1)
    d_lon = np.radians(lon2 - lon1)
    a = np.sin(d_lat/2)**2 + np.cos(np.radians(lat1))*np.cos(np.radians(lat2))*np.sin(d_lon/2)**2
    #c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    c = 2*np.arcsin(np.sqrt(a))
    d = r_tierra * c # distance
    return d


def distancia_pitagoras(lat_input, lon_input, lats, lons):
    """
    Calculate the distance between two points using the Pythagorean theorem.
    
    Parameters:
    -----------
    lat_input: float
        Latitude of the input point in degrees.
    lon_input: float
        Longitude of the input point in degrees.
    lats: numpy.ndarray
        Array of latitudes of the points to calculate the distance to, in degrees.
    lons: numpy.ndarray
        Array of longitudes of the points to calculate the distance to, in degrees.
    
    Returns:
    --------
    d: numpy.ndarray
        Array of distances in kilometers between the input point and each of the points specified by `lats` and `lons`.

    """
    
    R_t = 6371 # Radius of the Earth in kilometers
    dy = (np.pi/180)*(lat_input-lats)*R_t
    dx = (np.pi/180)*(lon_input-lons)*R_t*np.cos(np.radians(lats))
    d = np.sqrt(dx**2+dy**2) # distance
    #d = (2*np.pi*6371/360)*np.sqrt((lats-lat_input)**2 + (lons-lon_input)**2)
    return d


def distancia_min(lat_input, lon_input, lats, lons):
    """
    Calculates the minimum distance between an input point and a set of points defined by latitudes and longitudes.

    Parameters
    ----------
    lat_input : float
        Latitude of the input point in degrees.
    lon_input : float
        Longitude of the input point in degrees.
    lats : array-like
        Array of latitudes for the set of points in degrees.
    lons : array-like
        Array of longitudes for the set of points in degrees.

    Returns
    -------
    index_min: int
        Index of the station closest to the reference point
    dist_min : float
        Approximate distance in kilometers between the reference point and the closest station.
        
    If the function does not find a valid minimum distance (for example, if the reference point is too far away),
    the distance and index obtained using the spherical distance method are returned.
    """
    
    dist_esf = distancia_esferica(lat_input, lon_input, lats, lons)
    dist_pit = distancia_pitagoras(lat_input, lon_input, lats, lons)
    # Find the minimum distance and the corresponding index to identify the closest point
    indice_min_esf = np.argmin(dist_esf)
    indice_min_pit = np.argmin(dist_pit)
    distancia_min_esf = dist_esf[indice_min_esf]
    distancia_min_pit = dist_pit[indice_min_pit]
    # If the minimum distances calculated by both methods are approximately the same and the indices of the points with the minimum distances are the same,
    #then the function returns the index and the average of the two distances
    if (indice_min_esf==indice_min_pit) and np.round(distancia_min_esf,0)==np.round(distancia_min_pit,0):
        index_min = indice_min_esf
        dist_min = np.round((distancia_min_esf+distancia_min_pit)/2,2)
        return (index_min,dist_min)
    else:
        print('The spill is too far away from the observation stations.')
        return (indice_min_esf,distancia_min_esf) # returns the index and distance calculated using the spherical distance method
    

def estaciones_proximas_rad_solar(Lat,Lon):
    """
    Given a latitude and longitude coordinate, this function finds the nearest weather observation station for 
    solar radiation data, based on a pre-defined dataset of station locations. 
   
    Parameters
    ----------
    Lat : float
        Latitude of the point of interest.
    Lon : float
        Longitude of the point of interest.

    Returns
    -------
    index_rad: int
        Index of the closest weather station in the MeteoGz_Rad dataframe.
    dist_rad : float
        distance between the point of interest and the closest weather observation
        station, in kilometers.

    """
    
    # Path to the data file
    ruta = os.getcwd()+'\\common\\'
    datafile = 'Rias_UWWTP.xlsx'
    
    # Read in the data
    MeteoGz_Rad = pd.read_excel(ruta + datafile, sheet_name= 'MeteoGz_Rad')
    Data_Links = pd.read_excel(ruta + datafile, sheet_name= 'Data_Links')
    
    # Coordinates of the MeteoGz weather stations for solar radiation data
    MeteoGz_Rad_lats = MeteoGz_Rad.Lat.values
    MeteoGz_Rad_lons = MeteoGz_Rad.Lon.values
    
    # Find the nearest station and its distance
    index_rad,dist_rad = distancia_min(Lat, Lon, MeteoGz_Rad_lats, MeteoGz_Rad_lons)
    
    # Display relevant information
    print('You can check the following links to obtain solar radiation data:\n')
    print('- Solar Radiation Data from MeteoGalicia:\n', Data_Links.values[2][0])
    print("\nThe nearest weather observation station to ({}, {}) ".format(Lat,Lon)+
          "is {} ".format(MeteoGz_Rad.values[index_rad][0])+
          "located at ({}, {}) and approximately {} km away.".format(
              MeteoGz_Rad_lats[index_rad],MeteoGz_Rad_lons[index_rad],dist_rad))
    print('-'*40)
    return (index_rad, dist_rad)


def estaciones_proximas_TS(Lat,Lon):
    """
    Find the closest buoy station to a given latitude and longitude, and return its index and distance.

    Parameters
    ----------
    Lat : float
        Latitude of the point of interest.
    Lon : float
        Longitude of the point of interest.

    Returns
    -------
    index_TS: int
        Index of the closest buoy station in the MeteoGZ_Buoys dataframe.
    dist_TS : float
        Distance (in km) between the input coordinates and the closest buoy station.

    """
    
    # Path of the file of interest
    path = os.getcwd()+'\\common\\'
    datafile = 'Rias_UWWTP.xlsx'
    
    # Read the file
    MeteoGZ_Buoys = pd.read_excel(path + datafile, sheet_name= 'MeteoGZ_Buoys')
    Data_Links = pd.read_excel(path + datafile, sheet_name= 'Data_Links')

    # Coordinates of the MeteoGz stations where to consult salt+temp
    MeteoGZ_Buoys_lats = MeteoGZ_Buoys.Lat.values
    MeteoGZ_Buoys_lons = MeteoGZ_Buoys.Lon.values
    
    # Get the closest position and distance
    index_TS,dist_TS = distancia_min(Lat, Lon, MeteoGZ_Buoys_lats, MeteoGZ_Buoys_lons)
    
    # Show the relevant information
    print('You can check the following links to obtain the T and S data:\n')
    print('- T & S data from MeteoGalicia:\n', Data_Links.values[3][0])
    # Show the result
    print("\nThe nearest observation station to ({}, {}) ".format(Lat,Lon)+
          "is {} ".format(MeteoGZ_Buoys.values[index_TS][0])+
          "located at ({}, {}) at a distance of approximately {} km.".format(
              MeteoGZ_Buoys_lats[index_TS],MeteoGZ_Buoys_lons[index_TS],dist_TS))
    print('-'*40)
    
    print('\nYou can also check data from other links such as:')
    print('- T & S data from the CTD profiles of INTECMAR:', Data_Links.values[4][0]) 
    return(index_TS,dist_TS)


def tryint(s):
    """
    Attempts to convert the input argument `s` into an integer. If it is not possible, returns `s` unchanged.
    """
    
    try:
        return int(s)
    except ValueError:
        return s


def alphanum_key(s):
    """
    Turns a string `s` into a list of string and number chunks by splitting it at every number. Returns the resulting list.
    
    Example:
    >>> alphanum_key("z23a")
    ["z", 23, "a"]
    """
    
    return [tryint(c) for c in re.split('([0-9]+)', s)]


def human_sort(l):
    """
    Sorts a list `l` in the way that humans expect, by splitting each string element into chunks of strings and numbers
    and sorting them according to the values of the chunks. The sorting is done in place.
    """
    
    l.sort(key=alphanum_key)

    
def download_MOHID_files(fecha,key):
    """
    Downloads MOHID files for a given date and key

    Parameters
    ----------
    fecha : datetime.datetime
        Date for which to download files.
    key : list
        List with two elements, where the first element is the MOHID key and the second element is the file name.

    Returns
    -------
    str
        'Done' if all files were downloaded successfully, or None if there was an HTTP error during file download.
        
    """
    
    # Set start and end dates
    delta_time = timedelta(days=1)
    data = fecha - delta_time
    init_time = datetime(data.year, data.month, data.day)
    end_time = datetime(data.year, data.month, data.day)
    
    # Loop through the dates
    date = init_time
    print('-'*30)
    while (date <= end_time):
        print('Downloading file...')
        print('-'*30) 
        
        # Extract year, month, and day from date object
        year = date.strftime("%Y")
        month = date.strftime("%m")
        day = date.strftime("%d")
        
        # Define the file URLs and names
        file_URL = 'http://mandeo.meteogalicia.gal/thredds/fileServer/mohid_' + key[0] + '/fmrc/files/%s/MOHID_' + key[1] + '_%s_0000.nc4' \
                   %(year+month+day,year+month+day)
        file_name_out = "MOHID_" + key[1] + "_%s_0000.nc4" %(year+month+day)

        print(file_URL, file_name_out)
        
        try:
            # Download the file
            #urllib.request.urlretrieve(file_URL,file_name_out)
            with tqdm(unit='B', unit_scale=True, miniters=1, desc=file_name_out) as tqdm_instance:
                urllib.request.urlretrieve(file_URL, file_name_out, reporthook=lambda block_num, block_size, total_size: 
                tqdm_instance.update(block_num * block_size - tqdm_instance.n))            
        except urllib.error.HTTPError as e:
            # If there is an HTTP error, print the error and return None
            print(f' 02 Error downloading file: {e}')
            return None
        else:
            # If the file was downloaded successfully, print a success message
            print(file_name_out+' has been downloaded successfully.')
        finally:
            # Advance the date and print separator lines
            date = date + delta_time
            print('-'*30)
    
    # If all files were downloaded successfully, return 'Done'
    return 'Done'


def delete_nc4_files():
    """
    Deletes all .nc4 files in the current working directory.
    This function searches for all files in the current working directory 
    that have a .nc4 file extension and deletes them if they exist.

    Returns:
        None. 
        The function does not return anything, it simply deletes the files.
        
    """
    
    # Get the current working directory path and create a 'ruta' variable with a backslash at the end
    ruta = os.getcwd()+'\\'
    
    # Use the glob module to find all files with the ".nc4" extension in the 'ruta' directory
    all_files = glob.glob(os.path.join(ruta +'*.nc4'))
    
    # Iterate over the list of files and delete them using os.remove()
    for file_path in all_files:
        if os.path.exists(file_path):
            os.remove(file_path)
            print(f"The file {file_path} has been successfully deleted.")
        else:
            print(f"The file {file_path} does not exist.")

            
def extract_TS_nc_file(time,lat,lon,depth):
    """
    Extracts temperature and salinity values from NetCDF files in the current working directory.


    Parameters
    ----------
    time : str or datetime.datetime
        Time at which temperature and salinity values are desired.
    lat : float
        Latitude at which temperature and salinity values are desired.
    lon : float
        Longitude at which temperature and salinity values are desired.
    depth : float
        Depth at which temperature and salinity values are desired.

    Returns
    -------
    T: float
        Temperature value.
    S : float
        Salinity value.

    """
    
    # Get the current working directory and append a backslash
    ruta = os.getcwd()+'\\'
    # Find all files in the directory with a .nc4 extension
    all_files = glob.glob(os.path.join(ruta +'*.nc4'))
    # Sort the files in human-readable order
    human_sort(all_files)

    # Open the first file in the sorted list as an xarray dataset
    ds = xr.open_dataset(all_files[0])
    # Convert the time values in the dataset to pandas datetime objects
    datas = pd.to_datetime(ds.time.values)

    # Check if the provided time is in the dataset
    if str(time) in datas.astype(str):
        # If so, get the index of the corresponding time value
        time_index = np.where(datas == time)[0][0]
    else:
        # If not, print an error message
        print("The provided time is not in the dataset.")

    # Find the index of the latitude value closest to the provided latitude
    lat_index = np.abs(lat - ds.lat.values).argmin()
    # Find the index of the longitude value closest to the provided longitude
    lon_index = np.abs(lon - ds.lon.values).argmin()
    # Find the index of the depth value closest to the provided depth
    depth_index = np.abs(depth - ds.depth.values).argmin()

    # Get the temperature value at the specified time, latitude, longitude, and depth
    T = ds.temp.values[time_index][depth_index][lat_index][lon_index]
    # Get the salinity value at the specified time, latitude, longitude, and depth
    S = ds.salt.values[time_index][depth_index][lat_index][lon_index]

    # Return a tuple of the temperature and salinity values
    return (T, S)
  
def T90(T,S,iz):
    """
    Parameters
    ----------
    T : float
        Surrounding temperature in ºC.
    S : float
        Surrounding water salinity (psu).
    iz : float
        Light radiation (W/m2) at deph z (m).

    Returns
    -------
    T90 : float
        T90 is the time in which 90% of E.Coli population is no longer detectable
    """
    
    k = 2.533*(1.04**(T-20))*(1.012**S)+0.113*iz
    T90 = (2.303/k)*24 #T90 in hours
    return (T90)

def start_T90():
    """
    Calculates the 90% mortality time of E.Coli by taking input from the user for Temperature (ºC), Salinity (psu), and 
    Solar Radiation (W/m2) at a given depth z. The user needs to select the study Ria and provide the spill position 
    (Lat, Lon), and the date in the format YYYY-MM-DD HH. 
    
    Returns None if the user selects to exit the program or if there is an error in the inputs.
    
    """
    
    # Prints the heading of the function
    print(colored('\n - T90 CALCULATION: \n', attrs=['underline']))
    #print('\n - T90 CALCULATION: \n')
    
    # Prompts the user for input and describes the function
    print('This function aims to calculate the 90% mortality time of E.Coli.\n'
      'For this, we will need the values of Temperature (ºC), Salinity (psu) and\n'
      'Solar Radiation (W/m2) at a depth z.\n'
      'The first step is to select the study Ria, provide the spill position (Lat, Lon), and the date:')
    
    # List of keywords to select the Ria
    keyword = [['arousa','Arousa'],['vigo','Vigo'],['noiamuros','NoiaMuros','Noiamuros'],
               ['artabro','Artabro']]
    # Rias list
    Rias = ['Ria de Arousa', 'Ria de Pontevedra-Vigo', 'Ria de Noia-Muros', 'Ria do Artabro']
    
    nam = None
    while nam not in Rias:
        # Displays the available options for Rias and prompts the user to select one
        print('The available Rias for this study are:')
        for i in Rias:
            print(f'- {i}')
        # Ask user for Ria selection
        nam = input('Select one of the above options or type "exit" to exit: ')
        if nam == 'exit':
            return None
            break
        if nam not in Rias:
            print(f'Error: "{nam}" is not a valid option. Please try again.')

    # Select Ria based on user input
    key = keyword[Rias.index(nam)]
    # Initialize control variable
    continuar = True
    
    # Prompts the user for latitude
    while continuar:
        try:
            entrada = input("Enter the latitude in decimal format or type 'exit' to exit: ")
            if entrada == 'exit':
                continuar = False
                break
            else:
                lat = float(entrada)
                if -90 <= lat <= 90:
                    #if not 42.383 <= lat <= 42.680:
                        #print('La latitud no pertenece a la Ría de Arousa')
                    break
                else:
                    print("Latitude must be between -90 and 90 degrees.")
        except ValueError:
            print("Please enter a valid decimal number.")
    
    # Prompts the user for longitude
    while continuar:
        try:
            entrada = input("Enter the longitude in decimal format or type 'exit' to exit: ")
            if entrada == 'exit':
                continuar = False
                break
            else:
                lon = float(entrada)
                if -180 <= lon <= 180:
                    break
                else:
                    print("Longitude must be between -180 and 180 degrees.")
        except ValueError:
            print("Please enter a valid decimal number.")
    
    # Ask user for date
    while continuar:
        try:
            entrada = input("Enter the date in YYYY-MM-DD HH format or type 'exit' to exit: ")
            if entrada == 'exit':
                continuar = False
                break
            else:
                fecha = datetime.strptime(entrada, '%Y-%m-%d %H')
                break
        except ValueError:
            print("Please enter a valid date format (YYYY-MM-DD HH).")
    
    # Asks the user if they have the necessary data to calculate T90, or if they want to exit
    while continuar:
        try:
            entrada = input("Do you know the data to calculate T90? [y/n] or type 'exit' to exit: ")
            if entrada == 'exit':
                continuar = False
                break
            elif entrada.lower() in ['y', 'n']:
                answer1 = entrada.lower()
                
                # Ask for the value of z
                while answer1=='y':
                    z_input = input("Enter the value of the depth z (m) or type 'exit' to exit: ")
                    if z_input == 'exit':
                        answer1 = False
                        continuar = False
                        break
                    else:
                        try:
                            z = float(z_input)
                            break
                        except ValueError:
                            print("Please enter a valid number.")
                            
                # Ask for the value of T
                while answer1=='y':
                    T_input = input("Enter the value of temperature T (ºC) or type 'exit' to exit: ")
                    if T_input == 'exit':
                        answer1 = False
                        continuar = False
                        break
                    else:
                        try:
                            T = float(T_input)
                            break
                        except ValueError:
                            print("Please enter a valid number.")
                            
                # Ask for the value of S
                while answer1=='y':
                    S_input = input("Enter the value of salinity S (psu) or type 'exit' to exit: ")
                    if S_input == 'exit':
                        answer1 = False
                        continuar = False
                        break
                    else:
                        try:
                            S = float(S_input)
                            break
                        except ValueError:
                            print("Please enter a valid number.")
                            
                # Ask for the value of iz
                while answer1=='y':
                    iz_input = input("Enter the value of solar radiation iz (W/m2) or type 'exit' to exit: ")
                    if iz_input == 'exit':
                        answer1 = False
                        continuar = False
                        break
                    else:
                        try:
                            iz = float(iz_input)
                            break
                        except ValueError:
                            print("Please enter a valid number.")
                            
                while answer1=='n':
                    try:
                        entrada = input("Do you want the data to calculate T90 from observations[o] \nor numerical models[m]? [o/m] or type 'exit' to exit: ")
                        if entrada == 'exit':
                            answer1 = False
                            continuar = False
                            break
                        elif entrada.lower() in ['o', 'm']:
                            answer2 = entrada.lower()
                            if answer2 == 'o':
                                estaciones_proximas_rad_solar(lat,lon)
                                estaciones_proximas_TS(lat,lon)
                                
                                # Ask for the value of z
                                while answer2 == 'o':
                                    z_input = input("Enter the value of depth z (m) or type 'exit' to exit: ")
                                    if z_input == 'exit':
                                        answer1 = False; answer2 = False
                                        continuar = False
                                        break
                                    else:
                                        try:
                                            z = float(z_input)
                                            break
                                        except ValueError:
                                            print("Please enter a valid number.")
                                            
                                # Ask for the value of T
                                while answer2 == 'o':
                                    T_input = input("Enter the value of temperature T (ºC) or type 'exit' to exit: ")
                                    if T_input == 'exit':
                                        answer1 = False; answer2 = False
                                        continuar = False
                                        break
                                    else:
                                        try:
                                            T = float(T_input)
                                            break
                                        except ValueError:
                                            print("Please enter a valid number.")
                                            
                                # Ask for the value of S
                                while answer2 == 'o':
                                    S_input = input("Enter the value of salinity S (psu) or type 'exit' to exit: ")
                                    if S_input == 'exit':
                                        answer1 = False; answer2 = False
                                        continuar = False
                                        break
                                    else:
                                        try:
                                            S = float(S_input)
                                            break
                                        except ValueError:
                                            print("Please enter a valid number.")
                                            
                                # Ask for the value of iz
                                while answer2 == 'o':
                                    iz_input = input("Enter the value of solar radiation iz (W/m2) or type 'exit' to exit: ")
                                    if iz_input == 'exit':
                                        answer1 = False; answer2 = False
                                        continuar = False
                                        break
                                    else:
                                        try:
                                            iz = float(iz_input)
                                            break
                                        except ValueError:
                                            print("Please enter a valid number.")
                                            
                            elif answer2 == 'm':
                                print(colored("\nMeteoGalicia database:\n", attrs=['underline']))
                                print("\nWe'll use MOHID model data from MeteoGalicia for T and S.")
                                print("You can find the data at: http://mandeo.meteogalicia.gal/thredds/catalog.html")
                                print("You'll have to look up the Solar Rad data.")
                                estaciones_proximas_rad_solar(lat,lon)
                                
                                # Ask for the value of iz
                                while answer2 == 'm':
                                    iz_input = input("Enter the value of solar radiation iz (W/m2) or type 'exit' to exit: ")
                                    if iz_input == 'exit':
                                        answer1 = False; answer2 = False
                                        continuar = False
                                        break
                                    else:
                                        try:
                                            iz = float(iz_input)
                                            break
                                        except ValueError:
                                            print("Please enter a valid number.")
                                            
                                # Ask for the value of z
                                while answer2 == 'm':
                                    z_input = input("Enter the value of depth z (m) or type 'exit' to exit " +
                                                    "(If you don't know it, enter 0 to calculate on the surface): ")
                                    if z_input == 'exit':
                                        answer1 = False; answer2 = False
                                        continuar = False
                                        break
                                    else:
                                        try:
                                            z = float(z_input)
                                            break
                                        except ValueError:
                                            print("Please enter a valid number.")
                                            
                                # Now, we download and perform automatic reading of the MOHID files.
                                aux = download_MOHID_files(fecha,key) 
                                if aux == 'Done':
                                    # We extract the values of T and S
                                    T,S = extract_TS_nc_file(fecha,lat,lon,z)
                                    # We check if T or S are NaN
                                    if m.isnan(T) or m.isnan(S):
                                        print("Error: one of the variables is 'nan'.")
                                        print("Try to enter values further offshore or shallower.")
                                        answer1 = False; answer2 = False
                                        continuar = False
                                    else:
                                        print("T and S values obtained correctly.")
                                else:
                                    break
                                 
                            break
                        else:
                            print("Please enter 'o' or 'm'.")
                    except ValueError:
                        print("Please enter a valid response.") 
                break
            else:
                print("Please enter 'y' or 'n'.")
        except ValueError:
            print("Please enter a valid response.") 
            

    # Print user input values
    if continuar:
        print(colored("\nUser input values:\n", attrs=['underline']))
        print("- Latitude:", lat)
        print("- Longitude:", lon)
        print("- Date:", fecha)
        if answer1=='y' or answer2 == 'o' or aux == 'Done':
            print("- T (ºC):", T)
            print("- S (psu):", S)
            print("- z (m):", z)
            print("- iz (W/m2):", iz)
            print(colored("\nCalculate T90:\n", attrs=['underline']))
            # Calculate T90 and print the result
            time_90 = T90(T,S,iz)
            print("- The value of T90 is:", time_90, 'h')
            
    return(lat,lon,fecha,time_90,key)

#-------------------LCS PART-----------------------------------
def download_LCS_files(fecha,key):
    """
    This function downloads MYCOASTLCS netCDF files for a given date range and location.

    Parameters
    ----------
    fecha : datetime
        The date to start downloading files (datetime object).
    key : list
        The location key, containing two strings.

    Returns
    -------
    str
        A message indicating the function has finished ('Done') or None if there was an error.
    """
    
    # Calculate start and end dates
    delta_time = timedelta(days=1)
    init_time = datetime(fecha.year, fecha.month, fecha.day) - timedelta(days=2)
    end_time = datetime(fecha.year, fecha.month, fecha.day)

    # Loop through each day in the date range
    date = init_time
    while (date <= end_time):
        print('03:  Downloading file...')
        print('-' * 30)
        
        # Extract date info into separate variables
        year = date.strftime("%Y")
        month = date.strftime("%m")
        day = date.strftime("%d")
        
        # Remove location key if it is 'NoiaMuros'
        if key[1] == 'NoiaMuros':
            key.pop(1)
            
        # Define URL and file name based on date and location key
        #http: // thredds - gfnl.usc.es / thredds / fileServer / MYCOASTLCS / MYCOASTLCS_Vigo_20230318.nc
        file_URL = 'http://thredds-gfnl.usc.es/thredds/fileServer/MYCOASTLCS/MYCOASTLCS_' + key[1] + '_%s.nc' \
                   % (year + month + day)
        file_name_out = "MYCOASTLCS_" + key[1] + "_%s.nc" % (year + month + day)

        print(file_URL, file_name_out)

        try:
            # Download the file
            with tqdm(unit='B', unit_scale=True, miniters=1, desc=file_name_out) as tqdm_instance:
                urllib.request.urlretrieve(file_URL, file_name_out, reporthook=lambda block_num, block_size, total_size:
                tqdm_instance.update(block_num * block_size - tqdm_instance.n))
        except urllib.error.HTTPError as e:
            print(f'Error downloading file: {e}')
            return None
        else:
            print(file_name_out + ' downloaded successfully.')
        finally:
            # Advance the date
            date = date + delta_time
            print('-' * 30)
            
    return 'Done'


def delete_nc_files():
    """
    This function deletes all netCDF files in the current working directory.

    Returns:
        None
    """
    
    # Get current working directory and find all netCDF files
    ruta = os.getcwd()+'\\'
    all_files = glob.glob(os.path.join(ruta +'*.nc'))

    # Loop through files and delete them if they exist
    for file_path in all_files:
        if os.path.exists(file_path):
            os.remove(file_path)
            print(f"File {file_path} deleted successfully.")
        else:
            print(f"File {file_path} does not exist.")
            

def delete_png_files(ruta):
    """
    This function deletes all PNG files in the specified directory.

    Parameters
    ----------
    ruta : str
        The path to the directory containing the PNG files.

    Returns
    -------
    None.

    """
    
    # Find all PNG files in directory
    all_files = glob.glob(os.path.join(ruta +'*.png'))

    # Loop through files and delete them if they exist
    for file_path in all_files:
        if os.path.exists(file_path):
            os.remove(file_path)
            print(f"File {file_path} deleted successfully.")
        else:
            print(f"File {file_path} does not exist.")
            
def delete_gif_files():
    """
    This function deletes all GIF files in the current working directory.

    Returns:
        None
    """
    
    # Get current working directory and find all GIF files
    ruta = os.getcwd()+'\\'
    all_files = glob.glob(os.path.join(ruta +'*.gif'))

    # Loop through files and delete them if they exist
    for file_path in all_files:
        if os.path.exists(file_path):
            os.remove(file_path)
            print(f"File {file_path} deleted successfully.")
        else:
            print(f"File {file_path} does not exist.")
            
            
def extract_time_LCS_file(time):
    """
    This function extracts the time index and file index of the nearest time in LCS netCDF files.

    Parameters
    ----------
    time : str
        A string representing the time to search for, in format YYYY-MM-DDTHH:MM:SS.

    Returns
    -------
    list or None
        A list containing the time index and file index if the time is found, or None if not found.

    """
    
    # Get current working directory and find all netCDF files
    ruta = os.getcwd()+'\\'
    all_files = glob.glob(os.path.join(ruta +'*.nc'))
    human_sort(all_files)

    # Open the first two files and get their time values
    ds1 = xr.open_dataset(all_files[0])
    ds2 = xr.open_dataset(all_files[1])
    datas1 = pd.to_datetime(ds1.time.values)
    datas2 = pd.to_datetime(ds2.time.values)

    # Check if the time is in either of the files
    if str(time) in datas1.astype(str):
        time_index = np.where(datas1 == time)[0][0]
        file_index = 0
        return [time_index,file_index]
    elif str(time) in datas2.astype(str):
        time_index = np.where(datas2 == time)[0][0]
        file_index = 1
        return [time_index,file_index]
    else:
        print("The date is not in the files.")
        return None
    

def read_inputs(input_file):
    """Read keywords for options"""
    input_keys = ['path_in',
                  'file_in',
                  'path_out'
                  'file_out',
                  'nx',
                  'ny',
                  'resolution',
                  'scale',
                  'n_time',
                  'n_level',
                  'title',
                  'style',
                  'limits',
                  'vector',
                  'scalar',
                  'wms_url',
                  'wms_layers']
    return read_input(input_file, input_keys)

def crear_gift(ruta_maps, ruta_gift, vel_maps, loop):
    """
    Creates an animated GIF file from a directory of PNG images.

    Parameters
    ----------
    ruta_maps : str
        Path to the directory containing the PNG images.
    ruta_gift : str
        Path to the output GIF file.
    vel_maps : float
        Speed of the animation in frames per second.
    loop : int
        Number of loops for the animation (0 for infinite loop).

    Returns
    -------
    None.

    """
    
    # Set the path to the directory of images
    image_directory = ruta_maps
    
    # Get the list of file names in the directory
    image_files = os.listdir(image_directory)
    
    # Sort the list of file names alphabetically
    image_files = sorted(image_files)
    
    # Create a list of Image objects from Pillow
    image_list = []
    for filename in image_files:
        if filename.endswith(".png"):
            image_path = os.path.join(image_directory, filename)
            image = Image.open(image_path)
            image_list.append(image)
    
    # Save the list of images as an animated GIF file
    gif_path = ruta_gift
    frame_duration = 1e3 / vel_maps # Duration of each frame in milliseconds
    image_list[0].save(gif_path, save_all=True, append_images=image_list[1:],
                       duration=frame_duration, loop=loop)


def mostrar_gift(ruta_gift):
    """
    Displays an animated GIF file in a new window.

    Parameters
    ----------
    ruta_gift : str
        The path to the GIF file to be displayed.

    Returns
    -------
    None.

    """
    
    # Path of the GIF file to be displayed
    gif_path = ruta_gift
    
    # Open the GIF file using Pillow
    gif_im = Image.open(gif_path)
    
    # Create a new window using Tkinter
    root = tk.Tk()
    canvas = tk.Canvas(root, width=gif_im.width, height=gif_im.height)
    canvas.pack()
    
    # Convert each frame of the animation into a PhotoImage object
    frames = []
    for frame in range(0, gif_im.n_frames):
        gif_im.seek(frame)
        frame_im = ImageTk.PhotoImage(gif_im)
        frames.append(frame_im)
    
    
    def update_frame(frame_number):
        """
        Update the canvas with each frame of the animation.

        Parameters
        ----------
        frame_number : int
            The current frame number.

        Returns
        -------
        None.

        """
        # Delete all items on the canvas
        canvas.delete("all")
        # Add the current frame to the canvas
        print(frame_number, frames[frame_number])
        canvas.create_image(0, 0, image=frames[frame_number], anchor="nw")
        # Call the update_frame function again after a delay of 500 milliseconds
        root.after(500, update_frame, (frame_number + 1) % len(frames))
    
    # Start the animation
    update_frame(0)
    
    # Run the tkinter loop
    root.mainloop()
    
def main_drawmap(inputs,start,end,ref_lon,ref_lat,key,draw_scale):
    """
    Creates maps based on user inputs.

    Parameters
    ----------
    inputs : dict
        A dictionary containing user-defined inputs.
    start : int
        The start date index of the data to be plotted.
    end : int
        The end date index of the data to be plotted.
    ref_lon : float
        A reference longitude for the plot.
    ref_lat : float
        A reference latitude for the plot.
    key : list
        The location key, containing two strings.
    draw_scale : str
        Parameter to change the color scale of the map. If you enter 'jet' use the jet scale. 
        Otherwise it will use the viridis scale by default.

    Returns
    -------
    None.

    """
    
    # Start
    print(colored("\nCreating maps:\n", attrs=['underline']))
    
    if inputs['n_time'] == 'all':
        # Draw map for start-end period
        draw_map_24(inputs,start,end,ref_lon, ref_lat,key,draw_scale)
    else:
        if inputs['vector']:
            # Draw vector map
            draw_map_vector(inputs, inputs['n_time'])
        elif inputs['scalar']:
            # Draw scalar map
            draw_map_scalar(inputs, inputs['n_time'])
        else:
            # Print error message and quit
            print('You must choose between scalar and vector maps in the drawmap.json file.\nWords: vector, scalar')
            quit(1)
            
def delete_all():
    """
    Deletes all generated files.

    Returns
    -------
    None.

    """
    
    # Deletes all NetCDF files
    delete_nc_files()
    
    # Deletes all NetCDF4 files
    delete_nc4_files()
    
    # Deletes all PNG files in the maps directory
    delete_png_files(os.getcwd()+'\\maps\\')
    
    # Deletes all PNG files in the current directory
    delete_png_files(os.getcwd()+'\\')
    
    # Deletes all GIF files
    delete_gif_files()
    
    return None

def jet_map(start_index,file_index,T9,nam,key,ref_lon, ref_lat):
    """
    Creates maps based on user inputs.

    Parameters
    ----------
    start_index : int
        The start date index of the data to be plotted.
    file_index : int
        The index of the start file of the data to be plotted.
    T9 : int
        Time duration to be plotted.
    nam : str
        Name of the variable to represent
    key : list
        The location key, containing two strings.
    ref_lon : float
        A reference longitude for the plot.
    ref_lat : float
        A reference latitude for the plot.

    Returns
    -------
    None.

    """
    
    # Get all netCDF files in the current directory
    all_files = glob.glob(os.path.join('*.nc'))

    # Set the path for the drawmap.json file
    json_path = os.getcwd() + '\\common\\'

    # Read the input values from the drawmap.json file
    inputs = read_inputs(json_path + 'drawmap.json')
    inputs['path_in'] = os.getcwd()
    inputs['path_out'] = os.getcwd()
    inputs['file_out'] = "MYCOAST_" + nam + "_" + key[1] + "_during_T90"
    inputs['scalar_magnitude'] = nam
    inputs['title'] = "MYCOAST-" + nam + " Ria " + key[1] + ":"

    # Open the netCDF file for the specified time range
    time_nc = xr.open_dataset(all_files[file_index]).isel(time=slice(start_index, (start_index+1)))

    # Open the netCDF file for the entire day
    ds_aux = xr.open_dataset(all_files[file_index])

    # Calculate the sum of the variable values over the specified time range
    if (start_index+T9)>=24:
            end = 24
            aux = (ds_aux[nam]).isel(time=slice(start_index, end)).sum(dim='time', skipna=False).values
            T9 = T9-(end-start_index)
            if T9 == 0:
                time_nc[nam][0] = aux
                time_nc.to_netcdf('condensed.nc')
            elif  0<T9<=24:
                ds_aux2 = xr.open_dataset(all_files[file_index+1])
                aux2 = (ds_aux2[nam]).isel(time=slice(0, T9)).sum(dim='time', skipna=False).values
                time_nc[nam][0] = aux+aux2
                time_nc.to_netcdf('condensed.nc')
            elif T9>24:
                start_index = 0; end = 24
                ds_aux2 = xr.open_dataset(all_files[file_index+1])
                aux2 = (ds_aux2[nam]).isel(time=slice(start_index, end)).sum(dim='time', skipna=False).values
                T9 = T9-(end-start_index)
                ds_aux3 = xr.open_dataset(all_files[file_index+2])
                aux3 = (ds_aux3[nam]).isel(time=slice(0, T9)).sum(dim='time', skipna=False).values
                time_nc[nam][0] = aux+aux2+aux3
                time_nc.to_netcdf('condensed.nc')
    else:
            end = start_index+T9
            aux = (ds_aux[nam]).isel(time=slice(start_index, end)).sum(dim='time', skipna=False).values
            time_nc[nam][0] = aux
            time_nc.to_netcdf('condensed.nc')  
    
    # Set the input file name for main_drawmap function
    inputs['file_in'] = 'condensed.nc'
    # Set the scale for the map to 'jet'
    draw_scale = 'jet'

    # Create the jet map
    main_drawmap(inputs,0,1,ref_lon, ref_lat,key,draw_scale)

    # Directory where the images are located
    directory = os.getcwd()
    # Get the list of files in the directory
    files = os.listdir(directory)
    # Filter to get the name of the file that meets the conditions
    file_name = next((file for file in files if file.startswith(inputs['file_out'])), None)
    if file_name is not None:
        # Full path of the file
        file_path = os.path.join(directory, file_name)
        # Open the image
        image = Image.open(file_path)
        # Show the image
        image.show()
    else:
        print("No image found with the specified name.")


def start_MYCOAST(fecha,T90,ref_lon, ref_lat,key):
    """
    Downloads MYCOAST data for a given date and key, and creates an animation
    of the data for a given time window using the `main_drawmap` function.

    Parameters
    ----------
    fecha : str
        A string representing the date in the format YYYYMMDD.
    T90 : float
        A float representing the time window (in hours) for which to create the animation.
    ref_lon : float
        A reference longitude for the plot.
    ref_lat : float
        A reference latitude for the plot.
    key : list
        The location key, containing two strings.

    Returns
    -------
    None.

    """
    
    print(colored("\nMYCOAST_LCS Analysis:\n", attrs=['underline']))

    # Download the MYCOAST data for the specified date and variable
    aux = download_LCS_files(fecha,key)
    if aux == 'Done':
        # Get all netCDF files in the current directory
        all_files = glob.glob(os.path.join('*.nc'))

        # Extract the start index and file index for the specified time range
        T9 = int(np.round(T90,0))
        index = extract_time_LCS_file(fecha)
        start = index[0]
        file_index = index[1]

        # Get the list of variables in the netCDF file
        dt = xr.open_dataset(all_files[0])
        long_names = []
        for i in dt.data_vars:
            long_names.append(i)
        
        nam = None
        while nam not in long_names:
            # Prompt the user to select the variable to plot
            print('The available variables are:')
            for ln in long_names:
                print(f'- {ln}')
            nam = input('Enter the long_name of the variable to plot or type "exit" to quit: ')
            if nam == 'exit':
                break
            if nam not in long_names:
                print(f'Error: "{nam}" is not a valid long_name. Please try again.')
        
        # Create a plot for the specified time range
        jet_map(start,file_index,T9,nam,key,ref_lon, ref_lat)

        # Set the input parameters for the plot animation
        json_path = os.getcwd() + '\\common\\'
        inputs = read_inputs(json_path + 'drawmap.json')
        inputs['path_in'] = os.getcwd()
        inputs['path_out'] = os.getcwd()+'\\maps'
        inputs['file_out'] = "MYCOAST_" + nam + "_" + key[1]
        inputs['scalar_magnitude'] = nam
        inputs['title'] = "MYCOAST-" + nam + " Ria " + key[1] + ":"
        
        # Create a plot for each time window within the specified time range
        for i in all_files:
            inputs['file_in'] = i
          
            if (start+T9)>=24:
                end = 24
            else:
                end = start+T9
            
            # Set the scale for the map to 'viridis' (default)
            draw_scale = 'viridis'
            main_drawmap(inputs,start,end,ref_lon, ref_lat,key,draw_scale)
            
            T9 = T9-(end-start)
            if T9 <= 0:
                break
            else:
                start = 0
            
        # Create a GIF animation of the plots
        # Path where the maps are stored
        ruta_maps = os.getcwd()+'\\maps'
        # Path to store the generated GIF animation
        ruta_gift = os.getcwd() + '\\animacion.gif'
        # Create the gift with a duration of 0.5 seconds each frame and in an infinite loop
        crear_gift(ruta_maps, ruta_gift, 2, 0)
        # Display the generated GIF animation
        mostrar_gift(ruta_gift)
             
if __name__ == '__main__':
    # Delete any existing files from previous executions
    delete_all()
    # Start the T90 calculation and user input process
    lat,lon,fecha,time_90,key = start_T90()
    # Start the MYCOAST analysis with the user inputs from the previous step
    start_MYCOAST(fecha,time_90,lon,lat,key)
