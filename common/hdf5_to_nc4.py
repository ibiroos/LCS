# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 12:31:18 2023

@author: josmo
"""


import h5py
import xarray as xr
import numpy as np
import urllib
import pandas as pd
# Para las fechas
from datetime import datetime,timedelta
import os
from tqdm import tqdm
import glob

# Rutas no hdf5

# Abre o ficheiro lagranxiano e le


def MOHID_hdf5_to_nc(file_in):
   
    rootGroup = 'Results/'
    nameTimesGroup = '/Time/'
    nameTimesNameGroup = '/Time/Time_'
   
    f = h5py.File(file_in, 'r')
   
    # Le as latitudes e lonxitudes e busca os nodos do cadrado
   
    latIn = f['/Grid/Latitude']
    lonIn = f['/Grid/Longitude']
   
    rank = len(latIn.shape)
   
    if rank == 2:
        lat = latIn[0,]
        lon = lonIn[:, 0]
   
    if rank == 1:
        lat = latIn
        lon = lonIn
   
    # le os tempos pois o bucle vai ser por aqui
    timesGroup = f[nameTimesGroup]
    tempo = []
    u_t = []
    v_t = []
   
    for nameTime in timesGroup:
       
        numName = nameTime.split('_')[1]
        valNumName = int(numName)
        rootNameTime = nameTimesGroup + nameTime
        time = f[rootNameTime]
       
        dataIn = datetime(int(time[0]), int(time[1]), int(time[2]), int(time[3]), int(time[4]), int(time[5]))
        # nameData = '%04d/%02d/%02d %02d:%02d:%02d' % (int(time[0]),int(time[1]),int(time[2]),int(time[3]),int(time[4]),int(time[5]))
       
        tempo.append(dataIn)
        nameU = rootGroup + "velocity U/velocity U_" + numName
        nameV = rootGroup + "velocity V/velocity V_" + numName
       
        u = f[nameU]
        v = f[nameV]
        u_t.append(u[-1,:,:])
        v_t.append(v[-1,:,:])
   
   
    lon_c = (lon[1:]+lon[0:-1])/2.
    lat_c = (lat[1:]+lat[0:-1])/2.
    U = np.stack(u_t)
    V = np.stack(v_t)
    coords = {'lon':(['lon'],lon_c),'lat':(['lat'],lat_c),'time':(['time'],tempo)}    
    ds = xr.Dataset({'u':(['time','lon','lat'],U),'v':(['time','lon','lat'],V)},coords=coords)
   
    ds['u']=ds.u.where(ds.u != 0,-9.8999995E15)
    ds['v']=ds.v.where(ds.v != 0,-9.8999995E15)
    
    lon_attributtes = {'long_name': 'longitude',
    'standard_name': 'longitude',\
    'units': 'degrees_east',\
    'valid_min': -180.0,\
    'valid_max': 180.0}
    
    lat_attributtes = {'long_name': 'latitude',
    'standard_name': 'latitude',\
    'units': 'degrees_north',\
    'valid_min': -90.0,\
    'valid_max': 90.0}
    time_atts = {'long_name':'time'}
    ds.time.encoding['units'] = 'seconds since 1950-01-01 00:00:00'


#    depth_attributtes = {'long_name': 'depth',
#    'standard_name': 'depth',\
#    'units': 'meters',\
#    '_FillValue': -9.8999995E15,\
#    'valid_min': np.min(self.grid['depth']),\
#    'valid_max': np.max(self.grid['depth'])}
   
    ds.lon.attrs = lon_attributtes
    ds.lat.attrs = lat_attributtes
    #ds.depth.attrs = depth_attributtes
    ds.time.attrs = time_atts
   
    encoding = {'u': {'_FillValue': -9.8999995E15},'v': {'_FillValue': -9.8999995E15}}
    ds.to_netcdf(file_in[0:-5]+'.nc4',encoding=encoding)
    ds.close()
    return


def download_HDF5_files(fecha):
    # Fechas de inicio e fin
    delta_time=timedelta(days=1)
    data = fecha-delta_time
    init_time=datetime(data.year,data.month,data.day)
    end_time=datetime(data.year,data.month,data.day)
    # Bucle que recorre las fechas.
    date=init_time
    print('-'*30)
    while (date <= end_time):
        print('Descargando archivo...')
        print('-'*30) 
        # Extraemos la info del objeto fecha en diferentes variables caracter
        year=date.strftime("%Y")
        month=date.strftime("%m")
        day=date.strftime("%d")
        # Definimos los nombres de los ficheros de entrada y salida
        file_URL='https://mandeo.meteogalicia.es/thredds/fileServer/modelos/mohid/history/arousa/MOHID_Hydrodynamic_Arousa_%s_0000.hdf5' \
        %(year+month+day)
        file_name_out="MOHID_Arousa_%s_0000.hdf5" %(year+month+day)
        try:
            # Descargamos el fichero
            with tqdm(unit='B', unit_scale=True, miniters=1, desc=file_name_out) as tqdm_instance:
                urllib.request.urlretrieve(file_URL, file_name_out, reporthook=lambda block_num, block_size, total_size: 
                tqdm_instance.update(block_num * block_size - tqdm_instance.n))
        except urllib.error.HTTPError as e:
            print(f'Error al descargar el archivo: {e}')
            return None
        else:
            print(file_name_out)
        finally:
            # Avanzar la fecha
            date=date+delta_time
            print('-'*30)
    return 'Done'

if __name__ == '__main__':
    fecha = datetime(2021,12,7,7)
    download_HDF5_files(fecha)
    ruta = os.getcwd()+'\\'
    all_files = glob.glob(os.path.join(ruta +'*.hdf5'))
    MOHID_hdf5_to_nc(all_files[0])
    
    
    
    
    
    