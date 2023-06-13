
from datetime import datetime
import numpy as np
import h5py
from .reader import Reader


class ReaderHDF(Reader):

    names = {
        'northward_velocity': '/Results/velocity V/velocity V_',
        'eastward_velocity': '/Results/velocity U/velocity U_'
    }

    def open(self):
        self.dataset = h5py.File(self.file)
        return self.dataset

    def close(self):
        self.dataset.close()

    def get_latitudes(self):
        lat_in = self.dataset['/Grid/Latitude']
        if len(lat_in.shape) == 1:
            self.n_latitudes = lat_in.shape[0]
            return lat_in
        elif len(lat_in.shape) == 2:
            self.n_latitudes = lat_in.shape[1]
            return lat_in[0, ]

    def get_longitudes(self):
        lon_in = self.dataset['/Grid/Longitude']
        if len(lon_in.shape) == 1:
            self.n_longitudes = lon_in.shape[0]
            return lon_in
        elif len(lon_in.shape) == 2:
            self.n_longitudes = lon_in.shape[0]
            return lon_in[:, 1]

    def get_dates(self):
        return self.dataset['/Time']

    def get_date(self, n_time):

        date_in = self.dataset['/Time/Time_' + str(n_time).zfill(5)]
        return datetime(year=int(date_in[0]), month=int(date_in[1]), day=int(date_in[2]),
                        hour=int(date_in[3]), minute=int(date_in[4]), second=int(date_in[5]))

    def get_variable(self, name_var, n_time):
        path = self.names[name_var]
        variable = self.dataset[path + str(n_time).zfill(5)]
        if len(variable.shape) == 2:
            variable = np.transpose(variable)
        elif len(variable.shape) == 3:
            variable = np.transpose(variable, (0, 2, 1))
        return variable

    def get_ini_ntime(self):
        return 1




