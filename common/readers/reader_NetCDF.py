
from datetime import datetime
import netCDF4
from .reader import Reader


class ReaderNetCDF(Reader):

    def open(self):
        dataset = netCDF4.Dataset(self.file)
        self.variables = dataset.variables
        return dataset

    def close(self):
        self.dataset.close()

    def get_latitudes(self):
        lat_in = self.get_var('latitude')
        if len(lat_in.shape) == 1:
            self.n_latitudes = lat_in.shape[0]
        elif len(lat_in.shape) == 2:
            self.n_latitudes = lat_in.shape[1]
        return lat_in

    def get_longitudes(self):
        lon_in = self.get_var('longitude')
        if len(lon_in.shape) == 1:
            self.n_longitudes = lon_in.shape[0]
        elif len(lon_in.shape) == 2:
            self.n_longitudes = lon_in.shape[0]
        return lon_in

    def get_dates(self):
        times_in = self.get_var('time')
        return netCDF4.num2date(times_in[:], units=times_in.units)

    def get_date(self, n_time):
        return self.get_dates()[n_time]

    def get_variable(self, var_name, n_time):
        return self.get_var(var_name)[n_time, ]

    def get_var(self, var_name):
        """Return values using the CF standard name, long_name or the simple name of a variable in a netCDF file."""
        for var in self.variables:
            for atributo in (self.variables[var].ncattrs()):
                if atributo == 'standard_name':
                    nome_atributo = (getattr(self.variables[var], 'standard_name'))
                    if nome_atributo == var_name:
                        return self.variables[var]
                elif atributo == 'long_name':
                    nome_atributo = (getattr(self.variables[var], 'long_name'))
                    if nome_atributo == var_name:
                        return self.variables[var]
                elif var == var_name:
                    return self.variables[var]

    def get_ini_ntime(self):
        return 0
