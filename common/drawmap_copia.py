#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Para el tratamiento de datos
import numpy as np
import pandas as pd
# Para las representaciones
import matplotlib.pyplot as plt
from matplotlib.cm import jet
# Para los directorios
import os
# Importamos funciones de otros programas
from common.readers.reader_factory import read_factory
from common import read_input
from common.boundarybox import BoundaryBox
from common.drawcurrents import drawcurrents
from common.map import Map

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


def draw_map_vector(inputs, n):
    """draw 1 maps of a day"""
    draw_map = Vector(inputs)
    draw_map.read_head()
    draw_map.create_title(n)
    draw_map.reader_by_time()
    draw_map.draw()


def draw_map_scalar(inputs, n):
    draw_map = Scalar(inputs)
    draw_map.read_head()
    draw_map.create_title(n)

    if n is None:
        draw_map.reader_no_time()
    else:
        draw_map.reader_by_time()
    draw_map.draw()


def draw_map_24(inputs,start,end,ref_lon, ref_lat):
    """draw 24+1 maps of a day"""

    if inputs['scalar']:
        draw_map = Scalar(inputs)
    elif inputs['vector']:
        draw_map = Vector(inputs)
    
    draw_map.read_head()
    draw_map.options.time = 0 # creo que esto sobra

    for n in range(draw_map.reader.ini_ntime+start, draw_map.reader.ini_ntime+end):
        draw_map.create_title(n)
        draw_map.options.time = n
        draw_map.reader_by_time()
        draw_map.draw(ref_lon, ref_lat)


class OptionsMap:

    def __init__(self, inputs):

        self.file_path_in = inputs['path_in']
        self.file_path_out = inputs['path_out']
        self.file_in = inputs['file_in']
        self.file_name = os.path.join(self.file_path_in, self.file_in)
        self.file_hdf_out = inputs['file_out']
        self.file_out = os.path.join(self.file_path_out, self.file_hdf_out)

        self.nx = inputs['nx']
        self.ny = inputs['ny']
        self.scale = inputs['scale']
        self.resolution = inputs['resolution']
        self.style = inputs['style']
        self.title = inputs['title']

        self.level = inputs['n_level']
        self.time = inputs['n_time']
        limits = inputs['limits']
        self.boundary_box = BoundaryBox(limits[0], limits[1], limits[2], limits[3])

        self.vector_bool = inputs['vector']
        self.scalar_bool = inputs['scalar']

        self.wms_url = inputs['wms_url']
        self.wms_layers = inputs['wms_layers']


class DrawMap:

    def __init__(self, inputs):
        self.options = OptionsMap(inputs)
        self.reader = self.get_reader(self.options.file_name)

    def read_head(self):

        with self.reader.open():

            lat = self.reader.latitudes
            lon = self.reader.longitudes

            if self.reader.coordinates_rank == 1:
                self.lats = lat[0:self.reader.n_latitudes - 2]  # ATENCIóN:Esto y lo de abajo era -1, revisar
                self.lons = lon[0:self.reader.n_longitudes - 2]
            elif self.reader.coordinates_rank == 2:
                self.lats = lat[0:self.reader.n_longitudes - 2, 0:self.reader.n_latitudes - 2]
                self.lons = lon[0:self.reader.n_longitudes - 2, 0:self.reader.n_latitudes - 2]

    def create_title(self, n_time):
        if n_time is None:
            n_time = 0
        with self.reader.open():
            data = self.reader.get_date(n_time)
            data_str = data.strftime("%Y-%m-%d %H:%M UTC")
            data_comp = data.strftime("%Y%m%d%H%M")
            self.title_full = self.options.title + " " + data_str
            self.file_out_full = self.options.file_out + '_' + data_comp + '.png'

    def get_reader(self, file_in):
        print('Opening: {0}'.format(file_in))
        factory = read_factory(file_in)
        return factory.get_reader()


class Vector(DrawMap):

    def __init__(self, inputs):

        super().__init__(inputs)

        self.u_name = inputs['u']
        self.v_name = inputs['v']
        self.us = None
        self.vs = None
        self.modules = None

    def reader_by_time(self):

        with self.reader.open():
            u = self.reader.get_variable(self.u_name, self.options.time)
            v = self.reader.get_variable(self.v_name, self.options.time)

            if len(u.shape) == 3:
                self.us = u[self.options.level, :-1, :- 1]
                self.vs = v[self.options.level, :-1, :-1]

            elif len(u.shape) == 2:
                self.us = u[:-1, :-1]
                self.vs = v[:-1, :-1]

            self.modules = pow((pow(self.us, 2) + pow(self.vs, 2)), .5)

    def draw(self):
        drawcurrents(self.reader.coordinates_rank, self.options.nx, self.options.ny,
                     self.options.scale, self.options.resolution,
                     self.options.level, self.options.time,
                     self.lats, self.lons, self.us, self.vs, self.modules,
                     self.file_out_full, self.title_full, self.options.style,
                     self.options.boundary_box)


class Scalar(DrawMap):

    def __init__(self, inputs):

        super().__init__(inputs)

        self.scalar_name = inputs['scalar_magnitude']
        self.scalars = None

    def reader_by_time(self):

        with self.reader.open():
            scalar = self.reader.get_variable(self.scalar_name, self.options.time)

            if len(scalar.shape) == 3:
                self.scalars = scalar[self.options.level, :-1, :- 1]

            elif len(scalar.shape) == 2:
                self.scalars = scalar[:-1, :-1]

    def reader_no_time(self):
        with self.reader.open():
            scalar = self.reader.get_var(self.scalar_name)

            if len(scalar.shape) == 3:
                self.scalars = scalar[self.options.level, :-1, :- 1]

            elif len(scalar.shape) == 2:
                self.scalars = scalar[:-1, :-1]

    def draw(self,ref_lon, ref_lat):
        """
        :return:
        """
        # a Map object is created to deal with
        mapa = Map()
        mapa.set_extension(self.options.boundary_box.lat_min, self.options.boundary_box.lat_max,
                           self.options.boundary_box.lon_min, self.options.boundary_box.lon_max)
        fig = mapa.add_figure()
        mapa.add_wms(self.options.wms_url, self.options.wms_layers)
        #mapa.add_wms("https://wms.mapama.gob.es/sig/Agua/PVERT/2021/wms.aspx?", 
                       #["Puntos de vertido de depuradoras urbanas (Q2021. Dir 91/271/CEE)"])
        

        #  Draw the scalars
        if self.reader.coordinates_rank == 1:
            longitudes, latitudes = np.meshgrid(self.lons, self.lats)
        else:
            longitudes, latitudes = self.lons, self.lats

        im = mapa.ax.pcolormesh(longitudes, latitudes, self.scalars[:-1, :-1], shading='auto')  # since lats and lons are n-1xn-1 matrix
        #im.set_cmap(jet)
        # add title
        mapa.ax.set_title(self.title_full)
        # add point
        path_discharge = os.getcwd() +'\\common\\Rias_UWWTP.xlsx'

        key = ['arousa', 'Arousa']
        discharge_lats = (pd.read_excel(path_discharge, sheet_name = key[1])['Lat'].values)
        discharge_lons = (pd.read_excel(path_discharge, sheet_name = key[1])['Lon'].values)
        mapa.ax.scatter(discharge_lons, discharge_lats, marker='.', s=180, color='darkorange',
                        edgecolors = "black", linewidths = 1.)
        mapa.ax.scatter(ref_lon, ref_lat, marker='*', s=100, color='red')
        fig.colorbar(im, ax=mapa.ax)
        #
        fig.savefig(self.file_out_full, dpi=100, facecolor='w', edgecolor='w', format='png',
                    transparent=False, bbox_inches='tight', pad_inches=0.1)
        plt.clf()
        plt.close(fig)
        return
