#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
**reader.py**

* *Purpose:* Abstract class of reader object. Reader object deals with reading HDF, NetCDF son classes

* *python version:* 3.9
* *author:* Pedro Montero
* *license:* INTECMAR
* *requires:*
* *date:* 2021/10/13
* *version:* 1.0.0
* *date version* 2021/10/13


"""

from abc import ABC, abstractmethod


class Reader(ABC):
    """ Abstract class of Reader object"""
    def __init__(self, file):
        """ Open a file and get longitudes and latitudes and coordinates rank"""
        self.file = file
        self.dataset = self.open()
        self.n_longitudes = None
        self.n_latitudes = None
        self.longitudes = self.get_longitudes()
        self.latitudes = self.get_latitudes()
        self.coordinates_rank = self.get_rank(self.longitudes)
        self.ini_ntime = self.get_ini_ntime()
        self.close()

    @abstractmethod
    def open(self):
        pass

    @abstractmethod
    def close(self):
        pass

    @abstractmethod
    def get_latitudes(self):
        pass

    @abstractmethod
    def get_longitudes(self):
        pass

    @abstractmethod
    def get_dates(self):
        pass

    @abstractmethod
    def get_date(self, n_time):
        pass

    @abstractmethod
    def get_ini_ntime(self):
        pass

    @staticmethod
    def get_rank(array) -> int:
        """ Number of dimensions of an array"""
        return len(array.shape)




