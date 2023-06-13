import sys
from abc import ABC, abstractmethod
from .reader import Reader
from .reader_HDF import ReaderHDF
from .reader_NetCDF import ReaderNetCDF


class ReaderFactory(ABC):
    """Basic class of factory readers"""

    def __init__(self, file_in):
        self.file_in = file_in

    @abstractmethod
    def get_reader(self) -> Reader:
        """return a reader class"""


class ReaderHDFFactory(ReaderFactory):
    """Factory for Reader of HDF files"""

    def get_reader(self) -> Reader:
        return ReaderHDF(self.file_in)


class ReaderNetCDFFactory(ReaderFactory):
    """ Factoyr for Reader of NetCDF files"""

    def get_reader(self) -> Reader:
        return ReaderNetCDF(self.file_in)


def read_factory(file_in) -> ReaderFactory:
    """ Construct a reader factory based on the extension of file"""

    factories = {
        "nc": ReaderNetCDFFactory,
        "nc4": ReaderNetCDFFactory,
        "hdf": ReaderHDFFactory,
        "hdf5": ReaderHDFFactory
    }

    extension = file_in.split('.')[-1]
    if extension in factories:
        return factories[extension](file_in)
    print(f'Unknown file extension {extension}')
    sys.exit(1)





