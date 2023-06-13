
from readers.reader_factory import read_factory


def main():

    factory = read_factory('../datos/drawmap/MOHID_Hydrodynamic_Vigo_20180711_0000.hdf5')
    reader = factory.get_reader()
    print(reader.longitudes)
    print(reader.latitudes)
    print(reader.get_date(1))
    print(reader.get_date(10))

    factory = read_factory('../datos/drawmap/wrf_arw_det1km_history_d05_20180711_0000.nc4')
    reader = factory.get_reader()

    print(reader.latitudes)
    print(reader.longitudes)
    print(reader.get_date(1))
    print(reader.get_date(10))
    print(f'number of longitudes = {reader.n_longitudes} and latitudes = {reader.n_latitudes}')
    print(reader.get_variable('wind_module', 1))
    reader.close()


if __name__ == "__main__":
    main()
