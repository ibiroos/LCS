import cartopy.crs as ccrs
import matplotlib.pyplot as plt


class Map:

    def __init__(self, projection=ccrs.PlateCarree()) -> None:
        self.wms = None
        self.layers = None
        self.projection = projection
        self.extension = None
        self.ax = plt.axes(projection=self.projection)

    def add_figure(self, figure_size=(14, 8)):
        fig = plt.figure(figsize=figure_size)
        self.ax = fig.add_subplot(1, 1, 1, projection=self.projection)
        return fig

    def set_extension(self, lat_south: float, lat_north: float, lon_west: float, lon_east: float) -> None:
        self.extension = [lon_west, lon_east, lat_south, lat_north]
        self.ax.set_extent(self.extension)
    
    def add_wms(self, wms: str, layers) -> None:
        self.wms = wms
        self.layers = layers
        if type(layers) is list:
            self.ax.add_wms(wms, layers)
        elif type(layers) is str:
            self.ax.add_wmts(wms, layers)


    @staticmethod
    def show():
        plt.show()
    

def main():

    arousa = Map()   
    arousa.set_extension(41.7, 43.4, -9.5, -7.7)
   
    #arousa.add_wms('http://www.ign.es/wms-inspire/pnoa-ma?',
    #['fondo', 'OI.OrthoimageCoverage'])

    #arousa.add_wms('https://ideihm.covam.es/ihm-inspire/wms-regionesmarinas?', ['SR.Shoreline'])



    #arousa.add_wms('https://wms.mapama.gob.es/sig/Agua/RiosCompPfafs/wms.aspx?', ['HY.PhysicalWaters.Waterbodies'])

    arousa.add_wms('https://services.arcgisonline.com/arcgis/rest/services/Canvas/World_Light_Gray_Base/MapServer/WMTS?'
                   , 'Canvas_World_Light_Gray_Base')
    #arousa.add_wms('https://ows.emodnet-bathymetry.eu/wms?', ['emodnet:mean_atlas_land', 'coastlines'])
    arousa.add_wms('https://ows.emodnet-bathymetry.eu/wms?', ['coastlines'])
    #arousa.add_wms('https://map1c.vis.earthdata.nasa.gov/wmts-geo/wmts.cgi', 'VIIRS_CityLights_2012')
    #arousa.add_wms('https://services.arcgisonline.com/arcgis/rest/services/Canvas/World_Dark_Gray_Base/MapServer/WMTS?',
    #'Canvas_World_Dark_Gray_Base')
    #arousa.add_wms("https://wms.mapama.gob.es/sig/EvaluacionAmbiental/Residuos/ParcelasAgric/wms.aspx?",
                   #["Parcelas agrícolas de aplicación de lodos de EDAR"])
    
    arousa.add_wms("https://wms.mapama.gob.es/sig/Agua/EDAR/2021/wms.aspx?",
                   ["Depuradoras de aguas residuales (Q2021. Dir 91/271/CEE)"])
    arousa.add_wms("https://wms.mapama.gob.es/sig/Agua/PVERT/2021/wms.aspx?", 
                   ["Puntos de vertido de depuradoras urbanas (Q2021. Dir 91/271/CEE)"])
    arousa.add_wms("https://wms.mapama.gob.es/sig/Agua/CNV/wms.aspx?", 
                   ["Censo Nacional de Vertidos"])
    arousa.add_wms("https://wms.mapama.gob.es/sig/Agua/RedHidrograficaMDT/wms.aspx?", 
                   ["Red hidrográfica básica MDT 100x100"])
    
    
    #manager = plt.get_current_fig_manager()
    #manager.window.showMaximized()
    arousa.show()



if __name__ == '__main__':
    main()

