class BoundaryBox:
    lat_min: float
    lat_max: float
    lon_min: float
    lon_max: float

    def __init__(self, lat_min=-90., lat_max=90., lon_min=-180., lon_max=180.):
        self.lat_min = lat_min
        self.lat_max = lat_max
        self.lon_min = lon_min
        self.lon_max = lon_max

    def middle_lat(self):
        x = self.lat_min + 0.5 * (self.lat_max - self.lat_min)
        return x

    def middle_lon(self):
        return self.lon_min + 0.5 * (self.lon_max - self.lon_min)


