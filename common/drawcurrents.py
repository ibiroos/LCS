

def drawcurrents(coordinates_rank, nx, ny, scale, resolution,
                 level, time, lat, lon, ust, vst, mod,
                 file_out, title, style, boundary_box):
    """

    :return:
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

    dx = 0.01
    middle_lon = boundary_box.middle_lon()
    middle_lat = boundary_box.middle_lat()
    m = Basemap(llcrnrlon=boundary_box.lon_min - dx,
                llcrnrlat=boundary_box.lat_min - dx,
                urcrnrlon=boundary_box.lon_max + dx,
                urcrnrlat=boundary_box.lat_max + dx,
                resolution=resolution, projection='tmerc', lon_0=middle_lon, lat_0=middle_lat)

    m.drawcoastlines()
    m.fillcontinents(color='grey', lake_color='aqua')


    #m.drawmapboundary(fill_color='aqua')

    if coordinates_rank == 1:
        lon, lat = np.meshgrid(lon, lat)

    x, y = m(lon, lat)

    # draw filled contours.
    clevs = [0, 0.1, 0.2, 0.3, 0.4, 0.5]
    # clevs = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
    if style == 'contour':
        cs = m.contourf(x, y, mod, clevs, cmap=plt.cm.jet)
        plt.quiver(x[::nx, ::ny], y[::nx, ::ny], ust[::nx, ::ny], vst[::nx, ::ny], scale=scale)

    if style == 'cvector':
        clim = [0., 0.5]
        cs = plt.quiver(x[::nx, ::ny], y[::nx, ::ny], ust[::nx, ::ny], vst[::nx, ::ny], mod[::nx, ::ny], clim=clim, scale=scale)

    if style == 'windbarbs':
        clim = [0., 10]
        cs = plt.barbs(x[::nx, ::ny], y[::nx, ::ny], ust[::nx, ::ny], vst[::nx, ::ny], mod[::nx, ::ny], clim=clim,cmap=plt.cm.jet)

    if style == 'stream':
        print(lons.shape)
        print(x[::, 1])
        print(ust.shape)
        xx = x[1::, ::]
        xxx = xx[0]
        yyy = y[::, 1::][1]

        u = ust[::, ::]
        v = vst[::, ::]
        print(len(xxx))
        print(len(yyy))
        print(len(u))
        print(len(v))

        cs = plt.streamplot(lon, lat, u, v)

    cbar = m.colorbar(cs, location='bottom', pad="5%")
    cbar.ax.tick_params(labelsize=9)
    cbar.set_label('m/s')

    # add title
    plt.title(title)

    fig.savefig(file_out, dpi=100, facecolor='w', edgecolor='w', format='png',
                transparent=False,  bbox_inches='tight', pad_inches=0.1)
    plt.clf()
    plt.close(fig)

    return