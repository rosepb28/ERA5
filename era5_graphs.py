import os
from netCDF4 import Dataset
from netCDF4 import num2date
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader

import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib import cm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap

def c_cmap(file,sheet,clevs,inv=False):
    import pandas as pd
    from matplotlib.colors import ListedColormap, BoundaryNorm
    
    df = pd.read_excel(file,sheet) # lee el archivo excel y divide los valores entre 255
    df = df[['R','G','B']]/255
    if inv:
        df=df.iloc[::-1]   # Si es especificado invierte la barra de colores    
    cc=df.values.tolist()
    cmap = ListedColormap(cc[1:-1])
    cmap.set_under(cc[0])
    cmap.set_over(cc[-1])
    norm = BoundaryNorm(clevs,cmap.N)
    return cmap,norm

def met_variables(file):
    data = Dataset(file)
    for v in data.variables.keys():
        print('%-15s'%v,data.variables[v].long_name)

def coords_times_levels(file):
    data = Dataset(file)    
    lt = data.variables['latitude'][:]
    ln = data.variables['longitude'][:]
    ln,lt = np.meshgrid(ln,lt)
    times=data.variables['time']
    vtimes = num2date(times[:], times.units)
    levels = data['level'][:].tolist()
    return data, ln, lt, vtimes, levels

def save_img(name, level, vtimes, fig):
    '''Esta función guarda la imagen generada, siempre que se le haya indicado con un 'yes'
    a la función externa.

    Argumentos:
    * nombre: nombre del archivo en string y sin extensión.
    * nivel: nivel definido en la función externa.
    * vtimes: tiempos definido en la función externa.
    * fig: grafico definido en la función externa.
    '''
    year = vtimes[0].strftime("%Y")
    month = vtimes[0].strftime("%b")
    folder='Img/{}/{}/'.format(year, month)
    os.makedirs(folder) if not os.path.exists(folder) else 'No'
    fig.savefig(folder + '{}{}.png'.format(name,str(level)), dpi=fig.dpi, bbox_inches= 'tight')
    return

sud = [[-120,-30,-50,15], (16,10.5)]; peru = [[-85,-65,-22,5], (14,8.5)] #extensión, size
def projection_grids(ext, tpe = False, proj = ccrs.PlateCarree()):
    '''Esta función genera la proyecció, extensión, grillas y detalles del mapa
        para las extensiones de Sudamérica o Perú, y una extensión especial para los gráficos TPE.

    Argumentos:
    * ext : puede tomar los string 'sudamerica' o 'peru'.
    * tpe : por defecto es False. Activar con True, solo para graficar TPE.
    * proj : proyección PlateCarree por defecto.
    '''
    size = sud[1] if ext.lower() == 'sudamerica' else peru[1]
    extension = sud[0] if ext.lower() == 'sudamerica' else peru[0]
    prj = proj
    fig = plt.figure(figsize=size, facecolor='white')
    ax = plt.axes(projection = prj)

    ax.add_feature(cfeature.BORDERS,edgecolor='k')
    ax.add_feature(cfeature.COASTLINE,edgecolor='k')
    gl = ax.gridlines(draw_labels=True,linestyle='--')
    gl.top_labels = gl.right_labels = False
    ax.set_extent(extension, crs=prj)
    if ext.lower() == 'peru':
        dir_shp = '../SHP/DEPARTAMENTOS.shp'
        shp = list(shpreader.Reader(dir_shp).geometries())
        ax.add_geometries(shp,prj,edgecolor='black',facecolor='None')
        if tpe:
            ax.set_extent([-90,-65,-20,5], crs=prj)    
    return prj, fig, ax

def lcorr(file, level = 250, time = 0, save_png = 'no'):
    '''Esta función genera un gráfico de líneas de corriente.
    
    Argumentos:
        * archivo : ruta del archivo NetCDF.
        * nivel : nivel en hectopascales, por defecto es 250.
        * tiempo : tiempo a graficar, por defecto es 0.
        * save_png : opción para descargar gráfico, puede tomar los string 'Yes' o 'No'.
    '''
    data, ln, lt, vtimes, levels = coords_times_levels(file)
    u = data.variables['u'][time, levels.index(level), :, :]
    v = data.variables['v'][time, levels.index(level), :, :]

    levs = [0,5,10,15,20,25,30,35,40,50,60,65]
    cmap, norm = c_cmap('Colores/lcorr.xlsx', '200', levs)
    prj, fig, ax = projection_grids('sudamerica')
    sp = np.sqrt(u**2 + v**2)
    lw = 5*sp / sp.max()
    strm = ax.streamplot(ln, lt, u, v, color = sp, cmap=cmap,norm = norm, linewidth=0.8, density=6)
    fig.colorbar(strm.lines, shrink = 0.8, extend='both')

    save_png = save_png.lower()
    if save_png == 'yes':
        save_img(name = 'lcorr', vtimes = vtimes, level = level, fig = fig)

    return plt.show()

def divergencia(file, level = 250, time = 0, save_png = 'No'):
    '''Esta función genera un gráfico de divergencia.

    Argumentos:
    * archivo : ruta del archivo NetCDF.
    * nivel : nivel en hectopascales, por defecto es 250.
    * tiempo : tiempo a graficar, por defecto es 0.
    * save_png : opción para descargar gráfico, puede tomar los string 'Yes' o 'No'.
    '''
    data, ln, lt, vtimes, levels = coords_times_levels(file)
    div = data['d'][time,levels.index(level),:,:]*1e06

    levs = np.arange(-16,18,2)
    prj, fig, ax = projection_grids(ext = 'peru')    
    cx = ax.contourf(ln,lt,div,levs,transform=prj,cmap = 'RdBu_r', extend='both')
    cc=plt.colorbar(cx,shrink=0.9,drawedges=True)

    save_png = save_png.lower()
    if save_png == 'yes':
        save_img(name = 'div',vtimes = vtimes, level = level, fig = fig)

    return plt.show()

def hrp600_200(file, level = 600, time = 0, save_png = 'No'):
    '''Esta función genera un gráfico de humedad relativa promedio.
    Argumentos:
    * archivo : ruta del archivo NetCDF.
    * nivel : nivel inferior en hectopascales, por defecto es 600.
    * tiempo : tiempo a graficar, por defecto es 0.
    * save_png : opción para descargar gráfico, puede tomar los string 'Yes' o 'No'.
    '''
    data, ln, lt, vtimes, levels = coords_times_levels(file)
    r_sum = 0; u_sum = 0; v_sum = 0
    for i in range(levels.index(level)+1):
        r = data['r'][time,i,:,:]; u = data['u'][time,i,:,:]; v = data['v'][time,i,:,:]
        r_sum = r_sum + r; u_sum = u_sum + u; v_sum = v_sum + v
    hrp = (r_sum/len(levels[:levels.index(level)]))
    up = (u_sum/len(levels[:levels.index(level)]))
    vp = (v_sum/len(levels[:levels.index(level)]))

    levs = np.arange(0,110,10)
    cmap, norm = c_cmap('Colores/hrp2.xlsx', '550', levs)
    prj, fig, ax = projection_grids(ext = 'peru')  
    cx = ax.contourf(ln,lt,hrp,levs, transform=prj, cmap=cmap, norm = norm, extend='both')
    ax.quiver(ln[::4,::4],lt[::4,::4],up[::4,::4],vp[::4,::4], color = 'w', scale = 80)
    plt.colorbar(cx)

    save_png = save_png.lower()
    if save_png == 'yes':
        save_img(name = 'hrp',vtimes = vtimes, level = level, fig = fig)

    return plt.show()

def sp_hum(file, level = 550, time = 0, save_png = 'No'):
    '''Esta función genera un gráfico de humedad específica.

    Argumentos:
    * archivo : ruta del archivo NetCDF.
    * nivel : nivel en hectopascales, por defecto es 550.
    * tiempo : tiempo a graficar, por defecto es 0.
    * save_png : opción para descargar gráfico, puede tomar los string 'Yes' o 'No'.
    '''
    data, ln, lt, vtimes, levels = coords_times_levels(file)
    q = data['q'][time, levels.index(level),:,:]*1000
    u = data.variables['u'][time, levels.index(level), :, :]
    v = data.variables['v'][time, levels.index(level), :, :]

    levs = np.arange(3,6,0.5)
    cmap, norm = c_cmap('Colores/rmix.xlsx', '550', levs)
    prj, fig, ax = projection_grids(ext = 'peru')  
    sn = 5
    cx = ax.contourf(ln,lt,q,levs, transform=prj, cmap=cmap, norm= norm, extend='both')
    ax.quiver(ln[::sn,::sn],lt[::sn,::sn],u[::sn,::sn],v[::sn,::sn], color = 'k', scale = 50, headwidth=5, headlength=6)
    plt.colorbar(cx)

    save_png = save_png.lower()
    if save_png == 'yes':
        save_img(name = 'q',vtimes = vtimes, level = level, fig = fig)
    
    return plt.show()

def tpe(file, level = 950, time = 0, save_png = 'No'):
    """
    Esta función genera un gráfico de temperatura potencial equivalente para la costa del Perú.

    Argumentos:
    * archivo : ruta del archivo NetCDF.
    * nivel : nivel en hectopascales, por defecto es 950.
    * tiempo : tiempo a graficar, por defecto es 0.
    * save_png : opción para descargar gráfico, puede tomar los string 'Yes' o 'No'.
    """
    
    data, ln, lt, vtimes, levels = coords_times_levels(file)
    ln = ln[100:201,120:221]; lt = lt[100:201,120:221]
    t  = data.variables['t'][time,levels.index(level),100:201,120:221]
    r = data['r'][0,levels.index(level),100:201,120:221]
    u  = data.variables['u'][time,levels.index(level),100:201,120:221]
    v  = data.variables['v'][time,levels.index(level),100:201,120:221]
    es=6.11*pow(10,((7.5*(t-273))/(237.5+(t-273)))); e=es*r/100
    rm=0.622*e/(level-e); te=t*(1+(2529.8804*rm/t)); tpe=te*pow(1000/950,0.286)

    levs = np.arange(306,332,2)
    cmap, norm = c_cmap('Colores/tpe.xlsx', '950', levs)
    prj, fig, ax = projection_grids(ext = 'peru', tpe = True)  
    sn = 7
    cx = ax.contourf(ln,lt,tpe.round(2),levs, transform=prj, cmap=cmap, norm= norm, extend='both')
    cs = ax.contour(ln,lt,tpe.round(2),levs,transform=prj, colors='w',linewidths=0.8)
    ax.clabel(cs, cs.levels[6:], colors='k', inline=1, fontsize=10, inline_spacing=1)
    ax.quiver(ln[::sn,::sn],lt[::sn,::sn],u[::sn,::sn],v[::sn,::sn], color = 'k', scale = 200,
             headwidth=3, headlength=3)

    plt.colorbar(cx)
    save_png = save_png.lower()
    if save_png == 'yes':
        save_img(name = 'tpe',vtimes = vtimes, level = level, fig = fig)

    return plt.show()