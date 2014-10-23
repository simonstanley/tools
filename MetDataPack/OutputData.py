import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import iris.plot as iplt
import cartopy.crs as ccrs
import numpy
import scipy.stats
import math
import iris
from collections import OrderedDict
from Utils.iris_functions import get_coordinate_slice_dimensions, \
                                 cube_time_converter, \
                                 get_cubelist_data_range, get_xy_coords
from Utils.plotting_functions import get_sensible_contour_levels


def scale_axis_limits(amin, amax, scale):
    """
    Extend the min and max values by the given scale.
    
    """
    mean = (amax + amin) / 2.
    ax_range = abs(amax - mean)
    amin = mean - (ax_range * scale)
    amax = mean + (ax_range * scale)
    
    return amin, amax

def constrain_bounds(ax, lon, lat):
    
    xmax = numpy.max(lon.points)
    xmin = numpy.min(lon.points)
    ymax = numpy.max(lat.points)
    ymin = numpy.min(lat.points)
    
    ax.axis([xmin, xmax, ymin, ymax])

def regrid_cube(cube, grid_cube):
    return iris.analysis.interpolate.regrid(cube, grid_cube)

class Plot(object):
    """
    Class with plotting techniques for iris cubes.
    
    Kwargs:
    
    * plot_type: iris plotting method
        Default is iplt.contourf (colour filled contour plot).
    
    * projection: cartopy projection
        Default is ccrs.PlateCarree().
    
    """
    def __init__(self, plot_type=iplt.contourf, projection=ccrs.PlateCarree()):
        self.plot_type  = plot_type
        self.projection = projection
        self.configure_coord_names()
        
    def _plot(self, cube, *args, **kwargs):
        """
        Perform the plot.
        
        """
        return self.plot_type(cube, *args, **kwargs)
    
    def _arrow_plot(self, x_cube, y_cube, resolution, xy_coords, **kwargs):
        """
        Perform the arrow plot.
        
        """
        lon_dim = x_cube.coord_dims(xy_coords[0].name())[0]
        lat_dim = x_cube.coord_dims(xy_coords[1].name())[0]
        
        shape = x_cube.data.shape
        max_axis = max(shape)
        scaling_val = max_axis / float(resolution)
        shape = tuple([int(x/scaling_val) for x in shape])
        
        arr = numpy.zeros(shape)
        grid_cube = iris.cube.Cube(arr)
        
        lon = xy_coords[0]
        lat = xy_coords[1]
        
        new_lon = iris.coords.DimCoord(numpy.linspace(numpy.min(lon.points), 
                                                      numpy.max(lon.points), 
                                                      shape[lon_dim]), 
                                       standard_name='longitude', 
                                       units=lon.units,
                                       coord_system=lon.coord_system)
        new_lat = iris.coords.DimCoord(numpy.linspace(numpy.min(lat.points),   
                                                      numpy.max(lat.points), 
                                                      shape[lat_dim]), 
                                       standard_name='latitude', 
                                       units=lat.units,
                                       coord_system=lat.coord_system)

        
        grid_cube.add_dim_coord(new_lon, lon_dim)
        grid_cube.add_dim_coord(new_lat, lat_dim)
        
        x_cube = regrid_cube(x_cube, grid_cube)
        y_cube = regrid_cube(y_cube, grid_cube)

        U = x_cube.data.transpose([lat_dim, lon_dim])
        V = y_cube.data.transpose([lat_dim, lon_dim])
        X = x_cube.coord('longitude').points
        Y = y_cube.coord('latitude').points
         
        plt.quiver(X, Y, U, V, units='xy', headwidth=2,
                   transform=x_cube.coord('latitude').coord_system.\
                   as_cartopy_projection(), **kwargs)
        

    def _add_colourbar(self, fig, plot, units):
        """
        Add colourbar to the figure and label with units.
        
        """
        colorbar_axes = fig.add_axes([0.1, 0.1, 0.82, 0.03])
        cbar = plt.colorbar(plot, colorbar_axes, orientation='horizontal')
        if units and str(units) != 'unknown':
            cbar.set_label(units)
        else:
            cbar.set_label('')
        
    def _postage_title(self, cube, label_mems, label_ref_dates):
        """
        Build a title using realization numbers and forecast reference dates.
        
        """
        title = ''
        if label_mems:
            try:
                title += '%s: %s' % (self.realization.title(),
                                     cube.coord(self.realization).points[0])
            except:
                pass
        if label_ref_dates:
            try:
                time_unit = cube.coord(self.time_coord).units
                fcast_ref = cube_time_converter(
                            cube.coord(self.forecast_ref_time).points[0],
                            time_unit)
                title += '\nInit. date: %s' % fcast_ref.isoformat()[:10]
            except:
                pass
        return title
        
    def map_plot(self, cube, title=None, filename=None, show=False, 
                  num_of_levels=15, figsize=(9,9), coastlines=True, 
                  colourbar=True, frame=True, **kwargs):
        """
        Plot a map.
        
        Args:
        
        * cube: iris cube
            Cube to be plotted, must contain a single x/y array.
        
        Kwargs:
        
        * title: string
            Add a title to the plot.
        
        * filename: string
            Specify a full filename and path for the plot to be saved to. If no
            filename is given, no save will be attempted.
        
        * show: boolean
            Set to True to have the plot displayed on screen.
        
        * num_of_levels: integer
            Set the number of data value levels. This is overwritten if 
            specific level values are given in key word arguments.
        
        * figsize: tuple
            Specify the x and y lengths of the figure, default is (9,9).
        
        * coastlines: boolean
            Plot the coastlines on the map, default True
        
        * colourbar: boolean
            Add the colourbar, default True.
        
        * frame: boolean
            Plot a frame around the map, default True.
        
        * kwargs:
            Optional key word arguments for matplotlib.pyplot plot type.
        
        """
        fig = plt.figure(figsize=figsize)
        ax = plt.axes(projection=self.projection)
        ax.outline_patch.set_visible(frame)
        mapplt = self._plot(cube, num_of_levels, **kwargs)
        
        if title:
            plt.title(title, size='large')
        if coastlines:
            plt.gca().coastlines(resolution='50m')
        if colourbar:
            self._add_colourbar(fig, mapplt, cube.units)
        if filename is not None:
            plt.savefig(filename, bbox_inches="tight", pad_inches=0)
        if show:
            plt.show()
            
    def postage_plot(self, cubes, title='', filename=None, show=False, 
                       num_of_levels=15, label_mems=True, label_ref_dates=True, 
                       labels=None, rows_and_cols=None, figsize=(12,10), 
                       coastlines=True, colourbar=True, **kwargs):
        """
        Plot multiple maps.
        
        Args:
        
        * cubes: iris cube or cubelist        # We now have actual spatial data and a spatial "forecast" using the 
        # climatology mapping method. Get the spatial difference.
            Cubes to be plotted. If a single cube is provided, it must contain 
            one or more x/y array realizations. If a cubelist is provided each
            cube must contain a single x/y array.
        
        Kwargs:
        
        * title: string
            Add a title to the plot.
        
        * filename: string
            Specify a full filename and path for the plot to be saved to. If no
            filename is given, no save will be attempted.
        
        * show: boolean
            Set to True to have the plot displayed on screen.
        
        * num_of_levels: integer
            Set the number of data value levels. This is overwritten if 
            specific level values are given in key word arguments.
        
        * label_mems: boolean
            Attempt to label each plot with a realization number, default True.
        
        * label_ref_dates: boolean
            Attempt to label each plot with a initialisation date, default 
            True.
        
        * rows_and_cols: list
            Manually specify the number of rows and columns. Must be a list of 
            two integers.
        
        * figsize: tuple
            Specify the x and y lengths of the figure, default is (12,10).
        
        * coastlines: boolean
            Plot the coastlines on the map, default True
        
        * colourbar: boolean
            Add the colourbar, default True.
        
        * kwargs:
            Optional key word arguments for matplotlib.pyplot plot type.
        
        """
        if type(cubes) == iris.cube.Cube:
            real_slice = get_coordinate_slice_dimensions(
                                              cubes, 
                                              [self.forecast_ref_time,
                                               self.realization],
                                              ignore_missing_coords=True)
            cubes = [xyslice for xyslice in cubes.slices(real_slice)]
        
        # If contour levels are not specified, calculate them from the data 
        # range. Only when contour levels are specified for all sub-plots can a
        # consistent colour bar be produced.
        if kwargs.get('levels') is None:
            dmin, dmax = get_cubelist_data_range(cubes)
            kwargs['levels'] = get_sensible_contour_levels(dmin, dmax, 
                                                           num_of_levels)

        
        if not rows_and_cols:
            # Work out an appropriate number of rows and columns depending on 
            # the number of members.
            length = len(cubes)
            cols   = math.ceil((math.sqrt(length)))
            rows   = math.ceil(length / cols)
        else:
            assert len(rows_and_cols) == 2, '2 values must be specified for '\
                                            'rows_and_cols argument. % '\
                                            'provided.' % len(rows_and_cols)
            rows = rows_and_cols[0]
            cols = rows_and_cols[1]
        
        fig = plt.figure(figsize=figsize)
        for i, cube in enumerate(cubes):
            plt.subplot(rows, cols, i+1, projection=self.projection)
            mapplt = self._plot(cube, **kwargs)
            if labels is not None:
                title_str = labels[i]
            else:
                title_str = self._postage_title(cube, label_mems, 
                                                label_ref_dates)
            plt.title(title_str, size='small')
            if coastlines:
                plt.gca().coastlines(resolution='50m')
        plt.suptitle(title, size='large')
        
        if colourbar:
            self._add_colourbar(fig, mapplt, cube.units)
        if filename is not None:
            plt.savefig(filename, bbox_inches="tight", pad_inches=0.05)
        if show:
            plt.show()

    def plume_plot(self, cube, title='', filename=None, show=False,
                    colours=None, yaxis_scale_increase=1.7, zero_line=False, 
                    **kwargs):
        """
        Plot a time line of all members within a cube.
        
        Args:
        
        * cube: iris cube
            Cube to be plotted, must contain one or more time step and 
            realization and also represent a spatial average.
        
        Kwargs:
        
        * title: string
            Add a title to the plot.
        
        * filename: string
            Specify a full filename and path for the plot to be saved to. If no
            filename is given, no save will be attempted.
        
        * show: boolean
            Set to True to have the plot displayed on screen.
        
        * colours: list
            Specify a list of colours to plot the lines with. Colours are 
            repeated if not enough are given, e.g. if ['green', 'red'] is given
            and four colours are needed, ['green', 'red', 'green', 'red'] is 
            used. Default uses all black.
        
        * yaxis_scale_increase: float
            Extend the y axis by this scalar so plumes are not squashed.
        
        * zero_line: boolean
            If set True, a dashed line is added to y = 0.
            
        * kwargs:
            Optional key word arguments for matplotlib.pyplot plot type.
        
        """
        time_unit = cube.coord('time').units        
        myFmt = mdates.DateFormatter('%d %b')
        
        fig = plt.figure(figsize=(11,8))
        ax  = plt.axes()

        fcst_refs = sorted(set(cube.coord(self.forecast_ref_time).points),
                           reverse=True)
        fcst_refs = [cube_time_converter(time, time_unit) 
                     for time in fcst_refs]
        
        if colours:
            # Compare the number of colours with the number of initialisation
            # dates to see how many time the given colours need to be looped 
            # through.
            repeat_colours = math.ceil(float(len(fcst_refs)) / \
                                       float(len(colours)))
            colours = colours * int(repeat_colours)
        else:
            colours = ['k'] * len(fcst_refs)
        labels  = [date.strftime('%d/%m/%Y') for date in fcst_refs]
        
        
        real_slice = get_coordinate_slice_dimensions(cube, 
                                                     self.realization)
        for i, fcst_ref in enumerate(fcst_refs):
            fcst_ref_cube = cube.extract(
                            iris.Constraint(
                            **{self.forecast_ref_time:fcst_ref}))
                        
            for realization in fcst_ref_cube.slices(real_slice):
                iplt.plot(realization, color=colours[i], label=labels[i])
        
        if colours:
            # Add a legend.
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = OrderedDict(zip(labels, handles))
            ax.legend(by_label.values(), by_label.keys(), loc='lower left', 
                      fontsize='small', frameon=False, title='Init. Dates')
        
        ax.xaxis.set_major_formatter(myFmt)
        plt.xticks(rotation=40)
        
        # Stretch the y axis limits by given scale.
        ymin, ymax = plt.ylim()
        new_min, new_max = scale_axis_limits(ymin, ymax, yaxis_scale_increase)
        plt.ylim( (new_min, new_max) )
        
        plt.title(title, size='small')
        plt.ylabel(cube.units)
        plt.grid(which='major', axis='y')
        
        if zero_line:
            ax.axhline(y=0, ls='--', c='k')
        if filename:
            fig.tight_layout()
            plt.savefig(filename)
        if show:
            plt.show()

    def pdf_plot(self, cube, title='', filename=None, show=False,
                 xaxis_scale_increase=2, yaxis_scale_increase=1.1, 
                 data_range=None):
        """
        Plot the probability density of the members for the time spatial 
        point.

        Args:
        
        * cube: iris cube
            Cube to be plotted, must contain one time step and one or more
            realizations and also represent a spatial average.
        
        Kwargs:
        
        * title: string
            Add a title to the plot.
        
        * filename: string
            Specify a full filename and path for the plot to be saved to. If no
            filename is given, no save will be attempted.
        
        * show: boolean
            Set to True to have the plot displayed on screen.
        
        * xaxis_scale_increase: float
            Increase the width of the x axis by given factor. Default 2.
        
        * yaxis_scale_increase: float
            Increase the height of the y axis by given factor. Default 1.1.
        
        * data_range: list of 2 floats or None
            This sets the x axis limits according to the given min and max data 
            values (xaxis_scale_increase is still taken into account). This is 
            used for consistent x axes when comparing multiple plots. Default
            is None.

        """                
        num_bins = max(4, len(cube.data)//4)
        
        fig = plt.figure()
        # Get return data so bins can be used to set x axis ticks.
        plt.hist(cube.data, num_bins, normed=True, facecolor='grey', 
                 alpha=0.2, hatch='x', lw=2)
        
        if data_range:
            if len(data_range) != 2:
                raise UserWarning('data_range must have 2 values, [xmin, xmax'\
                                  ']')
            # This situation is for consistent x limits across multiple plots. 
            # data_min and data_max should be from the whole data across all
            # plots.
            new_xmin, new_xmax = scale_axis_limits(data_range[0], 
                                                   data_range[1], 
                                                   xaxis_scale_increase)
        else:
            xmin, xmax = plt.xlim()
            new_xmin, new_xmax = scale_axis_limits(xmin, xmax, 
                                                   xaxis_scale_increase)
            
        ymin, ymax = plt.ylim()
        _, new_ymax = scale_axis_limits(ymin, ymax, 
                                        yaxis_scale_increase)
        # Keep min as 0.
        new_ymin = 0
        plt.xlim( (new_xmin, new_xmax) )
        plt.ylim( (new_ymin, new_ymax) )
        # Generate kernel density estimate (PDF)
        kde = scipy.stats.gaussian_kde(cube.data.flatten(),
                                       bw_method="silverman")
        pdf_xval = numpy.linspace(new_xmin, new_xmax, 50)
        pdf = kde(pdf_xval)

        plt.plot(pdf_xval, pdf, 'r--', lw=2)
        plt.title(title, size='small')
        plt.xlabel(cube.units)
        plt.grid(which='major', axis='y')

        locs, _ = plt.yticks()
        labels = [round(x, 2) for x in locs]
        plt.yticks(locs, labels)
        
        if filename:
            fig.tight_layout()
            plt.savefig(filename)
        if show:
            plt.show()

    def arrow_map_plot(self, x_cube, y_cube, title='', filename=None, 
                        show=False, resolution=40, figsize=(9,9), 
                        background_map=True, frame=True, **kwargs):
        """
        Plot an arrow map.
        
        Args:
        
        * x_cube: iris cube
            Cube with x direction data, must contain a single x/y array.
        
        * y_cube: iris cube
            Cube with y direction data, must contain a single x/y array.
        
        Kwargs:
        
        * title: string
            Add a title to the plot.
        
        * filename: string
            Specify a full filename and path for the plot to be saved to. If no
            filename is given, no save will be attempted.
        
        * show: boolean
            Set to True to have the plot displayed on screen.
        
        * resolution: integer > 1
            A scaling number for the density of arrows within the plot. The 
            higher the number, the more (and smaller) arrows there are. Default
            is 40.
                
        * figsize: tuple
            Specify the x and y lengths of the figure, default is (9,9).
        
        * background_map: boolean
            Plot the earth (land, sea, coastlines etc.) under the arrows, 
            default True.
        
        * frame: boolean
            Plot a frame around the map, default True.
        
        * kwargs:
            Optional key word arguments for matplotlib.pyplot plot type.
        
        """
        xy_coords = get_xy_coords(x_cube)
        
        plt.figure(figsize=figsize)
        ax = plt.axes(projection=self.projection)
        ax.outline_patch.set_visible(frame)
        self._arrow_plot(x_cube, y_cube, resolution, xy_coords, **kwargs)
        if background_map:
            ax.coastlines()
            ax.gridlines()
            ax.stock_img()
        
        constrain_bounds(ax, 
                         xy_coords[0], 
                         xy_coords[1])
        
        plt.title(title, size='large')
                
        if filename:
            plt.savefig(filename, bbox_inches="tight", pad_inches=0.05)
        if show:
            plt.show()

    def arrow_postage_plot(self, x_cubes, y_cubes, title='', filename=None, 
                             show=False, resolution=20, label_mems=True, 
                             label_ref_dates=True, rows_and_cols=None, 
                             figsize=(12,10), background_map=True, **kwargs):
        """
        Plot multiple arrow maps. x cubes and y cubes must match in size and 
        order.
        
        Args:
        
        * x_cubes: iris cube or cubelist
            Cubes with x direction data. If a single cube is provided, it must 
            contain one or more x/y array realizations. If a cubelist is 
            provided each cube must contain a single x/y array.
        
        * y_cubes: iris cube or cubelist
            Cubes with y direction data. If a single cube is provided, it must 
            contain one or more x/y array realizations. If a cubelist is 
            provided each cube must contain a single x/y array.
        
        Kwargs:
        
        * title: string
            Add a title to the plot.
        
        * filename: string
            Specify a full filename and path for the plot to be saved to. If no
            filename is given, no save will be attempted.
        
        * show: boolean
            Set to True to have the plot displayed on screen.
        
        * resolution: integer > 1
            A scaling number for the density of arrows within the plot. The 
            higher the number, the more (and smaller) arrows there are. Default
            is 20.
        
        * label_mems: boolean
            Attempt to label each plot with a realization number, default True.
        
        * label_ref_dates: boolean
            Attempt to label each plot with a initialisation date, default 
            True.
        
        * rows_and_cols: list
            Manually specify the number of rows and columns. Must be a list of 
            two integers.
        
        * figsize: tuple
            Specify the x and y lengths of the figure, default is (9,9).
        
        * background_map: boolean
            Plot the earth (land, sea, coastlines etc.) under the arrows, 
            default True.
                
        * kwargs:
            Optional key word arguments for matplotlib.pyplot plot type.
        
        """
        diff_types_err = 'x_cubes and y_cubes are differing types.'
        if type(x_cubes) == iris.cube.Cube:
            if type(y_cubes) != iris.cube.Cube:
                raise UserWarning(diff_types_err)
                    
            real_slice = get_coordinate_slice_dimensions(
                                                  x_cubes, 
                                                  [self.forecast_ref_time,
                                                   self.realization],
                                                  ignore_missing_coords=True)
            x_cubes = [xyslice for xyslice in x_cubes.slices(real_slice)]
            y_cubes = [xyslice for xyslice in y_cubes.slices(real_slice)]
        
        else:
            if type(y_cubes) == iris.cube.Cube:
                raise UserWarning(diff_types_err)
        
        assert len(x_cubes) == len(y_cubes), 'x_cubes and y_cubes have '\
                                             'differing lengths.'
                                             
        xy_coords = get_xy_coords(x_cubes[0])
        
        if not rows_and_cols:
            # Work out an appropriate number of rows and columns depending on 
            # the number of members.
            length = len(x_cubes)
            cols   = math.ceil((math.sqrt(length)))
            rows   = math.ceil(length / cols)
        else:
            assert len(rows_and_cols) == 2, '2 values must be specified for '\
                                            'rows_and_cols argument. % '\
                                            'provided.' % len(rows_and_cols)
            rows = rows_and_cols[0]
            cols = rows_and_cols[1]
                
        plt.figure(figsize=figsize)
        for i, (x_cube, y_cube) in enumerate(zip(x_cubes, y_cubes)):
            ax = plt.subplot(rows, cols, i+1, projection=self.projection)
            self._arrow_plot(x_cube, y_cube, resolution, xy_coords, **kwargs)
            if background_map:
                ax.coastlines()
                ax.gridlines()
                ax.stock_img()
            constrain_bounds(ax, 
                             xy_coords[0], 
                             xy_coords[1])
            title_str = self._postage_title(x_cube, label_mems, 
                                            label_ref_dates)
            plt.title(title_str, size='x-small')
        plt.suptitle(title, size='large')
                
        if filename:
            plt.savefig(filename, bbox_inches="tight", pad_inches=0.05)
        if show:
            plt.show()

    def configure_coord_names(self, member_name='realization', 
                     period_name='time', 
                     initialistion_time_name='forecast_reference_time'):
        """
        Parse the coordinate names and the values to which they are to be set.
        
        """
        self.realization = member_name
        self.time_coord = period_name
        self.forecast_ref_time = initialistion_time_name
    