"""
Module for loading meteorological data and performing basic processing.

Notes:
All string format searches comply with Python regular expression formats.

Author: S Stanley
Date: March 2014
"""

from Utils.iris_functions import cube_time_converter, \
                                 get_time_bound_constraints, \
                                 get_list_of_time_bound_constraints, \
                                 get_coordinate_slice_dimensions, \
                                 get_xy_coords, \
                                 AreaBounds, \
                                 pad_coords
from Utils.useful_functions import search_directory, \
                                   get_all_dates_between_time_bounds, \
                                   is_leap_year
import iris
import iris.analysis
import numpy
import datetime
import math
import multiprocessing

# This allows date constraining on datetime objects.
iris.FUTURE.cell_datetime_objects = True

def ConfigureDefaults(area_bounds=None, 
                        area_bounds_format=['x_min','y_min','x_max','y_max'], 
                        area_bounds_range=None, years_are_bounds=False,
                        dates_are_bounds=False, init_date_str_format='%y%m%d',
                        member_name='realization', period_name='time', 
                        initialistion_time_name='forecast_reference_time'):
    """
    Configure default settings for all classes within module. Values can all be
    re-adjusted for each class separately, within each class. This class is for
    constant values over all classes only.
    
    """        
    global default_area_bounds
    global default_area_bounds_format
    global default_area_bounds_range
    global default_years_are_bounds
    global default_dates_are_bounds
    global default_init_date_str_format
    global default_member_name
    global default_period_name
    global default_initialistion_time_name
    
    default_area_bounds = area_bounds
    default_area_bounds_format = area_bounds_format
    default_area_bounds_range = area_bounds_range
    default_years_are_bounds = years_are_bounds
    default_dates_are_bounds = dates_are_bounds
    default_init_date_str_format = init_date_str_format
    default_member_name = member_name
    default_period_name = period_name
    default_initialistion_time_name = initialistion_time_name

ConfigureDefaults()

def print_defaults():
    """
    Print all the current default values.
    
    """
    print 'area_bounds :', default_area_bounds
    print 'area_bounds_format :', default_area_bounds_format
    print 'area_bounds_range :', default_area_bounds_range
    print 'years_bounds :', default_years_are_bounds
    print 'dates_are_bounds :', default_dates_are_bounds
    print 'init_date_str_format :', default_init_date_str_format
    print 'member_name :', default_member_name
    print 'period_name :', default_period_name
    print 'initialistion_time_name :', default_initialistion_time_name

def remove_attributes(cube, field, filename):
    """
    Callback function to remove attributes.
    
    """
    cube.attributes = None

def remove_coords(cube, unwanted_coords):
    """
    Remove all given coordinates from cube.
    
    Args:
    
    * cube: iris cube
    
    * unwanted_coords: list of strings
        List of coordinate names.
        
    Returns:
        iris cube
    
    """
    if type(unwanted_coords) != list:
        unwanted_coords = [unwanted_coords]
    for coord in unwanted_coords:
        try:
            cube.remove_coord(coord)
        except iris.exceptions.CoordinateNotFoundError:
            continue
    
    return cube

def check_dates(dates):
    """
    Check that all dates in a list are datetime objects.
    
    Args:
    
    * dates: list
    
    """
    for date in dates:
        if type(date) != datetime.datetime:
            raise TypeError('Input date, %s, not datetime object' % date)

def change_zeroth_hour(dates):
    """
    Make sure date is not set at 0 hours of the day, this causes extraction 
    problems because of bounds from previous day including the 0th hour of
    the next day.
    
    Args:
    
    * dates: list of datetimes
    
    Returns:
        list of datetimes
    
    """
    updated_dates = []
    for date in dates:
        if date.hour == 0:
            date += datetime.timedelta(hours=12)
        updated_dates.append(date)
    return updated_dates

def run_load_file(args):
    """
    
    """
    return _DataHandler._load_file(*args)
    

class _DataHandler(object):
    """
    Parent class containing all common methods.
    
    """            
    @staticmethod
    def _find_files(directory, dirs_to_look_in, files_to_search_for, 
                     current_dir, see_files):
        """
        Find all matches for the search criteria.

        Args:
        
        * directory: string
            Top directory to start search.

        * dirs_to_look_in: list of strings or string
            Specify any sub directories of the main directory (given during the 
            instantiation of the class) of which to search in. Default is to
            search all files and directories beneath the main directory.
        
        * files_to_search_for: list of strings or string
            The default loading technique is to find filenames which contain 
            string versions of the initialisation dates. If this is not 
            appropriate, specify here the string matches to search for in the 
            filenames. Use see_files to view the search results before loading.
        
        * current_dir: boolean:
            Set True to only search the files within the given directory.
        
        * see_files: boolean
            Set to True to return search results without loading the files. 
            Default is False.
            
        Returns:
            list
        
        """
        full_name = True
        if see_files:
            full_name = False
        files_to_load = search_directory(directory, 
                                         look_in=dirs_to_look_in,
                                         search_for=files_to_search_for,
                                         file_type='files',
                                         current_dir=current_dir,
                                         full_name=full_name)
        if not files_to_load:
            raise UserWarning('No files were found matching the search for %s'\
                              ' in the directory(s) %s%s' \
                              % (files_to_search_for, directory, 
                                 dirs_to_look_in))
        return files_to_load

    def _load_file(self, load_file):
        """
        Load file with iris and return all cubes within that file.

        Args:
        
        * files_to_load: list
                
        Returns:
            list of iris cubes
        
        """
        try:
            loaded_cubelist = iris.load(load_file, self.variable, 
                                        self._callback)
        except ValueError:
            return None
            
        if not loaded_cubelist:
            return None
        return [remove_coords(cube, self.unwanted_coords) 
                for cube in loaded_cubelist]
        
    def _load_all_cubes(self, files_to_load):
        """
        Load data, gather information and remove unwanted coordinates.
        
        Args:
        
        * files_to_load: list
        
        Returns:
            iris cubelist
        
        """
        if self.process_workers:
            arguments = [[self, load_file] for load_file in files_to_load]
            pool = multiprocessing.Pool(processes=self.process_workers)
            try:
                all_cubelists = pool.map(run_load_file, arguments)
                pool.close()
                pool.join()
            except KeyboardInterrupt:
                pool.terminate()
        else:
            all_cubelists = []
            for load_file in files_to_load:
                cubelist = self._load_file(load_file)
                if cubelist:
                    all_cubelists.append(cubelist)
        
        all_cubes = []
        for cubelist in all_cubelists:
            for cube in cubelist:
                all_cubes.append(cube)

        if len(all_cubes) == 0:
            raise UserWarning('No data loaded.')
        
        # Gather universal information from the first cube.
        if self.xy_coords is None:
            self.xy_coords = [coord.name() 
                              for coord in get_xy_coords(
                                           all_cubes[0])]
        if self._area_inst.bounds_range is None:
            self._area_inst.bounds_range = self._area_inst.\
                                           get_cube_area_bounds(all_cubes[0],
                                                                self.xy_coords)
        if self.area_bounds is None:
            self.area_bounds = self._area_inst.get_cube_area_bounds(
                                               all_cubes[0],
                                               self.xy_coords)
        self.time_unit = all_cubes[0].coord(self.time_coord).units
                
        return iris.cube.CubeList(all_cubes)

    def _sort_cubelist(self, cubelist):
        """
        Sort all data in cubelist into a single cube where all realisations are
        unique, missing time steps are masked and area bounds are extracted.
        
        Args:
        
        * cubelist: iris cubelist
            cubelist to be sorted and merged.
        
        Returns:
            iris cube
            
        """
        sorted_cubelist = []
        realization_num = 1
        cubelist = cubelist.merge(unique=False)
        for cube in cubelist:
            # If time is a scalar coordinate, promote it to a dimension 
            # coordinate, this is because all cubes must have the same number 
            # of dimensions to be compared.
            if len(cube.coord(self.time_coord).points) == 1:
                cube = iris.util.new_axis(cube, scalar_coord=self.time_coord)
            
            # Chop cubes into individual realizations for relabelling.
            member_slices = get_coordinate_slice_dimensions(
                            cube, [self.realization,self.forecast_ref_time],
                            ignore_missing_coords=True)
            for member_slice in cube.slices(member_slices):
                
                if self.realization in [coord.name() 
                                     for coord in member_slice.coords()]:
                    member_slice.coord(
                           self.realization).points = [realization_num]
                else:
                    realization_coord = iris.coords.AuxCoord([realization_num],
                                                             self.realization)
                    member_slice.add_aux_coord(realization_coord)
                
                member_slice.cell_methods = None
                sorted_cubelist.append(member_slice)
                realization_num += 1
        
        sorted_cubelist = iris.cube.CubeList(sorted_cubelist)
        # Mask missing time steps so merging can be done.
        sorted_cubelist = pad_coords(sorted_cubelist, self.time_coord)
        cube = sorted_cubelist.merge_cube()
        # Check x-y coordinates match the specified range.
        cube = self._area_inst.check_cube_area_bounds(cube, self.xy_coords, 
                                                     self.area_bounds)
        cube = self.extract_area_bounds(cubes=cube)
        
        if cube.coord_dims(cube.coord(self.realization)) == \
           cube.coord_dims(cube.coord(self.forecast_ref_time)):
            # Re order realizations in initialisation date order.
            ordered_inits = sorted(cube.coord('forecast_reference_time').points)
            ordered_mems  = range(1, len(cube.coord('realization').points)+1)
            ordered_cubes = []
            for member_slice in cube.slices(member_slices):
                mem_index = ordered_inits.index(
                            member_slice.coord(self.forecast_ref_time).points[0])
                member_slice.coord('realization').points = ordered_mems[mem_index]
                del ordered_inits[mem_index]
                del ordered_mems[mem_index]
                ordered_cubes.append(member_slice)
            cube = iris.cube.CubeList(ordered_cubes).merge_cube()
        
        return cube
    
    def _get_model_metadata(self, cube):
        """
        Retreive all relevant metadata.
        
        Args:
        
        * cube: iris cube
        
        Returns:
            dictionary
        
        """
        init_coord = cube.coord(self.forecast_ref_time)
        init_dates = [cube_time_converter(time, init_coord.units) 
                      for time in set(init_coord.points)]
        
        time_coord = cube.coord(self.time_coord)
        fcst_dates = [cube_time_converter(time, time_coord.units) 
                      for time in time_coord.points]
        
        area_bounds = self._area_inst.get_cube_area_bounds(cube, 
                                                           self.xy_coords)
        x_bounds = [area_bounds[self._area_inst.x_min], 
                    area_bounds[self._area_inst.x_max]]
        y_bounds = [area_bounds[self._area_inst.y_min], 
                    area_bounds[self._area_inst.y_max]]
        
        metadata = {}    
        metadata['VARIABLE'] = cube.name()
        metadata['UNITS'] = str(cube.units)
        metadata['INITIALISATION_DATES'] = init_dates
        metadata['MEMBERS'] = len(cube.coord(self.realization).points)
        metadata['FORECAST_DATES'] = fcst_dates
        metadata[self.xy_coords[0].upper()+'_BOUNDS'] = x_bounds
        metadata[self.xy_coords[-1].upper()+'_BOUNDS'] = y_bounds
        
        # Find additional coordinates in cube and add them to metadata.
        for coord in cube.coords():
            if coord.name() not in self.unwanted_coords and \
               coord.name() not in self._required_coords and \
               coord.name() not in self.xy_coords:
                metadata[coord.name().upper()] = coord.points
        
        return metadata

    def _get_obs_metadata(self, cube):
        """
        Retreive all relevant metadata.
        
        Args:
        
        * cube: iris cube
        
        Returns:
            dictionary
        
        """
        time_coord = cube.coord(self.time_coord)
        dates = [cube_time_converter(time, time_coord.units) 
                 for time in time_coord.points]
        
        area_bounds = self._area_inst.get_cube_area_bounds(cube, 
                                                           self.xy_coords)
        x_bounds = [area_bounds[self._area_inst.x_min], 
                    area_bounds[self._area_inst.x_max]]
        y_bounds = [area_bounds[self._area_inst.y_min], 
                    area_bounds[self._area_inst.y_max]]
        
        metadata = {}    
        metadata['VARIABLE'] = cube.name()
        metadata['UNITS'] = str(cube.units)
        metadata['DATES'] = dates
        metadata[self.xy_coords[0].upper()+'_BOUNDS'] = x_bounds
        metadata[self.xy_coords[-1].upper()+'_BOUNDS'] = y_bounds
        
        # Find additional coordinates in cube and add them to metadata.
        for coord in cube.coords():
            if coord.name() not in self.unwanted_coords and \
               coord.name() not in self._required_coords and \
               coord.name() not in self.xy_coords:
                metadata[coord.name().upper()] = coord.points
        
        return metadata
         
    def cube_ensemble_mean(self, cube):
        """
        Calculate the mean over all realizations.

        Args:
        
        * cube: iris cube
            Provide a cube to pass through the function. Default uses 
            self.cube.

        Returns:
            iris cube
            
        """
        success = False
        try:
            if len(cube.coord(self.realization).points) > 1 or \
               cube.coord(self.realization) in cube.coords(dim_coords=True):
                cube = cube.collapsed(self.realization, iris.analysis.MEAN)
            success = True
        except iris.exceptions.CoordinateNotFoundError:
            pass
        try:
            if len(cube.coord(self.forecast_ref_time).points) > 1  or \
               cube.coord(self.forecast_ref_time) in \
               cube.coords(dim_coords=True):
                cube = cube.collapsed(self.forecast_ref_time, 
                                      iris.analysis.MEAN)
            success = True
        except iris.exceptions.CoordinateNotFoundError:
            pass
        if not success:
            raise iris.exceptions.CoordinateNotFoundError(
                                  'No ensemble coordinates found.')
        return cube

    def cube_area_analysis(self, cube, method='MEAN'):
        """
        Calculate the spatial analysis.

        Args:
        
        * cube: iris cube
            Provide a cube to pass through the function.
        
        * method: string
            Specify the name of the iris.analysis method. Default is MEAN.

        Returns:
            iris cube
            
        """        
        try:
            if not cube.coord(self.xy_coords[0]).has_bounds():
                cube.coord(self.xy_coords[0]).guess_bounds()
                cube.coord(self.xy_coords[-1]).guess_bounds()
        except iris.exceptions.CoordinateNotFoundError:
            # If xycoords have been changed since load.
            self.xy_coords = get_xy_coords(cube)
            if not cube.coord(self.xy_coords[0]).has_bounds():
                cube.coord(self.xy_coords[0]).guess_bounds()
                cube.coord(self.xy_coords[-1]).guess_bounds()
        except ValueError:
            pass
        if 'longitude' in self.xy_coords[0] and \
           'latitude' in self.xy_coords[-1]:
            grid_areas = iris.analysis.cartography.area_weights(cube)
        else:
            grid_areas = None
        cube = cube.collapsed([self.xy_coords[0], self.xy_coords[-1]], 
                               getattr(iris.analysis, method), 
                               weights=grid_areas)
        return cube

    def cube_coordinate_analysis(self, cube, coordinate, method='MEAN'):
        """
        Calculate the analysis over coordinate points.
        
        Args:
        
        * cube: iris cube
            Provide a cube to pass through the function.
        
        * coordinate: string
            Name of coordinate.
            
        * method: string
            Specify the name of the iris.analysis method.

        Returns:
            iris cube
            
        """
        if len(cube.coord(coordinate).points) > 1 or \
           cube.coord(coordinate) in cube.coords(dim_coords=True):
            cube = cube.collapsed(coordinate, getattr(iris.analysis, method))
        return cube

    def advanced_search_pattern(self, prepend_search_string='', 
                                   append_search_string='',
                                   see_search_patterns=False):
        """
        Default is to use only the str_init_dates values to find
        files. This function can be used to prepend and append strings to 
        narrow the search.
        
        Args:
        
        * prepend_search_string: string
            String to prepend lag string.
        
        * append_search_string: string
            String to append lag string.
                
        """
        self.str_init_dates = [prepend_search_string + \
                               init_date + \
                               append_search_string
                               for init_date in self.str_init_dates]
        if see_search_patterns:
            print prepend_search_string + \
                  self.init_date_str_format + \
                  append_search_string

    def extract_area_bounds(self, area_bounds=None, cubes=None):
        """
        Extract area bounds. Note, instance data is already constrained with 
        instance area_bounds upon load so running this method without 
        specifying new data or bounds will result in no change.
        
        Args:

        * area_bounds: list of 4 floats
            Area bounds to extract, must be given in the format specified 
            during instantiation.

        * cubes: iris cube or cubelist
            Extract bounds (see area_bounds) from given cube(s). Default is 
            self.cube.
                    
        Returns:
            iris cube
            
        """
        if cubes is None:
            try:
                assert self.cube, ('Data has not been loaded, '\
                                   'use load method.')
                data = self.cube
            except AttributeError:
                assert self.cubelist, ('Data has not been loaded, '\
                                       'use load method.')
                data = self.cubelist
            if not area_bounds:
                print 'Area bounds have already been extracted on load.'
                return self.cube
            
            xy_coords = self.xy_coords
        else:
            if cubes:
                data = cubes
                if type(data) == iris.cube.CubeList:
                    ex_cube = data[0]
                else:
                    ex_cube = data
                xy_coords = [coord.name() for coord in get_xy_coords(ex_cube)]
            else:
                print 'No data in given cube(s).'
                return cubes

        if area_bounds:
            area_bounds = self._area_inst.check_area_bounds(area_bounds)
        else:
            area_bounds = self.area_bounds
                        
        area_constraint = self._area_inst.get_area_bound_constraints(
                                          area_bounds, 
                                          xy_coords)
        return data.extract(area_constraint)

    def extract_dates(self, dates=None, cubes=None, bounds=None, 
                        time_coord=None):
        """
        Extract a chosen set of dates from the data. Note, instance data is 
        already constrained with instance self.dates upon load so running this 
        method without specifying new data or dates will result in no change.
        
        Args:

        * dates: list of datetime objects
            Specify dates to extract. Default uses self.dates
            
        * cubes: iris cube or cubelist
            Extract dates from given cube(s). Default is self.cube.
                    
        * bounds: boolean
            If True, all dates between the min and max dates are extracted. 
            False means only specified dates are extracted. Default is None
            which uses the boolean set by self._dates_are_bounds.
            
        * time_coord: string
            The name of a time based coordinate. Default is self.time_coord.

        Returns:
            iris cube
            
        """
        if time_coord is None:
            time_coord = self.time_coord
        if cubes is None:
            try:
                assert self.cube, ('Data has not been loaded, use load '\
                                   'method.')
                single_cube = True
                data = iris.cube.CubeList([self.cube])
            except AttributeError:
                assert self.cubelist, ('Data has not been loaded, use load '\
                                   'method.')
                single_cube = False
                data = self.cubelist
            if dates is None:
                print 'Dates have already been extracted'
                return self.cube
            
            time_unit = self.time_unit
        else:
            if cubes:
                if type(cubes) == iris.cube.CubeList:
                    single_cube = False
                    data = cubes
                else:
                    single_cube = True
                    data = iris.cube.CubeList([cubes])
                time_unit = data[0].coord(time_coord).units
            else:
                print 'No data in given cube(s).'
                return cubes 
                
        if dates:
            check_dates(dates)
        else:
            dates = self.dates
        if bounds is None:
            dates_are_bounds = self._dates_are_bounds
        else:
            dates_are_bounds = bounds
        if dates_are_bounds:
            min_date = min(dates)
            min_date = datetime.datetime(min_date.year, min_date.month, 
                                         min_date.day, 0)
            max_date = max(dates)
            max_date = datetime.datetime(max_date.year, max_date.month, 
                                         max_date.day, 23, 59)
            time_constraint = get_time_bound_constraints(min_date, 
                                                         max_date,
                                                         time_coord)
        else:
            # Create daily bounds for the given dates.
            bound_dates = []
            for date in dates:
                bound_dates.append([datetime.datetime(date.year,
                                                      date.month,
                                                      date.day,
                                                      0),
                                    datetime.datetime(date.year,
                                                      date.month,
                                                      date.day,
                                                      23,
                                                      59)])            
            time_constraint = get_list_of_time_bound_constraints(bound_dates, 
                                                                 time_coord)
        extracted_cubes = data.extract(time_constraint)
        if not extracted_cubes:
            raise UserWarning('Data does not contain the given "%s" dates, %s' 
                              % (time_coord,
                                [date.strftime('%d/%m/%Y') for date in dates]))
            
        # Check if any of the requested dates were not in the data.
        # First compare the number of requested dates with the number of dates
        # extracted.
        extracted_dates = []
        for cube in extracted_cubes:
            for coord_point in cube.coord(time_coord).points:
                extracted_dates.append(coord_point)
        # If bounds is True, it is likely this will not be equal.
        if len(set(extracted_dates)) != len(dates):
            # If this is the case, check the dates manually.
            dates_not_in_data = []
            for date in dates:
                date_value = cube_time_converter(date, time_unit)
                found = False
                for cube in extracted_cubes:
                    # If time coordinate has bounds use these.
                    if cube.coord(time_coord).bounds is not None:
                        for bounds in cube.coord(time_coord).bounds:
                            if date_value >= min(bounds) and \
                               date_value <= max(bounds):
                                found = True
                                break
                    # Else use coordinate points.
                    else:
                        min_date = datetime.datetime(date.year, date.month, 
                                                     date.day, 0)
                        min_date_val = cube_time_converter(min_date, time_unit)
                        max_date = datetime.datetime(date.year, date.month, 
                                                     date.day, 23, 59)
                        max_date_val = cube_time_converter(max_date, time_unit)
                        for point in cube.coord(time_coord).points:
                            if min_date_val <= point <= max_date_val:
                                found = True
                                break
                if not found:
                    dates_not_in_data.append(date.strftime('%d/%m/%Y'))

            if dates_not_in_data:
                print 'Data does not contain the following "%s" dates, %s' \
                       % (time_coord, dates_not_in_data)
        
        if single_cube:
            extracted_cubes = extracted_cubes[0]
            
        return extracted_cubes

    def extract_coordinate_levels(self, coordinate, levels, cubes=None):
        """
        Extract level(s) from the given coordinate.
        
        Args:
        
        * coordinate: string
            Name of the coordinate from which to extract.
            
        * levels: integer, float or list of integers or floats
            Specify the level(s) to extract.
        
        * cubes: iris cube or cubelist
            cube(s) from which to extract level(s). Default is self.cube.
        
        Returns:
            iris cube
            
        """        
        if cubes is None:
            try:
                assert self.cube, ('Data has not been loaded, '\
                                   'use load method.')
                single_cube = True
                data = iris.cube.CubeList([self.cube])
            except AttributeError:
                assert self.cubelist, ('Data has not been loaded, '\
                                       'use load method.')
                single_cube = False
                data = self.cubelist
        else:
            if cubes:
                if type(cubes) == iris.cube.CubeList:
                    single_cube = False
                    data = cubes
                else:
                    single_cube = True
                    data = iris.cube.CubeList([cubes])
            else:
                print 'No data in given cube(s).'
                return cubes

        # Check coordinate is in all data.
        for cube in data:
            if coordinate not in [coord.name() for coord in cube.coords()]:
                raise UserWarning('There is no %s coordinate in the data.' 
                                  % coordinate)

        if type(levels) != list:
            levels = [levels]
        
        level_func = lambda l: l in levels
        level_constraint = iris.Constraint(**{coordinate : level_func})
        extracted_cubes = data.extract(level_constraint)
        if not extracted_cubes:
            valid_levels = []
            for cube in data:
                for coord_point in cube.coord(coordinate).points:
                    valid_levels.append(coord_point)
            raise UserWarning('%s has no specified level(s), %s. Valid '\
                              'level(s): %s' 
                              % (coordinate, levels,
                                 valid_levels))
        # Check if any of the requested levels were not in the data.
        extracted_levels = []
        for cube in extracted_cubes:
            for coord_point in cube.coord(coordinate).points:
                extracted_levels.append(coord_point)
        if len(set(extracted_levels)) != len(levels):
            coord_points = []
            for cube in extracted_cubes:
                for coord_point in cube.coord(coordinate).points:
                    coord_points.append(coord_point)
            levels_not_in_data = [level for level in levels
                                  if level not in coord_points]
            if levels_not_in_data:
                print 'The data did not contain the follows levels, %s' \
                       % levels_not_in_data
                   
        if single_cube:
            extracted_cubes = extracted_cubes[0]
            
        return extracted_cubes

    def configure_coord_names(self, member_name=default_member_name, 
                     period_name=default_period_name, 
                     initialistion_time_name=default_initialistion_time_name):
        """
        Parse the coordinate names and the values to which they are to be set.
        
        """
        self.realization = member_name
        self.time_coord = period_name
        self.forecast_ref_time = initialistion_time_name
                
        
    class MetaData(object):
        """
        Class for representing a metadata dictionary in a readable manner.
        
        Args:
        
        * meta_dict: dictionary
        
        * bound_names: list
            Specify the names of metadata refering to bounds.
        
        """
        def __init__(self, meta_dict, bound_names=None):
            assert type(meta_dict) == dict, 'MetaData requires a dictionary, '\
                                            'type %s given.' % type(meta_dict)
            self.meta_dict = meta_dict
            self.bound_names = bound_names
            self.keys = self.meta_dict.keys()
            self.meta_dict_rep = '--Metadata--\n\n---\n'
        
        def __str__(self):
            def add_line(key, remove=True):
                val = self.meta_dict.get(key)
                if val is not None:
                    self.meta_dict_rep += '{key}: {val}\n---\n'.format(
                                                                  key=key, 
                                                                  val=val)
                    if remove:
                        self.keys.remove(key)
            
            # Priority order.
            add_line('DATA_TYPE')
            add_line('VARIABLE')
            add_line('UNITS')
            add_line('INITIALISATION_DATES')
            add_line('FORECAST_DATES')
            add_line('DATES')
            add_line('YEARS')
            add_line('MEMBERS')
            add_line('TOTAL_MEMBERS')
            if self.bound_names:
                for bound_name in self.bound_names:
                    add_line(bound_name, remove=False)
                
            for key in self.keys:
                if key not in self.bound_names:
                    add_line(key, remove=False)
            
            return self.meta_dict_rep


class CubeData(_DataHandler):
    """
    Class for iris Cubes.
    
    """
    def __init__(self, cube, 
                  area_bounds_format=default_area_bounds_format, 
                  area_bounds_range=default_area_bounds_range,
                  unwanted_coords=[]):
        self.configure_coord_names()
        self.cube = cube
        self.xy_coords = [coord.name() 
                          for coord in get_xy_coords(self.cube)]
        self.time_unit = self.cube.coord(self.time_coord).units
        self._area_inst = AreaBounds(area_bounds_format, 
                                    area_bounds_range)
        self._required_coords = []
        self.unwanted_coords = []
        self.processes = []
        
    def time_analysis(self, method='MEAN'):
        """
        Calculate the analysis over time.

        Args:

        * method: string
            Specify the name of the iris.analysis method.

        Returns:
            iris cube
        
        """
        self.cube = self.cube_coordinate_analysis(self.cube, self.time_coord, 
                                                  method)
        self.processes.append('time_analysis')
        return self.cube

    def ensemble_mean(self):
        """
        Calculate the mean over all realizations.

        Returns:
            iris cube
            
        """
        self.cube = self.cube_ensemble_mean(self.cube)
        self.processes.append('ensemble_mean')
        return self.cube

    def area_analysis(self, method='MEAN'):
        """
        Calculate the spatial analysis.
        
        Args:

        * method: string
            Specify the name of the iris.analysis method.

        Returns:
            iris cube
            
        """
        self.cube = self.cube_area_analysis(self.cube, method)
        self.processes.append('area_analysis')
        return self.cube

    def coordinate_analysis(self, coordinate, method='MEAN'):
        """
        Calculate the analysis over coordinate points.
        
        Args:
        
        * coordinate: string
            Name of coordinate to be meaned.

        * method: string
            Specify the name of the iris.analysis method.

        Returns:
            iris cube
        
        """
        self.cube = self.cube_coordinate_analysis(self.cube, coordinate, method)
        self.processes.append('%s_analysis' % coordinate)
        return self.cube

class CubeListData(_DataHandler):
    """
    Class for iris CubeLists
    
    """
    def __init__(self, cubelist, 
                  area_bounds_format=default_area_bounds_format, 
                  area_bounds_range=default_area_bounds_range,
                  unwanted_coords=[]):
        self.configure_coord_names()
        self.cubelist = cubelist
        self.xy_coords = [coord.name() 
                          for coord in get_xy_coords(self.cubelist[0])]
        self.time_unit = self.cubelist[0].coord(self.time_coord).units
        self._area_inst = AreaBounds(area_bounds_format, 
                                    area_bounds_range)
        self._required_coords = []
        self.unwanted_coords = []
        self.processes = []

    def time_analysis(self, method='MEAN'):
        """
        Calculate the mean over time.

        * method: string
            Specify the name of the iris.analysis method.

        Returns:
            iris cubelist
        
        """
        new_cubelist = []
        for cube in self.cubelist:
            new_cubelist.append(self.cube_coordinate_analysis(cube, 
                                                              self.time_coord, 
                                                              method))
        self.cubelist = iris.cube.CubeList(new_cubelist)
        self.processes.append('time_analysis')
        return self.cubelist
    
    def ensemble_mean(self):
        """
        Calculate the mean over all realizations.
        
        Returns:
            iris cubelist
            
        """
        new_cubelist = []
        for cube in self.cubelist:
            new_cubelist.append(self.cube_ensemble_mean(cube))
        self.cubelist = iris.cube.CubeList(new_cubelist)
        self.processes.append('ensemble_mean')
        return self.cubelist

    def area_analysis(self, method='MEAN'):
        """
        Calculate the spatial mean.

        * method: string
            Specify the name of the iris.analysis method.

        Returns:
            iris cubelist
            
        """
        new_cubelist = []
        for cube in self.cubelist:
            new_cubelist.append(self.cube_area_analysis(cube, method))
        self.cubelist = iris.cube.CubeList(new_cubelist)
        self.processes.append('area_analysis')
        return self.cubelist

    def coordinate_analysis(self, coordinate, method='MEAN'):
        """
        Calculate the analysis over coordinate points.
        
        Args:
        
        * coordinate: string
            Name of coordinate to be meaned.

        * method: string
            Specify the name of the iris.analysis method.

        Returns:
            iris cubelist
        
        """
        new_cubelist = []
        for cube in self.cubelist:
            new_cubelist.append(self.cube_coordinate_analysis(cube, coordinate, method))
        self.cubelist = iris.cube.CubeList(new_cubelist)
        self.processes.append('%s_analysis' % coordinate)
        return self.cubelist

class ForecastData(CubeData):
    """
    Class for loading and handling forecast data.
    
    Args:
    
    * directory: string
        Specify the directory of the forecast data to be loaded. Additional sub
        directories can given in the load method if required.
        
    * variable: string
        Specify the meteorological variable name. This must be the name as 
        given in the data, e.g. precipitation_flux.
    
    * init_dates: list of datetime objects
        Specify the initialisation dates of the forecast runs. The default 
        loading technique is to find filenames which contain string versions of
        these dates.
        
    * fcst_dates: list of datetime objects
        Specify the forecast dates to be loaded. If dates_are_bounds (see 
        below) is True, all dates in between the max and min dates are loaded.
    
    * area_bounds: list of 4 floats or integers
        Specify the upper and lower bounds in the x and y direction. See 
        area_bounds_format for formatting. Default is [-180,-90,180,90].
        
    * dates_are_bounds: boolean
        State whether the fcst_dates given are a list of dates or a set of 
        bounds. Default is False.
    
    * area_bounds_format: list of 4 strings
        Provide a list with 'x_min','y_min','x_max' and 'y_max' in the 
        required order to specify the format. Default, ['x_min','y_min',
        'x_max','y_max'].
    
    * area_bounds_range: list of 4 floats or integers
        Specify the range for the x and y coordinates to fall into. Values must
        represent the full cycle of each coordinate system, e.g. 0 to 360 or 
        -180 to 180. This argument is used for shifting coordinate systems, NOT
        for extracting bounds. Values must be given in the same format 
        specified by area_bounds_format. Default is None which uses x range 
        -180 to 180 and y range -90 to 90 and automatically adjusts to the 
        given format.
        
    * init_date_str_format: string
        The default loading technique is to find filenames which contain string
        versions of the initialisation dates. The string format used is 
        specified here, using datetime string formating. Default is YYMMDD.
        
    * unwanted_coords: list of strings
        Specify any coordinate names that are to be removed. Note, removing or
        not removing coordinates can effect how data is merged together. If 
        coordinate name does not exist it is ignored and no exception is 
        raised. Default is to remove forecast_period and member_number.
    
    * multiprocess_workers: integer
        How many processes to run when loading data.
            
    """
    def __init__(self, directory, variable, init_dates, fcst_dates,
                  area_bounds=default_area_bounds, 
                  dates_are_bounds=default_dates_are_bounds,
                  area_bounds_format=default_area_bounds_format, 
                  area_bounds_range=default_area_bounds_range, 
                  init_date_str_format=default_init_date_str_format,
                  unwanted_coords=['forecast_period', 'member_number'],
                  multiprocess_workers=1):
        
        self.configure_coord_names()
        self._required_coords = [self.realization, self.time_coord, 
                                self.forecast_ref_time]
        self.directory   = directory
        self.variable    = variable
        # Create AreaBounds instance for use of methods later.
        self._area_inst   = AreaBounds(area_bounds_format, 
                                       area_bounds_range)
        self.area_bounds  = self._area_inst.check_area_bounds(area_bounds)
        self._dates_are_bounds = dates_are_bounds
        self.init_dates, self.dates = self.period_check(init_dates, 
                                                        fcst_dates)
        # Get init dates in string format to reference filenames.
        self.init_date_str_format = init_date_str_format
        self.str_init_dates = [date.strftime(init_date_str_format) 
                               for date in self.init_dates]
        if set(unwanted_coords) & set(self._required_coords):
            diff = set(unwanted_coords) & set(self._required_coords)
            raise ValueError('The coordinate(s) %s are required and can not '\
                             'be removed.' % list(diff))
        self.unwanted_coords = unwanted_coords
        self.xy_coords       = None
        self.processes       = []
        self.cube            = None
        self._callback       = remove_attributes
        self.process_workers = multiprocess_workers
        
    @staticmethod
    def period_check(init_dates, fcst_dates):
        """
        Check that the initialisation dates and forecast dates are valid and do
        not overlap.
        
        Args:
        
        * init_dates: list of datetime objects
            List of initialisation dates.
        
        * fcst_dates: list of datetime objects
            List of forecast dates.
        
        Returns
            Tuple containing two valid date lists.
        
        """
        check_dates(init_dates)
        check_dates(fcst_dates)
        
        if max(init_dates) > min(fcst_dates):
            raise ValueError('Forecast date, %s comes before initialisation '\
                             'date, %s.' % (min(fcst_dates), max(init_dates)))

        fcst_dates = change_zeroth_hour(fcst_dates)
        return init_dates, fcst_dates

    def _get_metadata(self):
        """
        Retrieve the metadata from the cube.
        
        Returns:
            class MetaData
        
        """    
        metadata = {'DATA_TYPE':'Forecast Data'}    
                
        cube_metadata = self._get_model_metadata(self.cube)
        
        self.cube_init_dates = cube_metadata['INITIALISATION_DATES']
        del cube_metadata['INITIALISATION_DATES']
        
        self.cube_dates = cube_metadata['FORECAST_DATES']
        del cube_metadata['FORECAST_DATES']
        
        for key, val in cube_metadata.items():
            # Find unique metadata which has not already been added by 
            # previous cubes. Years are the common one.
            current_vals = metadata.get(key)
            if current_vals is not None:
                for this_val in current_vals:
                    if hasattr(this_val, '__iter__'):     
                        try:    
                            if numpy.array_equal(this_val, val):
                                break
                        except AttributeError:
                            # If the array type is not comparable for 
                            # example array of strings.
                            equal = True
                            for this_item, item in zip(this_val, val):
                                if this_item != item:
                                    equal = False
                                    break
                            if equal:
                                break
                    else:
                        if this_val == val:
                            break
                    metadata[key].append(val)
            else:
                metadata[key] = [val]
        
        bound_names = []
        # Tidy up list of length 1.
        for key, val in metadata.items():
            if type(val) == list and len(val) == 1:
                metadata[key] = val[0]
            # Retrieve the exact bound names.
            if key[-7:] == '_BOUNDS':
                bound_names.append(key)
                
        metadata['INITIALISATION_DATES'] = [date.strftime('%d/%m/%Y') 
                                            for date in 
                                            self.cube_init_dates]
        metadata['FORECAST_DATES'] = [date.strftime('%d/%m/%Y') 
                                      for date in self.cube_dates]

        return self.MetaData(metadata, bound_names)

    def load(self, files_to_search_for=None, dirs_to_look_in='', 
             current_dir=False, see_files=False, raw_data=False):
        """
        Load data from file.
        
        Args:

        * files_to_search_for: list of strings or string
            The default loading technique is to find filenames which contain 
            string versions of the initialisation dates. If this is not 
            appropriate, specify here the string matches to search for in the 
            filenames. Use see_files to view the search results before loading.

        * dirs_to_look_in: list of strings or string
            Specify any sub directories of the main directory (given during the 
            instantiation of the class) of which to search in. Default is to
            search all files and directories beneath the main directory.
        
        * current_dir: boolean
            Only search the given directory for matching files.
            
        * see_files: boolean
            Set to True to return search results without loading the files. 
            Default is False.
            
        * raw_data: boolean
            Once the data is loaded it is sorted and constrained to the 
            parameters given during instantiation. If the raw data i.e. all the
            data within the loaded files is required set to True. Default is 
            False.

        Returns:
            iris cube
            
        """
        if not files_to_search_for:
            files_to_search_for = self.str_init_dates
        files_to_load = self._find_files(self.directory, dirs_to_look_in, 
                                         files_to_search_for, current_dir, 
                                         see_files)

        if see_files:
            print '\n'.join(files_to_load)
            return None
        if raw_data:
            return iris.load(files_to_load, self.variable)
        all_cubes = self._load_all_cubes(set(files_to_load))
        all_cubes = self.extract_dates(self.dates, all_cubes)
        all_cubes = self.extract_dates(self.init_dates, all_cubes, 
                                       bounds=False, 
                                       time_coord=self.forecast_ref_time)
        self.cube = self._sort_cubelist(all_cubes)
        self.metadata = self._get_metadata()
        return self.cube
    

class HindcastData(CubeListData):
    """
    Class for loading and handling hindcast data.
    
    Args:

    * directory: string
        Specify the directory of the hindcast data to be loaded. Additional sub
        directories can given in the load method if required.
    
    * variable: string
        Specify the meteorological variable name. This must be the name as 
        given in the data, e.g. precipitation_flux.
    
    * init_dates: list of datetime objects
        Specify the initialisation dates of the hindcast runs. The default 
        loading technique is to find filenames which contain string versions of
        these dates.
        
    * fcst_dates: list of datetime objects
        Specify the forecast dates to be loaded. If dates_are_bounds (see 
        below) is True, all dates in between the min and max dates are loaded.
        Note, the year within the dates is ignored (years are given by below
        argument) so only specify a single period e.g. 5/3/2014 - 10/4/2014, 
        and this will be replicated for all given years (meaning 2014 is not
        necessarily loaded unless specified in years argument).
        
    * years: list of integers
        Specify the hindcast years to be loaded. If years_bounds (see below) is
        True, all years in between the max and min years are loaded.
    
    * area_bounds: list of 4 floats or integers
        Specify the upper and lower bounds in the x and y direction. See 
        area_bounds_format for formatting. Default is [-180,-90,180,90].
    
    * dates_are_bounds: boolean
        State whether the fcst_dates given are a list of dates or a set of 
        bounds. Default is False.
    
    * years_bounds: boolean
        State whether the years given are a list of years or a set of bounds. 
        Default is False.
    
    * area_bounds_format: list of 4 strings
        Provide a list with 'x_min','y_min','x_max' and 'y_max' in the 
        required order to specify the format. Default, ['x_min','y_min',
        'x_max','y_max'].
    
    * area_bounds_range: list of 4 floats or integers
        Specify the range for the x and y coordinates to fall into. Values must
        represent the full cycle of each coordinate system, e.g. 0 to 360 or 
        -180 to 180. This argument is used for shifting coordinate systems, NOT
        for extracting bounds. Values must be given in the same format 
        specified by area_bounds_format. Default is None which uses x range 
        -180 to 180 and y range -90 to 90 and automatically adjusts to the 
        given format.
        
    * init_date_str_format: string
        The default loading technique is to find filenames which contain string
        versions of the initialisation dates. The string format used is 
        specified here, using datetime string formating. Default is YYMMDD.
        
    * unwanted_coords: list of strings
        Specify any coordinate names that are to be removed. Note, removing or
        not removing coordinates can effect how data is merged together. If 
        coordinate name does not exist it is ignored and no exception is 
        raised. Default is to remove forecast_period and member_number.

    * multiprocess_workers: integer
        How many processes to run when loading data.
        
    """
    def __init__(self, directory, variable, init_dates, fcst_dates, years,
                  area_bounds=default_area_bounds, 
                  dates_are_bounds=default_dates_are_bounds, 
                  years_are_bounds=default_years_are_bounds,
                  area_bounds_format=default_area_bounds_format, 
                  area_bounds_range=default_area_bounds_range, 
                  init_date_str_format=default_init_date_str_format,
                  unwanted_coords=['forecast_period', 'member_number'],
                  multiprocess_workers=1):

        self.configure_coord_names()
        self._required_coords = [self.realization, self.time_coord, 
                                self.forecast_ref_time]
        self.directory   = directory
        self.variable    = variable
        # Create AreaBounds instance for use of methods later.
        self._area_inst  = AreaBounds(area_bounds_format, 
                                      area_bounds_range)
        self.area_bounds = self._area_inst.check_area_bounds(area_bounds)
        self._dates_are_bounds = dates_are_bounds
        if years_are_bounds:
            self.years = range(min(years), max(years)+1)
        else:
            self.years = years
        self.init_dates, self.dates = self._sort_hcst_dates(init_dates,
                                                            fcst_dates,
                                                            self.years)
        # Get init dates in string format to reference filenames.
        self.init_date_str_format = init_date_str_format
        self.str_init_dates = [date.strftime(init_date_str_format) 
                               for date in self.init_dates]
        if set(unwanted_coords) & set(self._required_coords):
            diff = set(unwanted_coords) & set(self._required_coords)
            raise ValueError('The coordinate(s) %s are required and can not '\
                             'be removed.' % list(diff))
        self.unwanted_coords = unwanted_coords
        self.xy_coords       = None
        self.processes       = []
        self.cubelist        = None
        self._callback       = remove_attributes
        self.process_workers = multiprocess_workers

    def _sort_hcst_dates(self, init_dates, fcst_dates, years):
        """
        Sort the given dates into yearly lists of initialisation and forecast
        dates.
        
        Args:
        
        * init_dates: list of datetime objects
            List of initialisation dates.
        
        * fcst_dates: list of datetime objects
            List of forecast dates.
        
        * years: list of integers
        
        Returns
            Tuple containing two valid date lists.
        
        """
        check_dates(init_dates)
        check_dates(fcst_dates)

        fcst_dates = change_zeroth_hour(fcst_dates)

        sorted_init_dates = []
        sorted_fcst_dates = []
        # Get the lowest year of the forecast dates so it can be seen when 
        # years cross, e.g. Dec 2000, Jan 2001.
        fcst_start_year = min(fcst_dates).year
        for year in years:
            for init_date in init_dates:
                extra_years = init_date.year - fcst_start_year
                sorted_init_dates.append(datetime.datetime(year + extra_years, 
                                                           init_date.month,
                                                           init_date.day,
                                                           init_date.hour))
            # Create a separate list for each year.
            if self._dates_are_bounds:
                extra_years = max(fcst_dates).year - fcst_start_year
                yearly_fcst_dates = [datetime.datetime(year, 
                                                       min(fcst_dates).month,
                                                       min(fcst_dates).day,
                                                       min(fcst_dates).hour),
                                     datetime.datetime(year + extra_years, 
                                                       max(fcst_dates).month,
                                                       max(fcst_dates).day,
                                                       max(fcst_dates).hour)]
            else:
                yearly_fcst_dates = []
                for fcst_date in fcst_dates:
                    extra_years = fcst_date.year - fcst_start_year
                    yearly_fcst_dates.append(datetime.datetime(year + \
                                                               extra_years,
                                                               fcst_date.month,
                                                               fcst_date.day,
                                                               fcst_date.hour))
            sorted_fcst_dates.append(yearly_fcst_dates)
        return sorted_init_dates, sorted_fcst_dates

    def _sort_data(self, cubelist):
        """
        Sort all data in cubelist into yearly cubes.
        
        Args:
        
        * cubelist: iris cubelist
            cubelist of all data to be sorted into yearly cubes.
        
        Returns:
            iris cubelist
            
        """
        sorted_cubelist = []
        for dates in self.dates:
            year_cubelist = self.extract_dates(dates, cubelist)
            sorted_cubelist.append(self._sort_cubelist(year_cubelist))
        return iris.cube.CubeList(sorted_cubelist)

    def _get_metadata(self):
        """
        Retrieve the metadata from the cube.
        
        Returns:
            class MetaData
        
        """    
        def add_dates(date_list, dates):
            """
            Append dates to date_list which are not already within date_list.
                        
            """
            for date in dates:
                if date.strftime('%d-%b') not in date_list:
                    date_list.append(date.strftime('%d-%b'))
            return date_list
            
        metadata = {'DATA_TYPE':'Hindcast Data'}    
        
        members = 0
        self.cube_init_dates = []
        self.cube_dates = []
        years = []
        
        for cube in self.cubelist:
            cube_metadata = self._get_model_metadata(cube)
            members += cube_metadata['MEMBERS']
            del cube_metadata['MEMBERS']
            
            self.cube_init_dates = add_dates(self.cube_init_dates, 
                                   cube_metadata['INITIALISATION_DATES'])
            del cube_metadata['INITIALISATION_DATES']
            
            self.cube_dates = add_dates(self.cube_dates, 
                                             cube_metadata['FORECAST_DATES'])
            # Years are based on the earliest forecast date.
            years.append(min(cube_metadata['FORECAST_DATES']).year)
            del cube_metadata['FORECAST_DATES']
            
            for key, val in cube_metadata.items():
                # Find unique metadata which has not already been added by 
                # previous cubes. Years are the common one.
                current_vals = metadata.get(key)
                if current_vals is not None:
                    for this_val in current_vals:
                        if hasattr(this_val, '__iter__'):     
                            try:    
                                if numpy.array_equal(this_val, val):
                                    break
                            except AttributeError:
                                # If the array type is not comparable for 
                                # example array of strings.
                                equal = True
                                for this_item, item in zip(this_val, val):
                                    if this_item != item:
                                        equal = False
                                        break
                                if equal:
                                    break
                        else:
                            if this_val == val:
                                break
                        metadata[key].append(val)
                else:
                    metadata[key] = [val]
        
        bound_names = []
        # Tidy up lists of length 1.
        for key, val in metadata.items():
            if type(val) == list and len(val) == 1:
                metadata[key] = val[0]
            # Retrieve the exact bound names.
            if key[-7:] == '_BOUNDS':
                bound_names.append(key)
                
                
        metadata['INITIALISATION_DATES'] = self.cube_init_dates
        metadata['YEARS'] = sorted(list(set(years)))
        metadata['TOTAL_MEMBERS'] = members
        metadata['FORECAST_DATES'] = self.cube_dates

        return self.MetaData(metadata, bound_names)

    def load(self, files_to_search_for=None, dirs_to_look_in='', 
             current_dir=False, see_files=False, raw_data=False):
        """
        Load data from file.
        
        Args:

        * files_to_search_for: list of strings or string
            The default loading technique is to find filenames which contain 
            string versions of the initialisation dates. If this is not 
            appropriate, specify here the string matches to search for in the 
            filenames. Use see_files to view the search results before loading.

        * dirs_to_look_in: list of strings or string
            Specify any sub directories of the main directory (given during the 
            instantiation of the class) of which to search in. Default is to
            search all files and directories beneath the main directory.
        
        * current_dir: boolean
            Only search the given directory for matching files.
        
        * see_files: boolean
            Set to True to return search results without loading the files. 
            Default is False.
            
        * raw_data: boolean
            Once the data is loaded it is sorted and constrained to the 
            parameters given during instantiation. If the raw data i.e. all the
            data within the loaded files is required set to True. Default is 
            False.

        Returns:
            iris cubelist
            
        """
        if not files_to_search_for:
            files_to_search_for = self.str_init_dates
        files_to_load = self._find_files(self.directory, dirs_to_look_in, 
                                         files_to_search_for, current_dir,
                                         see_files)
        if see_files:
            print '\n'.join(files_to_load)
            return files_to_load
        if raw_data:
            return iris.load(files_to_load, self.variable)
        all_cubes = self._load_all_cubes(set(files_to_load))
        all_cubes = self.extract_dates(self.init_dates, all_cubes, 
                                       bounds=False, 
                                       time_coord=self.forecast_ref_time)
        self.cubelist = self._sort_data(all_cubes)
        self.metadata = self._get_metadata()
        return self.cubelist

    def create_climatology(self):
        """
        Mean all realizations from all years together to give a model 
        climatology. This creates a self.clim_cube attribute.
        
        Returns:
            iris cube
        """
        # Because data from all years are merged, the time coordinate must be 
        # made consistent with basic integer values. Monthly dates values are 
        # added to the attributes.
        time_points     = self.cubelist[0].coord(self.time_coord).points
        new_time_points = range(1, len(time_points) + 1)
        new_time_atts   = {'dates' : self.cube_dates}
        new_time_coord  = iris.coords.DimCoord(new_time_points,
                                               standard_name=self.time_coord,
                                               attributes=new_time_atts)
                
        new_cubelist = []
        realization_num = 1
        for cube in self.cubelist:
            if len(cube.coord(self.realization).points) > 1:
                cube = self.cube_ensemble_mean(cube)
            # Make sure all realization points are unique.
            cube.coord(self.realization).points = [realization_num]
            # Replace time dimension.
            time_dim = cube.coord_dims(cube.coord(self.time_coord))
            cube.remove_coord(self.time_coord)
            if time_dim:
                cube.add_dim_coord(new_time_coord, time_dim)
            else:
                # If no time_dim, coordinate is auxiliary or scalar.
                cube.add_aux_coord(new_time_coord)
                        
            new_cubelist.append(cube)
            realization_num += 1
            
        new_cube  = iris.cube.CubeList(new_cubelist).merge_cube()
        clim_cube = self.cube_ensemble_mean(new_cube)
        
        # The initialisation data is now a mean of all years, so like with time
        # replace the coordinate with a single point and monthly initialisation
        # dates added to the attributes.
        init_points     = clim_cube.coord(self.forecast_ref_time).points
        new_init_points = range(1, len(init_points) + 1)
        new_init_atts   = {'dates' : self.cube_init_dates}
        new_init_coord  = iris.coords.AuxCoord(new_init_points,
                                      standard_name=self.forecast_ref_time,
                                      attributes=new_init_atts)
        clim_cube.remove_coord(self.forecast_ref_time)
        clim_cube.add_aux_coord(new_init_coord)
        
        self.clim_cube = clim_cube
        return self.clim_cube
            

class ObservationData(CubeListData):
    """
    Class for loading and handling observation data.
    
     Args:

    * directory: string
        Specify the directory of the observation data to be loaded. Additional 
        sub directories can given in the load method if required.

    * variable: string
        Specify the meteorological variable name. This must be the name as 
        given in the data, e.g. precipitation.
            
    * dates: list of datetime objects
        Specify the dates to be loaded. If dates_are_bounds (see below) is 
        True, all dates in between the min and max dates are loaded.
        Note, If years are given (see below), the year of the dates is ignored. 
        In this case only specify a single period e.g. 5/3/2014 - 10/4/2014, 
        and this will be replicated for all given years (meaning 2014 is not
        necessarily loaded unless specified in years argument). If years are 
        not given then the exact dates are loaded.
        
    * years: list of integers
        Specify the years to be loaded. If years_bounds (see below) is True, 
        all years in between the max and min years are loaded. Default is None
        meaning only the dates given (see above) are loaded.
    
    * area_bounds: list of 4 floats or integers
        Specify the upper and lower bounds in the x and y direction. See 
        area_bounds_format for formatting. Default is None.
        
    * dates_are_bounds: boolean
        State whether the dates given are a list of dates or a set of 
        bounds. Default is False.
    
    * years_are_bounds: boolean
        State whether the years given are a list of years or a set of bounds. 
        Default is False.
    
    * area_bounds_format: list of 4 strings
        Provide a list with 'x_min','y_min','x_max' and 'y_max' in the 
        required order to specify the format. Default, ['x_min','y_min',
        'x_max','y_max'].
    
    * area_bounds_range: list of 4 floats or integers
        Specify the range for the x and y coordinates to fall into. Values must
        represent the full cycle of each coordinate system, e.g. 0 to 360 or 
        -180 to 180. This argument is used for shifting coordinate systems, NOT
        for extracting bounds. Values must be given in the same format 
        specified by area_bounds_format. Default is None which uses x range 
        -180 to 180 and y range -90 to 90 and automatically adjusts to the 
        given format.
                
    * unwanted_coords: list of strings
        Specify any coordinate names that are to be removed. Note, removing or
        not removing coordinates can effect how data is merged together. If 
        coordinate name does not exist it is ignored and no exception is 
        raised.
    
    * multiprocess_workers: integer
        How many processes to run when loading data.
        
    """
    def __init__(self, directory, variable, dates, years=None,
                  area_bounds=default_area_bounds, 
                  dates_are_bounds=default_dates_are_bounds,
                  years_are_bounds=default_years_are_bounds,
                  area_bounds_format=default_area_bounds_format, 
                  area_bounds_range=default_area_bounds_range, 
                  unwanted_coords=[], multiprocess_workers=1):
        
        self.configure_coord_names()
        self._required_coords = [self.time_coord]
        self.directory   = directory
        self.variable    = variable
        # Create AreaBounds instance for use of methods later.
        self._area_inst  = AreaBounds(area_bounds_format, 
                                      area_bounds_range)
        self.area_bounds = self._area_inst.check_area_bounds(area_bounds)
        self._dates_are_bounds = dates_are_bounds
        self._years_are_bounds = years_are_bounds
        self.dates = self._sort_dates(dates, years)
        if set(unwanted_coords) & set(self._required_coords):
            diff = set(unwanted_coords) & set(self._required_coords)
            raise ValueError('The coordinate(s) %s are required and can not '\
                             'be removed.' % list(diff))
        self.unwanted_coords = unwanted_coords
        self.xy_coords       = None
        self.processes       = []
        self.cubelist        = None
        self._callback       = None
        self.process_workers = multiprocess_workers

    def _sort_dates(self, dates, years):
        """
        Sort the given dates into yearly lists.
        
        Args:
        
        * dates: list of datetime objects
        
        * years: list of integers
        
        Returns
            list
            
        """
        check_dates(dates)
        dates = change_zeroth_hour(dates)
        
        sorted_dates = []
        if years is None:
            sorted_dates.append(dates)
        
        else:
            if self._years_are_bounds:
                years = range(min(years), max(years)+1)
            # Get the lowest year of the dates so it can be seen when 
            # years cross, e.g. Dec 2000, Jan 2001.
            start_year = min(dates).year
            for year in years:
                if self._dates_are_bounds:
                    extra_years = max(dates).year - start_year
                    yearly_dates = [datetime.datetime(year, 
                                                      min(dates).month,
                                                      min(dates).day,
                                                      min(dates).hour),
                                    datetime.datetime(year + extra_years, 
                                                      max(dates).month,
                                                      max(dates).day,
                                                      max(dates).hour)]
                else:
                    yearly_dates = []
                    for date in dates:
                        extra_years = date.year - start_year
                        yearly_dates.append(datetime.datetime(year + \
                                                              extra_years,
                                                              date.month,
                                                              date.day,
                                                              date.hour))
                sorted_dates.append(yearly_dates)
        return sorted_dates

    def _sort_data(self, cubelist):
        """
        Sort all data in cubelist into yearly cubes.
        
        Args:
        
        * cubelist: iris cubelist
            cubelist of all data to be sorted into yearly cubes.
        
        Returns:
            iris cubelist
            
        """
        sorted_cubelist = []
        for dates in self.dates:
            year_cubelist = self.extract_dates(dates, cubelist)
            for cube in year_cubelist.merge():
                # Check x-y coordinates match the specified range.
                cube = self._area_inst.check_cube_area_bounds(cube, 
                                                             self.xy_coords, 
                                                             self.area_bounds)
                cube = self.extract_area_bounds(cubes=cube)
                sorted_cubelist.append(cube)
        return iris.cube.CubeList(sorted_cubelist)

    def _get_metadata(self):
        """
        Retrieve the metadata from the cube.
        
        Returns:
            class MetaData
        
        """    
        def add_dates(date_list, dates):
            """
            Append dates to date_list which are not already within date_list.
                        
            """
            for date in dates:
                if date.strftime('%d-%b') not in date_list:
                    date_list.append(date.strftime('%d-%b'))
            return date_list
            
        metadata = {'DATA_TYPE':'Observation Data'}    
        
        self.cube_dates = []
        years = []
        
        for cube in self.cubelist:
            cube_metadata = self._get_obs_metadata(cube)
                        
            self.cube_dates = add_dates(self.cube_dates, 
                                        cube_metadata['DATES'])
            # Years are based on the earliest date.
            years.append(min(cube_metadata['DATES']).year)
            del cube_metadata['DATES']
            
            for key, val in cube_metadata.items():
                # Find unique metadata which has not already been added by 
                # previous cubes. Years are the common one.
                current_vals = metadata.get(key)
                if current_vals is not None:
                    for this_val in current_vals:
                        if hasattr(this_val, '__iter__'):
                            try:    
                                if numpy.array_equal(this_val, val):
                                    break
                            except AttributeError:
                                # If the array type is not comparable for 
                                # example array of strings.
                                equal = True
                                for this_item, item in zip(this_val, val):
                                    if this_item != item:
                                        equal = False
                                        break
                                if equal:
                                    break
                        else:
                            if this_val == val:
                                break
                        metadata[key].append(val)
                else:
                    metadata[key] = [val]
        
        bound_names = []
        # Tidy up lists of length 1.
        for key, val in metadata.items():
            if type(val) == list and len(val) == 1:
                metadata[key] = val[0]
            # Retrieve the exact bound names.
            if key[-7:] == '_BOUNDS':
                bound_names.append(key)
                
        metadata['YEARS'] = sorted(list(set(years)))
        metadata['DATES'] = self.cube_dates
        
        return self.MetaData(metadata, bound_names)

    def load(self, files_to_search_for, dirs_to_look_in='',
             current_dir=False, see_files=False, raw_data=False):
        """
        Load data from file.
        
        Args:

        * files_to_search_for: list of strings or string
            Specify here the string matches to search for in the filenames. Use
            see_files to view the search results before loading.

        * dirs_to_look_in: list of strings or string
            Specify any sub directories of the main directory (given during the 
            instantiation of the class) of which to search in. Default is to
            search all files and directories beneath the main directory.
        
        * current_dir: boolean
            Only search the given directory for matching files.
            
        * see_files: boolean
            Set to True to return search results without loading the files. 
            Default is False.
            
        * raw_data: boolean
            Once the data is loaded it is sorted and constrained to the 
            parameters given during instantiation. If the raw data i.e. all the
            data within the loaded files is required set to True. Default is 
            False.

        Returns:
            iris cube
            
        """
        files_to_load = self._find_files(self.directory, dirs_to_look_in, 
                                         files_to_search_for, current_dir,
                                         see_files)

        if see_files:
            print '\n'.join(files_to_load)
            return None
        if raw_data:
            return iris.load(files_to_load, self.variable)
        all_cubes = self._load_all_cubes(set(files_to_load))
        self.cubelist = self._sort_data(all_cubes)
        self.metadata = self._get_metadata()
        
        return self.cubelist


    def create_climatology(self):
        """
        Mean all data from all years together to give a climatology.
        
        Returns:
            iris cube
        """
        # Because data from all years are merged, the time coordinate must be 
        # made consistent with basic integer values. Monthly dates values are 
        # added to the attributes.
        time_points     = self.cubelist[0].coord(self.time_coord).points
        new_time_points = range(1, len(time_points) + 1)
        new_time_atts   = {'dates' : self.cube_dates}
        new_time_coord  = iris.coords.DimCoord(new_time_points,
                                               standard_name=self.time_coord,
                                               attributes=new_time_atts)
        
        # Create a years coordinate so they can be meaned.
        years_coord = iris.coords.AuxCoord(1, long_name='years')
                
        new_cubelist = []
        for this_cube in self.cubelist:
            cube = this_cube.copy()
            # Add year of earliest date.
            cube.add_aux_coord(years_coord)
            cube.coord('years').points = cube_time_converter(
                                         cube.coord('time').points[0], 
                                         self.time_unit).year
            # Replace time dimension.
            time_dim = cube.coord_dims(cube.coord(self.time_coord))
            cube.remove_coord(self.time_coord)
            if time_dim:
                cube.add_dim_coord(new_time_coord, time_dim)
            else:
                # If no time_dim, coordinate is auxiliary or scalar.
                cube.add_aux_coord(new_time_coord)
                        
            new_cubelist.append(cube.copy())
            
        new_cube  = iris.cube.CubeList(new_cubelist).merge_cube()
        if len(new_cube.coord('years').points) > 1:
            new_cube = new_cube.collapsed('years', iris.analysis.MEAN)
                
        self.clim_cube = new_cube
        return self.clim_cube


class FourierClimatology(ForecastData):
    """
    Class for loading climatologies from fourier coefficients.
    
    Args:
    
    * directory: string
        Specify the directory of the climatology data to be loaded. Additional 
        sub directories can given in the load method if required.
        
    * variable: string
        Specify the meteorological variable name. This must be the name as 
        given in the data, e.g. precipitation_flux.
    
    * init_dates: list of datetime objects
        Specify the initialisation dates of the required climatology. See note
        in dates.
        
    * dates: list of datetime objects
        Specify the dates of the required climatology. If dates_are_bounds (see 
        below) is True, all dates in between the max and min dates are loaded.
        Note, the specified year is irrelevant to the climatology, however, it
        is needed to establish and varify the exact order of initialisation 
        dates and climatology dates. E.g. dates covering Nov, Dec, Jan are 
        better specified if years are given, Nov 2000, Dec 2000, Jan 2001. This
        also allows easy matching of forecast dates.
    
    * area_bounds: list of 4 floats or integers
        Specify the upper and lower bounds in the x and y direction. See 
        area_bounds_format for formatting. Default is [-180,-90,180,90].
        
    * dates_are_bounds: boolean
        State whether the dates given are a list of dates or a set of 
        bounds. Default is False.
    
    * area_bounds_format: list of 4 strings
        Provide a list with 'x_min','y_min','x_max' and 'y_max' in the 
        required order to specify the format. Default, ['x_min','y_min',
        'x_max','y_max'].
    
    * area_bounds_range: list of 4 floats or integers
        Specify the range for the x and y coordinates to fall into. Values must
        represent the full cycle of each coordinate system, e.g. 0 to 360 or 
        -180 to 180. This argument is used for shifting coordinate systems, NOT
        for extracting bounds. Values must be given in the same format 
        specified by area_bounds_format. Default is None which uses x range 
        -180 to 180 and y range -90 to 90 and automatically adjusts to the 
        given format.
        
    * lag_str_format: string
        The default loading technique is to find filenames which contain string
        versions of the lags (days between initialisation date and 'forecast' 
        date). The string format used is specified here. Default is '%03d' 
        which corresponds to a number with 3 digits, e.g. 001 or 020.
        
    * unwanted_coords: list of strings
        Specify any coordinate names that are to be removed. Note, removing or
        not removing coordinates can effect how data is merged together. If 
        coordinate name does not exist it is ignored and no exception is 
        raised. Default is to remove forecast_period.
            
    """
    def __init__(self, directory, variable, init_dates, dates,
                  area_bounds=default_area_bounds, 
                  dates_are_bounds=default_dates_are_bounds,
                  area_bounds_format=default_area_bounds_format, 
                  area_bounds_range=default_area_bounds_range, 
                  lag_str_format='%03d',
                  unwanted_coords=['forecast_period']):
        
        self.configure_coord_names()
        self._required_coords = [self.forecast_ref_time, self.time_coord]
        self.directory   = directory
        self.variable    = variable
        # Create AreaBounds instance for use of methods later.
        self._area_inst   = AreaBounds(area_bounds_format, 
                                      area_bounds_range)
        self.area_bounds  = self._area_inst.check_area_bounds(area_bounds)
        self._dates_are_bounds = dates_are_bounds
        self.init_dates, self.dates = self.period_check(init_dates, 
                                                        dates)
        if self._dates_are_bounds:
            self.dates = get_all_dates_between_time_bounds(min(self.dates), 
                                                           max(self.dates))
        self.lag_str_format = lag_str_format
        self._lags_dict = self._get_lags_dict()
        if set(unwanted_coords) & set(self._required_coords):
            diff = set(unwanted_coords) & set(self._required_coords)
            raise ValueError('The coordinate(s) %s are required and can not '\
                             'be removed.' % list(diff))
        self.advanced_search_pattern()
        self.unwanted_coords = unwanted_coords
        self.xy_coords       = None
        self.processes       = []
        self.cube            = None
        
    @staticmethod
    def calculate_day_of_year(date):
        """
        Calculate the day number, 1st Jan is 1. Ignore leap years.
        
        Args:
        
        * date: datetime object
        
        Returns:
            Integer
        
        """
        day_of_year = date.timetuple().tm_yday
        if (is_leap_year(date.year) and \
            date > datetime.datetime(date.year, 2, 28)):
            day_of_year -= 1
            
        return day_of_year

    @staticmethod
    def _calculate_date(day_of_year):
        """
        Calculate the month and day of a day of year number. Day of year is 
        assumed to be from a non-leap year.
        
        Args:
        
        * day_of_year: integer
        
        Returns:
            String
        
        """
        date = datetime.datetime.strptime(str(day_of_year), '%j')
        return date.strftime('%d-%b')

    def _get_lags_dict(self):
        """
        Create a dictionary of all lag times and the forecast days of year that
        correspond to that lag.
        
        Returns:
            Dictionary
        
        """
        lags_dict = {}
        for fcst_date in self.dates:
            day_of_year = self.calculate_day_of_year(fcst_date)
            for init_date in self.init_dates:
                lag = day_of_year - self.calculate_day_of_year(init_date)
                days_of_year = lags_dict.get(lag)
                if days_of_year:
                    days_of_year.append(day_of_year)
                else:
                    lags_dict[lag] = [day_of_year]
                
        return lags_dict

    @staticmethod
    def _calculate_clim_data(day_of_year, fourier_coeffs):
        """
        Given a cubelist of fourier coefficients, return the array containing 
        the climatology for given day of year.
        
        Args:
        
        * day_of_year: integer
        
        * fourier_coeffs: iris cubelist
            Cubelist containing 9 cubes of fourier coeffients.
        
        Returns:
            numpy array
        
        """
        # Convert day_of_year into radians.
        rads = (day_of_year * 2 * math.pi) / 365
        clim_data = fourier_coeffs[0].data + (2 * \
                                  (fourier_coeffs[1].data*math.cos(rads)   + \
                                   fourier_coeffs[2].data*math.cos(2*rads) + \
                                   fourier_coeffs[3].data*math.cos(3*rads) + \
                                   fourier_coeffs[4].data*math.cos(4*rads) - \
                                   fourier_coeffs[5].data*math.sin(rads)   - \
                                   fourier_coeffs[6].data*math.sin(2*rads) - \
                                   fourier_coeffs[7].data*math.sin(3*rads) - \
                                   fourier_coeffs[8].data*math.sin(4*rads)))
            
        return clim_data

    def _get_metadata(self):
        """
        Retrieve the metadata from the cube.
                
        Returns:
            class MetaData
        
        """    
        metadata = {'DATA_TYPE':'Fourier Climatology'}    
        
        area_bounds = self._area_inst.get_cube_area_bounds(self.cube, 
                                                          self.xy_coords)
        x_bounds = [area_bounds[self._area_inst.x_min], 
                    area_bounds[self._area_inst.x_max]]
        y_bounds = [area_bounds[self._area_inst.y_min], 
                    area_bounds[self._area_inst.y_max]]
          
        metadata['VARIABLE'] = self.cube.name()
        metadata['UNITS'] = str(self.cube.units)
        metadata['INITIALISATION_DATES'] = self.cube_init_dates
        metadata['DATES'] = self.cube_dates
        metadata[self.xy_coords[0].upper()+'_BOUNDS'] = x_bounds
        metadata[self.xy_coords[-1].upper()+'_BOUNDS'] = y_bounds
        
        # Find additional coordinates in cube and add them to metadata.
        for coord in self.cube.coords():
            if coord.name() not in self.unwanted_coords and \
               coord.name() not in self._required_coords and \
               coord.name() not in self.xy_coords:
                metadata[coord.name().upper()] = coord.points
                
        bound_names = [self.xy_coords[0].upper()+'_BOUNDS',
                       self.xy_coords[-1].upper()+'_BOUNDS']
                
        return self.MetaData(metadata, bound_names)

    def advanced_search_pattern(self, prepend_lag_string='', 
                                   append_lag_string='',
                                   see_search_patterns=False):
        """
        Default is to use only the lag values and self.lag_str_format to find
        files. This function can be used to prepend and append strings to 
        narrow the search.
        
        Args:
        
        * prepend_lag_string: string
            String to prepend lag string.
        
        * append_lag_string: string
            String to append lag string.
                
        """
        self.prepend_lag_string = prepend_lag_string
        self.append_lag_string  = append_lag_string
        if see_search_patterns:
            print self.prepend_lag_string + \
                  '"lag_str_format"' + \
                  self.append_lag_string
        
    def load(self, dirs_to_look_in='', see_files=False, raw_data=False, 
             current_dir=False, pre_calculated=False):
        """
        Load data from file.
        
        Args:

        * dirs_to_look_in: list of strings or string
            Specify any sub directories of the main directory (given during the 
            instantiation of the class) of which to search in. Default is to
            search all files and directories beneath the main directory.
                
        * see_files: boolean
            Set to True to return search results without loading the files. 
            Default is False.
            
        * raw_data: boolean
            Once the data is loaded it is sorted and constrained to the 
            parameters given during instantiation. If the raw data i.e. all the
            data within the loaded files is required set to True. Default is 
            False.
        
        * current_dir: boolean
            Only search the given directory for matching files.
        
        * pre_calculated: boolean
            State whether the climatologies on file are still Fourier 
            coefficients or pre-calculated.
    
        Returns:
            iris cube
            
        """
        all_cubes = []
        get_info = True
        # List used when just looking at filenames or loading raw data.
        see_files_to_load = []
        
        for lag in self._lags_dict.keys():
            search_string = self.prepend_lag_string + \
                            self.lag_str_format % lag + \
                            self.append_lag_string 
            files_to_load = self._find_files(self.directory, dirs_to_look_in, 
                                             search_string, current_dir,
                                             see_files)
            for load_file in files_to_load:
                if see_files or raw_data:
                    see_files_to_load.append(load_file)
                    continue
                fourier_cubelist = iris.load(load_file, self.variable)
                if not fourier_cubelist:
                    continue
                if get_info:
                    # Gather universal information from the first data to be 
                    # loaded.
                    if self.xy_coords is None:
                        self.xy_coords = [coord.name() 
                                          for coord in get_xy_coords(
                                                       fourier_cubelist[0])]
                    self.time_unit = fourier_cubelist[0].coord(
                                     self.time_coord).units
                    # Take a copy of a cube.
                    skeleton_cube = fourier_cubelist[0].copy()
                    get_info = False
                
                for day_of_year in self._lags_dict[lag]:
                    if not pre_calculated:
                        clim_data = self._calculate_clim_data(day_of_year, 
                                                              fourier_cubelist)
                        clim_cube = skeleton_cube.copy()
                        clim_cube.data = clim_data
                    else:
                        assert len(fourier_cubelist) == 1, 'Not a valid pre '\
                        'calculated climatology.'
                        clim_cube = fourier_cubelist[0]
                        assert clim_cube.shape[0] == 365, 'Not a valid pre '\
                        'calculated climatology.'
                        clim_cube = clim_cube[day_of_year - 1]
                    # Remove time and forecast ref coords as they are to be 
                    # replaced.
                    clim_cube = remove_coords(clim_cube, 
                                              [self.time_coord,
                                               self.forecast_ref_time]+\
                                               self.unwanted_coords)
                    # Replace time coord with day of year number.
                    clim_cube.add_aux_coord(iris.coords.AuxCoord(
                                            day_of_year,
                                            standard_name=self.time_coord))
                    # Replace forecast ref coord with initialisation 'day of 
                    # year' number.
                    clim_cube.add_aux_coord(iris.coords.AuxCoord(
                                        day_of_year - lag,
                                        standard_name=self.forecast_ref_time))
                    all_cubes.append(clim_cube)
        
        if see_files:
            print '\n'.join(see_files_to_load)
            return
        if raw_data:
            return iris.load(see_files_to_load, self.variable)
        
        all_cubes = iris.cube.CubeList(all_cubes)
        cube = all_cubes.merge_cube()
        cube = self._area_inst.check_cube_area_bounds(cube, self.xy_coords, 
                                                     self.area_bounds)
        cube = self.extract_area_bounds(cubes=cube)
        
        # Get day and month dates from day of year numbers to be added to 
        # atrributes.
        self.cube_init_dates = [self._calculate_date(day_of_year) 
                                for day_of_year in 
                                cube.coord(self.forecast_ref_time).points]
        self.cube_dates = [self._calculate_date(day_of_year) 
                           for day_of_year in 
                           cube.coord(self.time_coord).points]
        cube.coord(self.forecast_ref_time).attributes = {
                                                  'dates':self.cube_init_dates}
        cube.coord(self.time_coord).attributes = {'dates':self.cube_dates}
        self.cube = cube
        self.metadata = self._get_metadata()
        return self.cube
