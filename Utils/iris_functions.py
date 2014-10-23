import iris
import numpy
import datetime
import re
from iris.experimental.regrid import regrid_bilinear_rectilinear_src_and_grid

def cube_time_converter(time, time_unit):
    """
    Convert time between datetime object and number value (standard format in
    an iris cube).

    Args:
    
    * time: datetime, or float, int
        Which ever type is given, it is converted to the other.
    
    * time_unit: iris.unit.Unit
        This describes the number so it can be converted. For example, the time
        unit may be 'hours since 1970'.

    Returns:
        datetime or float

    """
    assert type(time_unit) == iris.unit.Unit, 'time_unit must be iris.unit.Unit, not %s' % type(time_unit)
    if type(time) == datetime.datetime:
        converted_time = time_unit.date2num(time)
    else:
        converted_time = time_unit.num2date(time)
    return converted_time

def get_time_bound_constraints(time_start, time_end, time_coord='time'):
    """
    Return an iris constraint for the given time bounds.
    
    Args:
    
    * time_start: float, integer or datetime
    
    * time_end: float, integer or datetime

    Kwargs:
    
    * time_coord: string
        The name of a time based coordinate

    Returns:
        iris constraint
    """
    def time_func(time):
        if type(time_start) == datetime.datetime:
            time = time.point
        return time_start <= time <= time_end
    
    assert time_start <= time_end, 'Start time must come before end time'
    if type(time_start) == datetime.datetime or \
       type(time_end) == datetime.datetime:
        assert type(time_start) == type(time_end), 'Only one datetime object'\
                                                   ' given.'
    time_constraint = iris.Constraint(**{time_coord : time_func})
    return time_constraint

def get_list_of_time_bound_constraints(times, time_coord='time'):
    """
    Return an iris constraint for a list of times.
    
    Args:
    
    * times: list
        Time points

    Kwargs:
    
    * time_coord: string
        The name of a time based coordinate

    Returns:
        iris constraint
    
    """
    def time_func(time):
        within = False
        for t in times:
            time_start = t[0]
            time_end = t[1]
            if type(time_start) == datetime.datetime:
                if type(time) != datetime.datetime:
                    time = time.point
            if time_start <= time <= time_end:
                within = True
        return within
        
    time_constraint = iris.Constraint(**{time_coord : time_func})
    return time_constraint

def get_time_bounds_from_cube(cube):
    """
    Return the time bounds given by the cube.
    
    Args:
    
    * cube: iris cube
    
    Returns:
        Time bounds
    """
    try:
        cube_time_unit = cube.coord('time').units
        start_time = cube_time_converter(cube.coord('time').bounds[0][0], 
                                         cube_time_unit)
        end_time   = cube_time_converter(cube.coord('time').bounds[-1][-1],
                                         cube_time_unit)
    # If the cube has no time bounds
    except TypeError:
        start_time = cube_time_converter(cube.coord('time').points[0], 
                                         cube_time_unit)
        end_time   = cube_time_converter(cube.coord('time').points[-1],
                                         cube_time_unit)
    except:
        return None
    
    bounds = [start_time, end_time]
            
    return bounds

def get_coordinate_slice_dimensions(cube, coordinates, 
                                        ignore_missing_coords=False):
    """
    Return all dimensions which are not the specified coordinate(s). Used for
    cube slicing.
    
    Args:
    
    * cube: iris cube
    
    * coordinates: string or list of strings
        The name of the coordinate(s).

    Kwargs:
    
    * ignore_missing_coords: boolean
        If True, no errors are raised if a coordinate is not in cube. 
    
    Returns:
        List of coordinates
    """
    if type(coordinates) != list:
        coordinates = [coordinates]
        
    dim_coords = {}
    for coord in cube.coords(dim_coords=True):
        dim_coords[cube.coord_dims(coord)[0]] = coord.name()
    
    req_coords = []
    for coordinate in coordinates:
        if coordinate in [coord.name() for coord in cube.coords()]:
            associated_dims = cube.coord_dims(cube.coord(coordinate))
            for dim in associated_dims:
                dim_coord = dim_coords.get(dim)
                if dim_coord:
                    req_coords.append(dim_coord)
        else:
            if not ignore_missing_coords:
                raise ValueError('"%s" coordinate not in cube' % coordinate)
    
    return [coord for coord in dim_coords.values()
            if coord not in req_coords]

    
def get_cubelist_data_range(cubelist):
    """
    Return minimum and maximum data values over all the cubes.
    
    Args:
    
    * cubelist: iris cubelist
    
    Returns:
        Mininum and maximum values
    """
    dmin = None
    dmax = None
    for cube in cubelist:
        if dmin:
            dmin = min(dmin, cube.data.min())
        else:
            dmin = cube.data.min()
        if dmax:
            dmax = max(dmax, cube.data.max())
        else:
            dmax = cube.data.max()
    
    return dmin, dmax

def get_xy_coords(cube):
    """
    Return the spatial xy dimension coordinates, typically this would be 
    latitude and longitude.
    
    Args:
    
    * cube: iris cube
    
    Returns:
        List containing xy coordinates
    """
    return [cube.coord(axis='X', dim_coords=True), 
            cube.coord(axis='Y', dim_coords=True)]

def order_cubelist(cubelist, coord='time', value='first', reverse=False):
    """
    Order a list of cubes by a given coordinate, default is time.
    
    Args:
    
    * cubelist: list of iris cubes.
    
    Kwargs:
    
    * coord: string
        The name of the coordinate to order the cubelist by, default is time.
    
    * value: 'first' or 'last'
        Whether to order the cubes by the first or last date in each cube. 
        Default is first.
    
    * reverse: boolean
        Set True to return the reverse order cubelist.
    
    Returns:
        iris cubelist
    
    """
    if value == 'first':
        index = 0
    elif value == 'last':
        index = -1
    else:
        raise UserWarning("Invalid value argument %s. Use 'first' or 'last'" 
                          % value)
    coord_dict = {}
    for i, cube in enumerate(cubelist):
        coord_val = cube.coord(coord).points[index]
        coord_dict[coord_val] = i
    
    new_cubes = []  
    for coord_val in sorted(coord_dict.keys()):
        new_cubes.append(cubelist[coord_dict[coord_val]])
    if reverse:
        new_cubes.reverse()
    
    return iris.cube.CubeList(new_cubes)

class AreaBounds(object):
    """
    Class containing methods for area bounds.
    
    Kwargs:
    
    * bounds_format: list, tuple
        Provide a list with 'x_min','y_min','x_max' and 'y_max' in the 
        required order to specify the format. Default, ['x_min','y_min',
        'x_max','y_max'].
    
    * bounds_range: list, tuple
        Specify the min and max values of the x and y coordinates. The values
        are allocated in the order of the specified format.
    """
    def __init__(self, bounds_format=['x_min','y_min','x_max','y_max'],
                  bounds_range=None):
        
        assert len(bounds_format) == 4, '%s format value(s) provided, 4 required.' \
                                         % len(bounds_format)
        self.bounds_format = bounds_format
        # Index values
        self.x_min = bounds_format.index('x_min')
        self.y_min = bounds_format.index('y_min')
        self.x_max = bounds_format.index('x_max')
        self.y_max = bounds_format.index('y_max')
        
        self.bounds_range = bounds_range

    def check_area_bounds(self, bounds, bounds_range=None):
        """
        Check the bounds are in specified format and fall within the given range.
        
        Args:
        
        * bounds: list, tuple
            Provide the bounds to be checked.

        Kwargs:
    
        * bounds_range: 
            Specify alternative bounds_range to self.bounds_range. The values 
            are taken to be in the order of self.bounds_format. Default is to 
            use self.bounds_range.
            
        Returns:
            Valid area bounds.
        """
        if bounds is None:
            return self.bounds_range
        if bounds_range:
            assert len(bounds_range) == 4, '%s range value(s) provided, 4 required.' \
                                            % len(bounds_range)
        else:
            bounds_range = self.bounds_range
        
        if type(bounds) != list and type(bounds) != tuple:
            raise ValueError('Bounds must be given in a list or tuple')
        assert len(bounds) == 4, '%s bound(s) provided, 4 required.' \
                                  % len(bounds)
        if bounds[self.x_min] > bounds[self.x_max]:
            raise ValueError('Minimum x value, %s, is greater than the maximum, %s.' % (bounds[self.x_min],
                                                                                        bounds[self.x_max]))
        
        if bounds[self.y_min] > bounds[self.y_max]:
            raise ValueError('Minimum y value, %s, is greater than the maximum, %s.' % (bounds[self.y_min],
                                                                                        bounds[self.y_max]))
        if bounds_range:
            out_of_bounds = {}
            if round(bounds[self.x_min], 2) < bounds_range[self.x_min]:
                out_of_bounds['x_min'] = '%s > %s' % (bounds_range[self.x_min], bounds[self.x_min])
            if round(bounds[self.y_min], 2) < bounds_range[self.y_min]:
                out_of_bounds['y_min'] = '%s > %s' % (bounds_range[self.y_min],  bounds[self.y_min])
            if round(bounds[self.x_max], 2) > bounds_range[self.x_max]:
                out_of_bounds['x_max'] = '%s < %s' % (bounds_range[self.x_max], bounds[self.x_max])
            if round(bounds[self.y_max], 2) > bounds_range[self.y_max]:
                out_of_bounds['y_max'] = '%s < %s' % (bounds_range[self.y_max], bounds[self.y_max])
            
            if out_of_bounds:
                message = 'The following coordinates are out of bounds:\n'
                for key, val in out_of_bounds.items():
                    message += '{key} value {val}\n'.format(key=key,val=val)
                raise ValueError(message)
    
        return bounds
    
    def get_area_bound_constraints(self, bounds, xy_coords=None):
        """
        Return an iris constraint for the given area bounds. Default format,
        [x-min, y-min, x-max, y-max]
        
        Args:
        
        * bounds: list, tuple

        Kwargs:
        
        * xy_coords: list
            The names of the x and y coordinates in that order.
            
        Returns:
            iris constraint
        """
        assert len(bounds) == 4, '%s bound(s) provided, 4 required.' \
                                  % len(bounds)
        if not xy_coords:
            xy_coords = ['longitude', 'latitude']
        
        x_contraint_func = lambda x: bounds[self.x_min] <= x <= bounds[self.x_max]
        y_contraint_func = lambda y: bounds[self.y_min] <= y <= bounds[self.y_max]
        
        x_constraint = iris.Constraint(**{xy_coords[0] : x_contraint_func})
        y_constraint = iris.Constraint(**{xy_coords[-1]: y_contraint_func})
        
        return x_constraint & y_constraint
    
    def get_cube_area_bounds(self, cube, xy_coords=None):
        """
        Return the area bounds given by the cube, in the required format,
        
        Args:
        
        * cube: iris cube

        Kwargs:
        
        * xy_coords: list
            The x and y coordinates in that order.
        
        Returns:
            Area bounds
        """
        if not xy_coords:
            xy_coords = get_xy_coords(cube)
        else:
            xy_coords = [cube.coord(xy_coords[0]),
                         cube.coord(xy_coords[-1])]
        x = xy_coords[0].points
        y = xy_coords[-1].points
        
        bounds_dict = {str(self.x_min) : float(numpy.min(x)),
                       str(self.y_min) : float(numpy.min(y)),
                       str(self.x_max) : float(numpy.max(x)),
                       str(self.y_max) : float(numpy.max(y))}
        bounds = []
        for i in range(4):
            bounds.append(bounds_dict[str(i)])
            
        return bounds
    
    def check_cube_area_bounds(self, cube, xy_coords=None, bounds=None):
        """
        Check the self.bounds_range or given bounds are within the range of the
        cube area bounds. If not, shift the coordinate system to fit. E.g. if
        the cube longitude bounds are 0-360 but bounds_range specifies -180-180
        then the longitude coordinate is shifted.
        
        Args:
        
        * cube: iris cube
            cube to be checked.

        Kwargs:
        
        * xy_coords: list
            The names of the x and y coordinates in that order.
        
        * bounds: list, tuple
            Provide alternative bounds to check if they fall within the cube
            area bounds. This can be used if bounds which will be extracted 
            from the cube already fall within the cubes range and so a shift
            is not required. E.g. cube longitude bounds = [0,360], 
            self.bounds_range = [-180,180] but the bounds to be extracted are 
            [20,50], in which case no shift is required.
            
        Returns:
            iris cube
        
        """
        cube_bounds = self.get_cube_area_bounds(cube, xy_coords)
        if not xy_coords:
            xy_coords = [coord.name() for coord in get_xy_coords(cube)]
        if bounds:
            bounds = self.check_area_bounds(bounds)
        else:
            bounds = self.bounds_range
        try:
            # Test if bounds fall within the cube bounds.
            self.check_area_bounds(bounds, cube_bounds)
        except ValueError as err_message:
            # Read error message to identify coordinate to be shifted.
            if re.findall('x_', str(err_message)):
                new_range = [self.bounds_range[self.x_min],
                             self.bounds_range[self.x_max]]
                cube = self.coordinate_shift(cube, xy_coords[0], new_range)
            if re.findall('y_', str(err_message)):
                new_range = [self.bounds_range[self.y_min],
                             self.bounds_range[self.y_max]]
                cube = self.coordinate_shift(cube, xy_coords[-1], new_range)
        
        return cube

    @staticmethod
    def coordinate_shift(cube, coordinate, new_range, transpose_back=False):
        """
        Shift the given coordinate and associated data to the new range.
        
        Args:
        
        * cube: iris cube
            cube to be shifted.
        
        * coordinate: string
            Name of the coordinate within the cube which is to be shifted.
        
        * new_range: list of 2 floats or integers
            Specify the new range which the coordinate is to be shifted to. The
            range must be of the format [min, max]. Note, the new_range must 
            specify a full cycle of the coordinate system, e.g. [0,360] or 
            [-180,180] where there are 360 degrees. If a non-full cycle e.g. 
            [-100,100] for a 360 degree coordinate, is given, the function 
            will assume a 200 point system and the resulting data is likely to 
            contain errors.

        Kwargs:
        
        * transpose_back: boolean
            This function can change the dimension order of the cube. If 
            transpose_back is set to True, the original dimension order is 
            restored. Note, this can add on significant processing time for 
            large data sets.
            
        Returns:
            iris cube

        """        
        coord_points = cube.coord(coordinate).points
        if max(coord_points) > max(new_range) and \
           min(coord_points) < min(new_range):
            raise UserWarning('Coordinate shift cannot be performed for a '\
                              'new range, %s, which lies within the cubes '\
                              'existing range, %s.' % ([min(new_range),
                                                        max(new_range)],
                                                       [min(coord_points),
                                                        max(coord_points)]))
        # Get all data from cube which is outside required range.
        to_change_func = lambda cell: cell < new_range[0] or \
                                      cell > new_range[-1]
        to_change_constraint = iris.Constraint(**{coordinate : to_change_func})
        to_change_cube = cube.extract(to_change_constraint)
        # If there is none return cube.
        if not to_change_cube:
            return cube
        
        # Calculate the shift.
        # Direction of shift.
        min_diff = min(new_range) - min(cube.coord(coordinate).points)
        direction = min_diff / abs(min_diff)
        # Magnitude is given by the new range
        magnitude = max(new_range) - min(new_range)
        print 'The coordinate, %s, is being shifted and is taken to be a %s '\
              'point system.' % (coordinate, magnitude)
        shift = magnitude * direction
        
        # Shift the coordinate points
        coord = to_change_cube.coord(coordinate)
        coord.points = coord.points + shift
        
        # Extract all the data which does not need changing. Function is exact 
        # opposite of to_change_func.
        to_remain_func = lambda cell: new_range[0] <= cell <= new_range[-1]
        to_remain_constraint = iris.Constraint(**{coordinate : to_remain_func})
        to_remain_cube = cube.extract(to_remain_constraint)
        
        # To merge, both changed and unchanged cubes are sliced by the 
        # coordinate.
        coord_slices = get_coordinate_slice_dimensions(to_remain_cube, 
                                                       coordinate)
        cubelist = []
        for coord_slice in to_remain_cube.slices(coord_slices):
            cubelist.append(coord_slice)
        for coord_slice in to_change_cube.slices(coord_slices):
            cubelist.append(coord_slice)
        
        cubelist = iris.cube.CubeList(cubelist).merge()
        if len(cubelist) != 1:
            raise UserWarning('Coordinate shift unsuccessful.')
        shifted_cube = cubelist[0]
        
        if transpose_back:
            # Place the coordinate dimensions back in the same order.
            old_coord_dim_indx = cube.coord_dims(cube.coord(coordinate))[0]
            new_coord_dim_indx = shifted_cube.coord_dims(
                                 cube.coord(coordinate))[0]
            dim_order = range(len(shifted_cube.coords(dim_coords=True)))
            dim_order.pop(new_coord_dim_indx)
            dim_order.insert(old_coord_dim_indx, new_coord_dim_indx)        
            shifted_cube.transpose(dim_order)
                
        return shifted_cube

def pad_coords(cubelist, dim_coord_names):
    """
    Returns a :class:`CubeList` of cubes which have
    cubes which cover all coord values that appear
    in any of the source cubes, with masked data values
    inserted to pad the cubes to the correct shape.

    """

    if type(dim_coord_names) is str:
        dim_coord_names = [dim_coord_names]

    new_cubes = iris.cube.CubeList([c.copy() for c in cubelist])

    for dim_coord_name in dim_coord_names:
        # newer_cubes is updated at the end of
        # the for loop i.e after performing
        # operation of one dimension
        newer_cubes = iris.cube.CubeList([])

        a_cube = new_cubes[0]
        a_cube_coord_dim, = a_cube.coord_dims(dim_coord_name)
        a_cube_dim_coord, = a_cube.coords(
            dimensions=a_cube_coord_dim, dim_coords=True
        )
        a_cube_aux_coords = a_cube.coords(
            dimensions=a_cube_coord_dim, dim_coords=False
        )

        coord_is_compatible = (
            a_cube_dim_coord.is_compatible(
                other.coord(a_cube_dim_coord.name())
            )
            for other in cubelist
        )
        if not all(coord_is_compatible):
            raise ValueError("%s dim_coords are not compatible"
                             % a_cube_dim_coord.name())

        coord_is_compatible = (
            a_cube_aux_coord.is_compatible(
                other.coord(a_cube_aux_coord.name())
            )
            for other in cubelist
            for a_cube_aux_coord in a_cube_aux_coords
        )
        if not all(coord_is_compatible):
            raise ValueError("aux_coords are not compatible")

        all_pts = (pt for c in cubelist
                   for pt in c.coord(a_cube_dim_coord._as_defn()).points)
        new_dim_coord_points = sorted(list(set(all_pts)))
        new_dim_coord = iris.coords.DimCoord(
            new_dim_coord_points,
            **a_cube_dim_coord._as_defn()._asdict()
        )
        if hasattr(a_cube_dim_coord, 'bounds') \
                and a_cube_dim_coord.bounds is not None:
            all_pts = (tuple(pt)
                       for c in cubelist for pt in
                       c.coord(a_cube_dim_coord._as_defn()).bounds)
            new_dim_coord.bounds = sorted(list(set(all_pts)))

        new_aux_coords = []
        if len(a_cube_aux_coords) > 0:
            for a_cube_aux_coord in a_cube_aux_coords:
                new_aux_coord_points = \
                    numpy.array([numpy.nan]*len(new_dim_coord.points))
                new_aux_coord = iris.coords.AuxCoord(
                    new_aux_coord_points,
                    **a_cube_aux_coord._as_defn()._asdict()
                )
                new_aux_coords.append(new_aux_coord)

            for cube in cubelist:
                ind_data = numpy.in1d(new_dim_coord.points,
                                   cube.coord(dim_coord_name).points)
                for new_aux_coord in new_aux_coords:
                    cubes_aux_coord = cube.coord(new_aux_coord._as_defn())
                    # check if aux coords are different
                    # for same dim coord values
                    # in different cubes. i.e. if any
                    # replacement coord is different
                    # and present coord is not a nan
                    if any(
                        (cubes_aux_coord.points
                         != new_aux_coord.points[ind_data])
                            & (~numpy.isnan(new_aux_coord.points[ind_data]))):
                        raise ValueError(
                            "%s aux_coord values are different "
                            "for the same dim_coord "
                            "values on different cubes"
                            % new_aux_coord.name()
                        )
                    new_aux_coord.points[ind_data] = cubes_aux_coord.points
                    if (hasattr(cubes_aux_coord, 'bounds')
                            and cubes_aux_coord.bounds is not None):
                        all_pts = (
                            tuple(pt)
                            for c in cubelist
                            for pt in cubes_aux_coord.bounds
                        )
                        new_aux_coord.bounds = sorted(list(set(all_pts)))

        for cube in new_cubes:
            old_dim_coord = cube.coord(dim_coord_name)
            new_shape = list(cube.shape)
            new_shape[numpy.array(a_cube_coord_dim)] \
                = len(new_dim_coord.points)
            new_shape = tuple(new_shape)

            new_data = numpy.ma.empty(new_shape)
            new_data.mask = True
            insert_real_data_index = [slice(None)]*len(cube.shape)
            insert_real_data_index[numpy.array(a_cube_coord_dim)] \
                = numpy.in1d(new_dim_coord.points, old_dim_coord.points)

            new_data[insert_real_data_index] = cube.data
            try:
                new_data.mask[insert_real_data_index] = cube.data.mask
            except AttributeError:
                new_data.mask[insert_real_data_index] = False

            dim_coords_and_dims = [[dc, dim]
                                   for (dc, dim)
                                   in cube._dim_coords_and_dims
                                   if dim != a_cube_coord_dim]
            dim_coords_and_dims.append([new_dim_coord, (a_cube_coord_dim)])
            aux_coords_and_dims = [[ac, dim] for (ac, dim)
                                   in cube._aux_coords_and_dims
                                   if dim != (a_cube_coord_dim,)]
            aux_coords_and_dims.extend(
                [[new_aux_coord, (a_cube_coord_dim)]
                 for new_aux_coord in new_aux_coords]
            )

            new_cube = iris.cube.Cube(
                new_data,
                dim_coords_and_dims=dim_coords_and_dims,
                aux_coords_and_dims=aux_coords_and_dims
            )

            new_cube.metadata = cube.metadata

            newer_cubes.append(new_cube)
        new_cubes = newer_cubes[:]

    return new_cubes


#-----To be tested/developed-----#
def load_mask(mask_path):
    """
    Load mask file and make sure it has GeogCS coordinate system.
    
    """
    mask = iris.load_cube(mask_path)
    # Irrespective of the metadata in the files, we know both masks
    # should be defined on WGS84.
    geog = iris.coord_systems.GeogCS(6378137, 6356752)
    mask.coord('latitude').coord_system = geog
    mask.coord('longitude').coord_system = geog

    # Ensure we have grid cell bounds.
    mask.coord('latitude').guess_bounds()
    mask.coord('longitude').guess_bounds()

    return mask

def apply_mask(cube, maskfile, keep_existing_mask=True):
    """
    Apply mask in maskfile to cube.
    
    """
    mask = load_mask(maskfile)
    xy_coords = get_xy_coords(cube)
    xy_coord_names = [coord.name() for coord in xy_coords]
    
    masked_xy_slices = []
    for xy_slice in cube.slices(xy_coord_names):
        xy_slice = regrid_bilinear_rectilinear_src_and_grid(xy_slice, mask)
        if xy_slice.data.shape != mask.data.shape:
            mask.transpose([1,0])
        # Make sure data is a masked array.
        xy_slice.data = numpy.ma.array(xy_slice.data)
        if keep_existing_mask:
            xy_slice.data.mask = numpy.logical_or(xy_slice.data.mask, 
                                                  mask.data.mask)
        else:
            xy_slice.data.mask = xy_slice.data.mask
        masked_xy_slices.append(xy_slice)
    cube = iris.cube.CubeList(masked_xy_slices).merge_cube()
    return cube
