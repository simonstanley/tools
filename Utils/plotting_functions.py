import numpy
import math
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from useful_functions import round_data_limits, find_integer_divisor

def get_sensible_contour_levels(dmin, dmax, about_number_of_levels=21,
                                   centre_val=None, number_of_levels_within=3):
    """
    Given a data range and a requested number of levels, return a list of level
    values where the number of levels and data boundaries are adjusted to 
    compensate sensible step sizes. E.g. Data range -5.43 to 5.39 with about 15
    levels, will return a level value list of 16 steps from -6 to 6 in step 
    sizes of 0.8, as oppose to 15 steps from -5.43 to 5.39 in step sizes of 
    0.7729.
    
    Args:
    
    * dmin: float
        Minimum value
    
    * dmax: float
        Maximum value

    Kwargs:

    * about_number_of_levels: integer
        The requested number of levels. This is not necessarily the number of
        levels returned. Default of 21.
    
    * center_val: float
        Specify a value for the level values to be centred around.
    
    * number_of_levels_within: integer
        Specify the limit of how far either side of about_number_of_levels can
        be tested in order to find sensible contour levels. Default 3. Setting 
        this to be larger, results in preference of adjusting the number of 
        levels over adjustment of the data limits. Setting it smaller does the 
        opposite.
    
    Returns:
        list
    
    """
    if about_number_of_levels <= 1:
        about_number_of_levels = 2
    about_number_of_levels -= 1.
    
    if dmin == dmax:
        raise UserWarning('Maximum and minimum data values are equal.')
    if dmin > dmax:
        dmin, dmax = dmax, dmin
            
    if centre_val is not None:
        dmax_diff = abs(dmax - centre_val)
        dmin_diff = abs(dmin - centre_val)
        if dmax_diff > dmin_diff:
            dmin = centre_val - dmax_diff
        else:
            dmax = centre_val + dmin_diff
    
    # Using the current step size, establish the log10 scale to be worked on.
    # E.g. 0.01, 0.1, 1, 10. The scalar is then used to "normalise" the data
    # range.
    d_range     = dmax - dmin
    step_size   = d_range / about_number_of_levels
    log10_val   = math.log10(step_size)
    log10_scale = math.floor(log10_val)
    scalar      = 10 ** log10_scale
    
    # Round the data limits according to the scalar.
    adj_dmin, adj_dmax = round_data_limits(dmin, dmax, scalar)
    # Get the new data range and scale it.
    scaled_d_range   = (adj_dmax - adj_dmin) / scalar
    scaled_step_size = scaled_d_range / about_number_of_levels
    
    if (scaled_step_size) % 0.5 != 0:
        # If the scaled step size is not a multiple of 0.5, try values around 
        # the requested number of levels.
        result = find_integer_divisor(scaled_d_range, about_number_of_levels,
                                      multiple_of=0.5,
                                      search_range=number_of_levels_within)
        if result:
            scaled_step_size = scaled_d_range / result
        else:
            # Try rounding the data limits to a larger scale.
            re_adj_dmin, re_adj_dmax = round_data_limits(dmin, dmax, 
                                                         (scalar*10))
            scaled_d_range = (re_adj_dmax - re_adj_dmin) / scalar
            result = find_integer_divisor(scaled_d_range, 
                                          about_number_of_levels,
                                          multiple_of=0.5,
                                          search_range=number_of_levels_within)
            if result:
                # Set the re adjusted limits.
                adj_dmin, adj_dmax = re_adj_dmin, re_adj_dmax
                scaled_step_size = scaled_d_range / result
            else:
                # Round to the nearest 0.5 step size for even spacing. Note, 
                # this means the upper data limit may be shifted up more, 
                # leaving a slightly uncentred data range.
                scaled_step_size = round(scaled_step_size * 2.) / 2.
                    
    step_size = scaled_step_size * scalar
    
    contour_levels = []
    d_point = adj_dmin
    # Round everything to avoid floating point precision issue.
    while round(d_point, 8) < round((adj_dmax + step_size), 8):
        contour_levels.append(round(d_point, 8))
        d_point += step_size
    
    contour_levels = numpy.array(contour_levels)
            
    return contour_levels

def colour_map(colours, match_colour=None, match_value=None, dmin=None, 
                dmax=None, data=None, cmap_len=256, extend='neither'):
    """
    Return a matplotlib colour map from a list of colours. A single colour 
    within the list can be assigned a value by providing the data limits that 
    the colour map will be used on. Note, if using this functionality, 
    match_colour, match_value and data limits (or all the data) must be 
    provided.
    
    Args:
    
    * colours: list
        A list of matplotlib accepted colours. These include names (see 
        http://www.w3schools.com/html/html_colornames.asp), html hex colour 
        codes and RGB arrays.
    
    Kwargs:
    
    * match_colour: string or RBG array
        Specify one of the colours in the colour list (but not the first or 
        last) to be matched against a given value (see below).
    
    * match_value: float
        Specify a value to which a given colour is to be matched.
    
    * dmin: float
        Data minimum.
    
    * dmax: 
        Data maximum
        
    * data: array like
        Alternative to providing the limits. Limits are calculated using the 
        data within the function.
    
    * cmap_len: integer
        Total number of colours in the colour map.
    
    * extend: 'neither' or 'both'
        If 'both', the first and last colours are set to under and over data
        range colours.
    
    Returns:
        matplotlib.colors.Colormap
    
    """
    cmap = LinearSegmentedColormap.from_list('cmap', colours, N=cmap_len)
    if match_colour is not None:       
        assert match_value is not None, 'A value must be given with which to '\
        'match the colour.'
        colours = [colour.lower() for colour in colours]
        match_colour = match_colour.lower()
        assert match_colour in colours, 'The colour to match, %s, is not in'\
        ' the given colours list, %s.' % (match_colour, colours)
        if dmin is None or dmax is None:
            assert data is not None, 'To scale the colour map, data or data '\
            'minimum and maximum must be provided.'
            dmin = numpy.min(data)
            dmax = numpy.max(data)
        else:
            assert dmin is not None and dmax is not None, 'Both dmin and dmax'\
            ' must be provided.'
            assert dmin < dmax, 'dmin must be smaller than dmax.'
        
        assert dmin <= match_value <= dmax, 'match_value, %s, value must fall'\
        ' within the data range, %s & %s.' % (match_value, dmin, dmax) 
                                      
        colour_position = float(colours.index(match_colour)) / \
                          float(len(colours) - 1)
        if colour_position in [0., 1]:
            raise UserWarning('The colour to match the value cannot be a '\
                              'colour on the limits of the colour map.')
        value_position = float(match_value - dmin) / \
                               float(dmax - dmin)
                               
        if value_position > colour_position:
            # Cut off the top end of the colour map using equation...
            x = (colour_position * cmap.N) / value_position
            # Take colours from 0 to x (+1 for range to reach x value)
            colour_RGBs = cmap(range(int(round(x + 1))))
            cmap = ListedColormap(colour_RGBs)
        elif value_position < colour_position:
            # Cut off the bottom end of the colour map using equation...
            x = ((colour_position - value_position) * cmap.N) / \
                 (1. - value_position)
            # Take colours from x to end colour index (+1 for range to reach x 
            # value)
            colour_RGBs = cmap(range(int(round(x)), (cmap.N + 1)))
            cmap = ListedColormap(colour_RGBs)
    else:
        assert match_value is None, 'A value has been specified without a '\
        'colour to match it with.'
    if extend == 'both':
        over_colour = cmap(cmap.N)
        under_colour = cmap(0)
        colour_RGBs = cmap(range(1, cmap.N - 1))
        cmap = ListedColormap(colour_RGBs)
        cmap.set_over(over_colour)
        cmap.set_under(under_colour)
        
    return cmap
