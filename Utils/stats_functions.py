import numpy
import scipy.stats
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def rmse(x, y):
    """
    Root-mean-square error
    
    Args:
    
    * x: 1D array like
    
    * y: 1D array like
    
    Returns:
        float
    
    """
    return numpy.sqrt(((numpy.array(x) - numpy.array(y)) ** 2).mean())

def sum_of_squares(data):
    """
    Calculate the sum of squares for an array like.
    
    Args:
    
    * data: array like
    
    Returns:
        float
    
    """
    return numpy.sum([x**2 for x in data])

def spread_data(data, scale):
    """
    Spread all the data points (from the mean) by the given scale.
    
    Args:
    
    * data: array like
    
    * scale: float
    
    Returns:
        array like
    
    """
    mean = numpy.mean(data)
    spread_data = []
    for val in data:
        spread_data.append(((val - mean) * scale) + mean)
    return spread_data

def shift_data(data, shift):
    """
    Shift all data points by given amount.
    
    Args:
    
    * data: array like
    
    * shift: float
    
    Returns:
        array like
        
    """
    return numpy.array(data) + shift

def blend_data(blend_data, fixed_data, blend):
    """
    Blend data towards fixed data using some crazy maths.
    
    Args:
    
    * blend_data: array like
        Data to be blended.
    
    * fixed_data: array like
        Data for which the blend_data is blended towards.
    
    * blend: float
        Percentage value of blend.
    
    Returns:
        array like
    
    """
    fcst_mean = numpy.mean(blend_data)
    fcst_std  = numpy.std(blend_data)
    clim_mean = numpy.mean(fixed_data)
    xbar_of_blend = (((100. - blend) * fcst_mean) + (blend * clim_mean)) / 100.
    xbar_2n = (xbar_of_blend ** 2) * 100.
    sx_2f = ((sum_of_squares(blend_data) * (100. - blend)) / len(blend_data)) \
            + ((sum_of_squares(fixed_data) * blend) / len(fixed_data))
    stdv_of_blend = ((sx_2f - xbar_2n) / 100.) ** 0.5
    blended_data = []
    for val in blend_data:
        adjusted_val = (((val - fcst_mean) / fcst_std) * stdv_of_blend) + \
                       xbar_of_blend
        blended_data.append(adjusted_val)
    
    return blended_data

def gerrity_score(contingency_table, see_scoring_weights=False):
    """
    Calculate the Gerrity score for a given contingency table. The table 
    columns must refer to observations and rows to forecasts.
    
    Args:
    
    * contingency_table: array like
        For the array to represent a vaild contingency table it must have 2
        dimensions only and each must have the same length, e.g. 3x3 or 4x4.
    
    Kwargs:
    
    * see_scoring_weights: boolean
        Set True to print the calculated scoring weights for the table.
    
    Returns:
        float
    
    """
    contingency_table = numpy.array(contingency_table, dtype=numpy.float64)
    # Check the table shape and allocate the number of rows and columns.
    tab_shape = contingency_table.shape
    assert len(tab_shape) == 2 and tab_shape[0] == tab_shape[1], \
           'Invalid contingency table. Table array must have 2 dimensions of '\
           'the same size. Shape of this table: (%s, %s)' % (tab_shape[0], 
                                                             tab_shape[1])
    rows_cols = tab_shape[0]
    
    total = numpy.sum(contingency_table)
    probs_table = contingency_table / total
    obs_distribution = numpy.sum(probs_table, axis=0)
    
    odds_ratios = []
    for j in xrange(rows_cols - 1):
        prob_sum = sum([obs_distribution[i] for i in range(j+1)])
        odds_ratios.append((1 - prob_sum) / prob_sum)
    
    scoring_weights = numpy.zeros(tab_shape)
    for i in xrange(rows_cols):
        for j in xrange(rows_cols):
            if j >= i:
                comp1 = sum([1. / odds_ratios[x] for x in range(i)])
                comp2 = sum([odds_ratios[x] for x in range(j, rows_cols-1)])
                comp3 = j - i
                scoring_weights[i,j] = (comp1 + comp2 - comp3) / \
                                       (rows_cols - 1.)
            else:
                scoring_weights[i,j] = scoring_weights[j,i]
    if see_scoring_weights:
        print scoring_weights
    return numpy.sum(scoring_weights * probs_table)

def percentile_boundaries(data, num_of_categories):
    """
    Return the boundary values which split the given data set into the 
    requested number of categories. E.g. data = [1,2,3,4] split into 3 
    categories would return [2.0, 3.0] as the tercile boundaries.
    
    Args:
    
    * data: array like
    
    * num_of_categories: integer
        The number of categories wanted. Note, the function will always return
        num_of_categories - 1 values.
    
    Returns:
        list
    
    """
    percentiles = numpy.linspace(0, 100, num_of_categories+1)[1:-1]
    bounds = [round(numpy.percentile(data, percentile), 8)
              for percentile in percentiles]
    return bounds

def value_category(values, bounds, boundary_val_cat='outer', 
                     middle_val_cat='upper'):
    """
    Given a set of values and boundaries, return each value's category. 
    Categories are named numerically starting from 1. There are always 
    1 + number of bounds categories.
    
    Args:
    
    * values: float or list of floats
    
    * bounds: list
        A list of boundary values. These are automatically sorted into numeric 
        order.

    Kwargs:

    * boundary_val_cat:
        If a value equals a boundary value, specify whether it is placed in an 
        inner or outer category. Default is outer.
    
    * middle_val_cat:
        If a value equals the middle boundary value (only for odd number of 
        boundaries), specify whether it is placed in the upper or lower 
        category. Default is upper.
    
    Returns:
        list
    
    """
    if boundary_val_cat not in ['inner', 'outer']:
        raise ValueError('%s is not a valid input, use "inner" or "outer"' 
                         % boundary_val_cat)
    if middle_val_cat not in ['upper', 'lower']:
        raise ValueError('%s is not a valid input, use "upper" or "lower"' 
                         % middle_val_cat)
    if not hasattr(values, '__iter__'):
        values = [values]
    bounds.sort()
    num_of_bounds = len(bounds)
    middle_index  = float(num_of_bounds - 1) / 2.
    
    categories = []
    for value in values:
        category_found = False
        for index, bound in enumerate(bounds):
            if value > bound:
                continue
            elif value < bound:
                category_found = True
                categories.append(index + 1)
                break
            else:
                # When value equals a bound.
                if index < middle_index:
                    if boundary_val_cat == 'inner':
                        category = index + 2
                    else:
                        category = index + 1
                    category_found = True
                    categories.append(category)
                    break
                elif index > middle_index:
                    if boundary_val_cat == 'inner':
                        category = index + 1
                    else:
                        category = index + 2
                    category_found = True
                    categories.append(category)
                    break
                else:            
                    # When value equals the middle bound.
                    if middle_val_cat == 'lower':
                        category = index + 1
                    else:
                        category = index + 2
                    category_found = True
                    categories.append(category)
                    break
        if not category_found:
            # The value is above all boundaries.
            categories.append(index + 2)
    return categories

def category_probabilities(values, bounds, boundary_val_cat='outer', 
                              middle_val_cat='upper', return_counts=False):
    """
    Given a set of values and boundaries, return the associated probabilities 
    for each category. There are always 1 + number of bounds categories.
    
    Args:
    
    * values: list
        A list of values.
    
    * bounds: list
        A list of boundary values. These are automatically sorted into numeric 
        order.

    Kwargs:

    * boundary_val_cat:
        If a value equals a boundary value, specify whether it is placed in an 
        inner or outer category. Default is outer.
    
    * middle_val_cat:
        If a value equals the middle boundary value (only for odd number of 
        boundaries), specify whether it is placed in the upper or lower 
        category. Default is upper.

    Returns:
        list
    
    """
    category_counts = [0.] * (len(bounds) + 1)
    num_of_vals     = float(len(values))
    categories = value_category(values, bounds, boundary_val_cat, 
                                middle_val_cat)
    for category in categories:
        category_counts[category - 1] += 1
    
    if return_counts:
        return [int(val) for val in category_counts]
    else:
        category_probs = [val / num_of_vals for val in category_counts]
        return category_probs

def skill_score(accuracy_score, reference_score, perfect_score):
    """
    Calculate the skill score.
    
    Args:
    
    * accuracy_score: float
    
    * reference_score: float
    
    * perfect_score: float
    
    Returns:
        float
    
    """
    skill_score = (accuracy_score - reference_score) / \
                  (perfect_score  - reference_score)
    return skill_score

def pdf_probabilities(pdf, bounds):
    """
    Calculate the area of the PDF in between each bound, hence the probability.

    Args:
    
    * pdf: instance of scipy.stats.gaussian_kde
    
    * bounds: list
        A list of boundary values. These are automatically sorted into numeric 
        order.

    """
    bounds.sort()
    extended_boundaries = [-numpy.inf] + bounds
    extended_boundaries.append(numpy.inf)
    probs = []
    for i in range(len(extended_boundaries) - 1):
        probs.append(pdf.integrate_box_1d(extended_boundaries[i], extended_boundaries[i+1]))
    return probs

def pdf_percentile_boundaries(pdf, num_of_categories, accuracy_factor=50):
    """
    Estimate the boundary values when splitting a PDF in to equally sized 
    areas.
    
    Args:
    
    * pdf: instance of scipy.stats.gaussian_kde
    
    * num_of_categories: integer
        The number of equally sized areas the PDF is split into.
    
    Kwargs:
    
    * accuracy_factor: integer
        The estimation is calculated using iteration, this value specifies how
        many values to split the PDF into and iterate over. Therefore, the
        higher the factor, the longer the calculation takes but the more 
        accurate the returned values. Default is 50.
    
    Returns:
        list of bounds
    
    """
    dmin = numpy.min(pdf.dataset)
    dmax = numpy.max(pdf.dataset)
    x_vals = numpy.linspace(dmin, dmax, accuracy_factor)
    
    required_area_size = 1. / float(num_of_categories)
    bounds = []
    lower_bound = -numpy.inf
    for i, x_val in enumerate(x_vals):        
        this_area_size = pdf.integrate_box_1d(lower_bound, x_val)
        if this_area_size > required_area_size:
            upper_diff = this_area_size - required_area_size
            lower_diff = required_area_size - \
                         pdf.integrate_box_1d(lower_bound, x_vals[i-1])
            total_diff = upper_diff + lower_diff
            proportion_diff = upper_diff / total_diff
            
            val_diff = x_val - x_vals[i-1]
            proportion_val_diff = val_diff * proportion_diff
            adjusted_x_val = x_val - proportion_val_diff
            bounds.append(adjusted_x_val)
            if len(bounds) == num_of_categories - 1:
                break
            lower_bound = adjusted_x_val
            
    return bounds

def calculate_pdf_limits(pdf, levels=50, range_limiter=20):
    """
    Calculate the values where the PDF stops. The range_limiter determines the 
    value at which to cut the PDF outer limits. It is a proportional value not 
    an actual value. The larger the given value the further out the extremes 
    will be returned.

    Args:
    
    * pdf: instance of scipy.stats.gaussian_kde

    Kwargs:
    
    * levels : integer
        This determines the step size when calculating the limits.
    
    * range_limiter: scalar
        This value is used to calculate the range of the PDF. A PDF function 
        can take a while to converge to 0, so to calculate sensible stop and 
        start points, some proportional value above 0 is calculated. The given
        range_limiter value is used as factor to determine what that above 0 
        value is. Simply, the higher the given value the wider the PDF limits.
        See nested function calculate_pdf_limits for more details.

    """
    dmin = numpy.min(pdf.dataset)
    dmax = numpy.max(pdf.dataset)
    pdf_min = numpy.mean([pdf(dmin)[0], pdf(dmax)[0]]) / float(range_limiter)
    # First calculate the appropriate step size given the data range and number
    # of levels.
    step_size = (dmax - dmin) / float(levels)
    while pdf(dmin)[0] > pdf_min:
        dmin -= step_size
    while pdf(dmax)[0] > pdf_min:
        dmax += step_size
    return dmin, dmax

def generate_pdf_values(data, levels=50, range_limiter=20, 
                          bandwidth='silverman', return_pdf=False):
    """
    Calculate the PDF function and return a set of values and points along the 
    curve.
    
    Args:
    
    * data : 1D array like
        List of data values.
    
    Kwargs:
    
    * levels : integer
        This determines how many points are returned. If plotting, higher 
        values lead to smoother plots.
    
    * range_limiter: scalar
        This value is used to calculate the range of the PDF. A PDF function 
        can take a while to converge to 0, so to calculate sensible stop and 
        start points, some proportional value above 0 is calculated. The given
        range_limiter value is used as factor to determine what that above 0 
        value is. Simply, the higher the given value the wider the PDF limits.
        See nested function calculate_pdf_limits for more details.
    
    * bandwidth: string, scalar or callable
        The method used to calculate the estimator bandwidth. This can be 
        'scott', 'silverman', a scalar constant or a callable. If a scalar, 
        this will be used directly as kernel-density estimate (kde) factor. If 
        a callable, it should take a scipy.stats.gaussian_kde instance as only 
        parameter and return a scalar. Default is 'silverman'.
    
    * return_pdf: boolean
        If True, the callable scipy.stats.gaussian_kde instance is also 
        returned.
    
    Returns:
        PDF values, PDF points, PDF function (optional)
    
    """    
    # Generate kernel density estimate (PDF)
    pdf = scipy.stats.gaussian_kde(data, bw_method=bandwidth)
    dmin, dmax = calculate_pdf_limits(pdf, levels, range_limiter)
    pdf_points = numpy.linspace(dmin, dmax, levels)
    if return_pdf:
        return pdf(pdf_points), pdf_points, pdf
    else:
        return pdf(pdf_points), pdf_points

def array_correlation(x, y, method='pearson'):
    """
    Calculate the correlation between each matching element of the arrays in x 
    and y. Note, x and y must be identical in shape.
    
    Args:
    
    * x and y: List of arrays
    
    Kwargs:
    
    * method: 'pearson' or 'spearman'
        Pearson's correlation or Spearman's ranked correlation.
    
    Returns:
        Array of correlation values for corresponding elements.
    
    """
    assert method in ['pearson', 'spearman'], 'Invalid method %s.' % method
    x = numpy.ma.masked_array(x)
    y = numpy.ma.masked_array(y)
    assert x.shape[0] == y.shape[0], 'x and y must contain the same number of'\
    ' arrays.'
    assert x.shape[1:] == y.shape[1:], 'All arrays in x and y must be the '\
    'same shape.'
    
    if method == 'pearson':
        corr_method = scipy.stats.pearsonr
    elif method == 'spearman':
        corr_method = scipy.stats.spearmanr
    
    correlation_array = numpy.empty(x.shape[1:])
    for index in numpy.ndindex(x.shape[1:]):
        array_index = tuple( [slice(None)] + list(index) )
        element_corr = corr_method(x[array_index], y[array_index])[0]
        correlation_array[index] = element_corr
    
    return numpy.ma.masked_invalid(correlation_array)

class ProbabilityAccuracyScores(object):
    """
    Class containing accuracy score methods for probabilistic methods. Lists 
    for observations and probabilities must be given in the same order as each
    other. I.e. ob_categories[x] must refer to category_probs[x].
    
    Args:
    
    * ob_categories: integer or list of integers
        The numbers of the observed categories, 1 referring to the 1st. E.g. 
        for tercile categories, 1 = lower, 2 = middle, 3 = upper.
    
    * category_probs: list
        A list of the probabilities for each category. The sum must equal 1 
        when rounded to 5 decimal places, so irrational numbers need at least
        6 decimal places, e.g. 0.333333.
    
    """
    def __init__(self, ob_categories, category_probs):
        ob_categories  = numpy.array(ob_categories)
        category_probs = numpy.array(category_probs)
        if category_probs.dtype == object:
            raise ValueError('Not all sets of probabilities contain the same'\
                             ' number of categories: %s' % category_probs)
        if category_probs.ndim > 2:
            raise ValueError('Too many dimensions in category_probs.')
        if (ob_categories.ndim + 1) != category_probs.ndim:
            raise ValueError('Must be exactly 1 observation for each set of '\
                             'probabilities.')            
        if ob_categories.ndim == 0:
            ob_categories  = numpy.array([ob_categories])
            category_probs = numpy.array([category_probs])
        if ob_categories.shape[0] != category_probs.shape[0]:
            raise ValueError('Number of observations, %s, does not match the '\
                             'number of probability sets, %s' % (
                              ob_categories.shape[0], category_probs.shape[0]))
        if numpy.max(ob_categories) > category_probs.shape[1] or \
           numpy.min(ob_categories) < 1:
            raise ValueError('Observation category, %s, out of range.' % 
                             numpy.max(ob_categories))
        prob_totals = [round(val, 5) 
                       for val in numpy.sum(category_probs, axis=1)]
        if prob_totals != [1] * category_probs.shape[0]:
            raise ValueError('All probability sets must sum to 1.')
                
        self.ob_categories     = ob_categories
        self.category_probs    = category_probs
        self.categories        = range(1, category_probs.shape[1] + 1)
        self.num_of_categories = len(self.categories)

    def _ROC_plot(self, all_hit_rates, all_false_alarm_rates, ROC_scores,
                   categoriy_names, title, colours, save):
        """
        Plot the ROC curves.
        
        """
        cmap = LinearSegmentedColormap.from_list('cmap', colours)
        colour_index = [int(round(val)) 
                        for val in numpy.linspace(0, 256, 
                                                  len(categoriy_names))]
        plt.figure()
        legend_labels = []
        for i, (hit_rates, fa_rates, ROC_score, category) in enumerate(
                                                         zip(
                                                         all_hit_rates, 
                                                         all_false_alarm_rates,
                                                         ROC_scores,
                                                         categoriy_names)):
            if hit_rates is not None and fa_rates is not None:
                plt.plot(fa_rates, hit_rates, 'o-', 
                         color=cmap(colour_index[i]))
                label= '{cat} = {score}'.format(cat=category,
                                                score=round(ROC_score,
                                                3))
                legend_labels.append(label)
            
        plt.xlabel('False Alarm Rate')
        plt.ylabel('Hit Rate')
        if title is not None:
            plt.title(title)
        else:
            plt.title('ROC Curves')
        plt.legend(tuple(legend_labels),'lower right', fontsize='small', 
                   title='ROC Scores')
        plt.plot([0,1], [0,1], 'k--')
        plt.grid()
        if save is not None:
            plt.savefig(save)
        else:
            plt.show()

    def _brier_score(self, ob_categories, category_probs):
        """
        Method to do the Brier score calculations.
        
        Returns:
            float
        
        """
        num_of_trails = len(category_probs)
        cumalative_score = 0.
        for ob_category, probability_set in zip(ob_categories, 
                                                category_probs):
            for category, category_prob in zip(self.categories, 
                                               probability_set):
                ob_value = 0.
                if category == ob_category:
                    ob_value = 1.
                cumalative_score += (category_prob - ob_value) ** 2
        brier_score = (1. / (2. * num_of_trails)) \
                      * cumalative_score
        return brier_score

    def _ranked_probability_score(self, ob_categories, category_probs):
        """
        Method to do the RPS calculations.
        
        Returns:
            float
        
        """
        num_of_trails = len(category_probs)
        cumalative_score = 0.
        for ob_category, probability_set in zip(ob_categories, 
                                                category_probs):
            ob_cumalative    = 0.
            prob_cumalative  = 0.
            for category, category_prob in zip(self.categories, 
                                               probability_set):
                if category == ob_category:
                    ob_cumalative += 1.
                prob_cumalative  += category_prob
                cumalative_score += (prob_cumalative - ob_cumalative) ** 2
        ranked_prob_score = 1. - (1. / 
                                 ((self.num_of_categories - 1.) \
                                  * num_of_trails) \
                                  * cumalative_score)
        return ranked_prob_score

    def _calculate_score(self, score_method, split_categories):
        """
        Run the given method, either for each category or for all combined.
        
        """
        if split_categories:
            scores = []
            for category in self.categories:
                category_indx = numpy.where(self.ob_categories == category)
                if len(category_indx[0]) > 0:
                    scores.append(round(score_method(
                                  self.ob_categories[category_indx],
                                  self.category_probs[category_indx]), 8))
                else:
                    scores.append(None)
            return scores
        else:
            return round(score_method(self.ob_categories, self.category_probs),
                         8)

    def brier_score(self, split_categories=True):
        """
        Used for unordered categories, e.g. rain or no rain.
        
        Kwargs:
        
        * split_categories: boolean
            If True, the score is calculated separately for each category. 
            Otherwise a single value is returned derived from all categories.
            Default is True.
                
        Returns:
            list or float
        
        """
        return self._calculate_score(self._brier_score, split_categories)

    def probability_score(self, split_categories=True):
        """
        The negative orientation of the Brier score, i.e. 1 = perfect instead 
        of 0.
        
        Kwargs:
        
        * split_categories: boolean
            If True, the score is calculated separately for each category.
            Otherwise a single value is returned derived from all categories.
            Default is True.
            
        Returns:
            list or float
        
        """
        if split_categories:
            return [1 - score for score in self.brier_score(split_categories)]
        else:
            return 1. - self.brier_score(split_categories)

    def ranked_probability_score(self, split_categories=True):
        """
        Used for ordered categories, e.g. tercile categories.

        
        Kwargs:
        
        * split_categories: boolean
            If True, the score is calculated separately for each category.
            Otherwise a single value is returned derived from all categories.
            Default is True.
            
        Returns:
            list or float
        
        """
        return self._calculate_score(self._ranked_probability_score, 
                                     split_categories)

    def ROC_score(self, num_of_thresholds=6, outer_categories_only=False, 
                   plot=False, title=None, category_names=None, 
                   colours=['blue','green','red'], save=None):
        """
        Calculate the Relative operating characteristic score for each 
        category.
        
        Kwargs:
        
        * num_of_thresholds: integer
            The number of thresholds between 0 and 1 at which to compare the 
            probabilities. Adjusting the number of thresholds is used to smooth
            the ROC curves. Higher numbers make less difference for smaller 
            data sets and add on computing time. Adjusting can also affect the 
            ROC scores, so it is recommended to use a consistent number when 
            comparing results. Default is 6 i.e. [0, .2, .4, .6, .8, 1.]
        
        * outer_categories_only: boolean
            If True, only the ROC scores for the lowest and highest categories
            are calculated (and plotted if required).
                
        * plot: boolean
            Set True to plot and show the ROC curves. Default is False.

        Note, all kwargs from here are only used if plot=True

        * title: string
            Specify a title for the ROC plot.

        * category_names: list
            Provide a list of category names. There must be the same number of
            names as there are categories and they categories are labelled from
            lowest to highest.

        * colours: list
            List of the colours from which to create a colour map for each
            category plot.

        * save: string
            If a string is given, the plot is saved (and not shown).

        Returns:
            list of ROC scores for each category.
        
        """
        thresholds = numpy.linspace(0., 1., num_of_thresholds)
        if outer_categories_only:
            categories = [self.categories[0], self.categories[-1]]
        else:
            categories = self.categories
        all_hit_rates         = []
        all_false_alarm_rates = []
        ROC_scores            = []
        for category in categories:
            occurances = (self.ob_categories == category).sum()
            if occurances == 0:
                all_hit_rates.append(None)
                all_false_alarm_rates.append(None)
                ROC_scores.append(None)
                continue
            non_occurances = len(self.ob_categories) - occurances
            if non_occurances == 0:
                raise UserWarning('ROC scores cannot be calculated when all '\
                                  'observations are in the same category.')
            hit_rates         = []
            false_alarm_rates = []
            for threshold in thresholds:
                hits         = 0.
                false_alarms = 0.
                for ob_category, probability_set in zip(self.ob_categories, 
                                                        self.category_probs):
                    # Look at the relevant category only.
                    category_prob = probability_set[category - 1]
                    if category_prob >= threshold:
                        if ob_category == category:
                            hits += 1.
                        else:
                            false_alarms += 1.
                hit_rates.append(hits / float(occurances))
                false_alarm_rates.append(false_alarms / float(non_occurances))
            ROC_score = 0.
            for i in xrange(num_of_thresholds - 1):
                # Calculate the area under each bit of curve using trapezium 
                # rule.
                ROC_score += ((hit_rates[i] + hit_rates[i+1]) / 2.) * \
                             (false_alarm_rates[i] - false_alarm_rates[i+1])
            all_hit_rates.append(hit_rates)
            all_false_alarm_rates.append(false_alarm_rates)
            ROC_scores.append(ROC_score)
        if plot:
            if category_names:
                assert len(category_names) == len(categories), '%s '\
                'category_names must be provided.' % len(categories)
            else:
                category_names = ['Category %s' % cat for cat in categories]
            self._ROC_plot(all_hit_rates, all_false_alarm_rates, ROC_scores,
                           category_names, title, colours, save)
        return ROC_scores


class ArrayRegression(object):
    """
    Given a series of values (predictor) and a series of arrays (predictand) of
    the same length, calculate the regression equations for each point.
    
    Args:
    
    * predictor: 1D array or list
    
    * predictand: array like
        The first dimension must be the same length as the predictor.
    
    Kwargs:
    
    * num_of_polynomials: integer
        Specify the number of polynomials in the calculated regression 
        equations.
        
    """
    def __init__(self, predictor, predictand, num_of_polynomials=1):
        predictor  = numpy.array(predictor)
        predictand = numpy.ma.masked_array(predictand)
        if predictor.ndim != 1:
            raise UserWarning('Predictor must be 1 dimensional.')
        if len(predictor) <= 1:
            raise UserWarning('Predictor must contain at least 2 values.')
        if predictor.shape[0] != predictand.shape[0]:
            raise UserWarning('Predictor and predictand do not have the same '\
                              'first dimension size. %s %s' % (
                                                            predictor.shape,
                                                            predictand.shape))
        self.predictor  = predictor
        self.predictand = predictand
        # For numpy.polyfit to work with masked arrays, all masked values must 
        # be replaced with numpy.nan.
        if self.predictand.mask.any() == True:
            self.predictand[self.predictand.mask] = numpy.nan
                
        self.regression_array = numpy.empty(self.predictand.shape[1:],
                                            dtype=numpy.object)
        # Break up the predictand array into 2 dimensional arrays which numpy
        # can handle.
        if self.predictand.ndim > 2:
            for latter_indicies in numpy.ndindex(self.predictand.shape[2:]):
                index = tuple( [slice(None),slice(None)] + 
                               list(latter_indicies) )
                
                regress_coeffs = numpy.polyfit(self.predictor,
                                                  self.predictand[index],
                                                  num_of_polynomials)
                for i, coeffs in enumerate(regress_coeffs.transpose()):
                    regression_equation = numpy.poly1d(coeffs)
                    self.regression_array[index[1:]][i] = regression_equation
        else:
            regress_coeffs = numpy.polyfit(self.predictor, self.predictand,
                                              num_of_polynomials)
            for i, coeffs in enumerate(regress_coeffs.transpose()):
                regression_equation = numpy.poly1d(coeffs)
                self.regression_array[i] = regression_equation
                
    def __repr__(self):
        def repr_func(poly):
            poly_str = ''
            for coeff, order in zip(poly.coeffs, range(poly.order, -1, -1)):
                if poly.order != order:
                    if coeff < 0:
                        poly_str += ' - '
                    else:
                        poly_str += ' + '
                else:
                    if coeff < 0:
                        poly_str += ' -'
                poly_str += '%s' % round(abs(coeff),2)
                if order > 0:
                    poly_str += 'x'
                    if order >= 2:
                        poly_str += '^%s' % order
            return poly_str
        vfunc = numpy.vectorize(repr_func)
        repr_arr = vfunc(self.regression_array)
        return '%s' % repr_arr
    
    def __call__(self, value):
        """
        Generate an array of values by evaluating each regression equation with
        a specified value.
        
        Args:
        
        * value: float
            Value for which to evaluate each regression equation.
        
        Returns
            numpy array
        
        """
        def evaluate_regression(reg_equation, value):
            return reg_equation(value)
        vfunc = numpy.vectorize(evaluate_regression)
        result = numpy.ma.masked_array(vfunc(self.regression_array, value))
        result = numpy.ma.masked_invalid(result)
        return result
    
class ArrayCategoryProbabilities(object):
    """
    Class for calculating probability values across arrays.
    
    Args:
    
    * value_array: list or array of arrays
        A set of arrays of values
    
    * bounds_arrays: list or array of arrays
        A set of arrays of bounds
    
    Kwargs:
    
    * boundary_val_cat:
        If a value equals a boundary value, specify whether it is placed in an 
        inner or outer category. Default is outer.
    
    * middle_val_cat:
        If a value equals the middle boundary value (only for odd number of 
        boundaries), specify whether it is placed in the upper or lower 
        category. Default is upper.
    
    """
    def __init__(self, value_array, bounds_arrays, boundary_val_cat='outer', 
                  middle_val_cat='upper'):
        bounds_arrays = [numpy.array(arr) for arr in bounds_arrays]
        array_shape   = bounds_arrays[0].shape
        value_array   = numpy.array(value_array)
        array_shape_err_message = 'Array shapes do not match.'
        for array in bounds_arrays:
            assert array.shape == array_shape, array_shape_err_message
        assert value_array.shape[1:] == array_shape, array_shape_err_message
        
        # For indexing the value array during calculation. 
        self.indices = list(numpy.ndindex(array_shape))
        # numpy.vectorize runs the function once to check so start index on -1.
        self.indices_index = -1
        
        self.value_array   = value_array
        self.bounds_arrays = bounds_arrays
        self.boundary_val_cat = boundary_val_cat
        self.middle_val_cat   = middle_val_cat
    
    def caluculate(self):
        """
        Run the calculation of the probabilities for each array point.
        
        Returns:
            tuple of probability arrays, starting with lowest category
        
        """
        def calculate_probs(*bounds):
            index = tuple( [slice(None)] + \
                           list(self.indices[self.indices_index]) )
            self.indices_index += 1
            return tuple(category_probabilities(self.value_array[index], 
                                                list(bounds),
                                                self.boundary_val_cat,
                                                self.middle_val_cat))
        vfunc = numpy.vectorize(calculate_probs)
        return vfunc(*self.bounds_arrays)

def array_category_probabilities(value_array, bounds_arrays,  
                                    boundary_val_cat='outer', 
                                    middle_val_cat='upper'):
    """
    Given a set of value arrays and a set of bounds arrays, calculate the 
    probabilities for each category (defined by the bounds) at each point by 
    counting the number of values within each of the categories. E.g. Given
    10 realisations of an x-y grid and 2 sets of bounds on the same x-y grid, 
    return 3 x-y grids of probabilities, these being probabilities of, beneath 
    the lower bounds, between the bounds, and above the upper bounds. For each 
    array point the probabilities will add up to 1.
    
    Args:
    
    * value_array: list or array of arrays
        A set of arrays of values
    
    * bounds_arrays: list or array of arrays
        A set of arrays of bounds
    
    Kwargs:
    
    * boundary_val_cat:
        If a value equals a boundary value, specify whether it is placed in an 
        inner or outer category. Default is outer.
    
    * middle_val_cat:
        If a value equals the middle boundary value (only for odd number of 
        boundaries), specify whether it is placed in the upper or lower 
        category. Default is upper.
    
    Returns:
        tuple of probability arrays, starting with lowest category
    
    """
    return ArrayCategoryProbabilities(value_array, bounds_arrays, 
                                      boundary_val_cat, 
                                      middle_val_cat).caluculate()