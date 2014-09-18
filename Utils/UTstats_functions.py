import unittest
from stats_functions import *

class Test_rmse(unittest.TestCase):
    def test_rmse(self):
        x = [3,7,6,8,4,6,5,7,9,3]
        y = [4,3,4,3,2,3,5,4,7,5]
        result = rmse(x, y)
        self.assertEqual(result, 2.7568097504180442)

class Test_gerrity_score(unittest.TestCase):
    def test_gerrity_score(self):
        cont_table = [[50, 91, 71],[47,2364,170],[54,205,3288]]
        g_score = gerrity_score(cont_table)
        self.assertAlmostEqual(g_score, 0.572, 3)

class Test_percentile_boundaries(unittest.TestCase):
    def test_boundaries(self):
        bounds = percentile_boundaries([1,2,3,4],3)
        self.assertEqual(bounds, [2,3])

class Test_value_category(unittest.TestCase):
    def setUp(self):
        self.values = [1,2,3,4,5,6,7,8]
        self.bounds = [3,5,7]
        
    def test_argument_error(self):
        self.assertRaises(ValueError, value_category, self.values, self.bounds, 'bad_input')
        self.assertRaises(ValueError, value_category, self.values, self.bounds, 'inner', 'bad_input')
        
    def test_inner_outer_catgories(self):
        outer_cats = value_category(self.values, self.bounds, 'outer', 'upper')
        self.assertEqual(outer_cats, [1,1,1,2,3,3,4,4])
        inner_cats = value_category(self.values, self.bounds, 'inner', 'upper')
        self.assertEqual(inner_cats, [1,1,2,2,3,3,3,4])
        lower_mid_cats = value_category(self.values, self.bounds, 'inner', 'lower')
        self.assertEqual(lower_mid_cats, [1,1,2,2,2,3,3,4])
    

class Test_category_probabilities(Test_value_category):
    def setUp(self):
        super(Test_category_probabilities, self).setUp()
    
    def test_probabilities(self):
        probs = category_probabilities(self.values, self.bounds)
        self.assertEqual(probs, [0.375, 0.125, 0.25, 0.25])
    
    
class Test_skill_score(unittest.TestCase):
    def test_skill_score(self):
        score = skill_score(0.6, 0.2, 1)
        self.assertAlmostEqual(score, 0.5)


class Test_generate_pdf_values(unittest.TestCase):
    def test_pdf_vals(self):
        data = [5,6,7,8,9]
        pdf_vals, pdf_points = generate_pdf_values(data, levels=3)
        self.assertEqual(list(pdf_vals), [0.0066132171305349121, 0.19319522907583744, 0.0066132171305348904])
        self.assertEqual(list(pdf_points), [2.3333333333333339, 7.0000000000000009, 11.666666666666668])


class Test_ProbabilityAccuracyScores(unittest.TestCase):
    def test_bad_input(self):
        # Bad array shape
        self.assertRaises(ValueError, ProbabilityAccuracyScores, [3,2], [[0.2,0.4,0.1,0.2],[0.5,0.3,0.2]])
        self.assertRaises(ValueError, ProbabilityAccuracyScores, [3,2], [[0.2,0.4,0.4],[0.5,0.3,0.2],[1,0,0]])
        # Bad number of observations
        self.assertRaises(ValueError, ProbabilityAccuracyScores, [3,2,1], [[0.2,0.4,0.4],[0.5,0.3,0.2]])
        # Bad obs category value
        self.assertRaises(ValueError, ProbabilityAccuracyScores, [3,4], [[0.2,0.4,0.4],[0.5,0.3,0.2]])
        # Probability set not equal to 1.
        self.assertRaises(ValueError, ProbabilityAccuracyScores, [3,2], [[0.4,0.4,0.4],[0.5,0.3,0.2]])


class Test_brier_score(unittest.TestCase):
    def test_scores(self):
        bad_score = ProbabilityAccuracyScores(1, [0,0,1]).brier_score(split_categories=False)
        self.assertEqual(bad_score, 1)
        good_score = ProbabilityAccuracyScores(1, [1,0,0]).brier_score(split_categories=False)
        self.assertEqual(good_score, 0)
        av_score = ProbabilityAccuracyScores(3, [0.2,0.2,0.2,0.2,0.2]).brier_score(split_categories=False)
        self.assertAlmostEqual(av_score, 0.4)


class Test_probability_score(unittest.TestCase):
    def test_score(self):
        bad_score = ProbabilityAccuracyScores(1, [0,0,1]).probability_score(split_categories=False)
        self.assertEqual(bad_score, 0)
        good_score = ProbabilityAccuracyScores(1, [1,0,0]).probability_score(split_categories=False)
        self.assertEqual(good_score, 1)
        av_score = ProbabilityAccuracyScores(3, [0.2,0.2,0.2,0.2,0.2]).probability_score(split_categories=False)
        self.assertAlmostEqual(av_score, 0.6)


class Test_ranked_probability_score(unittest.TestCase):
    def test_scores(self):
        bad_score = ProbabilityAccuracyScores(1, [0,0,1]).ranked_probability_score(split_categories=False)
        self.assertEqual(bad_score, 0)
        good_score = ProbabilityAccuracyScores(1, [1,0,0]).ranked_probability_score(split_categories=False)
        self.assertEqual(good_score, 1)
        av_score = ProbabilityAccuracyScores(3, [0.2,0.2,0.2,0.2,0.2]).ranked_probability_score(split_categories=False)
        self.assertEqual(av_score, 0.9)


class Test_ROC_score(unittest.TestCase):
    def test_scores(self):
        bad_score = ProbabilityAccuracyScores([1,3,3], [[0,0,1],[1,0,0],[1,0,0]]).ROC_score()
        self.assertEqual(bad_score, [0,None,0])
        good_score = ProbabilityAccuracyScores([1,2,3], [[1,0,0],[0,1,0],[0,0,1]]).ROC_score()
        self.assertEqual(good_score, [1,1,1])
        av_score = ProbabilityAccuracyScores([1,2,3,4], [[0.25,0.25,0.25,0.25],[0.25,0.25,0.25,0.25],[0.25,0.25,0.25,0.25],[0.25,0.25,0.25,0.25]]).ROC_score()
        self.assertEqual(av_score, [0.5,0.5,0.5,0.5])
    
    def test_error(self):
        bad_input = ProbabilityAccuracyScores([1,1,1], [[0,0,1],[1,0,0],[1,0,0]])
        self.assertRaises(UserWarning, bad_input.ROC_score)


class Test_ArrayRegression(unittest.TestCase):
    def test_point_regression(self):
        predictor  = numpy.random.rand(20)
        predictand = numpy.random.rand(20, 10, 10)
        reg_arr = ArrayRegression(predictor, predictand, 3)
        point_coeffs = numpy.polyfit(predictor, predictand[:,4,5], 3)
        point_eq = numpy.poly1d(point_coeffs)
        self.assertAlmostEqual(reg_arr.regression_array[4,5](4.3), point_eq(4.3))
        self.assertAlmostEqual(reg_arr(4.3)[4,5], point_eq(4.3))


class Test_array_category_probabilities(unittest.TestCase):
    def test_probs(self):
        value_array = [[3,6,8], [4,5,5], [2,4,8], [2,6,9]]
        bounds_arrays = [[3.5,5.5,7.5]]
        prob_arrays = array_category_probabilities(value_array, bounds_arrays)
        self.assertEqual(list(prob_arrays[0]), [0.75,  0.5,  0.25])
        
        
        