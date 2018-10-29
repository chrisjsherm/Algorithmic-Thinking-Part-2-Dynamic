"""
Test suite of the matrix utility module.
"""
import unittest
import utility as matrix_utility


class TestMatrixUtility(unittest.TestCase):
    """
    Unit tests. Each test case must begin with "test_".
    """

    def setUp(self):
        """
        Run before each test.
        """
        pass

    def test_build_scoring_matrix(self):
        alphabet = frozenset(['A', ''])
