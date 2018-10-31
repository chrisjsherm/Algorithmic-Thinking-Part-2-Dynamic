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
        self.maxDiff = None
        self._diag_score = 10
        self._off_diag_score = 4
        self._dash_score = -4

    def test_build_scoring_matrix(self):
        alphabet = set(['A', 'C', 'G', 'T'])
        self.assertDictEqual(matrix_utility.build_scoring_matrix(
            alphabet, self._diag_score, self._off_diag_score, self._dash_score
        ), {
            'A': {
                'A': 10,
                'C': 4,
                'G': 4,
                'T': 4,
                '-': -4
            },
            'C': {
                'A': 4,
                'C': 10,
                'G': 4,
                'T': 4,
                '-': -4
            },
            'G': {
                'A': 4,
                'C': 4,
                'G': 10,
                'T': 4,
                '-': -4
            },
            'T': {
                'A': 4,
                'C': 4,
                'G': 4,
                'T': 10,
                '-': -4
            },
            '-': {
                'A': -4,
                'C': -4,
                'G': -4,
                'T': -4,
                '-': 10
            },
        })

    def test_compute_alignment_matrix(self):
        seq_x = ['A', 'A', 'T']
        seq_y = ['A', 'G', 'C', 'T']

        alphabet = set(seq_x)
        alphabet |= set(seq_y)

        scoring_matrix = matrix_utility.build_scoring_matrix(
            alphabet, self._diag_score, self._off_diag_score, self._dash_score
        )

        self.assertDictEqual(matrix_utility.compute_alignment_matrix(
            seq_x, seq_y, scoring_matrix, True
        ), {
            0: {0: 0, 1: -4, 2: -8, 3: -12, 4: -16},
            1: {0: -4, 1: 10, 2: 6, 3: 2, 4: -2},
            2: {0: -8, 1: 6, 2: 14, 3: 10, 4: 6},
            3: {0: -12, 1: 2, 2: 10, 3: 18, 4: 20}
        })

        self.assertDictEqual(matrix_utility.compute_alignment_matrix(
            seq_x, seq_y, scoring_matrix, False
        ), {
            0: {0: 0, 1: 0, 2: 0, 3: 0, 4: 0},
            1: {0: 0, 1: 10, 2: 6, 3: 4, 4: 4},
            2: {0: 0, 1: 10, 2: 14, 3: 10, 4: 8},
            3: {0: 0, 1: 6, 2: 14, 3: 18, 4: 20}
        })

    def test_determine_alignment_score(self):
        self.assertEqual(matrix_utility.determine_alignment_score(0, True), 0)
        self.assertEqual(matrix_utility.determine_alignment_score(0, False), 0)
        self.assertEqual(matrix_utility.determine_alignment_score(-5, True), 0)
        self.assertEqual(matrix_utility.determine_alignment_score(5, True), 5)
        self.assertEqual(
            matrix_utility.determine_alignment_score(-5, False), -5)
        self.assertEqual(matrix_utility.determine_alignment_score(5, False), 5)


test_suite = unittest.TestLoader().loadTestsFromTestCase(
    TestMatrixUtility
)
unittest.TextTestRunner(verbosity=2).run(test_suite)
