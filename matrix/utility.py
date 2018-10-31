"""
Utility module for matrix functions.
"""
DASH_CHAR = '-'


def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
    """
    Creates a scoring matrix based on an alphabet and scoring values.

    @param Set characters that make up the matrix rows and columns.
    @param int scoring value for a cell with equal row and column values.
    @param int scoring value for a cell with unequal row and column values,
      neither of which contain a dash.
    @param int scoring value for a cell in which either the row or the column
      contains a dash.

    @return dict{dict} indexed by pairs of characters in alphabet plus '-'. The
      outer dictionary keys represent the rows, while the inner dictionary keys
      represent the columns. The value of the inner dictionary is the score.
    """
    try:
        diag_score = int(diag_score)
        off_diag_score = int(off_diag_score)
        dash_score = int(dash_score)
    except ValueError:
        print('Scores must be integer values.')

    scoring_dict = {}

    for i in alphabet:
        scoring_dict[i] = {}
        for j in alphabet:
            if i == j:
                scoring_dict[i][j] = diag_score
            else:
                scoring_dict[i][j] = off_diag_score
        # Add dash score.
        scoring_dict[i][DASH_CHAR] = dash_score

    # Add dash row and column.
    scoring_dict[DASH_CHAR] = {letter: dash_score for letter in alphabet}
    scoring_dict[DASH_CHAR][DASH_CHAR] = diag_score

    return scoring_dict


def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
    """
    Computes the alignment matrix for seq_x and seq_y.

    @param list Sequence to align.
    @param list Sequence to align.
    @param dict row/column dictionary of scores.
    @param boolean If False, provide the local alignment sequence.

    @return dict Alignment matrix.
    """
    m = len(seq_x)
    n = len(seq_y)

    s = {
        0: {
            0: 0
        }
    }

    for i in range(1, m + 1):
        value = determine_alignment_score(
            s[i - 1][0] + scoring_matrix[seq_x[i - 1]][DASH_CHAR],
            not global_flag
        )
        s[i] = {0: value}

    for j in range(1, n + 1):
        value = determine_alignment_score(
            s[0][j - 1] + scoring_matrix[DASH_CHAR][seq_y[j - 1]],
            not global_flag
        )
        s[0][j] = value

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            max_value = max(
                s[i - 1][j - 1] + scoring_matrix[seq_x[i - 1]][seq_y[j - 1]],
                s[i - 1][j] + scoring_matrix[seq_x[i - 1]][DASH_CHAR],
                s[i][j - 1] + scoring_matrix[DASH_CHAR][seq_y[j - 1]]
            )

            s[i][j] = determine_alignment_score(max_value, not global_flag)

    return s


def determine_alignment_score(value, is_local_alignment):
    """
    Replaces negative values with zero when calculating local alignment scores.
    """
    if is_local_alignment and value < 0:
        return 0

    return value
