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

    for i_idx in alphabet:
        scoring_dict[i_idx] = {}
        for j_idx in alphabet:
            if i_idx == j_idx:
                scoring_dict[i_idx][j_idx] = diag_score
            else:
                scoring_dict[i_idx][j_idx] = off_diag_score
        # Add dash score.
        scoring_dict[i_idx][DASH_CHAR] = dash_score

    # Add dash row and column.
    scoring_dict[DASH_CHAR] = {letter: dash_score for letter in alphabet}
    scoring_dict[DASH_CHAR][DASH_CHAR] = dash_score

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
    m_var = len(seq_x)
    n_var = len(seq_y)

    s_var = {
        0: {
            0: 0
        }
    }

    for i_idx in range(1, m_var + 1):
        value = determine_alignment_score(
            s_var[i_idx - 1][0] + scoring_matrix[seq_x[i_idx - 1]][DASH_CHAR],
            not global_flag
        )
        s_var[i_idx] = {0: value}

    for j_idx in range(1, n_var + 1):
        value = determine_alignment_score(
            s_var[0][j_idx - 1] + scoring_matrix[DASH_CHAR][seq_y[j_idx - 1]],
            not global_flag
        )
        s_var[0][j_idx] = value

    for i_idx in range(1, m_var + 1):
        for j_idx in range(1, n_var + 1):
            max_value = max(
                s_var[i_idx - 1][j_idx - 1] +
                scoring_matrix[seq_x[i_idx - 1]][seq_y[j_idx - 1]],
                s_var[i_idx - 1][j_idx] +
                scoring_matrix[seq_x[i_idx - 1]][DASH_CHAR],
                s_var[i_idx][j_idx - 1] +
                scoring_matrix[DASH_CHAR][seq_y[j_idx - 1]]
            )

            s_var[i_idx][j_idx] = determine_alignment_score(
                max_value, not global_flag)

    alignment_matrix = []
    for row in range(len(s_var)):
        alignment_matrix.append([value for value in s_var[row].itervalues()])

    return alignment_matrix


def determine_alignment_score(value, is_local_alignment):
    """
    Replace negative values with zero when calculating local alignment scores.

    @param int Value to consider for the alignment matrix.
    @param boolean Indicate whether the alignment being calculated is local or
        global.

    return int Value to insert into the alignment matrix.
    """
    if is_local_alignment and value < 0:
        return 0

    return value


def compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    Takes two sequences whose elements share a common alphabet with the scoring
    matrix. Computes a global alignment of the sequences using the global 
    alignment matrix.

    @param list Sequence to align.
    @param list Sequence to align.
    @param dict matrix of alphabet pairing scores.
    @param dict matrix of alignment scores.

    @return tuple (score, alignment of x, alignment of y) Tuple containing the
        score and corresponding alignment of sequences x and y.
    """
    i_idx = len(seq_x)
    j_idx = len(seq_y)

    x_prime = ''
    y_prime = ''

    while i_idx != 0 and j_idx != 0:
        if alignment_matrix[i_idx][j_idx] == alignment_matrix[i_idx - 1][j_idx - 1] + \
                scoring_matrix[seq_x[i_idx - 1]][seq_y[j_idx - 1]]:
            x_prime = seq_x[i_idx - 1] + x_prime
            y_prime = seq_y[j_idx - 1] + y_prime
            i_idx -= 1
            j_idx -= 1
        elif alignment_matrix[i_idx][j_idx] == alignment_matrix[i_idx - 1][j_idx] + \
                scoring_matrix[seq_x[i_idx - 1]][DASH_CHAR]:
            x_prime = seq_x[i_idx - 1] + x_prime
            y_prime = DASH_CHAR + y_prime
            i_idx -= 1
        else:
            x_prime = DASH_CHAR + x_prime
            y_prime = seq_y[j_idx - 1] + y_prime
            j_idx -= 1

    while i_idx != 0:
        x_prime = seq_x[i_idx - 1] + x_prime
        y_prime = DASH_CHAR + y_prime
        i_idx -= 1

    while j_idx != 0:
        x_prime = DASH_CHAR + x_prime
        y_prime = seq_y[j_idx - 1] + y_prime
        j_idx -= 1

    if len(x_prime) != len(y_prime):
        raise ValueError("The sequences are not of the same length.")

    score = 0
    for idx in range(len(x_prime)):
        score += scoring_matrix[x_prime[idx]][y_prime[idx]]

    return (score, x_prime, y_prime)


def compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    Takes two sequences whose elements share a common alphabet with the scoring
    matrix. Computes a global alignment of the sequences using the local 
    alignment matrix.

    @param list Sequence to align.
    @param list Sequence to align.
    @param dict matrix of alphabet pairing scores.
    @param dict matrix of alignment scores.

    @return tuple (score, alignment of x, alignment of y) Tuple containing the
        score and corresponding alignment of sequences x and y.
    """
    i_idx = len(seq_x)
    j_idx = len(seq_y)

    x_prime = ''
    y_prime = ''

    max_value = 0
    for row_idx in range(len(alignment_matrix)):
        for col_idx in range(len(alignment_matrix[row_idx])):
            if alignment_matrix[row_idx][col_idx] > max_value:
                max_value = alignment_matrix[row_idx][col_idx]
                i_idx = row_idx
                j_idx = col_idx    

    while i_idx != 0 and j_idx != 0 and alignment_matrix[i_idx][j_idx] != 0:        
        if alignment_matrix[i_idx][j_idx] == alignment_matrix[i_idx - 1][j_idx - 1] + \
                scoring_matrix[seq_x[i_idx - 1]][seq_y[j_idx - 1]]:
            x_prime = seq_x[i_idx - 1] + x_prime
            y_prime = seq_y[j_idx - 1] + y_prime
            i_idx -= 1
            j_idx -= 1
        elif alignment_matrix[i_idx][j_idx] == alignment_matrix[i_idx - 1][j_idx] + \
                scoring_matrix[seq_x[i_idx - 1]][DASH_CHAR]:
            x_prime = seq_x[i_idx - 1] + x_prime
            y_prime = DASH_CHAR + y_prime
            i_idx -= 1
        else:
            x_prime = DASH_CHAR + x_prime
            y_prime = seq_y[j_idx - 1] + y_prime
            j_idx -= 1

    if len(x_prime) != len(y_prime):
        raise ValueError("The sequences are not of the same length.")

    score = 0
    for idx in range(len(x_prime)):
        score += scoring_matrix[x_prime[idx]][y_prime[idx]]

    return (score, x_prime, y_prime)
