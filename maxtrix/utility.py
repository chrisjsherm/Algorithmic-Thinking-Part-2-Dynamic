"""
Utility module for matrix functions.
"""

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
  pass