# hungarian algorithm
class HungarianError(Exception):
  pass

try:
  import numpy as np
except:
  raise HungarianError("NumPy is not installed.")

class Hungarian:
  def __init__(self, input_matrix=None, is_profit_matrix=False):
    '''
    input matrix is a list of Lists
    input_matrix is assumed to be a cost matrix unless is_profit_matrix is True
    '''
    if input_matrix is not None:
      my_matrix = np.array(input_matrix)
      
