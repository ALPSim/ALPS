import copy
import numpy as np

class ResultProperties:
    def __deepcopy__(self,memo):
        ret = ResultProperties()
        ret.props = copy.deepcopy(self.props,memo)
        return ret
    
    def __init__(self):
        self.props = {}

class DataSet(ResultProperties):
    def __deepcopy__(self,memo):
        ret = DataSet()
        ret.props = copy.deepcopy(self.props,memo)
        ret.x = copy.deepcopy(self.x,memo)
        ret.y = copy.deepcopy(self.y,memo)
        ret.dx = copy.deepcopy(self.dx,memo)
        ret.dy = copy.deepcopy(self.dy,memo)
        return ret

    def __init__(self):
        ResultProperties.__init__(self)
        
        self.x = np.array([])
        self.y = np.array([])
        
class ResultFile(ResultProperties):
    def __deepcopy__(self,memo):
        ret = ResultFile
        ret.props = copy.deepcopy(self.props,memo)
        return ret
    
    def __init__(self):
        ResultProperties.__init__(self)
        
        