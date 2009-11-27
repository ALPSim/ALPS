import copy
import numpy as np

class DataSet:
	def __deepcopy__(self,memo):
		ret = DataSet()
		ret.props = copy.deepcopy(self.props,memo)
		ret.x = copy.deepcopy(self.x,memo)
		ret.y = copy.deepcopy(self.y,memo)
		ret.dx = copy.deepcopy(self.dx,memo)
		ret.dy = copy.deepcopy(self.dy,memo)
		return ret

	def __init__(self):
		self.props = {}
		self.x = np.array([])
		self.y = np.array([])

		self.dx = np.array([])
		self.dy = np.array([])