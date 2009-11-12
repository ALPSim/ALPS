class InvalidInput(Exception):
	def __init__(self,what):
		self.what = what
	
	def __str__(self):
		return self.what

class EmptyInputPort(Exception):
	def __init__(self,which):
		self.which = which
	
	def __str__(self):
		return 'Missing input on ' + self.which
