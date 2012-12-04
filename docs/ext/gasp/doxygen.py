from docutils import nodes

class DoxygenNode(nodes.General, nodes.Element):    
    def __init__(self,  **kwargs):
        super(DoxygenNode,self).__init__('')
        self.__dict__.update(kwargs)

    def render(self):
        return "The render function of %s has not been overridden" % self.__class__.__name__
