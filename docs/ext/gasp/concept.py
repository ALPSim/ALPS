from doxygen import DoxygenNode
from sphinx.util.compat import Directive

class concept(DoxygenNode):
    def __init__(self, name, inherits, **kwargs):
        super(concept,self).__init__(**kwargs)
        self.inherits = inherits
        self.name = name

    def render(self): 
#        template = self.environment.get_template(self.filename+'.html')
#        return template.render(**self.context)
        return "<a name=\"concept_%s\"></a>"% self.name.replace(" ","_")


class ConceptDirective(Directive):
    has_content = True
    def run(self):
        env = self.state.document.settings.env
        try:
            name = self.content[0]
        except:
            raise BaseException("Please specify a unique name for the concept.") 

        try:
            inherits = self.content[1]
        except:
            inherits = None

        return [concept(name, inherits)]


class link_concept(DoxygenNode):
    def __init__(self, variable, scope, name, **kwargs):
        super(link_concept,self).__init__(**kwargs)
        self.variable = variable
        self.scope = scope
        self.name = name

    def render(self): 
#        self.ref_register.set_current_doc(self.doc)
#        template = self.environment.get_template(self.filename+'.html')
        return "" #template.render(**self.context)


class LinkConceptDirective(Directive):

    has_content = True
    def run(self):
        env = self.state.document.settings.env
        try:
            variable = self.content[0]
        except:
            raise BaseException("Please specify a variable/function name.") 

        try:
            scope = self.content[1]
        except:
            raise BaseException("Please specify a scope.") 

        try:
            name = self.content[2]
        except:
            raise BaseException("Please specify a concept name.") 

        return [link_concept(variable, scope, name)]
