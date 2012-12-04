import re
from docutils import nodes
from doxygen import DoxygenNode

_ = lambda x: x
from jinja2 import Environment, PackageLoader


def list_to_dict(lst):
    nl = [q[1:].split(":",1) for q in lst if ":" in q ]
    return dict([(q[0],q[1].strip()) for q in nl])


import copy
class reference_object:
    def __init__(self, **kwargs):
        self.as_string = ""
        self.__dict__.update(kwargs)

    def __str__(self):
        return self.as_string

class doxygen(DoxygenNode):
    def __init__(self, docname, filename, which, **kwargs):
        self.docname = docname
        self.filename = filename
        self.which = which

        if not "exclude" in kwargs or kwargs["exclude"].strip() == "":
            kwargs["exclude"]=[]
        else:
            exclude = [x.strip() for x in kwargs["exclude"].split(",")]
            kwargs["exclude"]=exclude

        for (a,b) in kwargs.iteritems():
            if hasattr(b,"strip") and b.strip() == "":
                kwargs[a] = True

        self.kwargs = kwargs

        super(doxygen,self).__init__(**kwargs)

    def create_references(self, objs):
        refs = []
        for o in objs:
            no = reference_object(**o.attributes)
            setattr(no,"as_string", str(o))
            setattr(no,"as_dict", o.attributes)

            if len(o.children) > 0: setattr(no,"children", self.create_references(o.children))
            setattr(no,"docname", self.docname)
            refs.append(no)
        return refs

    def prepare_node(self,**ctx):
        self.path = self.doxygen.scope_to_path(self.which)
        nodes = self.doxygen.get(self.path)
        if len(nodes) == 0:
            self.context = {}        
            print "WARNING: could not find '%s'"%self.which

            return ""

        if len(self.exclude) == 0:
            objs = nodes
        else:

            objs = []
            for n in nodes:
                if not "name" in n.attributes or not n.attributes["name"] in self.exclude:
                    objs.append(n)
        self.references=self.create_references(objs)
        self.context = copy.copy(self.kwargs)
        self.context["objects"] = self.references
        self.context.update(ctx)

    def render(self): 
        try:
            self.ref_register.set_current_doc(self.doc)
            template = self.environment.get_template(self.filename+'.html')
            return template.render(**self.context)
        except:
            print "Filename: ", self.filename            
            print "XML Path: ", self.XML_PATH
            raise

def visit_doxygen_node(self, node):
    x= node.render()    
    self.body.append(x)

def depart_doxygen_node(self, node):
    pass

from sphinx.util.compat import Directive
class DoxygenDirective(Directive):

    has_content = True
    def run(self):
        env = self.state.document.settings.env
        try:
            filename = self.content[0]
        except:
            raise BaseException("Please specify how you want to render the information, i.e. class/synopsis.") 

        try:
            which = self.content[1]
        except:
            raise BaseException("Please specify a class memener to document.")

        ret =  doxygen( env.docname,filename, which, **list_to_dict(self.content[2:]) )
        return [ret]




def purge_doxygen(app, env, docname):
    pass

#    env.all_doxs = [todo for todo in env.todo_all_todos
#                          if todo['docname'] != docname]


class reference_register:
    def __init__(self):
        self.doc = None
        self.builder = None
        self._register = {}
        self._path_register = {}
        self._request = []
        self._production_mode = False

    def set_current_doc(self, doc):
        self.doc = doc

    def set_builder(self, builder):
        self.builder = builder

    def get_id(self, obj):
        i = None
        if hasattr(obj,"id"):
            i = obj.id
        else:
            print "Warning: object has no id."
        return i

    def production_mode(self):
        self._production_mode = True

    def create_reference(self,obj):
        i = self.get_id(obj)
        if i in self._register and not self._production_mode:
            print "Warning: label for '%s' already exists. Skipping label creation." %i
        self._register[i] = obj.docname
        if hasattr(self,"path"): self._path_register[obj.path] = i
        return "<a name=\"%s\"></a>"% i

    def get_reference(self,obj, text = None):
        i = self.get_id(obj)
        if text is None: 
            text = obj.name
        if not i in self._request: self._request.append(i)
        if not self._production_mode:
            return "<a href=\"javascript:void(0);\">%s</a>"%text
        elif not i in self._register:
            return text
        else:
            uri = self.builder.get_relative_uri(self.doc, self._register[i]) 
            return "<a href=\"%s#%s\">%s</a>"%(uri,i, text)

    def link_scoped_references(self,text, pat="([\w\d]+\<.*\>::)+[\w\d]+\<.+\>"):
        pattern = re.compile(pat)
        newtext = ""
        for ret in pattern.finditer(text):
            pass # TODO
        return text

    def link_doxygen_references(self,text, pat="\[\@ref:(?P<id>[\w\d\_\-]+)=(?P<name>[\w\d\_\-]+)\]"):
        pattern = re.compile(pat)
        newtext = ""
        for ret in pattern.finditer(text):
            pass # TODO
        return text


    def get_or_create_reference(self,obj):
        pass


from doxtools import DoxygenDocumentation
import glob
DOXYGEN_DOC = None
JINJA_ENVIRONMENT = None
REFERENCE_REGISTER = None
CONCEPT_REGISTER = None

from concept import concept,ConceptDirective, link_concept, LinkConceptDirective

class concept_register:

    def __init__(self):
        self.scopes = {}

    def register_type(self, concept, name, type):
        if not concept in self.scopes:
            self.scopes[concept] ={}
        self.scopes[concept][name] = type

    def get_instance_name(self,obj, ifnone = "a"):
        # TODO: implement
        return ifnone

    def simplify_repr(self,name, obj):
        obj = self.get_instance_name(obj)
        ## TODO: generalise
        remove = [("operator()", obj, True),  ("operator+=", obj+"+=",False)]
        has_paran = True
        for bef,aft,par in remove:
            if bef in name:
                name = name.replace(bef, aft)
                has_paran = par
                break
            
        return {"name":name, "has_paranthesis": has_paran}

def process_doxygen(app, doctree, fromdocname):
#    print "XXXX; ", fromdocname
    global DOXYGEN_DOC, JINJA_ENVIRONMENT, REFERENCE_REGISTER, CONCEPT_REGISTER
    if app.config.doxygen_xml is None:
        raise BaseException("Please specify the path to the Doxygen XML in the conf.py using the variable 'doxygen_xml'.")

    env = app.builder.env
    XML_PATH = app.config.doxygen_xml
    if not XML_PATH.endswith("/"): XML_PATH+="/"
    testfiles = glob.glob("%s*.xml"%XML_PATH)

    if DOXYGEN_DOC is None:
        print "Loading Doxygen documentation ... "
        DOXYGEN_DOC = DoxygenDocumentation(testfiles)
        print "Generating Doxygen output"
    dox = DOXYGEN_DOC

    if JINJA_ENVIRONMENT is None:
        JINJA_ENVIRONMENT =Environment(loader=PackageLoader('gasp', 'templates'))
    jenv = JINJA_ENVIRONMENT

    if REFERENCE_REGISTER is None:
        REFERENCE_REGISTER = reference_register()

    if CONCEPT_REGISTER is None:
        CONCEPT_REGISTER = concept_register()

        
    reg = REFERENCE_REGISTER
    reg.set_builder(app.builder)


    def lst(obj):
        s = ""
        for o in dir(obj):
            s+="<li>%s</li>" % o
        return "<ul>%s</ul>" % s

    jenv.filters['get_abstract_object_instance'] = CONCEPT_REGISTER.get_instance_name
    jenv.globals['simplify_member_representation'] = CONCEPT_REGISTER.simplify_repr
    jenv.filters['list_contents'] = lst
    jenv.filters['ref'] = reg.get_reference
    jenv.filters['label'] = reg.create_reference
    jenv.filters['link_refs'] = reg.link_doxygen_references

 #   uri = app.builder.get_relative_uri()
    for node in doctree.traverse(doxygen):
        node.XML_PATH = XML_PATH
        node.ref_register = reg
        node.doxygen = dox
        node.environment = jenv
        node.doc = fromdocname
        node.prepare_node()
        node.render()
#        node.relative_uri = uri
    reg.production_mode()



def setup(app):
    app.add_config_value('doxygen_xml', None, True)

    app.add_node(doxygen,
                 html=(visit_doxygen_node, depart_doxygen_node))
    app.add_node(concept,
                 html=(visit_doxygen_node, depart_doxygen_node))
    app.add_node(link_concept,
                 html=(visit_doxygen_node, depart_doxygen_node))


    app.add_directive('doxygen', DoxygenDirective)
    app.add_directive('concept', ConceptDirective)
    app.add_directive('link-concept', LinkConceptDirective)


    app.connect('doctree-resolved', process_doxygen)
    app.connect('env-purge-doc', purge_doxygen)


