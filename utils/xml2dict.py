class dict2xml(object):
    #doc     = Document()

    def __init__(self, structure):
        self.doc     = Document()
        if len(structure) == 1:
            rootName    = str(structure.keys()[0])
            self.root   = self.doc.createElement(rootName)

            self.doc.appendChild(self.root)
            self.build(self.root, structure[rootName])

    def build(self, father, structure):
        if type(structure) == dict:
            for k in structure:
                #print(k)
                if k[0].islower() : 
                    #print(k)
                    #tag = self.doc.createElement(k)                   
                    #print(tag)
                    #tag.attributes[k]=structure[k]
                    father.attributes[k]=structure[k]
                    
                else:
                    tag = self.doc.createElement(k)                   
                    father.appendChild(tag)
                    self.build(tag, structure[k])

        elif type(structure) == list:
            #print(structure)
            grandFather = father.parentNode
            tagName     = father.tagName
            grandFather.removeChild(father)
            for l in structure:
                tag = self.doc.createElement(tagName)
                #self.doc.ATTRIBUTE_NODE()
                #print(l)
                for atr in l:
                    #grandFather.setAttribute(atr, l[atr])
                    #print(tagName)
                    #tag = self.doc.createElement(tagName)
                    
                    tag.attributes[atr]=l[atr]
                grandFather.appendChild(tag)
                    #tag.__setattr__(atr,l[atr])
                    #print(atr + "=" +l[atr])
                #self.build(tag, l)
                    #grandFather.
                    #grandFather.__setattr__(atr,l[atr])
                #grandFather.appendChild(tag)

        else:
            data    = str(structure)
            tag     = self.doc.createTextNode(data)
            father.appendChild(tag)

    def display(self):
        #print self.doc.toprettyxml(indent="  ")
        return self.doc.toprettyxml(indent="  ")




def etree_to_dict(t):
    d = {t.tag: {} if t.attrib else None}
    children = list(t)
    if children:
        dd = defaultdict(list)
        for dc in map(etree_to_dict, children):
            for k, v in dc.items():
                dd[k].append(v)
        d = {t.tag: {k:v[0] if len(v) == 1 else v for k, v in dd.items()}}
    if t.attrib:
        #d[t.tag].update(('@' + k, v) for k, v in t.attrib.items())
        d[t.tag].update((k, v) for k, v in t.attrib.items())
    if t.text:
        text = t.text.strip()
        if children:
            if text:
              d[t.tag]['#text'] = text
        elif t.attrib:
            d[t.tag]['#text'] = [text]
        else:
            d[t.tag] = text
    return d
