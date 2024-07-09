#Updated maximal faces to reduce

import re
import numpy as np
import sympy as sym
import copy
import itertools

def get_ordered_vertices(a):
    c=str(a).replace('**','^').replace('*','')
    factors=list(re.finditer('\(1,.\)',c)); f1,f2=factors

    letter1,letter2=[re.findall('['+''.join(Algebra.Symbols)+']',c[f.start():f.end()])[0] for f in factors]
    clist=[c[:f1.start()],c[f1.start():f1.end()],c[f1.end():f2.start()],c[f2.start():f2.end()],c[f2.end():]]
    v0=Cell(clist[0]+''+clist[2]+''+clist[4])
    v1=Cell(clist[0]+letter1+clist[2]+''+clist[4])
    v2=Cell(clist[0]+letter1+clist[2]+letter2+clist[4])
    v3=Cell(clist[0]+''+clist[2]+letter2+clist[4])
    return [str(v0),str(v1),str(v2),str(v3),str(v0)]

def recenter(W,pos,prev_bound):
    coords=[pos[k] for k in W]
    centroid=np.mean(coords,axis=0)
    r=np.sqrt(np.max(np.sum((np.array(coords)-centroid)**2,axis=1)))
    i=np.argmax(np.sum((np.array(coords)-centroid)**2,axis=1))
    pos.update({k:pos[k]-centroid+prev_bound+r+1 for k in W})
    return pos, pos[W[i]]

v1=np.array([1,0,-np.sqrt(2)/2])
v2=np.array([-1,0,-np.sqrt(2)/2])

v3=np.array([0,1,np.sqrt(2)/2])
v4=np.array([0,-1,np.sqrt(2)/2])

fixed_vertices=['1','A','C','G','T']+[letter+'**'+str(d) for letter in ['A','T','C','G'] for d in range(2,6)]
pos0={'1':(0,0,0), 'A':v1, 'C':v2,'G':v3,'T':v4}
pos0b={letter+'**'+str(d):d*pos0[letter] for letter in ['A','T','C','G'] for d in range(2,6)}
pos0.update(pos0b)


class Algebra:
    Symbols=['A','C','G','T']
    Sigma={}; Epsilon={}
    for s in Symbols:
        Sigma[s]= sym.Symbol(s,commutative=False)
        Epsilon[Sigma[s]]= sym.Symbol('(1,_)'.replace('_',s),commutative=False)
    One=sym.core.numbers.One()
    Zero=sym.core.numbers.Zero()

    @staticmethod
    def update_symbols(Symbols):
        Algebra.Symbols=Symbols
        Algebra.Sigma={s:sym.Symbol(s,commutative=False) for s in Symbols}
        Algebra.Epsilon={Algebra.Sigma[s]:sym.Symbol('(1,_)'.replace('_',s),commutative=False) for s in Symbols}

def convert_coeff(string):
    if string in ['','+','-']:
        string+='1'
    return int(string)

class Chain(Algebra):
    def __init__(self,expression=None,*,list_coeffs=None,list_cells=None,dict_cells=None):
        if not(dict_cells is None):
            self.dict_cells=dict_cells
            self.update_expression()
        else:
            if isinstance(expression,str):
                list_coeffs,list_cells=self.str_constructor(expression)
                list_cells=[Cell(c) for c in list_cells]

            #Create dictionary
            self.dict_cells={}
            for i,c in enumerate(list_cells):
                if not c in self.dict_cells:
                    self.dict_cells[c]=list_coeffs[i]
                else:
                    self.dict_cells[c]+=list_coeffs[i]
                    if self.dict_cells[c]==0:
                        del self.dict_cells[c]
            self.update_expression()
           
    def update_expression(self):
        #Create expression
        self.expression=self.Zero+sum([self.dict_cells[c]*c.expression for c in self.dict_cells])
 
    
    def str_constructor(self,expression):
        expression=re.sub('\s+','',expression)
        split_pattern=r'((?:^\+?\-?\d*)|(?:\+\d*)|(?:\-\d*))'
        terms_list=re.split(split_pattern,expression)
        terms_list.pop(0)
        list_coeffs=[convert_coeff(c) for c in terms_list[::2]]
        list_cells=terms_list[1::2]
        return list_coeffs, list_cells
    
    def to_chain(self):
        return self
    
    def _repr_latex_(self):
        return self.expression._repr_latex_()
    
    def __add__(self,other):
        new_dict=self.dict_cells
        other=other.to_chain()
        for c in other.dict_cells:
            if c in new_dict:
                new_dict[c]+=other.dict_cells[c]
                if new_dict[c]==0:
                    del new_dict[c]
            else:
                new_dict[c]=other.dict_cells[c]
        return Chain(dict_cells=new_dict)
    
    def __rmul__(self,coeff):
        if coeff!=0:
            new_dict={c:coeff*self.dict_cells[c] for c in self.dict_cells}
        else:
            new_dict={}
        return Chain(dict_cells=new_dict)
    
    def __sub__(self,other):
        return self+(-1)*other
    
    def boundary(self):
        new_dict={}
        for c in self.dict_cells:
            dict_boundary=c.boundary(c0=self.dict_cells[c]).dict_cells
            for new_c in dict_boundary:
                if new_c in new_dict:
                    new_dict[new_c]+=dict_boundary[new_c]
                else:
                    new_dict[new_c]=dict_boundary[new_c]
        #Delete zeros
        keys_new_dict=list(new_dict.keys())
        for c in keys_new_dict:
            if new_dict[c]==0:
                del new_dict[c]
        return Chain(dict_cells=new_dict)
    
    def get_cells(self):
        return list(self.to_chain().dict_cells.keys())       
    
    

class Cell(Chain):    
    def __init__(self,expression=None,*,alphas=[Algebra.Zero],edges=[],extended_word=None,complete_sequence=None, subindices=None):
        if isinstance(expression, str):
            if expression=='1':
                alphas=[self.One]; edges=[]
            else:
                alphas,edges=self.str_constructor(expression)
        self.list_constructor(alphas=alphas,edges=edges)
        self.extended_word=extended_word
        self.complete_sequence=complete_sequence
        self.subindices=subindices
        
    def list_constructor(self,alphas,edges):    
        if len(alphas)==0:
            self.alphas=[self.Zero]
            self.edges=[]
            self.dim=-1
            
        elif len(alphas)-len(edges)!=1:
            raise Exception('WrongLengths')
        
        #Commute Terms
        self.alphas=alphas
        self.edges=edges
        self.dim=len(edges)
        self.commute_terms()
        
        #Check is a valid cell
        if not self.is_valid_cell():
            self.alphas=[self.Zero]
            self.edges=[]
            self.dim=-1
        
        #Expression
        Edges=[self.Epsilon[a] for a in self.edges]+[self.One]
        self.expression=sym.prod([val for pair in zip(self.alphas, Edges) for val in pair])

    def commute_terms(self):
        for k in range(self.dim,0,-1):
            factors=list(self.alphas[k].as_coeff_mul()[1])
            try:
                first_factor=factors.pop(0) 
            except IndexError:
                first_factor=Algebra.One
            if first_factor.as_base_exp()[0]==self.edges[k-1]:
                self.alphas[k-1]=self.alphas[k-1]*first_factor
                self.alphas[k]=sym.Mul(sym.prod(factors))
    
    def is_valid_cell(self):
        if self.dim>0:
            position_ones=[i for i, s in enumerate(self.alphas) if s==self.One]
            for k in position_ones:
                if 0<k<self.dim and self.edges[k]==self.edges[k-1]:
                    return False
        elif self.dim==0:
            if self.alphas[0]==self.Zero:
                return False
        return True              
    
    def str_constructor(self,expression):
        expression=expression.replace('**',r'^').replace('*','')
        symbs='['+''.join(self.Symbols)+']'
        pattern=r'((?:\(1,'+symbs+'\)(?:\^\d+)?)|(?:'+symbs+'(?:\^\d+)?))'
        list_factors=re.split(pattern,expression)
        list_factors.pop(0)

        pattern_alpha=r'(?P<base>'+symbs+')(?:\^(?P<power>\d+))?'
        re_alpha=re.compile(pattern_alpha)

        pattern_edge=r'\(1,(?P<base>'+symbs+')\)(?:\^(?P<power>\d+))?'
        re_edge=re.compile(pattern_edge)

        list_factors=[f for f in list_factors if not(f=='')]

        degree_zero=True
        complete_list=[]
        current_factor=self.One

        alphas=[]
        edges=[]
        for factor in list_factors:
            if degree_zero:
                base=re_alpha.match(factor)
                if base is None:
                    alphas.append(current_factor)
                    current_factor=self.One
                    degree_zero=False
                else:
                    if base['power'] is None:
                        power=1
                    else:
                        power=int(base['power'])
                    current_factor*=self.Sigma[base['base']]**power
            if not(degree_zero):
                edge=re_edge.match(factor)
                if not(edge['power'] is None):
                    return [self.Zero],[]
                edges.append(self.Sigma[edge['base']])
                degree_zero=True
        alphas.append(current_factor)
        return alphas, edges
    
    def to_chain(self):
        if self.dim>=0:
            return Chain(dict_cells={self:1})
        else:
            return Chain(dict_cells={})
   
    def __hash__(self):
        return self.expression.__hash__()
    
    def __eq__(self,other):
        return self.expression==other.expression
    
    def __rmul__(self,coeff):
        return Chain(dict_cells={self:coeff})
    
    def __add__(self,other):
        return self.to_chain()+other.to_chain()   

    def boundary(self,c0=1):
        collapsing_inds=[]
        for i in range(1,self.dim-1):
            if self.edges[i-1]==self.edges[i+1]:
                middle_factor=self.edges[i-1]*self.alphas[i]*self.alphas[i+1]
                if len(middle_factor.as_coeff_mul()[1])<=1:
                    collapsing_inds.append(i)
        regular_inds=[i for i in range(self.dim) if not i in collapsing_inds]

        alphas=self.alphas.copy(); edges=self.edges.copy()
        result={}
        for i in regular_inds:
            coeff=c0*(-1)**(i)
            edges_d=(edges[:i]+edges[(i+1):])
            alphas_d0=(alphas[:i]+[alphas[i]*alphas[i+1]]+alphas[(i+2):])
            alphas_d1=(alphas[:i]+[alphas[i]*edges[i]*alphas[i+1]]+alphas[(i+2):])
            
            result[Cell(alphas=alphas_d1,edges=edges_d)]=coeff
            result[Cell(alphas=alphas_d0,edges=edges_d)]=-coeff

        for i in collapsing_inds:
            coeff=c0*(-1)**(i)
            edges_d=(edges[:i]+edges[(i+1):])
            alphas_d1=(alphas[:i]+[alphas[i]*edges[i]*alphas[i+1]]+alphas[(i+2):])

            result[Cell(alphas=alphas_d1,edges=edges_d)]=coeff
        return Chain(dict_cells=result)
      
    def __str__(self):
        return self.expression.__str__()
    
    def __repr__(self):
        return self.expression.__str__()
    
    def get_facets(self):
        return self.boundary().get_cells()
    
    def get_all_faces(self,include_self=False):
        if include_self:
            all_faces=set([self])
        else:
            all_faces=set()
        F=[self]
        for _ in range(self.dim):
            F=[f for c in F for f in c.get_facets()]
            all_faces=all_faces.union(F)
        all_faces=list(all_faces)
        all_faces.sort(key=lambda x:x.dim)
        return all_faces

def get_word_dimension(w):
    try:
        if w.expression==Algebra.One:
            return 0
    except:
        if w==Algebra.One:
            return 0
    return len(str(w).replace('**',r'^').split('*'))

def check_edge(w1,w1_str,w2_str):
    for i in range(len(w1_str)):
        if w2_str==w1_str[:i]+w1_str[i+1:]:
            exps=np.cumsum([b.as_base_exp()[1] for b in w1.expression.as_ordered_factors()])
            j=np.argwhere((exps-i-1)>=0).min()
            return True,j
    return False, None

class Filtration(Algebra):
    def __init__(self,name=None,id_exp=None,*,Symbols=Algebra.Symbols):
        self.id=id_exp
        self.name=name
        self.Symbols=Symbols
        self.Sigma={}; self.Epsilon={}; self.Doubles={}
        for s in Symbols:
            self.Sigma[s]= sym.Symbol(s,commutative=False)
            self.Epsilon[self.Sigma[s]]= sym.Symbol('(1,_)'.replace('_',s),commutative=False)
            self.Doubles[s+s]=self.Epsilon[self.Sigma[s]]
        self.One=sym.core.numbers.One()
        self.Zero=sym.core.numbers.Zero()
        self.dim=-1
        self.Complex_dict={}
        self.Boundary_dict={}
        self.filtration_values=[]
        self.indexed_ordered_pairs={}
        self.cells_indexes={}
        
        self.Incidence_dict={}
        self.root=None
    
    def set_root(self,root):
        self.root=Cell(root)
        
        
    def convert_word(self,w):
        return sym.prod([self.Sigma[s] for s in w],start=self.One)
    
    def maximal_faces_word(self,w,indices,dim):
        l=len(w.as_coeff_mul()[1])
        faces=set()
        for subindices in itertools.combinations(indices,dim):
            product=self.One
            for ind,pair in enumerate(w.as_coeff_mul()[1]):
                base,exp=pair.as_base_exp()
                if ind in subindices:
                    if exp==1:
                        product*=self.Epsilon[base]
                    else:
                        factor=base**(exp-1)*self.Epsilon[base]
                        product*=factor
                else:
                    product*=base**exp
            cell=Cell(str(product),complete_sequence=w,subindices=list(subindices))       
            faces.add(cell)
        result=list(faces)
        if Cell(alphas=[0]) in result:
            result.remove(Cell(alphas=[0]))
        return result
    
    def get_boundary_cell(self,cell):
        if cell in self.Boundary_dict:
            return self.Boundary_dict[cell].get_cells()
        else:
            B=cell.boundary()
            self.Boundary_dict[cell]=B
            return B.get_cells()
    
    def compute_zero_skeleton(self,list_words,list_weights,already_cells=False):
        #Assume they come as words
        if not(already_cells):
            list_zero_cells=[Cell(w,complete_sequence=self.convert_word(w),extended_word=w) for w in list_words]
        else:
            list_zero_cells=list_words
            
        self.Complex_dict[0]=dict(zip(list_zero_cells,list_weights))

        self.Incidence_dict[0]=dict()

        list_zero_cells.sort(key=lambda x: len(x.extended_word),reverse=True)
        for i,w1 in enumerate(list_zero_cells):
            w1_str=w1.extended_word
            l1=len(w1_str)
            for w2 in list_zero_cells[(i+1):]:
                w2_str=w2.extended_word
                l2=len(w2_str)
                if l1-l2==1:
                    is_edge,j=check_edge(w1,w1_str,w2_str)
                    if is_edge:                
                        if w1.complete_sequence in self.Incidence_dict[0]:
                            self.Incidence_dict[0][w1.complete_sequence].add(j)
                        else:
                            self.Incidence_dict[0][w1.complete_sequence]=set([j])
                elif l1-l2>1:
                    break
    
    def add_faces(self,cell):
        d=cell.dim
        if d in self.Complex_dict and cell in self.Complex_dict[d]:
            return True, self.Complex_dict[d][cell]
        if d==0:
            if (0 in self.Complex_dict and cell in self.Complex_dict[0]):
                return True, self.Complex_dict[0][cell]
            else:
                return False,0
     
        F=self.get_boundary_cell(cell)
        subfaces_added=[self.add_faces(f) for f in F]
        
        bool_values=[i for i,j in subfaces_added]
        weights=[j for i,j in subfaces_added]
        
        if all(bool_values):
            weight=max(weights)
            w=cell.complete_sequence
            
            if w in self.Incidence_dict[d]:
                indices= self.Incidence_dict[d][w]
                indices=indices.union(cell.subindices)
            else:
                indices=set(cell.subindices)
            self.Incidence_dict[d][w]=indices
            
            if d in self.Complex_dict:            
                self.Complex_dict[d][cell]=weight
            else:
                self.Complex_dict[d]={cell:weight}
            return True, weight
        else:
            return False, 0   
        
    def __repr__(self):
        self.Complex_dict
    
    
    def compute_d_skeleton(self,W,weights,max_dim=10,already_cells=False):
        self.filtration_values=list(weights)
        self.filtration_values.sort()
        self.compute_zero_skeleton(W,weights,already_cells)

        dim=0
        while (dim <max_dim) and (dim in self.Incidence_dict):
            dim+=1
            if not(dim in self.Incidence_dict):
                self.Incidence_dict[dim]=dict()
            for w in self.Incidence_dict[dim-1]:
                indices=self.Incidence_dict[dim-1][w]
                MF=self.maximal_faces_word(w,indices,dim)
                for M in MF:
                    self.add_faces(M)
        self.dim=max(self.Complex_dict)
        self.index_cells() 
    
    def index_cells(self):
        L=[c for d in self.Complex_dict for c in self.Complex_dict[d].items()]
        L.sort(key=lambda c:(c[1],c[0].dim))
        self.indexed_ordered_pairs=dict(enumerate(L))
        self.cells_indexes={pair[0]:i for i,pair in enumerate(L)}
      
    
    def filtrated_complex(self,height):
        Cells={}
        for d in self.Complex_dict:
            for c in self.Complex_dict[d]:
                if self.Complex_dict[d][c]<=height:
                    if d in Cells:
                        Cells[d].add(c)
                    else:
                        Cells[d]=set([c])
        return(Cells)
    
    def get_Maximal_Faces(self,height,reduced=False):
        Cells=self.filtrated_complex(height)
        Max_Faces=copy.deepcopy(Cells)
        D=max(Cells)
        for d in range(D, 0,-1):
            for f in Cells[d]:
                Max_Faces[d-1]=Max_Faces[d-1].difference(self.get_boundary_cell(f))
        if reduced:
            return self.reduce_maximal_faces(Max_Faces)
        else:    
            return Max_Faces
    
    def reduce_maximal_faces(K,max_dict):
        Vertices_dict={}
        for d in max_dict:
            for c in max_dict[d]:
                Vertices_dict[c]=set([v for v in c.get_all_faces(True) if v.dim==0])
        max_dim=max(max_dict)
        for d in range(max_dim):
            to_remove=set()
            for c1 in max_dict[d]:
                finish=False
                for k in range(d+1,max_dim+1):
                    for c2 in max_dict[k]:
                        if Vertices_dict[c1].issubset(Vertices_dict[c2]):
                            to_remove.add(c1)
                            finish=True
                            break
                    if finish:
                        break
            if len(list(to_remove))>0:
                print('Removed')
                max_dict[d]=max_dict[d].difference(to_remove)
        return max_dict
    
    
    def matlab_constructor_cell(self,id_cell):
        cell=self.indexed_ordered_pairs[id_cell][0]
        
        string_0='cell_{id_cell}=Cell();\n'+\
        'stream.addElement(cell_{id_cell},{val});'
        
        string_d='cell_{id_cell}=Cell({dim},{list_boundary});\n'+\
        'stream.addElement(cell_{id_cell},{val});'

        if cell.dim >0:
            boundary=self.Boundary_dict[cell].dict_cells
            n=len(boundary)
            interleaving_ones=[-1,1]*int(n/2)+[-1]*(n%2)

            cells=list(boundary.keys())
            coefs=list(boundary.values())
            id_terms=[self.cells_indexes[c] for c in cells]

            lista=[]
            for i in interleaving_ones:
                id_c=coefs.index(i)
                coefs.pop(id_c)
                lista.append(id_terms.pop(id_c))
            
            list_boundary=str(['cell_'+str(id_c) for id_c in lista]).replace("'","")
            
            return string_d.format(
                    id_cell=id_cell, dim=cell.dim,list_boundary=list_boundary,
                    val=self.indexed_ordered_pairs[id_cell][1])                 
        else:
            return string_0.format(id_cell=id_cell,val=self.indexed_ordered_pairs[id_cell][1])
       
    def save_matlab_code(self,file_code,file_results):
        with open(file_code,'w') as matlab_file:
            matlab_file.write('load_javaplex\n')
            matlab_file.write('import edu.stanford.math.plex4.*;\n')
            matlab_file.write('import edu.stanford.math.plex4.homology.chain_basis.*;\n')
            matlab_file.write('\nstream=api.Plex4.createExplicitCellStream(1.01);\n')
            matlab_file.write('\n\n %----------------------------------\n\n')

            matlab_file.write('\n'.join(
                [self.matlab_constructor_cell(id_cell) for id_cell in self.indexed_ordered_pairs]))
            
            matlab_file.write('stream.finalizeStream();')
            matlab_file.write('persistence=api.Plex4.getDefaultCellularAlgorithm(3);\n')
            matlab_file.write('intervals=persistence.computeIntervals(stream)\n')
            matlab_file.write('intervals0=fix_intervals(intervals.getIntervalsAtDimension(0));\n')
            matlab_file.write('intervals1=fix_intervals(intervals.getIntervalsAtDimension(1));\n')
            matlab_file.write('intervals2=fix_intervals(intervals.getIntervalsAtDimension(2));\n')

            matlab_file.write('save("'+file_results+'","intervals0","intervals1","intervals2");')
            
    def get_graph(self,show_labels=True,dim_3=True,values=None):
        from sage.all import Graph,text3d, polygon, animate, show
        import networkx as nx
        
        print(self.id,': ',self.name)
        
        if values is None:
            values=self.filtration_values

        #Initial pos
        Weights_by_dim=[1,1,2,2,4,4]
        Cells=self.filtrated_complex(values[-1])

        Edges0=[]
        Edges=[]
        max_dim=max(Cells)
        for dim in range(max_dim,0,-1):
            for cell in Cells[dim]:
                all_cells=cell.get_all_faces(True)
                one_cells=[f for f in all_cells if f.dim==1]
                for a in one_cells:
                    b=list(map(str,a.get_facets()))+[Weights_by_dim[dim]]
                    Edges.append(b)
                    Edges0.append(b)

                two_cells=[f for f in all_cells if f.dim==2]
                for a in two_cells:
                    vv=get_ordered_vertices(a)
                    b1=[str(vv[0]),str(vv[2]),0]
                    b2=[str(vv[1]),str(vv[3]),0]
                    Edges0.append(b1)
                    Edges0.append(b2)

        Grafo0=Graph([list(map(str,Cells[0])),Edges0],format='vertices_and_edges')
        Grafo=Graph([list(map(str,Cells[0])),Edges],format='vertices_and_edges')

        Grafo_nx=Grafo0.networkx_graph()
        pos1={v:pos0[v] for v in pos0 if v in Grafo0.vertices()}
        if len(pos1)>0:
            pos = nx.spring_layout(Grafo_nx, iterations=500,dim=3,scale=5,pos=pos1,
                                   fixed=[v for v in fixed_vertices if v in Grafo0.vertices()])
        else:
            pos = nx.spring_layout(Grafo_nx, iterations=500,dim=3,scale=5)
            

        #Make positions reasonable
        Components=Grafo.connected_components()
        #Most external point
        C0=Components[0]
        coords=[pos[k] for k in C0]
        centroid=np.mean(coords,axis=0)
        r=np.sqrt(np.max(np.sum((np.array(coords)-centroid)**2,axis=1)))
        i=np.argmax(np.sum((np.array(coords)-centroid)**2,axis=1))
        prev_bound=coords[i]
        for num in range(1,len(Components)):
            pos,prev_bound=recenter(Components[num],pos,prev_bound)

        PPs=[]

        for val in values:
            Cells=self.get_Maximal_Faces(val)
            Vv=list(map(str,self.filtrated_complex(val)[0]))
            Hrafo=Grafo.subgraph(Vv)
            max_dim=max(Cells)

            PP=Hrafo.plot3d(pos3d={v:pos[v] for v in Vv},vertex_size=.1,vertex_labels=False)
            if show_labels:
                for v in Vv:
                    PP+=text3d(v.replace('*',''),pos[v],size=4)

            Colors_by_dim=[0,0,'yellow','red','blue','purple']
            for dim in range(max_dim,1,-1):
                for cell in Cells[dim]:
                    two_cells=[f for f in cell.get_all_faces(True) if f.dim==2]
                    for a in two_cells:
                        vv=map(str,get_ordered_vertices(a))

                        PP+=polygon([pos[a] for a in vv],alpha=0.2,axes=False,
                                    color=Colors_by_dim[dim],vertex_size=.1).plot()
            #PP+=text3d("{:.2e}".format(val),(0,10,-10),fontsize=20)
            PPs.append(PP)
        plot = animate(PPs,).interactive(online=True, delay=100,loop=False)
        show(plot, delay=5, auto_play=False, projection='orthographic')