#%% [markdown]
import numpy as np
from math import comb
import copy
import itertools

#%%
class MetricVector:
    def __init__(self, w, W, p=0.5): 
        # w:available channels on link. 
        # W:maximum number of channels a link may have.
        # p:the success probability of generating an entanglement.
        self.W = W
        self.p = p
        self.calc_mv(w,W)

    def __calc_pk(self, w,k):
        if k>w:
            return 0
        return comb(w,k)*(self.p**k)*((1-self.p)**(w-k))
    
    def calc_mv(self, w, W):
        P = np.zeros(W+1)
        for k in range(1, W+1):
            P[k] = self.__calc_pk(w,k)
        P[0] = 1-np.sum(P[1:])
        self.P = P

    def ext(self):
        return np.sum(self.P*np.arange(self.P.shape[0]))

    def enswap(self, mv, q): # entanglement swapping
        if mv == None:
            print('mv None')
            return self
        assert self.W == mv.W, f'metric vector dimension is not same. {self.W}!={mv.W} ' 
        W = self.W
        x = MetricVector(0,W)
        P = np.zeros(W+1)
        for k in range(1, W+1):
            p_xz=0
            for i in range(k, W+1):
                p_xy = self.P[i]
                for j in range(k, W+1):
                    if j < i:
                        p_xz += p_xy * mv.P[j] * comb(j,k)*(q**k)*((1-q)**(j-k))
                    else:
                        p_xz += p_xy * mv.P[j] * comb(i,k)*(q**k)*((1-q)**(i-k))        
            P[k] = p_xz
        P[0] = 1-np.sum(P[1:])
        x.P = P
        return x


class Vertex:
    def __init__(self, w, W, name, p=0.5, q=0.8): 
        # W:maximum number of channels a link may have.
        # p:the success probability of generating an entanglement.
        # q:the success probability of swapping a pair of adjacent.
        assert name != '', 'name value is empty.' 
        self.name = name
        self.w = w
        self.p = p
        self.W = W
        self.q = q


class Graph:
    def __init__(self): 
        self.graph = {}

    def add_vertex(self, v):
        assert v.name not in self.graph.keys(), f'duplicate vertex name.' 
        self.graph[v.name] = v

    def calc_mv_by_name(self, x, y):
        x = self.graph[x]
        y = self.graph[y]
        assert x.W == y.W, f'metric vector dimension is not same. {x.W}!={y.W} ' 
        W = x.W
        if type(x.p) == dict:
            assert y.name in x.p.keys()
            p=x.p[y.name]
        else:
            p=x.p
        
        if type(x.w) == dict:
            if y.name in x.w.keys():
                mv = MetricVector(x.w[y.name], W, p)
            else:
                mv = MetricVector(0, x.W, p)
        else:
                mv = MetricVector(x.w, x.W, p)
        return mv

    def calc_mv_with_enswap(self, x, y ,z):
        # connecting nodes x and z, after performing entanglement swapping at node y
        mv_xy = self.calc_mv_by_name(x,y)
        mv_yz = self.calc_mv_by_name(y,z)
        mv_xz = mv_xy.enswap(mv_yz,self.graph[y].q)
        return mv_xz

    def calc_mv_with_order(self, order):
        g = copy.deepcopy(self.graph)
        linked = {}
        for o in order:
            l = list(g[o].w.keys())
            l0_linked = f'{o}_{l[0]}' in linked.keys() or f'{l[0]}_{o}' in linked.keys()
            l1_linked = f'{o}_{l[1]}' in linked.keys() or f'{l[1]}_{o}' in linked.keys()
            if l0_linked==False and l1_linked==False:
                mv = self.calc_mv_with_enswap(l[0], o, l[1])
                
            elif l0_linked==False and l1_linked==True:
                m = self.calc_mv_by_name(o, l[0])
                p_link = linked[f'{o}_{l[1]}']
                mv = p_link.enswap(m, g[o].q)
            elif l0_linked==True and l1_linked==False:
                m = self.calc_mv_by_name(o, l[1])
                p_link = linked[f'{o}_{l[0]}']
                mv = p_link.enswap(m, g[o].q)
            elif l0_linked==True and l1_linked==True:
                p_link0 = linked[f'{o}_{l[0]}']
                p_link1 = linked[f'{o}_{l[1]}']
                mv = p_link0.enswap(p_link1, g[o].q)

            
            linked[f'{l[0]}_{l[1]}'] = mv
            linked[f'{l[1]}_{l[0]}'] = mv
            g[l[0]].w[l[1]]=-1
            g[l[1]].w[l[0]]=-1
            del g[l[0]].w[o]
            del g[l[1]].w[o]
            del g[o]
        result_key = ""
        for c in list(g.keys()):
            if result_key=="":
                result_key = c
            else:
                result_key += f'_{c}'
        return linked[result_key]
    
#%%
print('\n')
print('   1     2     3     3   ')
print('s-----x-----y-----z-----t')
s = Vertex(w={'x':1}, W=3, name='s', p=0.5, q=0.8)
x = Vertex(w={'s':1, 'y':2}, W=3, name='x', p=0.5, q=0.8)
y = Vertex(w={'x':2, 'z':3}, W=3, name='y', p=0.5,q=0.8)
z = Vertex(w={'y':3, 't':3}, W=3, name='z', p=0.5, q=0.8)
t = Vertex(w={'z':3}, W=3, name='t', p=0.5, q=0.8)

graph =Graph()
graph.add_vertex(s)
graph.add_vertex(x)
graph.add_vertex(y)
graph.add_vertex(z)
graph.add_vertex(t)

order = ['x', 'y', 'z']
mv = graph.calc_mv_with_order(order)
print(f'metric vector: {np.round(mv.P, 3)}', f'EXT: {np.round(mv.ext(), 3)}')

order = ['x', 'z', 'y']
mv = graph.calc_mv_with_order(order)
print(f'metric vector: {np.round(mv.P, 3)}', f'EXT: {np.round(mv.ext(), 3)}')

order = ['y', 'x', 'z']
mv = graph.calc_mv_with_order(order)
print(f'metric vector: {np.round(mv.P, 3)}', f'EXT: {np.round(mv.ext(), 3)}')

order = ['y', 'z', 'x']
mv = graph.calc_mv_with_order(order)
print(f'metric vector: {np.round(mv.P, 3)}', f'EXT: {np.round(mv.ext(), 3)}')

order = ['z', 'y', 'x']
mv = graph.calc_mv_with_order(order)
print(f'metric vector: {np.round(mv.P, 3)}', f'EXT: {np.round(mv.ext(), 3)}')

#%%
#%%
ext = [0]
for order in list(itertools.permutations(['x', 'y', 'z'])):
    mv = graph.calc_mv_with_order(order)
    if mv.ext() > max(ext):
        opt_order = order 
    ext.append(mv.ext())

print(f'opimize_order: {opt_order}', f'EXT: {round(max(ext),3)}')


#%%
print('\n')
print('   3     4     5     5     5     5     4     3')
print('s-----a-----b-----c-----d-----e-----f-----g-----t')
s = Vertex(w={'a':3}, W=5, name='s', p=0.9, q=0.95)
a = Vertex(w={'s':3, 'b':4}, W=5, name='a', p=0.9, q=0.95)
b = Vertex(w={'a':4, 'c':5}, W=5, name='b', p=0.9, q=0.95)
c = Vertex(w={'b':5, 'd':5}, W=5, name='c', p=0.9, q=0.95)
d = Vertex(w={'c':5, 'e':5}, W=5, name='d', p=0.9, q=0.95)
e = Vertex(w={'d':5, 'f':5}, W=5, name='e', p=0.9, q=0.95)
f = Vertex(w={'e':5, 'g':4}, W=5, name='f', p=0.9, q=0.95)
g = Vertex(w={'f':4, 't':3}, W=5, name='g', p=0.9, q=0.95)
t = Vertex(w={'g':3}, W=5, name='t', p=0.9, q=0.95)

graph =Graph()
graph.add_vertex(s)
graph.add_vertex(a)
graph.add_vertex(b)
graph.add_vertex(c)
graph.add_vertex(d)
graph.add_vertex(e)
graph.add_vertex(f)
graph.add_vertex(g)
graph.add_vertex(t)

order = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
mv = graph.calc_mv_with_order(order)
print(f'metric vector: {np.round(mv.P, 3)}', f'EXT: {np.round(mv.ext(), 3)}')


order = ['a', 'c', 'e', 'g', 'b', 'f', 'd']
mv = graph.calc_mv_with_order(order)
print(f'metric vector: {np.round(mv.P, 3)}', f'EXT: {np.round(mv.ext(), 3)}')


order = ['c', 'b', 'e', 'f', 'd', 'g', 'a']
mv = graph.calc_mv_with_order(order)
print(f'metric vector: {np.round(mv.P, 3)}', f'EXT: {np.round(mv.ext(), 3)}')

#%%
ext = [0]
for order in list(itertools.permutations(['a', 'b', 'c', 'd', 'e', 'f', 'g'])):
    mv = graph.calc_mv_with_order(order)
    if mv.ext() > max(ext):
        opt_order = order 
    ext.append(mv.ext())

print(f'opimize_order: {opt_order}', f'EXT: {round(max(ext),3)}')
