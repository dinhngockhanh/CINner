import csv
import numpy as np
import pandas as pd
import scipy
import sys
import os
import pathos.multiprocessing as mp
from gurobipy import *

########################################################################
############################ WCND functions ############################
########################################################################

def compress_cn_mat(mat):
    n,m = mat.shape
    y = mat[:,1:] != mat[:,:m-1]
    x = np.any(y,axis=0)
    x = np.hstack((True,x))
    return mat[:,x]

def validateProfiles(u,v):
    return not np.logical_and(np.array(u)==0,np.array(v)>0).any()

def print_events_summaries(e,name):
    x = np.clip(np.hstack((e[0],e[1:]-e[:-1],0)),0,np.inf)
    y = np.clip(np.hstack((0,e[:-1]-e[1:],e[-1])),0,np.inf)
    print(name+' starts',list(zip(list(np.where(x>0)[0]),list(x[x>0]))))
    print(name+' ends',list(zip(list(np.where(y>0)[0]),list(y[y>0]))))

def directed_cnd_prefix(u,v):
    if len(v)!=len(u):
        raise ValueError("Length doesn't match")
    u,v=np.array(u),np.array(v)
    u,v=u[np.logical_not(np.logical_and(u==0,v==0))],v[np.logical_not(np.logical_and(u==0,v==0))]
    n = len(u)
    if n==0:
        return u,v,n,0
    M = max(max(v),max(u))
    return u,v,n,M

def calc_Q(u,v):
    n=len(u)
    u_m = 0
    prev = -1
    Q = {}
    for i in range(n):
        if v[i]==0:
            u_m = max(u_m,u[i])
        else:
            Q[i] = (u_m,prev)
            prev = i
            u_m = 0
    Q[n] = (u_m,prev)
    return Q

def CalcPmin(ui,vi):
    return max(ui-vi,0)

def CalcPmax(ui):
    return max(ui-1,0)

def CalcM(ui,vi,ui1,vi1):
    return (vi-ui)-(vi1-ui1)

def limit(x,down,up):
    if x<down:
        return down
    if x>up:
        return up
    return x

class MyFunc:
    """
    Representation of the pairwise linear function of used in the linear time algorithm for the unweighted copy number distance.
    """
    def __init__(self,ui,vi,a,b,base):
        self.ui = ui
        self.vi = vi
        self.pmin = CalcPmin(ui,vi)
        self.pmax = CalcPmax(ui)
        self.a = a
        self.b = b
        self.base = base
    def CalcP(self,p):
        if p<self.pmin or p>self.pmax:
            raise("Error - unvalid p")
        if p<=self.a:
            return self.base
        if p<=self.b:
            return self.base + p - self.a
        if p<=self.pmax:
            return self.base - self.b - self.a + 2*p
    def CalcNext(self,ui,vi,Qi):
        Mi = CalcM(ui,vi,self.ui,self.vi)
        pmin = CalcPmin(ui,vi)
        pmax = CalcPmax(ui)
        if Mi>=0:
            if Qi<=self.a:
                nextBase = self.base
                nextA = self.a-Mi
                nextB = self.b
            elif Qi<=self.b:
                nextBase = self.base + Qi - self.a
                nextA = Qi-Mi
                nextB = self.b
            elif Qi>=self.b:
                nextBase = self.base + Qi - self.a
                nextA = self.b-Mi
                nextB = Qi
        elif Mi<=0:
            if Qi<=self.a:
                nextBase = self.base
                nextA = self.a
                nextB = self.b-Mi
            elif Qi<=self.b:
                nextBase = self.base + Qi - self.a
                nextA = Qi
                nextB = self.b-Mi
            elif Qi>=self.b:
                nextBase = self.base + Qi - self.a
                nextA = min(self.b-Mi,Qi)
                nextB = max(Qi,self.b-Mi)
        if pmin > nextA and pmin<=nextB:
            nextBase = nextBase + pmin - nextA
        if pmin > nextB:
            nextBase = nextBase - nextB - nextA + 2*pmin
        nextA = limit(nextA,pmin,pmax)
        nextB = limit(nextB,nextA,pmax)
        return MyFunc(ui,vi,nextA,nextB,nextBase)

def DirectedCopyNumberDistanceLinear(u,v):
    """
    Calculates the unweighted Copy Number Distance with the linear time algorithm

    Parameters
    ----------
    u : list/np.array of integers
    v : list/np.array of integers
    Returns
    -------
    int
        The unweighted copy number distance from u to v, i.e., the number of a amplifications and deletions in a shortest transformation from u to v

    """
    u,v,n,M = directed_cnd_prefix(u,v)
    if len(v)!=len(u):
        raise ValueError("Length doesn't match")
    if n==0 or (u==v).all():
        return 0
    else:
        u,v=compress_cn_mat(np.array([u,v]))
        n=len(u)
    Q = calc_Q(u,v)
    prevFunc = MyFunc(M+1,M+1,0,0,0)
    for i in range(n):
        if v[i]>0:
            prevFunc = prevFunc.CalcNext(u[i],v[i],Q[i][0])
    u_m,prev = Q[n]
    if prev<0:
        return u_m
    p_min = CalcPmin(u[prev],v[prev])
    d = prevFunc.CalcP(p_min)+max(u_m-p_min,0)
    for p in range(int(p_min)+1,int(prevFunc.pmax)+1):
        d = min(d,prevFunc.CalcP(p)+max(u_m-p,0))
    return d

EuclideanDistance=lambda x,y: np.sqrt(np.sum((np.array(x)-np.array(y))**2))
LOneDistance=lambda x,y: np.sum(np.abs(np.array(x)-np.array(y)))

################
# WCND functions
################

def op_to_subkind(i,j,n):
    if i==0 and j==n:
        return 'whole'
    elif j-i==1:
        return 'small'
    elif i==0 or j==n:
        return 'arm'
    else:
        return 'segmental'

def weigh_ops_with_prob(stat):
    return lambda s,i,j,n: -np.log(stat.loc[s,op_to_subkind(i,j,n)])

def cn_breakpoints(vec):
    return [0]+list(np.where(vec[1:] != vec[:-1])[0]+1)+[len(vec)]

def WeightedDistance(u,v,debug=False,weight_op=lambda s,i,j,n: 1,wgd=0,min_ops=False):
    """
    Calculates the Weighted Copy Number Distance using a linear programming formulation

    Parameters
    ----------
    u : list/np.array of integers
    v : list/np.array of integers
    weight_op: a function for weiging different events.
                The fucntion should get 4 paramaters s,i,j,n:
                        s is the type of the event,
                        i is the start index of the event,
                        j is the end index of the event (excluding j),
                        n is the length of the profile
    min_ops: boolean - if False finds the minimum weighted transformation,
                        if True finds the minimum weighted transformation among the shortest transfomations
    wgd: int - forces the transformation to have at least wgd whole genome duplications
    debug: boolean
    Returns
    -------
    a tuple (d,ops)
        d - The weight of the minimum weighted transformation from u to v
        ops - a dictionary: maps each phase to a list of events,
                            each event is a tuple (i,j,amplitude) representing an event from index i to index j (excluding) with an amplitude
    """
    if not validateProfiles(u,v):
        return np.inf,{}
    u,v,n,M = np.array(u),np.array(v),len(u),max(max(v),max(u))
    if (u==v).all():
        return 0,{},pd.DataFrame(columns=['sign','start','end'])
    model = Model("WeightedDistance")
    d,x = {},{}
    signs = ['- before','+','- after']
    for i in range(n):
        for s in signs:
            d[s,i] = model.addVar(lb=0, ub=M, obj=0)
    breakpoints = sorted(list(set(cn_breakpoints(u)).union(cn_breakpoints(v))))
    pairs = [(i,j) for ind,i in enumerate(breakpoints[:-1]) for j in breakpoints[ind+1:]]
    for i,j in pairs:
        for s in signs:
            x[s,i,j] = model.addVar(lb=0, ub=M, obj=0)
    model.setObjective(quicksum(weight_op(s,i,j,n)*x[s,i,j] for i,j in pairs for s in signs),GRB.MINIMIZE)
    model.update()
    for k in range(n):
        for s in signs:
            model.addConstr(d[s,k]==quicksum(x[s,i,j] for i,j in pairs if i<=k and k<j))
    for i in range(n):
        if v[i]==0:
            model.addConstr(u[i],'<=',d['- before',i])
        else:
            model.addConstr(u[i]-d['- before',i]+d['+',i]-d['- after',i]==v[i])
            model.addConstr(d['- before',i]<=u[i]-1)
    if wgd>0:
        model.addConstr(x['+',0,n]>=wgd)
    if min_ops:
        min_len = DirectedCopyNumberDistanceLinear(u,v)
        model.addConstr(quicksum(x[s,i,j] for i,j in pairs for s in signs)<=min_len)
    model.update()
    model.setParam('OutputFlag', False )
    model.optimize()
    if model.status == GRB.status.INFEASIBLE:
        print(model.ModelName,"Infeasible!!!")
        return np.inf
    if debug:
        for s in signs:
            print('d',s,[(i,d[s,i].x) for i in range(n) if d[s,i].x>0])
            print('x',s,[(i,j,x[s,i,j].x) for i,j in pairs if x[s,i,j].x>0])
    ####################################################################
    df_focal = pd.DataFrame()
    df_focal['sign'] = []
    df_focal['start'] = []
    df_focal['end'] = []
    for s in signs:
        if s == '- before':
            s_ID = '-'
        elif s == '- after':
            s_ID = '-'
        else:
            s_ID = '+'
        for i,j in pairs:
            if x[s,i,j].x>0:
                df_focal_next = pd.DataFrame()
                df_focal_next['sign'] = [s_ID]
                df_focal_next['start'] = [i]
                df_focal_next['end'] = [j]
                df_focal = pd.concat([df_focal,df_focal_next])
    ####################################################################
    return model.objVal,{s:[(i,j,x[s,i,j].x) for i,j in pairs if x[s,i,j].x>0] for s in signs},df_focal

def WeightedDistance_pratial(args):
    return WeightedDistance(args[0],args[1],weight_op=args[2],wgd=args[3],min_ops=args[4])

def semi_directed_cnd(u,v,dist_func=DirectedCopyNumberDistanceLinear):
    if validateProfiles(u,v):
        return dist_func(u,v)
    u_new,v_new=np.copy(u),np.copy(v)
    x = np.where(np.logical_and(u_new==0,v_new>0))[0]
    u_new[x]=1
    dummy_v=np.ones(len(v_new))
    dummy_v[x]=0
    additional_d = dist_func(np.ones(len(u_new)),dummy_v)
    d = dist_func(u_new,v_new)
    return d+additional_d

def semi_symmetrized_cnd(u,v,dist_func=DirectedCopyNumberDistanceLinear):
    return (semi_directed_cnd(u,v,dist_func)+semi_directed_cnd(v,u,dist_func))/2

########################################################################
############################# MAIN ROUTINE #############################
########################################################################
def WCND_arm(model):
    print(sys.path)
    # model = sys.argv[1]
    filename = model + '/' + 'ID_list.csv'
    ID_list = pd.read_csv(filename)
    for row in range(ID_list.shape[0]):
    # for row in range(10):
        ID = ID_list.iloc[row,0]
        print('===============================================================')
        print(row)
        print(ID_list.shape[0])
        print(ID)
        filename = model + '/' + ID + '.csv'
        df = pd.read_csv(filename)
        CN_from = np.array(df.iloc[:,0])
        CN_to = np.array(df.iloc[:,1])
        a, b, df_focal = WeightedDistance(CN_from, CN_to)
        print('---------------------------------------------------------------')
        print(df_focal)
        print('---------------------------------------------------------------')
        filename = model + '/' + ID + '_focal_events.csv'
        print(filename)
        df_focal.to_csv(filename)
