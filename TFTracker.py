import numpy as np
import tensorflow as tf
from copy  import deepcopy
from array import array
import ROOT

class trdpyinterface():
    tfworker  = None
    def __init__(self,trackmodel):
        self.tfworker = tf.saved_model.load(trackmodel)

    def makedf(self,x,v,det,tag):
        # mask dead detectors
        #mask = np.isin(det,dead)
        #det  = det[~mask]
        #x    = x[~mask]
        #v    = v[~mask]
        #tag  = tag[~mask]

        sector =det//30
        vsector=sector
        stack  = (det%30)//6
        layer  = (det%6)

        dist = np.sqrt(v[:,0]*v[:,0] + v[:,1]*v[:,1])
        c    = v[:,0]/dist
        s    = v[:,1]/dist
        data = np.stack([vsector,sector,det,stack,layer,x[:,0],x[:,1],x[:,2],v[:,0],v[:,1],v[:,2],c,s,tag]).T
        return data

    def findtrklts(self,data):
        vsector=data[:,0]
        alltriplets=[]
        allquartets=[]
        allquintets=[]
        allsixlets=[]
        allpara3=[]
        allpara4=[]
        allpara5=[]
        allpara6=[]

        for vs in range(18):
            mask = vsector==vs
            dataTF = tf.constant(data[mask],dtype=tf.float32,name="data")
            triplets,quartets,quintets,sixlets,parameter3,parameter4,parameter5,parameter6 = self.tfworker.TF_tracker(dataTF,vs)
            alltriplets.append(triplets.numpy())
            allquartets.append(quartets.numpy())
            allquintets.append(quintets.numpy())
            allsixlets.append(sixlets.numpy())
            allpara3.append(parameter3)
            allpara4.append(parameter4)
            allpara5.append(parameter5)
            allpara6.append(parameter6)
        alltriplets=np.concatenate(alltriplets,axis=0)
        allquartets=np.concatenate(allquartets,axis=0)
        allquintets=np.concatenate(allquintets,axis=0)
        allsixlets=np.concatenate(allsixlets,axis=0)
        allpara   = [ allsixlets,allquintets,allquartets,alltriplets ,allpara6,allpara5,allpara4,allpara3 ]
        return allpara

    def worker(self):
        print("trdpyinterface loaded")

    def gettrackparameter(self,para):
        results=[]
        for idx,result in enumerate(para):
            #print("sector ",idx," found ",np.array(result["L"])," tracks")
            for index in range(result["L"]):
                #print("      ",idx,index)
                tmp = result["parameter"][index]
                results.append(np.array(tmp))
        #print(results)
        results=np.concatenate([results]).flatten()
        return results
        
        #print("x0,y0,z0,tgl,C,PHI0");

    def getdata(self,inputs):
        x = inputs[:,:3]
        v = inputs[:,3:6]
        det = inputs[:,6]
        tag = inputs[:,7]
        data = self.makedf(x,v,det,tag)
        allpara = self.findtrklts(data)
        p6 = allpara[4]
        cls6 = allpara[0]
        results=self.gettrackparameter(p6)
        #print(allpara[0])
        return cls6.flatten(), results


