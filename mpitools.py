#====================#
# Python tools - MPI #
#====================# 
#
# Author: Max Hantke
# Email: maxhantke@gmail.com

from pylab import *
import time
from mpi4py import MPI

def multinode_master(comm,Njobs,getwork,logres=None,logger=None):
    comm_size = comm.Get_size()
    t_start = time.time()

    Njobs_done = 0
    Njobs_started = 0
    results = []
    texec = zeros(shape=(comm_size,2))
    # send initial jobs
    for i in range(1,comm_size):
        comm.send(getwork(),dest=i,tag=1)
        texec[i,0] = time.time()
        Njobs_started += 1

    dummybuffer = array(1, dtype='i')        
    status = MPI.Status()
    Dt = {"wait":[],"send":[],"write":[],"read":[],"receive":[],"work":[]}
    while Njobs_done < Njobs:
        t0 = time.time()
        request = comm.Irecv(dummybuffer,MPI.ANY_SOURCE,0)
        MPI.Request.Wait(request, status)
        i_done = status.source
        texec[i_done,1] = t1 = time.time()
        Njobs_done += 1
        if logger != None:
            logger.info("Datarate %.1f Hz job %i/%i rank %i/%i",Njobs_done/(time.time()-t_start),Njobs_done,Njobs,i_done,comm_size)
        result = comm.recv(source=i_done,tag=1)
        t2 = time.time()
        if logres != None:
            logres(result)
        t3 = time.time()
        Dt["work"].append((texec[i_done,1]-texec[i_done,0])/(1.*(comm_size-1)))
        if Njobs_started < Njobs:
            work = getwork()
            t4 = time.time()
            comm.send(work,dest=i_done,tag=1)
            t5 = texec[i_done,0] = time.time()
            Njobs_started += 1
        else:
            t5 = t4 = time.time()
        Dt["wait"].append(t1-t0)
        Dt["receive"].append(t2-t1)
        Dt["write"].append(t3-t2)
        Dt["read"].append(t4-t3)
        Dt["send"].append(t5-t4)
        S = ""
        for k in Dt.keys():
            S += "%s %f " % (k,Dt[k][-1])
        if logger != None:
            logger.info("Speed: " + S)

    for i in range(1,comm_size):
        comm.send("fika",dest=i,tag=1)

    
    return Dt

def multinode_slave(comm,worker,logger=None):
    work_package = comm.recv(source=0,tag=1)
    dummydata = array(1, dtype='i')
    while work_package != "fika":
        result = worker(work_package)
        comm.Send([dummydata,MPI.INT],0,0)
        comm.send(result,dest=0,tag=1)
        work_package = comm.recv(source=0,tag=1)
