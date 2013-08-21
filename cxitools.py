#====================#
# Python tools - CXI #
#====================# 
#
# Author: Max Hantke
# Email: maxhantke@gmail.com


from pylab import *
import h5py

# CXI pixelmask bits
PIXEL_IS_PERFECT = 0
PIXEL_IS_INVALID = 1
PIXEL_IS_SATURATED = 2
PIXEL_IS_HOT = 4
PIXEL_IS_DEAD = 8
PIXEL_IS_SHADOWED = 16
PIXEL_IS_IN_PEAKMASK = 32
PIXEL_IS_TO_BE_IGNORED = 64
PIXEL_IS_BAD = 128
PIXEL_IS_OUT_OF_RESOLUTION_LIMITS = 256
PIXEL_IS_MISSING = 512
PIXEL_IS_IN_MASK = PIXEL_IS_INVALID |  PIXEL_IS_SATURATED | PIXEL_IS_HOT | PIXEL_IS_DEAD | PIXEL_IS_SHADOWED | PIXEL_IS_IN_PEAKMASK | PIXEL_IS_TO_BE_IGNORED | PIXEL_IS_BAD | PIXEL_IS_MISSING 

class CXIWriter:
    def __init__(self,filename,N):
        self.filename = filename
        self.f = h5py.File(filename,"w")
        self.N = N
    def write(self,d):
        for k in d.keys():
            if isinstance(d[k],dict):
                if k not in self.f:
                    self.f.create_group(k)
                for j in d[k].keys():
                    name = k+"/"+j
                    self.write_to_dataset(name,d[k][j],d["i"])
            else:
                self.write_to_dataset(k,d[k],d["i"])
    def write_to_dataset(self,name,data,i):
        if name not in self.f:
            if isscalar(data):
                s=tuple([self.N])
                t=dtype(type(data))
                if t == "S":
                    t = h5py.new_vlen(str)
                axes = "experiment_identifier:value"
            else:
                s = list(data.shape)
                s.insert(0,self.N)
                t=data.dtype
                axes = "experiment_identifier:y:x"
            self.f.create_dataset(name,s,t)
            self.f[name].attrs.modify("axes",axes)
        if isscalar(data):
            self.f[name][i] = data
        else:
            self.f[name][i,:] = data[:]
    def close(self):
        self.f.close()

class CXIReader:
    # location can be either a file or a directory
    def __init__(self,location,def_stack_ds="/i",**kwargs):
        
        [self.directories,self.filenames] = self._resolve_location(location)
        self.ifile = 0
        self.ifile_opened = None
        self.Nfiles = len(filenames)

        self.def_stack_ds = def_stack_ds
        self.Nevents_files = self.get_Nevents_files()
        self.Nevents_tot = 0
        for N in self.Nevents_files: self.Nevents_tot += N
        self.ievent_file = -1
        self.ievent_tot = -1

        nevents = kwargs.get('nevents',0)
        ifirst = kwargs.get('ifirst',0)
        if nevents < 0:
            print "ERROR: Events to read smaller 1. Change your configuration."
            return
        elif nevents+ifirst >= self.Nevents_tot:
            print "ERROR: Not enough events to read. Change your configuration."
            return
        
        self.is_event_to_process = self._pick_events(cg.get('pick_events','in_sequence'),ifirst,nevents)
        self.Nevents_process = 0
        for N in self.is_event_to_process: self.Nevents_process += N.sum()
        self.ievent_process = -1

    def get_next(self,dsnames):
        if self._next():
            return self._get(dsnames)
        else:
            return None

    def get_Nevents_files(self):
        N = []
        for i in range(len(self.filenames)):
            F = h5py.File(self.directories[i]+'/'+self.filenames[i],'r')
            N.append(F[self.def_stack_ds].shape[0])
            F.close()
        return N
        
    def close(self):
        if self.ifile_opened != None:
            self.F.close()

    def _resolve_location(self,location):
        if os.path.isdir(location):
            fs = filter(lambda x: x[:-4] == ".cxi",os.listdir(location))
            directories = []
            filenames = []
            for f in fs:
                directories.append(location+"/")
                filenames.append(f.split("/")[-1])
        else:
            filenames = [location.split("/")[-1]]
            directories = [location[:-len(filenames[0])]]
        return [directories,filenames]

    def _pick_events(self,mode,ifirst,nevents):
        if mode == 'in_sequence':
            temp = ones(self.Nevents_tot,dtype='bool')
            temp[:ifirst] = False
            if nevents != 0:
                if nevents+ifirst < self.Nevents_tot:
                    temp[ifirst+nevents:] = False
        elif mode == 'randomly':
            temp = zeros(self.Nevents_tot,dtype='bool')
            r = range(ifirst,self.Nevents_tot)
            for i in range(nevents):
                k = self.Nevents_tot-i
                if k >= 0:
                    j = r[randint(k)].pop()
                    temp[j] = True
        else:
            print "ERROR: No valid picking mode chosen: %s" % mode
            return
        L = []
        i = 0
        for N in self.Nevents_files:
            L.append(temp[i:i+N])
            i += N
        return L

    # move to next event that shall be processed
    def _next(self):
        # skip events that shall not be processed
        while True:
            self.ievent_file += 1
            # return none if end of file list is reached
            if self.ifile >= self.Nfiles:
                logger.info("Reached end of list of files.")
                self.F.close()
                return False
            if self.ievent_file >= self.Nevents_files[self.ifile]:
                logger.info("Reached end of file (%i) %s/%s.",self.ifile,self.directories[self.ifile],self.filenames[self.ifile])
                self.ifile += 1
                self.ievent_file = 0
            if self.is_event_to_process[self.ifile][self.ievent_file] == False:
                logger.info("Skip event %i in file %i.",self.ievent_file,self.ifile)
            else:
                self.ievent_process += 1
                break
        return True

    def _open_file(self):
        if self.ifile_opened != self.ifile:
            try:
                self.F.close()
            except:
                pass
            self.F = h5py.File(self.directories[self.ifile]+'/'+self.filenames[self.ifile],'r')
            self.ifile_opened = self.ifile
        
    def _get(self,dsnames):
        self._open_file()
        i = self.ievent_process
        D = {}
        D["i"] = i
        for (key,dsname) in dsnames.items():
            D[tag] = self.F[dsname][i]
        return D

