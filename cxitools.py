#====================#
# Python tools - CXI #
#====================# 
#
# Author: Max Hantke
# Email: maxhantke@gmail.com


from pylab import *
import h5py,os

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
    def __init__(self,filename,N,logger=None):
        self.filename = filename
        self.f = h5py.File(filename,"w")
        self.N = N
        self.logger = logger
    def write(self,d,prefix=""):
        for k in d.keys():
            name = prefix+"/"+k
            if isinstance(d[k],dict):
                if name not in self.f:
                    self.f.create_group(name)
                self.write(d[k],name)
            else:
                self.write_to_dataset(name,d[k],d["i"])
    def write_to_dataset(self,name,data,i):
        if name not in self.f:
            if isscalar(data):
                if i == -1:
                    s = [1]
                else:
                    s= [self.N]
                t=dtype(type(data))
                if t == "S":
                    t = h5py.new_vlen(str)
                axes = "experiment_identifier:value"
            else:
                s = list(data.shape)
                if i != -1:
                    s.insert(0,self.N)
                t=data.dtype
                axes = "experiment_identifier:y:x"
            self.f.create_dataset(name,s,t)
            self.f[name].attrs.modify("axes",axes)
        if i == -1:
            if isscalar(data):
                self.f[name][0] = data
            else:
                self.f[name][:] = data[:]
        else:
            if isscalar(data):
                self.f[name][i] = data
            else:
                self.f[name][i,:] = data[:]
    def close(self):
        self.f.close()

class CXIReader:
    # location can be either a file or a directory
    def __init__(self,location,dsnames={},**kwargs):
        self.logger = kwargs.get("logger",None)
        nevents = kwargs.get("nevents",0)
        ifirst = kwargs.get("ifirst",0)
        def_stack_ds = kwargs.get("def_stack_ds","/i")
        pick_events = kwargs.get("pick_events","in_sequence")
        event_filters = kwargs.get("event_filters",{})

        [self.directories,self.filenames] = self._resolve_location(location)
        self.ifile = 0
        self.ifile_opened = None
        self.Nfiles = len(self.filenames)

        self.def_stack_ds = def_stack_ds
        self.Nevents_files = self.get_Nevents_files()
        self.Nevents_tot = 0
        for N in self.Nevents_files: self.Nevents_tot += N
        self.ievent_file = -1
        self.ievent_tot = -1

        if nevents < 0:
            sys.exit("ERROR: Events to read smaller 1. Change your configuration.")
        elif nevents+ifirst > self.Nevents_tot:
            sys.exit("ERROR: Not enough events to read. Change your configuration.")
        
        to_process = []
        for N in self.Nevents_files: to_process.append(ones(N,dtype="bool"))
        to_process = self._filter_events(to_process,event_filters)
        to_process = self._pick_events(to_process,pick_events,ifirst,nevents)
        self.Nevents_process = 0
        for N in to_process: self.Nevents_process += N.sum()
        self.is_event_to_process = to_process
        self.ievent_process = -1

        self.dsnames = dsnames

    def get_next(self):
        if self._next():
            return self._get(self.dsnames)
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
            print "isdir"
            fs = filter(lambda x: x[-4:] == ".cxi",os.listdir(location))
            directories = []
            filenames = []
            for f in fs:
                directories.append(location+"/")
                filenames.append(f.split("/")[-1])
        else:
            print "is not dir"
            filenames = [location.split("/")[-1]]
            directories = [location[:-len(filenames[0])]]
        return [directories,filenames]        

    def _pick_events(self,to_process,mode,ifirst,nevents):
        temp = ones(self.Nevents_tot,dtype='bool')
        offset = 0
        for t in to_process:
            temp[offset:offset+len(t)] = t[:]        
        if mode == 'in_sequence':
            temp[:ifirst] = False
            if nevents != 0:
                if nevents <= temp.sum():
                    s = 0
                    for i in range(ifirst,self.Nevents_tot):
                        s += temp[i]
                        if s == nevents:
                            break
                    temp[i+1:] = False
        elif mode == 'randomly':
            to_pick_from = arange(self.Nevents_tot)
            to_pick_from = list(to_pick_from[temp])
            temp = zeros_like(temp)
            N = len(to_pick_from)
            for i in range(nevents):
                if N-i > 0:
                    j = to_pick_from[randint(N-i)].pop()
                    temp[j] = True
        else:
            print "ERROR: No valid picking mode chosen: %s" % mode
            return
        to_process_new = []
        i = 0
        for N in self.Nevents_files:
            to_process_new.append(temp[i:i+N])
            i += N
        return to_process_new

    def _filter_events(self,to_process,event_filters):
        for flt_name in event_filters.keys():
            flt = event_filters[flt_name]
            for (i,dty,fle) in zip(range(self.Nfiles),self.directories,self.filenames):
                f = h5py.File(dty+"/"+fle,"r")
                filter_ds = f[flt["dataset_name"]].value
                to_process[i] *= (filter_ds >= f["vmin"]) * (filter_ds <= f["vmax"])
        return to_process

    # move to next event that shall be processed
    def _next(self):
        # skip events that shall not be processed
        while True:
            self.ievent_file += 1
            # return none if end of file list is reached
            if self.ifile >= self.Nfiles:
                if self.logger != None:
                    self.logger.info("Reached end of list of files.")
                self.F.close()
                return False
            if self.ievent_file >= self.Nevents_files[self.ifile]:
                if self.logger != None:
                    self.logger.info("Reached end of file (%i) %s/%s.",self.ifile,self.directories[self.ifile],self.filenames[self.ifile])
                self.ifile += 1
                self.ievent_file = 0
            if self.is_event_to_process[self.ifile][self.ievent_file] == False:
                if self.logger != None:
                    self.logger.info("Skip event %i in file %i.",self.ievent_file,self.ifile)
            else:
                self.ievent_process += 1
                break
        return True

    def _open_file(self):
        if self.ifile_opened != self.ifile:
            print "open %s" % (self.directories[self.ifile]+'/'+self.filenames[self.ifile])
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
            D[key] = self.F[dsname][i]
        return D

