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
                self.write_to_dataset(name,d[k],d.get("i",-1))
    def write_to_dataset(self,name,data,i):
        #if self.logger != None:
        #    self.logger.info("Write dataset %s of event %i." % (name,i))
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

        if self.logger != None:
            for d,f,N in zip(self.directories,self.filenames,self.Nevents_files):
                self.logger.info("Found file %s/%s with %i events.",d,f,N)

        self.ievent_file = -1
        self.ievent_tot = -1

        if nevents < 0:
            sys.exit("ERROR: Events to read smaller 0. Change your configuration.")
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
            fs = filter(lambda x: x[-4:] == ".cxi",os.listdir(location))
            fs.sort()
            directories = []
            filenames = []
            for f in fs:
                directories.append(location+"/")
                filenames.append(f.split("/")[-1])
        else:
            filenames = [location.split("/")[-1]]
            directories = [location[:-len(filenames[0])]]
        return [directories,filenames]        

    def _pick_events(self,to_process,mode,ifirst,nevents):
        temp = ones(self.Nevents_tot,dtype='bool')
        offset = 0
        for t in to_process:
            temp[offset:offset+len(t)] = t[:]        
            offset += len(t)
        if mode == 'in_sequence':
            temp[:ifirst] = False
            if nevents != 0:
                if nevents < temp.sum():
                    s = 0
                    for i in range(ifirst,self.Nevents_tot):
                        if temp[i]:
                            s += 1
                        if s == nevents:
                            break
                    temp[i+1:] = False
        elif mode == 'random':
            to_pick_from = arange(self.Nevents_tot)
            to_pick_from = list(to_pick_from[temp])
            temp = zeros_like(temp)
            for i in range(nevents):
                N = len(to_pick_from)
                if N != 1:
                    j = to_pick_from.pop(randint(N))
                    temp[j] = True
                else:
                    temp[to_pick_from[0]] = True
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
        for (i,dty,fle) in zip(range(self.Nfiles),self.directories,self.filenames):
            if self.logger != None:
                self.logger.info("Filtering %s/%s",dty,fle)
            f = h5py.File(dty+"/"+fle,"r")
            for flt_name in event_filters.keys():
                flt = event_filters[flt_name]
                filter_ds = f[flt["dataset_name"]].value.flatten()
                if "vmin" in flt.keys() and "vmax" in flt.keys():
                    F = (filter_ds >= flt["vmin"]) * (filter_ds <= flt["vmax"])
                elif "indices" in flt.keys():
                    if i != 0:
                        if self.logger != None:
                            self.logger.warning("Filter indices are applied to every file!")
                    F = zeros_like(to_process[i])
                    for index in flt["indices"]:
                        F[filter_ds == index] = True
                else:
                    if self.logger != None:
                        self.logger.warning("No valid filter arguments given for filter %s!" % flt_name)
                    F = ones_like(to_process[i])
                to_process[i] *= F
                if self.logger != None:
                    self.logger.info("Filter %s - yield %.3f %% -> total yield %.3f %%",flt_name,100.*F.sum()/len(F),100.*to_process[i].sum()/len(F))
                    self.logger.info("Filter %s - First index: %i",flt_name,(arange(len(to_process[i]))[to_process[i]])[0])
        return to_process

    # move to next event that shall be processed
    def _next(self):
        # skip events that shall not be processed
        while True:
            if self.ievent_process == self.Nevents_process-1:
                if self.logger != None:
                    self.logger.info("Reached last event to process.")
                return False
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
                self.ievent_file = -1
            if self.is_event_to_process[self.ifile][self.ievent_file] == False:
                pass
                #if self.logger != None:
                #    self.logger.info("Skip event %i in file %i.",self.ievent_file,self.ifile)
            else:
                self.ievent_process += 1
                break
        return True

    def _open_file(self):
        if self.ifile_opened != self.ifile:
            if self.ifile_opened != None:
                if self.logger != None:
                    self.logger.info("Closing file: %s/%s",self.directories[self.ifile_opened],self.filenames[self.ifile_opened])            
                self.F.close()
            if self.logger != None:
                self.logger.info("Opening file: %s/%s",self.directories[self.ifile],self.filenames[self.ifile])
            self.F = h5py.File(self.directories[self.ifile]+'/'+self.filenames[self.ifile],'r')
            self.ifile_opened = self.ifile
        
    def _get(self,dsnames):
        self._open_file()
        D = {}
        D["i"] = self.ievent_process
        for (key,dsname) in dsnames.items():
            D[key] = self.F[dsname][self.ievent_file].copy()
        return D

def get_filters(C):
    filters = filter(lambda x: "filter" in x,C.keys())
    Cfilters = {}
    for f in filters:
        Cfilters[f] = C[f]
    return Cfilters


def cxi_to_spimage(filename,i,ds_img="/entry_1/image_2/data",ds_msk="/entry_1/image_2/mask"):
    import spimage
    f = h5py.File(filename,"r")
    img_cxi = f[ds_img][i,:,:]
    img_cxi[img_cxi<0] = 0
    msk_cxi = f[ds_msk][i,:,:]
    Nx = msk_cxi.shape[1]
    Ny = msk_cxi.shape[0]
    img = spimage.sp_image_alloc(Nx,Ny,1)
    img.image[:,:] = img_cxi[:,:]
    img.mask[:,:] = (msk_cxi[:,:] & PIXEL_IS_IN_MASK) == 0
    filename_new = filename[:-4] + ("_%i.h5" % i)
    spimage.sp_image_write(img,filename_new,0)
    spimage.sp_image_free(img)
    f.close()
