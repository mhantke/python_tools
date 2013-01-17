import os,re,sys,h5py,pylab,spimage,numpy,fnmatch,time,datetime,ConfigParser
sys.path.append("/home/hantke/pythonscripts/tools")
import imgtools, gentools, hawktools
from imgtools import *
from gentools import *

def estimateType(var):
    #first test bools
    if var == 'True':
            return True
    elif var == 'False':
            return False
    else:
        #int
        try:
            return int(var)
        except ValueError:
            pass
        #float
        try:
            return float(var)
        except ValueError:
            pass
        #string
        try:
            return str(var)
        except ValueError:
            raise NameError('Something Messed Up Autocasting var %s (%s)' 
                            % (var, type(var)))

def read_configfile(configfile):
    config = ConfigParser.ConfigParser()
    try:
        config.readfp(open(configfile))
    except IOError:
        print "ERROR: Can't read configuration-file."

    conf = []
    for section in config.sections(): 
        conf.extend(config.items(section))
    conf = dict(conf)
    for key in conf.keys(): 
        conf[key] = estimateType(conf[key])
        if "directory" in key:
            if conf[key][-1] == '/': conf[key] = conf[key][:-1]
    return conf

def get_filenames_to_process(c):
    # Find H5 files and *.conf to process
    # ------------------------------------------------
    fn = os.listdir(c['images_directory'])
    fimages = fnmatch.filter(fn,"*.h5")
    fimages.sort()
    if c['nimages'] != 0:
        fimages = fimages[:c['nimages']]
    fn = os.listdir(c['uwrapc_conf_directory'])
    fconf = fnmatch.filter(fn,"*.conf")
    fconf.sort()
    return [fimages,fconf]

def init_directory_tree(c):
    # Initialize output diretory tree
    # -----------------------------------------

    c['root_directory'] = os.path.abspath(os.path.curdir)
    l = os.listdir(c['output_directory'])
    t = datetime.date.today()
    if ("recrun%i%02i%02i" % (t.year,t.month,t.day)) in l: 
        os.system("rm %s/recrun%i%02i%02i -r" % (c['output_directory'],t.year,t.month,t.day))
    c['run_directory'] = c['output_directory']+"/recrun%i%02i%02i" % (t.year,t.month,t.day)
    os.system("mkdir %s" % c['run_directory'])
    c['reconstructions_png_directory'] = c['run_directory']+"/reconstructions_png"
    os.system("mkdir %s" % c['reconstructions_png_directory'])
    c['reconstructions_h5_directory'] = c['run_directory']+"/reconstructions_h5"
    os.system("mkdir %s" % c['reconstructions_h5_directory'])
    c['reconstructions_h5_interpolated_directory'] = c['run_directory']+"/reconstructions_h5_interpolated"
    os.system("mkdir %s" % c['reconstructions_h5_interpolated_directory'])
    c['prtfs_png_directory'] = c['run_directory']+"/prtfs_png"
    os.system("mkdir %s" % c['prtfs_png_directory'])
    c['intensities_png_directory'] = c['run_directory']+"/intensities_png"
    os.system("mkdir %s" % c['intensities_png_directory'])
    for nf in range(0,c['NImages']):
        c['image%.06i_directory' % nf] = c['run_directory']+"/image%.06i" % nf
        os.system("mkdir %s" % c['image%.06i_directory' % nf])
        for nc in range(0,c['NConfs']):
            c['image%.06iconf%.06i_directory' % (nf,nc)] = c['image%.06i_directory' % nf]+"/conf%.06i" % nc
            os.system("mkdir %s" % c['image%.06iconf%.06i_directory' % (nf,nc)])

def cp_images_png(c):
    # expecting that pngs already made... !?
    for i in range(len(c['images_filenames'])):
        #print c['images_filenames'][i][:-3]
        os.system('cp %s/%s.png %s/img%06i.png' % (c['images_directory'],
                                                   c['images_filenames'][i][:-3],
                                                   c['intensities_png_directory'],
                                                   i))
def cp_confs(c):
    os.chdir(c['uwrapc_conf_directory'])
    for nc in range(0,c['NConfs']):
        conf_filename = c['conf_filenames'][nc]
        for nf in range(0,c['NImages']):
            os.system("cp %s %s/uwrapc.conf" % (conf_filename,c['image%.06iconf%.06i_directory' % (nf,nc)]))
    os.chdir(c['root_directory'])


def cp_images(c):
    os.chdir(c['images_directory'])
    for nf in range(0,c['NImages']):
        image_filename = c['images_filenames'][nf]
        for nc in range(0,c['NConfs']):
            os.system("cp %s %s/img.h5" % (image_filename,c['image%.06iconf%.06i_directory' % (nf,nc)]))
    os.chdir(c['root_directory'])

def init_general_logfile(rund):
    #try:
    logf = open(rund+"/recrun.log","w")
    logf.write("# image filename\timage number stamp\tconf filename\tconf number stamp\titerations\tEreal\tEfourier\t<Fc/Fo>\n")
    return logf
    #except:
    #    return False

def repeat_reconstruction(c,nf,nc):
    os.chdir(c['image%.06iconf%.06i_directory' % (nf,nc)])
    os.system("repeat_reconstruction.pl %i" % c['nseeds'])
    os.chdir(c['root_directory'])

def get_confd(c,nf,nc):
    return c['image%.06iconf%.06i_directory' % (nf,nc)]

def log_reconstruction_parameters(c,nf,nc):
    logs = []
    EfourierEnds = pylab.zeros(c['nseeds'])
    confd = get_confd(c,nf,nc)
    for seed in range(0,c['nseeds']):
        params= hawktools.get_log_values("%s/seed%06i_uwrapc.log" % (confd,seed))
        It = params['It(out)']
        Ereal = params['Ereal']
        Efourier = params['Efourier']
        FcFo = params['<Fc/Fo>']
        ErealEnd = Ereal[-1]
        EfourierEnd = Efourier[-1]
        EfourierEnds[seed] = EfourierEnd
        FcFoEnd = FcFo[-1]
        ItEnd = It[-1]
        c['logfile'].write("%s\t%.06i\t%s\t%.06i\t%i\t%f\t%f\t%f\n" % (c['images_filenames'][nf],nf,c['conf_filenames'][nc],nc,ItEnd,EfourierEnd,ErealEnd,FcFoEnd))
        logs.append(params)
    print "Average Fourier error: %e" % EfourierEnds.mean()
    return logs

def clean_reconstruction_output(c,nf,nc):
    confd = c['image%.06iconf%.06i_directory' % (nf,nc)]
    for s in range(0,c['nseeds']):
        seedd = confd + "/%.06i" % s
        os.system("mv %s/real_space-%s.h5 %s/seed%.06i_real_space-%s.h5" % (seedd,c['lastlog_str'],confd,s,c['lastlog_str']))
        os.system("mv %s/fourier_space-%s.h5 %s/seed%.06i_fourier_space-%s.h5" % (seedd,c['lastlog_str'],confd,s,c['lastlog_str']))
        os.system("mv %s/uwrapc.log %s/seed%.06i_uwrapc.log" % (seedd,confd,s))
        os.system("rm %s/ -r" % seedd)

def get_best_reconstructions(logs,c):
    minp_str = "Efourier"
    pend = pylab.zeros(len(logs))
    for i in range(0,len(logs)):
        pend[i] = logs[i][minp_str][-1]
    median_pend = pylab.median(pend)
    Ns = len(logs)
    prtf_filenames = []
    for s in range(0,Ns):
        if pend[s] <= median_pend:
            prtf_filenames.append("seed%.06i_real_space-%s.h5" % (s,c['lastlog_str']))
    return prtf_filenames

def get_mean_pend(logs,pname):
    pend = pylab.zeros(len(logs))
    for i in range(0,len(logs)):
        pend[i] = logs[i][pname][-1]
    return pend.mean()

def perform_prtf(flist,c,nf,nc):
    flist_str = ""
    i = nf*c['NConfs']+nc
    confd = get_confd(c,nf,nc)
    for fp in flist:
        flist_str += " %s/%s" % (confd,fp)
    os.chdir(confd)
    prtf_filename = "prtf_image%.06i_conf%.06i" % (nf,nc)
    os.system("prtf %s %s" % (prtf_filename,flist_str))
    os.system("rm *.png")
    prtf = pylab.loadtxt(prtf_filename).T
    os.chdir(c['root_directory'])
    if 'prtfs' not in c.keys():
        plog = open(c['prtfslogfilePath'],"w")
        plog.write("# image filename\timage number stamp\tconf filename\tconf number stamp\t1/e-resolution\n")
        c['prtfs'] = pylab.zeros(shape=(1+c['NImages']*c['NConfs'],prtf.shape[1]))
        c['prtfs'][0,:] = prtf[0,:]
        ress = pylab.zeros(c['NImages']*c['NConfs'])
    else:
        plog = open(c['prtfslogfilePath'],"a")
    c['prtfs'][i+1,:] = prtf[1,:]
    res = pylab.find(prtf[1,:]<pylab.exp(-1.0))
    if len(res) > 0:
        res = pylab.find(prtf[1,:]<pylab.exp(-1.0))[0]
    else:
        res = 0
    plog.write("%s\t%.06i\t%s\t%.06i\t%i\n" % (c['images_filenames'][nf],nf,c['conf_filenames'][nc],nc,res))
    plog.close()
    return prtf

def plot_prtf(data,c,nf,nc):

    class Setup:
        def __init__(self,w,d,p,n,R,o):
            self.wavelength = w
            self.distance = d
            self.pixel_size = p
            self.number_of_images = n
            self.detector_center_edge_distance = R
            self.real_space_oversampling = o
            if self.wavelength > 0.0 and\
                    self.distance > 0.0 and\
                    self.pixel_size > 0.0:
                self.known = True
            else:
                self.known = False

        def convert(self,pixels,nm=False):
            if self.known:
                if nm:
                    return pixels*self.pixel_size/self.distance/(self.wavelength/1.E-09)
                else:
                    return pixels*self.pixel_size/self.distance/self.wavelength
            else:
                return pixels
            
        def rescale(self,value):
            if self.number_of_images != 0:
                return (value-1.0/pylab.sqrt(float(self.number_of_images)))/(1.0-1.0/pylab.sqrt(float(self.number_of_images)))
            else:
                return value
    
        def smooth(self,values):
            if self.sigma > 0:
                return abs(gentools.smooth(values,self.sigma))
            else:
                return values

    setup = Setup(c['wavelength'], c['detector_distance'],
                  c['pixel_size']*c['downsampling'], c['nseeds'],
                  c['detector_center_edge_distance'], c['magnification_factor'])

    # calc resolution
    i_min_res = 0
    while setup.rescale(data[1,i_min_res]) > pylab.exp(-1) and not i_min_res == data.shape[1]-1:
        i_min_res += 1
    
    y0 = pylab.exp(-1)
    y1 = setup.rescale(data[1,i_min_res-1])
    y2 = setup.rescale(data[1,i_min_res])
    x1 = setup.convert(data[0,i_min_res-1],True)
    x2 = setup.convert(data[0,i_min_res],True)
    x_res = x1 + (y0-y1)*(x2-x1)/(1.0*(y2-y1))
    
    cryst_resolution_nm = 1/x_res
    cryst_resolution = cryst_resolution_nm*1.E-09

    from matplotlib import rc
    rc('text', usetex=True)
    rc('font', family='serif')
    dpi = 100
    fsize = 30
    f_linewidth = 2
    p_linewidth = 2
    pylab.rcParams['axes.linewidth'] = f_linewidth
    rc('lines', linewidth=p_linewidth)
    pylab.rc(('xtick.major','xtick.minor','ytick.major','ytick.minor'), pad=15)
    fig = pylab.figure()
    fig.set_size_inches((8,5))
    ax = fig.add_axes([0.2,0.25,0.7,0.5])#,frameon=False, xticks=[], yticks=[])
    for label in ax.xaxis.get_ticklabels():
        label.set_fontsize(fsize)
    for label in ax.yaxis.get_ticklabels():
        label.set_fontsize(fsize)
    pylab.ylabel(r"$PRTF$",fontsize = fsize,labelpad=10)

    # plot data
    ax.plot(pylab.array([setup.convert(i,True) for i in data[0,:]]),
            pylab.array([setup.rescale(i) for i in data[1,:]]),lw=p_linewidth,color="#000099")
    # plot horizontal line at y=1/e
    ax.plot([setup.convert(data[0,0],True),setup.convert(data[0,-1],True)],
            [pylab.exp(-1),pylab.exp(-1)],linestyle="dashed",color="black",lw=p_linewidth)
    # plot vertical line at x(y=1/e)
    ax.plot([x_res,x_res],[0,y0],color="#990000")
    
    for tl in pylab.gca().get_xticklines() + pylab.gca().get_yticklines(): 
        tl.set_markeredgewidth(f_linewidth) 

    xmax = pylab.ceil(setup.convert(data[0,-1],True)/4.0/1.0E-3)*4*1.0E-03
    ax.set_ylim([0.0,1.0])
    ax.set_xlim([0.0,xmax])
    ax.set_xticks(pylab.arange(0,xmax,0.01))
    ax.set_yticks(pylab.arange(0,1.5,0.5))
    
    if setup.known and i_min_res != -1:
        ax.set_title(r"$R_{\mbox{crystall}} = %.1f \mbox{ nm}$" % (1.0/x_res),fontsize=fsize,ha="right",x=1.0,color="#990000")
        pylab.xlabel(r"$1/R$  [1/nm]",fontsize = fsize,labelpad=10)

    print "Resolution %i: %e" % (i,cryst_resolution)
    prtf_filename_png = "%s/prtf_image%06i_conf%06i.png" % (c['prtfs_png_directory'],nf,nc)
    reconstruction_filename_png = "%s/prtf_image%06i_conf%06i-avg_image_amplitudes_magn%i.png" % (c['reconstructions_png_directory'],nf,nc,c['magnification_factor'])
    intensities_filename_png = "%s/img%06i.png" % (c['intensities_png_directory'],nf)
    os.system("python_script_scalingbar -f %s -d %s -w %e -R %e -D %e -o %i -c black -r %e" % (reconstruction_filename_png,reconstruction_filename_png,setup.wavelength,setup.detector_center_edge_distance,setup.distance,setup.real_space_oversampling,cryst_resolution))
    pylab.savefig(prtf_filename_png,dpi=dpi)
    os.system("python_script_put_images_together -A %s -B %s -C %s -R" % (intensities_filename_png,reconstruction_filename_png,reconstruction_filename_png))
    os.system("python_script_put_images_together -A %s -B %s -C %s -R" % (reconstruction_filename_png,prtf_filename_png,prtf_filename_png))

def padding(c,nf,nc):
    os.chdir(get_confd(c,nf,nc))
    prtf_filename = "prtf_image%.06i_conf%.06i" % (nf,nc)
    os.system("python_script_pad_zeros -i %s-avg_image.h5 -f %i -H -P -p -a -c %i -C %s" % (prtf_filename,
                                                                                            c['magnification_factor'],
                                                                                            c['crop_reconstruction'],
                                                                                            c['colormap']))
    #os.system("cp %s-avg_image_interp.h5" % (c['root_directory']))
    os.chdir(c['root_directory'])
    
def move_averaged_images(c,nf,nc):
    confd = get_confd(c,nf,nc)
    prtf_filename = "prtf_image%.06i_conf%.06i" % (nf,nc)
    os.system("mv %s/%s-avg_image_phases_magn%i.png %s/" % (confd,prtf_filename,c['magnification_factor'],c['reconstructions_png_directory']))
    os.system("mv %s/%s-avg_image_amplitudes_magn%i.png %s/" % (confd,prtf_filename,c['magnification_factor'],c['reconstructions_png_directory']))
    os.system("cp %s/%s-avg_image.h5 %s/" % (confd,prtf_filename,c['reconstructions_h5_directory']))
    os.system("cp %s/%s-avg_image_magn%i.h5 %s/" % (confd,prtf_filename,c['magnification_factor'],c['reconstructions_h5_interpolated_directory']))
               
def generate_intensity_pngs(c):
    os.chdir(c['images_directory'])
    os.system('python_script_to_png Mask Jet Log')
    os.chdir(c['root_directory'])
