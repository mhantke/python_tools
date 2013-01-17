import os,re,sys,h5py,pylab,spimage,numpy,fnmatch,time,datetime
sys.path.append("/home/hantke/pythonscripts/tools")
import imgtools, gentools, hawktools
from imgtools import *
from gentools import *

RUNMODE_DEFAULT = 1
RUNMODE_OPTIMIZE_PARAMETERS = 2
RUNMODE_OPTIMIZE_SUPPORT_SIZE = 4
RUNMODE_AUTOCORRELATION_SUPPORT_SIZE = 8

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
    #if conf['directory'][-1] != '/': conf['directory'] += '/'
    return conf


def get_configuration_parameters(config):
    # Load configuration for given experiment
    # ---------------------------------------------------------
    # LCLS201101 517 eV
    if config == "LCLS201101_517":
        return {'magnification_factor' : 7,
                'crop_reconstruction' : 212,
                'wavelength' : 2.40E-09,
                'detector_distance' : 0.741,
                'binning' : 2,
                'downsampling' : 2,
                'pixel_size' : 75.0E-06,
                'N' : 1024,
                'auto_xmin' : 17,
                'auto_xmax' : 23,
                'auto_lowcut' : 6,
                'auto_highcut' : 50,
                'auto_interpolation' : 4.0,
                'auto_smoothing' : 15.0
                }
    # LCLS201101 1096 eV
    elif config == "LCLS201101_1096":
        return {'magnification_factor' : 7,
                'crop_reconstruction' : 449,
                'wavelength' : 1.13E-09,
                'detector_distance' : 0.741,
                'binning' : 1,
                'downsampling' : 2,
                'pixel_size' : 75.0E-06,
                'N' : 1024,
                'auto_xmin' : 34,
                'auto_xmax' : 46,
                'auto_lowcut' : 20,
                'auto_highcut' : 110,
                'auto_interpolation' : 4.0,
                'auto_smoothing' : 25.0
                }
    # LCLS201101 1096 eV
    elif config == "LCLS201207_1096_carbo":
        return {'magnification_factor' : 15,
                'crop_reconstruction' : 349,
                'wavelength' : 1.13E-09,
                'detector_distance' : 0.741,
                'binning' : 1,
                'downsampling' : 4,
                'pixel_size' : 75.0E-06,
                'N' : 1024,
                'auto_xmin' : 34,
                'auto_xmax' : 46,
                'auto_lowcut' : 20,
                'auto_highcut' : 110,
                'auto_interpolation' : 4.0,
                'auto_smoothing' : 25.0
                }
    # FLASH201009
    elif config == "FLASH201009":
        return {'magnification_factor' : 13,
                'crop_reconstruction' : 255,
                'wavelength' : 5.7E-09,
                'detector_distance' : 0.150,
                'binning' : 4,
                'downsampling' : 4,
                'pixel_size' : 15.0E-06,
                'N' : 1600,
                'auto_xmin' : 9,
                'auto_xmax' : 17,
                'auto_lowcut' : 6,
                'auto_highcut' : 23,
                'auto_interpolation' : 4.0,
                'auto_smoothing' : 15.0
                }
    # Simulation
    elif config == "simulation":
        return {'magnification_factor' : 7,
                'crop_reconstruction' : 200,
                'wavelength' : 5.7E-09,
                'detector_distance' : 0.150,
                'binning' : 4,
                'downsampling' : 2,
                'pixel_size' : 15.0E-06,
                'N' : 1600,
                'auto_xmin' : 9,
                'auto_xmax' : 17,
                'auto_lowcut' : 6,
                'auto_highcut' : 23
                }
    elif config == "simulation2":
        return {'magnification_factor' : 7,
                'crop_reconstruction' : 200,
                'wavelength' : 2.4E-09,
                'detector_distance' : 0.741,
                'binning' : 1,
                'downsampling' : 1,
                'pixel_size' : 75.0E-06,
                'N' : 512,
                }
    else:
        return {}

def get_filenames_to_process(rmode):
    # Find H5 files and *.conf to process
    # ------------------------------------------------
    try:
        fn = os.listdir(os.path.curdir)
        fimages = fnmatch.filter(fn,"*.h5")
        fimages.sort()
        fconf = fnmatch.filter(fn,"*.conf")
        fconf.sort()
    except:
        return []
    if rmode & RUNMODE_OPTIMIZE_SUPPORT_SIZE:
        if len(fconf) > 2:
            print "ERROR: Too many *.conf files."
            exit(0)
    return [fimages,fconf]

def init_directory_tree(Nf,Nc,rootd,rmode):
    # Initialize output diretory tree
    # -----------------------------------------

    dirs = {}

    l = os.listdir(rootd)
    l.sort()
    t = datetime.date.today()
    
    if ("recrun%i%02i%02i" % (t.year,t.month,t.day)) in l: 
        os.system("rm %s/recrun%i%02i%02i -r" % (rootd,t.year,t.month,t.day))
    dirs['run'] = rootd+"/recrun%i%02i%02i" % (t.year,t.month,t.day)
    os.system("mkdir %s" % dirs['run'])
    dirs['reconstructions_png'] = dirs['run']+"/reconstructions_png"
    os.system("mkdir %s" % dirs['reconstructions_png'])
    dirs['reconstructions_h5'] = dirs['run']+"/reconstructions_h5"
    os.system("mkdir %s" % dirs['reconstructions_h5'])
    dirs['prtfs_png'] = dirs['run']+"/prtfs_png"
    os.system("mkdir %s" % dirs['prtfs_png'])
    dirs['intensities_png'] = dirs['run']+"/intensities_png"
    os.system("mkdir %s" % dirs['intensities_png'])
    for nf in range(0,Nf):
        dirs['image%.06i' % nf] = dirs['run']+"/image%.06i" % nf
        os.system("mkdir %s" % dirs['image%.06i' % nf])
        for nc in range(0,Nc):
            dirs['image%.06iconf%.06i' % (nf,nc)] = dirs['image%.06i' % nf]+"/conf%.06i" % nc
            os.system("mkdir %s" % dirs['image%.06iconf%.06i' % (nf,nc)])
    dirs['root'] = rootd
    return dirs

def init_general_logfile(rund):
    try:
        logf = open(rund+"/recrun.log","w")
        logf.write("# image filename\timage number stamp\tconf filename\tconf number stamp\titerations\tEreal\tEfourier\t<Fc/Fo>\n")
        return logf
    except:
        return False

def repeat_reconstruction(f,rootd,confd,Ns):
    try:
        os.system("cp %s/%s %s/img.h5" %  (rootd,f,confd))
        os.chdir(confd)
        os.system("repeat_reconstruction.pl %i" % Ns)
        os.chdir(rootd)
        return 1
    except:
        return 0

def log_reconstruction_parameters(logf,confd,nf,nc,f,u,Ns):
    logs = []
    Efourier_sum = 0
    for seed in range(0,Ns):
        params= hawktools.get_log_values("%s/%.06i/uwrapc.log" % (confd,seed))
        print params
        It = params['It(out)']
        Ereal = params['Ereal']
        Efourier = params['Efourier']
        FcFo = params['<Fc/Fo>']
        ErealEnd = Ereal[-1]
        EfourierEnd = Efourier[-1]
        FcFoEnd = FcFo[-1]
        ItEnd = It[-1]
        logf.write("%s\t%.06i\t%s\t%.06i\t%i\t%f\t%f\t%f\n" % (f,nf,u,nc,ItEnd,EfourierEnd,ErealEnd,FcFoEnd))
        logs.append(params)
        Efourier_sum += Efourier
    print "Average Fourier error: %e" % (Efourier_sum/(1.0*Ns))
    return logs

def clean_reconstruction_output(confd,Ns,ll_str):
    for s in range(0,Ns):
        seedd = confd + "/%.06i" % s
        os.system("mv %s/real_space-%s.h5 %s/seed%.06i_real_space-%s.h5" % (seedd,ll_str,confd,s,ll_str))
        os.system("mv %s/fourier_space-%s.h5 %s/seed%.06i_fourier_space-%s.h5" % (seedd,ll_str,confd,s,ll_str))
        os.system("mv %s/uwrapc.log %s/seed%.06i_uwrapc.log" % (seedd,confd,s))
        os.system("rm %s/ -r" % seedd)

def get_best_reconstructions(logs,ll_str):
    minp_str = "Efourier"
    pend = pylab.zeros(len(logs))
    for i in range(0,len(logs)):
        pend[i] = logs[i][minp_str][-1]
    median_pend = pylab.median(pend)
    Ns = len(logs)
    prtf_filenames = []
    for s in range(0,Ns):
        if pend[s] <= median_pend:
            prtf_filenames.append("seed%.06i_real_space-%s.h5" % (s,ll_str))
    return prtf_filenames

def get_mean_pend(logs,pname):
    pend = pylab.zeros(len(logs))
    for i in range(0,len(logs)):
        pend[i] = logs[i][pname][-1]
    return pend.mean()

def perform_prtf(f,flist,nf,nc,Nf,Nc,confN,confd,rund,rootd,plogp,ps):
    flist_str = ""
    i = nf*Nc+nc
    for fp in flist:
        flist_str += " %s/%s" % (confd,fp)
    os.chdir(confd)
    prtf_filename = "prtf_image%.06i_conf%.06i" % (nf,nc)
    os.system("prtf %s %s" % (prtf_filename,flist_str))
    os.system("rm *.png")
    prtf = pylab.loadtxt(prtf_filename).T
    os.chdir(rootd)
    if ps == []:
        plog = open(plogp,"w")
        plog.write("# image filename\timage number stamp\tconf filename\tconf number stamp\t1/e-resolution\n")
        ps = pylab.zeros(shape=(1+Nf*Nc,prtf.shape[1]))
        ps[0,:] = prtf[0,:]
        ress = pylab.zeros(Nf*Nc)
    else:
        plog = open(plogp,"a")
    ps[i+1,:] = prtf[1,:]
    res = pylab.find(prtf[1,:]<pylab.exp(-1.0))
    if len(res) > 0:
        res = pylab.find(prtf[1,:]<pylab.exp(-1.0))[0]
    else:
        res = 0
    plog.write("%s\t%.06i\t%s\t%.06i\t%i\n" % (f,nf,confN,nc,res))
    plog.close()
    return ps

def padding(m,c,rootd,confd,nf,nc,colormap="WhiteJetShadow"):
    os.chdir(confd)
    prtf_filename = "prtf_image%.06i_conf%.06i" % (nf,nc)
    os.system("python_script_pad_zeros -i %s-avg_image.h5 -f %i -H -P -p -a  -c %i -C %s" % (prtf_filename,m,c,colormap))
    os.system("cp *.png %s" % (rootd))
    os.chdir(rootd)
    
#def scalingbar(rootd,confd,m,D,R, w,nf,nc):
#    os.chdir(confd)
#    filenames = ["prtf_image%.06i_conf%.06i-avg_image_amplitudes_magn%i.png" % (nf,nc,m),
#                 "prtf_image%.06i_conf%.06i-avg_image_phases_magn%i.png" % (nf,nc,m)]
#    colors = ["black","white"]
#    for i in range(0,len(filenames)):
#        os.system("python_script_scalingbar -f %s -R %f -w %f -D %f -o %f -c %s" % (filename[i],R,w,D,m,color[i]))
#    os.chdir(rootd)

def move_averaged_images(confd,recspngd,recsh5d,nf,nc,m):
    prtf_filename = "prtf_image%.06i_conf%.06i" % (nf,nc)
    os.system("mv %s/%s-avg_image_phases_magn%i.png %s/%s-avg_image_phases_magn%i.png" % (confd,prtf_filename,m,recspngd,prtf_filename,m))
    os.system("mv %s/%s-avg_image_amplitudes_magn%i.png %s/%s-avg_image_amplitudes_magn%i.png" % (confd,prtf_filename,m,recspngd,prtf_filename,m))
    os.system("cp %s/%s-avg_image.h5 %s/%s-avg_image.h5" % (confd,prtf_filename,recsh5d,prtf_filename))
               
def set_support_size_optimum(rlogs,confd,rootd,confN):
    optimalSupSizes = pylab.zeros(len(rlogs))
    for seed in range(0,len(rlogs)):
        rlog = rlogs[seed]
        SupSize = rlog['SupSize(%)']
        Efourier = rlog['Efourier']
        optimalSupSizes[seed] = SupSize[Efourier.argmin()]
    #print optimalSupSizes
    new_supportsize = optimalSupSizes.mean()/100.0
    write_new_conf(rootd,confN,new_supportsize)

def write_new_conf(rootd,confd,confN,new_supportsize):
    fconf_org = open("%s/%s" % (rootd,confN),"r")
    lines = fconf_org.readlines()
    for i_line in range(0,len(lines)):
        line = lines[i_line] 
        if line[:13] == "  object_area":
            break
    confNameAbsolute = "%s/uwrapc.conf" % confd
    fconf_new = open(confNameAbsolute,"w")
    for j_line in range(0,len(lines)):
        if j_line == i_line:
            line = lines[j_line][:-18] + "%.6e ] );\n" % new_supportsize
        else:
            line = lines[j_line]
        fconf_new.write(line)
    fconf_new.close()

def shrink_support(confd_opt,confd,confN,rootd,ll_str,nseeds):
    threshold_parameter = 0.1
    area_sum = 0
    for seed in range(0,nseeds):
        img = spimage.sp_image_read("%s/seed%.06i_real_space-%s.h5" % (confd_opt,seed,ll_str),0)
        i = abs(img.image)
        m = img.mask
        med = pylab.median(i[m==1])
        area_sum += len((i[i>med*threshold_parameter]).flatten())
        spimage.sp_image_free(img)
    new_supportsize = pylab.ceil(area_sum/(1.0*nseeds))/(1.0*i.shape[0]*i.shape[1])
    write_new_conf(rootd,confd,confN,new_supportsize)
    return new_supportsize


def make_filter(image_size,a):
    x_array = pylab.arange(image_size) - image_size/2 + 0.5
    y_array = pylab.arange(image_size) - image_size/2 + 0.5
    X_array, Y_array = pylab.meshgrid(x_array, y_array)
    r = pylab.sqrt(X_array**2 + Y_array**2)
    kernel = (r/2.0/a)**4*pylab.exp(2.0-r**2/2.0/a**2)
    kernel[r > 2.0*a] = 1
    return kernel, r

def get_autoradius(filename,sigma,outcut,xmin,xmax,interpolation=4.0,smoothing=15.0):
    img = spimage.sp_image_read(filename,0)
    Nx = img.image.shape[1]
    Ny = img.image.shape[0]
    kernel = make_filter(Nx,sigma)[0]
    img.image.real *= kernel
    img.image.real *= img.mask 
    auto = spimage.sp_image_fftw3(img)
    auto_cr = pylab.fftshift(abs(auto.image))
    auto_cr /= auto_cr.mean()
    auto_cr[(auto_cr.shape[0]-outcut)/2:(auto_cr.shape[0]+outcut)/2,:] = pylab.Inf
    auto_cr[:,(auto_cr.shape[1]-outcut)/2:(auto_cr.shape[1]+outcut)/2] = pylab.Inf
    pylab.imsave(filename[:-3]+"_auto.png",pylab.log10(auto_cr),cmap=pylab.get_cmap("hot"))
    spimage.sp_image_free(img)
    spimage.sp_image_free(auto)
    [r,A] = imgtools.radial_pixel_average(auto_cr)
    r = r[xmin:xmax]
    A = A[xmin:xmax]
    r = r[pylab.isfinite(A)]
    A = A[pylab.isfinite(A)]
    A /= A.mean()
    f = interpolation
    interpolated = imgtools.interpolate1d(A,f)
    r_interpolated = pylab.arange(0,len(r),1/f)+r.min()
    smoothed = abs(gentools.smooth(interpolated,smoothing))
    smoothed /= smoothed.mean()
    pylab.clf()
    pylab.plot(r_interpolated,smoothed,r_interpolated,interpolated,r,A)
    pylab.savefig(filename[:-3]+"_auto_radial.png")
    #x = r_interpolated[interpolated.argmax()]
    x = r_interpolated[smoothed.argmax()]
    return x
