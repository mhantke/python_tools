import os,re,sys,h5py,pylab,numpy,time,fileinput,imgtools,spimage

# final log values from uwrapc.log
def get_log_values(filename):
    variableNames = ["It(out)","lognumber?","Ereal","Efourier","<Fc/Fo>","SupSize(%)","Beta","Threshold","Algorithm","dEreal","dSupSize","Blur Radius","Int Cum Fluct","Object Area","Phase Relation Error","Correlation with solution","Phase Blur Radius","Delta Rho","Iterations/s"]
    #variableNames = []
    outDict = []
    atflag = False
    valueflag = False
    fh = open(filename,"r")
    lines = fh.readlines()
    L = len(lines)
    #print lines
    for j in range(0,len(lines)):
        line = lines[j]
        if line[0:3] == "@ s":
            atflag = True
        elif atflag == True:
            line = line.split("\t")        
            if valueflag == False:
                variables = pylab.zeros(shape=(L-j,len(variableNames)))
                valueflag = True
                i = 0
            line = pylab.array(line,dtype="float64")
            variables[i,:] = line[:]
            i += 1
    for i in range(0,len(variableNames)):
        outDict.append((variableNames[i],variables[:,i]))
    fh.close()
    return dict(outDict)

def get_log_value(filename,variableName):
    [VN,V] = get_log_values(filename)
    vindex = VN.index(variableName)
    return V[vindex]

def plot_log_values():
    try:
        l = os.popen('ls').readlines()
        directories = [f[:-1] for f in l if not re.search('reconstructed_h5',f) and not re.search('result_png',f)]
    except:
        print "Error finding folders"
        exit(1)
    pylab.figure()
    for directory in directories:
      	log_values = get_log_values(directory + "/000000/uwrapc.log")
        pylab.plot(log_values[1][:,0],log_values[1][:,3])
        pylab.draw()
    pylab.show()

def find_images(variable,maxvalue):	
    try:
        l = os.popen('ls').readlines()
        directories = [f[:-1] for f in l if not re.search('reconstructed_h5',f) and not re.search('result_png',f)]
    except:
        print "Error finding folders"
        exit(1)
    list_bigger = []
    list_smaller = []
    for directory in directories:
        lnum = os.popen('ls %s' % directory).readlines()
        directoriesnum = [f[:-1] for f in lnum if re.search('0000',f)]
        valueSum = 0
        valueN = 0
        for directorynum in directoriesnum:
            log_values = get_log_values(directory + "/" + directorynum + "/uwrapc.log")
            valueSum += log_values[1][-1,list(log_values[0]).index(variable)]
            valueN += 1
        valueAverage = valueSum/valueN
        if valueAverage > maxvalue:
            #list_bigger.append([directory,valueAverage])
            list_bigger.append(directory+"-avg_image.h5")
        else:
            list_smaller.append(directory+"-avg_image.h5")
            #list_smaller.append([directory,valueAverage])
    #print "Smaller than threshold:"
    #for i in list_smaller:
    #    print i
    #print "Bigger than threshold:"
    #for i in list_bigger:
    #    print i
    return list_smaller

def parameters(filenames=None):
    if not filenames:
        try:
            l = os.popen('ls').readlines()
            filenames = [f[:-1] for f in l if re.search('image.h5',f)]
        except:
            print "Error finding folders"
            exit(1)
    ellipticities = []
    areas = []
    mindiameters = []
    maxdiameters = []
    for filename in filenames:
        img = spimage.sp_image_read(filename,0)
        #img = spimage.sp_image_shift(img)
        A = abs(img.image)
        #B = imgtools.crop(A,50,None,4)
        C = A.copy()
        threshold = C.max()/100.0
        C[C<threshold] = -1
        C[C>=threshold] = 0
        C += 1
        area = float(C.sum())/float(len(C.flatten()))
        [mindiameter,maxdiameter] = diameter_extrema(C)
        ellipticity = imgtools.ellipticity(C)
        ellipticities.append(ellipticity)
        areas.append(area)
        mindiameters.append(mindiameter)
        maxdiameters.append(maxdiameter)
    return [areas,mindiameters,maxdiameters,ellipticities]

def make_colorbar(colorflag,minvalue=0.0,maxvalue=1.0,filename="colorbar.png"):
    Ny = 1000
    Nx = 40
    spbar = spimage.sp_image_alloc(Nx,Ny,1)
    a = pylab.arange(0.0,1.0,1.0/(1.0*Ny)).T
    A = pylab.ones(shape=(Nx,Ny))
    A *= a
    A = imgtools.turnccw(A)
    spbar.image[:,:] = A[:,:]
    spimage.sp_image_write(spbar,filename,colorflag)
