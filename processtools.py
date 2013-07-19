import os,re,sys,h5py,pylab,spimage,numpy,time
sys.path.append("/home/hantke/pythonscripts/tools")
import imgtools

def recenter(filename,region,crop):
    program = "/home/hantke/programs/tomas/find_center/find_center"
    import subprocess as sp
    p = sp.Popen(["%s" % program,"-i%s" % filename,"-s%d" % region,"-c%d" % crop], stdout = sp.PIPE)
    res = p.stdout.read()
    res = res.split()
    cx = float(res[0])
    cy = float(res[1])
    img = spimage.sp_image_read(filename,0)
    img.detector.image_center[0] = cx
    img.detector.image_center[0] = cy
    spimage.sp_image_write(img,filename,0)
    spimage.sp_image_free(img)
    #os.system("rm debug_mask.h5")
    return [cx,cy]

def stitch(pattern):
    dx = 22
    dy = 1
    Nx = pattern.shape[1]
    Ny = pattern.shape[0]
    patternStitched = pylab.zeros(shape=(Ny+dy,Nx+dx+1))
    patternStitched[dy:Ny/2+dy,0:Nx/2] = pattern[0:Ny/2,0:Nx/2]
    patternStitched[0:Ny/2,Nx/2+dx:Nx+dx] = pattern[0:Ny/2,Nx/2:Nx]
    patternStitched[Ny/2+dy:Ny+dy,1:Nx/2+1] = pattern[Ny/2:Ny,:Nx/2]
    patternStitched[Ny/2:Ny,Nx/2+dx+1:Nx+dx+1] = pattern[Ny/2:Ny,Nx/2:Nx]
    return patternStitched

def odd_column_maskout(mask):
    mode = 1
    maskout = pylab.ones_like(mask)
    maskout[mask==-3] = 0
    (Ny,Nx) = mask.shape
    X,Y = pylab.meshgrid(pylab.arange(mode,Nx+mode,1),pylab.arange(mode,Ny+mode,1))
    oscillator = X%2
    oscillator[1:-1,1:-1] *= (maskout[1:Ny-1,1:Nx-1])*(maskout[1:Ny-1,0:Nx-2])*(maskout[1:Ny-1,2:Nx])
    mask[oscillator!=0] = -2

def even_column_maskout(mask):
    mode = 0
    maskout = pylab.ones_like(mask)
    maskout[mask==-3] = 0
    (Ny,Nx) = mask.shape
    X,Y = pylab.meshgrid(pylab.arange(mode,Nx+mode,1),pylab.arange(mode,Ny+mode,1))
    oscillator = X%2
    oscillator[1:-1,1:-1] *= (maskout[1:Ny-1,1:Nx-1])*(maskout[1:Ny-1,0:Nx-2])*(maskout[1:Ny-1,2:Nx])
    mask[oscillator!=0] = -2

def column_interpolation(pattern,mask):
    (Ny,Nx) = mask.shape
    interpolationmask = pylab.zeros_like(mask)
    interpolationmask[mask==-2] = 1
    tempmask_big_c = pylab.zeros(shape=(Ny,Nx))
    tempmask_big_l = pylab.zeros(shape=(Ny,Nx))
    tempmask_big_r = pylab.zeros(shape=(Ny,Nx))
    tempmask_big_l[1:Ny-1,0:Nx-2] = interpolationmask[1:-1,1:-1]
    tempmask_big_c[1:Ny-1,1:Nx-1] = interpolationmask[1:-1,1:-1]
    tempmask_big_r[1:Ny-1,2:Nx] = interpolationmask[1:-1,1:-1]
    pattern[tempmask_big_c!=0] = (pattern[tempmask_big_l!=0]+pattern[tempmask_big_r!=0])/2.0

def even_column_interpolation(pattern,mask):
    (Ny,Nx) = mask.shape
    X,Y = pylab.meshgrid(pylab.arange(1,Nx-1,1),pylab.arange(1,Ny-1,1))
    oszillator = X%2
    oszillator -= 1
    oszillator = abs(oszillator)
    tempmask = mask[1:Ny-1,1:Nx-1]*mask[1:Ny-1,0:Nx-2]*mask[1:Ny-1,2:Nx]*oszillator
    tempmask_big_c = pylab.zeros(shape=(Ny,Nx))
    tempmask_big_l = pylab.zeros(shape=(Ny,Nx))
    tempmask_big_r = pylab.zeros(shape=(Ny,Nx))
    tempmask_big_l[1:Ny-1,0:Nx-2] = tempmask[:,:]
    tempmask_big_c[1:Ny-1,1:Nx-1] = tempmask[:,:]
    tempmask_big_r[1:Ny-1,2:Nx] = tempmask[:,:]
    pattern[tempmask_big_c!=0] = (pattern[tempmask_big_l!=0]+pattern[tempmask_big_r!=0])/2.0

def define_row_list(pattern,mode,value):
    y_list = []    
    if mode == "sumThreshold":
        for y in range(0,1024):
            if pattern[y,:512].sum() > value:
              y_list.append(pylab.array([y,0]))
        for y in range(0,1024):
            if pattern[y,512:].sum() > value:
                y_list.append(pylab.array([y,1]))
    elif mode == "sumThresholdStripe":
        left_min = -1
        right_min = -1
        left_max = -1
        right_max = -1
        for y in range(0,1024):
            if pattern[y,:512].sum() > value:
                if left_min == -1:
                    left_min = y
                left_max = y
            if pattern[y,512:].sum() > value:
                if right_min == -1:
                    right_min = y
                right_max = y
        y_list1 = pylab.arange(left_min,left_max+1,1)
        y_list2 = pylab.arange(right_min,right_max+1,1)
        y_list = pylab.ones(shape=(len(y_list1)+len(y_list2),2))
        y_list[:len(y_list1),0] = 0
        y_list[:len(y_list1),1] = y_list1.T
        y_list[len(y_list1):,1] = y_list2.T
        print y_list
    elif mode == "saturationThreshold":
        for y in range(0,1024):
            if pattern[y,:512].max() > value[0]:
                #print pattern[y,:512].sum()
                y_list.append(pylab.array([y,0]))
            if pattern[y,512:].max() > value[1]:
                y_list.append(pylab.array([y,1]))            
                #print pattern[y,512:].sum()
    elif mode == "staticStripe":
        for y in range(0,1024):
            if abs(y-512) < value/2:
                y_list.append(pylab.array([y,0]))            
                y_list.append(pylab.array([y,1]))            
    return pylab.array(y_list)

def offset_correction_border(pattern,mask):
    patternCorrected = pattern.copy()
    for y in range(0,pattern.shape[0]):
        if mask[y,:512].any() >= 0:
            (patternCorrected[y,:512])[mask[y,:512]>=0] -= patternCorrected[y,1:12].mean()
        if mask[y,:512].any() >= 0:
            (patternCorrected[y,512:])[mask[y,512:]>=0] -= patternCorrected[y,-12:-1].mean()
    return patternCorrected

# mask: 
# value=-3 -> masked out
# value=-2 -> pixels for horizontal interpolation
# value=-1 -> good pixel
# value>=0 -> member of chunk corresponding to value
def chunkmap(pattern,mask,mode,value):
    y_list = define_row_list(pattern,mode,value)
    #chunkR = pylab.zeros(len(y_list)*4)
    #chunks_ok = [[0,3],[1,0]]
    #chunks_ok = [[0,3],[1,0],[0,0],[1,3]]
    chunks_ok = []
    l = 0
    for iy in range(0,len(y_list)):
        [y,i] = y_list[iy]
        for chunk in range(0,4):
            #if [i,chunk] in chunks_ok:
            #    pass
            #else:
            (mask[y,128*(chunk+i*4):128*(chunk+1+i*4)])[mask[y,128*(chunk+i*4):128*(chunk+1+i*4)]==-1] = l
                #print len(mask[mask==l].flatten())
                #if i == 0:
                #    chunkR[l] = 3-chunk
                #elif i == 1:
                #    chunkR[l] = chunk
            l += 1
    #return chunkR

def get_chunkpixel_coordinates(chunkmap,maskR=[]):
    if maskR == []:
        calcR = False
    else:
        calcR = True
    N = int(chunkmap.max()+1)
    chunkmap_flat = chunkmap.flatten()
    chunkPixelCoordinates = pylab.zeros(shape=(N,128),dtype="int")
    if calcR: 
        maskR_flat = maskR.flatten()
        chunkPixelCoordinatesR = pylab.zeros(shape=(N,128),dtype="int")
    chunkPixelCoordinatesL = pylab.zeros(N,dtype="int")
    j = 0
    for i in range(0,N):
        coordinates = pylab.find(chunkmap_flat==i)
        L = len(coordinates)
        chunkPixelCoordinatesL[i] = L
        chunkPixelCoordinates[i,:L] = coordinates[:]
        if calcR:
            positionmap = pylab.zeros_like(chunkmap)
            positionmap = positionmap.flatten()
            positionmap[coordinates] = 1
            R = maskR_flat[positionmap==1]
            chunkPixelCoordinatesR[i,:L] = R[:]
    if calcR:
        return [chunkPixelCoordinatesL,chunkPixelCoordinates,chunkPixelCoordinatesR]
    else:
        return [chunkPixelCoordinatesL,chunkPixelCoordinates]

def get_scalingCorrectionR(patternReconstructed,patternRaw,maskR):
    Qpattern = patternReconstructed/patternRaw
    L = min([patternReconstructed.shape[0],patternReconstructed.shape[1]])/2
    scalingCorrectionR = pylab.zeros(L)
    for i in range(0,L):
        scalingCorrectionR[i] = Qpattern[maskR==i].mean()
    return scalingCorrectionR

def get_maskR(mask,data):
    Y,X = pylab.meshgrid(pylab.arange(0,mask.shape[0]),pylab.arange(0,mask.shape[1]))
    maskR = pylab.sqrt((X-mask.shape[1]/2)**2+(Y-mask.shape[0]/2)**2)
    maskR = pylab.array(maskR,'int16')
    #maskR = maskR/128
    maskR[mask!=-1] = -1
    maskR[data==0] = -1
    #print maskR
    return maskR

def get_scalingCorrectionR128(patternReconstructed,patternRaw,maskR):
    Qpattern = patternReconstructed/patternRaw
    scalingCorrectionR = pylab.zeros(4)
    for i in range(0,4):
        scalingCorrectionR[i] = Qpattern[maskR==i].mean()
    return scalingCorrectionR

def get_maskR128(mask,data):
    Y,X = pylab.meshgrid(pylab.arange(0,mask.shape[0]),pylab.arange(0,mask.shape[1]))
    maskR = pylab.sqrt((X-mask.shape[1]/2)**2+(Y-mask.shape[0]/2)**2)
    maskR = pylab.array(maskR,'int16')
    maskR = maskR/128
    maskR[mask!=-1] = -1
    maskR[data==0] = -1
    #print maskR
    return maskR


def set_scalers(pattern,mask,chunkmap,rawdata,darkdata,scalers=[],scalers0=[],maskR=[]):
    print "%i chunks have to be rescaled." % int(mask.max()+1)
    if scalers == []:
        calculateFirst = True
        scalers = pylab.ones(int(mask.max()+1))
        #TEMP CHANGE
        rowintensity = pylab.ones(int(mask.max()+1))
        
    else:
        calculateFirst = False
    Nx = pattern.shape[1]
    Ny = pattern.shape[0]
    rawtest = rawdata.copy()
    pattern[pattern.real<0] = 0
    #mask = pylab.array(mask,dtype='float64')
    rawdata = pylab.array(rawdata,dtype='float64')
    darkdata = pylab.array(darkdata,dtype='float64')
    pattern = pylab.array(pattern,dtype='float64')
    #print len(darkdata)
    #print len(pattern)
    if not calculateFirst:
        scalingCorrectionR = get_scalingCorrectionR(pattern+darkdata,rawdata,maskR)
        #print scalingCorrectionR
    mask = mask.flatten()
    pattern = abs(pattern).flatten()
    patternScaled = (rawdata-darkdata).flatten()
    patternScaled = pylab.array(patternScaled,dtype='float64')
    rawdata = rawdata.flatten()
    darkdata = darkdata.flatten()
    newScalers = 0
    for chunkindex in range(0,len(chunkmap[0])):
        L = (chunkmap[0])[chunkindex]
        if L == 0:
            pass
        else:
            C = (chunkmap[1])[chunkindex,:L]
            valuesA = (pattern+darkdata)[C]
            valuesB = rawdata[C]            
            #valuesB[rawdata[C]*scalers[chunkindex]-darkdata[C]<0] = darkdata[C]
            #valuesC = rawdata[C]*scalers[chunkindex]
            #error_before = abs(valuesA-valuesC).sum()
            (scalers[chunkindex],offset) = pylab.lstsq(pylab.vstack([valuesB,pylab.zeros_like(valuesB)]).T,valuesA)[0]
            #scalers[chunkindex] = pylab.mean(valuesA/valuesB)
            if not calculateFirst:
                #scalingCorrection = scalingCorrectionR[(chunkmap[2])[chunkindex,:L]]
                #scalingCorrection[scalingCorrection*scalers[chunkindex]<scalers0[chunkindex]] = scalers0[chunkindex]/scalers[chunkindex]            
                scalingCorrection = 1
#scalers[chunkindex] *= scalingCorrection
                if scalers[chunkindex] < scalers0[chunkindex]:
                    scalers[chunkindex] = scalers0[chunkindex]
                #                changed = False
                    print "Avoid change of scaler."
                else:
                    scalers[chunkindex] += scalers0[chunkindex]
                    scalers[chunkindex] /= 2.0
            #else:
                newScalers += 1
                changed = True
            else:
                changed = False
                scalingCorrection = 1
                #rowintensity[chunkindex] = 
            patternScaled[C] = scalers[chunkindex]*scalingCorrection*rawdata[C]-darkdata[C]
            #valuesA_scaled = (patternScaled+darkdata)[C]
            #error_after = abs(valuesA_scaled-valuesC).sum()
            #if error_after-error_before>=0 and changed:
            #    print error_after-error_before
    print "%i/%i chunks were rescaled." % (newScalers,int(mask.max()+1))
    if calculateFirst:
        scalers0 = scalers
    if (rawdata.reshape(Ny,Nx)-rawtest).sum() != 0:
        print "RAW CHANGED!!!"
        print (rawdata.reshape(Ny,Nx)-rawtest).sum()
    rawdata = rawdata.reshape(Ny,Nx)
    darkdata = darkdata.reshape(Ny,Nx)
    if calculateFirst:
        scalers0 = scalers.copy()
    return [scalers,scalers0,patternScaled.reshape(Ny,Nx)]

def saturation_maskout(pattern,mask,saturationThreshold1,saturationThreshold2):
    mask[:,:len(pattern)/2][pattern[:,:len(pattern)/2]>saturationThreshold1] = -3
    mask[:,len(pattern)/2:][pattern[:,len(pattern)/2:]>saturationThreshold2] = -3
    mask = pylab.array(mask,dtype='int32')
    mleft = mask[:,:-2]
    mright = mask[:,2:]
    mcenter = mask[:,1:-1]
    mask[:,1:-1][mask[:,1:-1]==-2][mleft[mask[:,1:-1]==-3]] = -3 
    mask[:,1:-1][mask[:,1:-1]==-2][mright[mask[:,1:-1]==-3]] = -3 

def zero(data,mask):
    data += pylab.absolute(data)
    data /= 2.0
    data[mask==-3] = 0

def support_update_threshold(data,support,support0,threshold):
    support[abs(data)<threshold] = 0
    support[abs(data)>=threshold] = 1
    support *= support0
    lowpass = imgtools.shift(imgtools.draw_circle(support.shape[1],support.shape[0],int(0.2*support.shape[0])))
    supportfourier = pylab.fft2(support)
    supportfourier *= lowpass
    support = pylab.round_(abs(pylab.ifft2(supportfourier)))
    support *= support0
    #pylab.figure()
    #pylab.imshow(imgtools.crop(imgtools.shift(abs(support)),200))
    #pylab.show()
    #exit(1)

def support_update_area(data,support,support0,area):
    sigma = 0.3
    support = pylab.zeros_like(support)
    max_val = pylab.sort(data[support0==1])[-area]
    support[data>=max_val] = 1
    lowpass = imgtools.shift(imgtools.draw_circle(support.shape[1],support.shape[0],int(sigma*support.shape[0])))
    supportfourier = pylab.fft2(support)
    supportfourier *= lowpass
    #support = pylab.round_(abs(pylab.ifft2(supportfourier)))
    support = abs(pylab.ifft2(supportfourier))
    support *= support0    
    #pylab.figure()
    #pylab.imshow(imgtools.crop(imgtools.shift(abs(support)),200))
    #pylab.show()
    #exit(1)
    return support

def support_update_positivity(data,support,support0,area):
    support *= 0
    support[data.real>0] = 1
    support *= support0    

def support_update_shrinking(support,supportDiameterStart,supportDiameterEnd,N,i):
    supportDiameter = int(supportDiameterStart-float(supportDiameterStart-supportDiameterEnd)/(N-1.0)*i)
    print "Current support diameter is %i" % supportDiameter
    supportUnshifted = imgtools.draw_circle(support.shape[1],support.shape[0],supportDiameter)
    return imgtools.shift(supportUnshifted)

def ER(scaleddata,rawdata,darkdata,scalers,scalers0,mask,supportDiameter,Niterations,startRowConstraint,outputinterval,center,chunkmap,maskR,plotflag=False,writeflag=False,filename=None):
    scaleddata_new = scaleddata.copy()
    if plotflag:
        pylab.close('all')
    beta0 = 0.8
    intervalRowConstraint = 1
    supportThreshold = 3.0E5
    supportArea = 550
    supportUpdateShrinking = False
    supportDiameterStart = 40
    supportDiameterEnd = 30
    supportDiameterN = 100
    startpositivity = 1000
    error = []
    amplitude = []
    supportUnshifted = imgtools.draw_circle(scaleddata.shape[1],scaleddata.shape[0],supportDiameter)
    support = imgtools.shift(supportUnshifted)
    #support0 = support.copy()
    if plotflag:
        fig = pylab.figure(figsize=(12,8))
        pylab.show()
    dataIt = scaleddata.copy()
    dataItsh = imgtools.shift_back(dataIt,center)
    masksh = imgtools.shift_back(mask,center)
    dtSum = 0
    simpledata = scaleddata.copy()
    simpledata[simpledata.real<0] = 0
    simpledata = abs(simpledata)
    simpledata[mask==-3] = 0
    simpledata[mask>-3] = pylab.log10(simpledata[mask>-3]+0.5)
    simpledata[mask==-3] = pylab.log10(simpledata[mask==-3])
    pylab.imsave(filename[:-4]+"_processed_simple.png",simpledata)
    for i in range(Niterations):
        print "%i/%i" % (i+1,Niterations)
        start = time.time()/60.0
        if i+1 >= startRowConstraint:
        #if i%intervalRowConstraint == intervalRowConstraint-1:
            rowConstraint = True
        else:
            rowConstraint = False
        #beta = beta0+(1-beta0)*(1-pylab.exp(-(i/7.0)**3))
        print "Calculate patterson"
        autosh = pylab.fft2(dataItsh)
        print "Apply support constraint"
        scaler = (autosh.real).sum()
        autosh *= support
        scaler /= (autosh.real).sum()
        print scaler
        #autosh *= scaler
        #autosh[autosh.real<0] = 0        
        #autosh = abs(autosh)
        auto = imgtools.shift(autosh,center)
        dataItsh = pylab.ifft2(autosh)
        dataIt = imgtools.shift(dataItsh,center)
        #scalingCorrection = pylab.mean(scaleddata[mask==-1].sum()/dataIt.real[mask==-1].sum())
        #print "Scaling: measured/reconstructed  %f" % scalingCorrection
        scalingCorrection = 1
        if (i%outputinterval == outputinterval-1 and plotflag) or i>=startRowConstraint:
            showimage1 = dataIt.copy()
        if rowConstraint:
            print "Apply intensity constraint to healthy rows"
            dataIt[mask==-1] = scaleddata[mask==-1]
            print "Apply intensity constraint to ill rows"
            [scalers,scalers0,scaleddata_new] = set_scalers(dataIt,mask,chunkmap,rawdata,darkdata,scalers,scalers0,maskR)
            dataIt[mask>=0] = scaleddata_new[mask>=0]
        else:
            print "Apply intensity constraint to all rows"
            dataIt[mask>=-1] = scaleddata[mask>=-1]
        error = list(error)
        error.append(abs(scaleddata_new-dataIt)[mask!=-3].sum())
        error = pylab.array(error)
        amplitude.append(scaler)
        print "Apply reality constraint to all rows"
        dataIt[dataIt.real<0] = 0
        dataIt = abs(dataIt)
        column_interpolation(dataIt,mask)
        dataItsh = imgtools.shift(dataIt,center)
        #if i>startpositivity:
        #    print "Positivity support update"
        #    support_update_positivity(autosh,support,support0,supportArea)
        if supportUpdateShrinking and i < supportDiameterN:
            print "Shrink support."
            support = support_update_shrinking(support,supportDiameterStart,supportDiameterEnd,supportDiameterN,i)
        if (i%outputinterval == outputinterval-1 or rowConstraint) and plotflag:
            fig.clf()
            ax1 = fig.add_subplot(234)
            ax1.imshow(abs(imgtools.crop(imgtools.shift(support*autosh,center),100,center)),interpolation="nearest")
            ax21 = fig.add_subplot(232)
            ax22 = fig.add_subplot(233)
            showimage2 = dataIt.copy()
            showimage2[dataIt.real<0] = 0
            #if rowConstraint:
            #    stripe_scaling_refinement(showimage2,mask)
            showimage1[showimage1.real<0] = 0
            #showimage2 = imgtools.gaussian_blur(showimage2,100)
            showimage2 = abs(showimage2)
            column_interpolation(showimage2,mask)
            if rowConstraint:
                calculate_row_contrasts(showimage2,mask,filename,i)
            minval = min([showimage1.min(),showimage2.min()])
            maxval = min([showimage1.max(),showimage2.max()])
            showimage1[0,0] = minval
            showimage1[0,1] = maxval
            showimage2[0,0] = minval
            showimage2[0,1] = maxval
            ax21.imshow(pylab.log10(showimage1),interpolation="nearest")
            ax22.imshow(pylab.log10(showimage2),interpolation="nearest")
            ax3 = fig.add_subplot(231)
            ax3.plot(pylab.arange(1,len(error)+1,1),pylab.array(error)/error.max(),pylab.arange(1,len(amplitude)+1,1),pylab.array(amplitude))
            pylab.ylim([0.0,1.0])
            ax41 = fig.add_subplot(235)
            ax42 = fig.add_subplot(236)
            yshow1 = center[0]-112
            yshow2 = center[0]-12
            showimage3 = showimage1.copy()
            showimage3[mask==-3] = 0
            showimage4 = showimage2.copy()
            showimage4[mask==-3] = 0
            showline11 = pylab.log10(showimage3[yshow1,:])
            showline21 = pylab.log10(showimage3[yshow2,:])
            showline12 = pylab.log10(showimage4[yshow1,:])
            showline22 = pylab.log10(showimage4[yshow2,:])
            X = pylab.arange(0,len(showline11),1)
            ax41.plot(X,showline11,'.',X,showline12,'.')
            ax42.plot(X,showline21,'.',X,showline22,'.')
            pylab.draw()
            pylab.show()
            if i+1 == outputinterval:
                pylab.draw()
                pylab.show()
            if writeflag and rowConstraint:
                showimage2[mask==-3] = 0
                showimage2[mask==-3] = pylab.log10(showimage2[mask==-3])
                showimage2[mask!=-3] = pylab.log10(showimage2[mask!=-3]+0.5)
                pylab.imsave(filename[:-4]+("_processing_%i.png" % (i+1)),showimage2)    
        stop = time.time()/60.0
        dt = stop-start
        dtSum += stop-start
        dtMean = (dtSum)/(i+1)
        tLeft = (Niterations-(i+1))*dtMean
        print "dt = %.2f min ; <dt> = %.2f ; tLeft = %.2f" % (dt,dtMean,tLeft)
    dataIt = abs(dataIt)
    column_interpolation(dataIt,mask)
    dataIt[dataIt==0] = (dataIt[dataIt!=0]).min()
    return dataIt

def get_h5filenamelist(h5Folder="./"):
    try:
        l = os.popen('ls %s' % h5Folder).readlines()
        filenames = [f for f in l if re.search('.h5$',f) and not re.search('processed',f)]
    except:
        print "Error finding h5-files in %s" % h5Folder
        exit(1)
    return filenames

def random_collect(filenames,N):
    while len(filenames) > N:
        filenames.remove(filenames[pylab.randint(len(filenames))])

def check_cassh5(filename):
    f = h5py.File(filename,"r")
    if len(f.values()[1].values()) == 3:
        newCassFlag = True
    else:
        newCassFlag = False
    return newCassFlag

def calculate_row_contrasts(pattern,mask,filename,i):
    chunk_limits = pylab.zeros(shape=(8,2))
    for chunk in range(0,8):
        if chunk < 4:
            positions = pylab.find(mask==chunk)
        else:
            positions = pylab.find(mask==mask.max()-7+chunk)
        chunk_limits[chunk,0] = positions[0]%mask.shape[1]
        chunk_limits[chunk,1] = positions[-1]%mask.shape[1]
    #print chunk_limits
    row_contrasts1 = pylab.zeros(shape=(mask.shape[0]-1,8))
    row_contrasts2 = pylab.zeros(shape=(mask.shape[0]-1,8))
    FIG1 = pylab.figure()
    FIG2 = pylab.figure()
    AX1 = FIG1.add_subplot(111)
    AX2 = FIG2.add_subplot(111)
    for chunk in range(0,8):
        for row in range(0,mask.shape[0]-1):
            xstart = chunk_limits[chunk,0]
            xstop = chunk_limits[chunk,1]
            values = pattern[row:row+2,xstart:xstop]
            maskvalues = mask[row:row+2,xstart:xstop].copy()
            maskvalues[maskvalues!=-3] = 1
            maskvalues[maskvalues==-3] = 0
            maskvalues = maskvalues[0,:]*maskvalues[1,:]
            row_contrasts1[row,chunk] = pylab.mean(abs(values[0,:][maskvalues==1]-values[1,:][maskvalues==1]))/pylab.mean(abs(values[0,:][maskvalues==1]+values[1,:][maskvalues==1]))
            row_contrasts2[row,chunk] = (pylab.mean(values[0,:][maskvalues==1])-pylab.mean(values[1,:][maskvalues==1]))/pylab.mean(values[0,:][maskvalues==1]+values[1,:][maskvalues==1])
        AX1.plot(pylab.arange(0,mask.shape[0]-1),row_contrasts1[:,chunk].T+chunk,'.')
        AX2.plot(pylab.arange(0,mask.shape[0]-1),row_contrasts2[:,chunk].T+chunk,'.')
        #pylab.draw()
    pylab.show()
    FIG1.savefig(filename[:-4]+("%i_rowContrasts1.png" % (i+1)))
    FIG2.savefig(filename[:-4]+("%i_rowContrasts2.png" % (i+1)))

def stripe_scaling_refinement(pattern,mask):
    n = 1
    top = pylab.zeros_like(mask)
    top[:top.shape[0],:] = 1
    bot = pylab.zeros_like(mask)
    bot[bot.shape[0]:,:] = 1
    mask_left = mask.copy()
    mask_right = mask.copy()
    mask_left[:,mask.shape[1]/2:] = -1
    mask_left[mask_left<0] = -1
    mask_right[:,:mask.shape[1]/2] = -1
    mask_right_first = mask_right[mask_right>=0].min()
    mask_right[mask_right>=0] -= mask_right_first + mask_right_first%4
    mask_right[mask_right<0] = -1
    mask_left[mask_left>=0] = mask_left[mask_left>=0]%4
    mask_right[mask_right>=0] = mask_right[mask_right>=0]%4
    chunk_limits = pylab.zeros(shape=(8,2))
    neighbour_rows = pylab.zeros(shape=(8,2))
    for chunk in range(0,8):
        if chunk < 4:
            mask_side = mask_left
        else:
            mask_side = mask_right
        positions = pylab.find(mask_side==chunk%4)
        if len(positions) == 0:
            break
        chunk_limits[chunk,0] = min(positions%mask.shape[1])
        chunk_limits[chunk,1] = max(positions%mask.shape[1])
        neighbour_rows[chunk,0] = positions[0]/mask.shape[1]
        neighbour_rows[chunk,1] = positions[-1]/mask.shape[1]+1
        x_min = chunk_limits[chunk,0]
        x_max = chunk_limits[chunk,1]
        #y2 = neighbour_rows[chunk,0]
        #y5 = neighbour_rows[chunk,1] 
        #y1 = y2-n
        #y3 = y2+n
        #y4 = y5-n
        #y6 = y5+n
        #upper_unscaled = (pattern[y1:y2,x_min:x_max+1])[mask[y1:y2,x_min:x_max+1]!=-3]
        #upper_scaled = (pattern[y2:y3,x_min:x_max+1])[mask[y2:y3,x_min:x_max+1]!=-3]
        #upper_scaler = pylab.mean(upper_unscaled)/ pylab.mean(upper_scaled)
        #lower_scaled = (pattern[y4:y5,x_min:x_max+1])[mask[y4:y5,x_min:x_max+1]!=-3]
        #lower_unscaled = (pattern[y5:y6,x_min:x_max+1])[mask[y5:y6,x_min:x_max+1]!=-3]
        #lower_scaler = pylab.mean(lower_unscaled)/pylab.mean(lower_scaled)
        values = pattern[neighbour_rows[chunk,0]:neighbour_rows[chunk,0]+2,x_min:x_max+1]
        maskvalues = mask[neighbour_rows[chunk,0]:neighbour_rows[chunk,0]+2,x_min:x_max+1].copy()
        maskvalues[maskvalues!=-3] = 1
        maskvalues[maskvalues==-3] = 0
        maskvalues[values==0] = 0
        maskvalues = maskvalues[0,:]*maskvalues[1,:]
        row_contrast1 = pylab.mean(values[0,:][maskvalues==1]/values[1,:][maskvalues==1])
        pattern[(mask_side+1)*top-1==chunk%4] *= row_contrast1
        values = pattern[neighbour_rows[chunk,1]-1:neighbour_rows[chunk,1]+1,x_min:x_max+1]
        maskvalues = mask[neighbour_rows[chunk,1]-1:neighbour_rows[chunk,1]+1,x_min:x_max+1].copy()
        maskvalues[maskvalues!=-3] = 1
        maskvalues[maskvalues==-3] = 0
        maskvalues[values==0] = 0
        maskvalues = maskvalues[0,:]*maskvalues[1,:]
        row_contrast2 = pylab.mean(values[1,:][maskvalues==1]/values[0,:][maskvalues==1])
        pattern[(mask_side+1)*bot-1==chunk%4] *= row_contrast2
        print row_contrast1
        print row_contrast2
        #print (lower_scaler+upper_scaler)/2.0
        #pattern[mask_side==chunk%4] *= (lower_scaler+upper_scaler)/2.0
        #pattern[mask_side==chunk%4] *= (row_contrast1+row_contrast2)/2.0
    print chunk_limits
    print neighbour_rows
