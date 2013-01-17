import numpy,pylab,os,re,csv,string,Image,spimage

# in_filter can be a string or a list of strings
def get_filenames(in_filter=None,path="./"):
    filenames = os.popen('ls %s' % path).readlines()
    if in_filter:
        if isinstance(in_filter,list):
            for filt in in_filter:
                filenames = [f for f in filenames if re.search(filt,f)]
        else:
            filenames = [f for f in filenames if re.search(in_filter,f)]
    filenames.sort()
    return filenames

def estimate_type(var):
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

def png_mask_to_h5(filename,Nx,Ny):
    I = Image.open(filename)
    D = pylab.array(I.getdata())[:,0]
    D=D.reshape((Nx,Ny))
    img = spimage.sp_image_alloc(Nx,Ny,1)
    img.mask[:,:] = D[:,:]/255
    spimage.sp_image_write(img,filename[:-4]+'.h5',0)
    spimage.sp_image_free(img)

def get_png_mask(filename):
    I = Image.open(filename)
    D = pylab.array(I.getdata())[:,0]
    D = D.reshape((I.size[1],I.size[0]))
    return D[:,:]/255.

def save_to_csv(filename,list_of_arrays,list_of_array_names=[]):
    """Save given array values to a new csv file."""
    f = open(filename, 'wb')
    writer = csv.writer(f)
    if list_of_array_names != []: writer.writerow(list_of_array_names)
    for i in range(0,len(list_of_arrays[0])):
        row = []
        for j in range(0,len(list_of_arrays)): 
            row.append(list_of_arrays[j][i])
        writer.writerow(row)
    f.close()

def load_from_csv(filename,has_header=False):
    """Loads csv file and gives out one array for each column and the header."""
    f = open(filename,'r')
    reader = csv.reader(f)
    rows = []
    for row in reader: rows.append(row)
    f.close()
    if has_header:
        header = rows[0]
        rows.remove(rows[0])
    else:
        header = []
        for i in range(0,len(row)):
            header.append(string.letters[i])
    print header
    out = []
    for i in range(0,len(row)): out.append([header[i],[]])
    out = dict(out)
    for i in range(0,len(rows)):
        for j in range(0,len(row)):
            out[header[j]].append(estimate_type(rows[i][j]))
    return out

def create_numberstring(number,number_of_letters):
    number_str = ""
    for j in -numpy.arange(-number_of_letters+1,1,1):
        number_str =  number_str + ("%i" % (number/pow(10,j)%10))
    return number_str


def smoothList(list,strippedXs=False,degree=10):  
    if strippedXs==True: return Xs[0:-(len(list)-(len(list)-degree+1))]  
    smoothed=[0]*(len(list)-degree+1)  
    for i in range(len(smoothed)):  
        smoothed[i]=sum(list[i:i+degree])/float(degree)  
    return smoothed  

def smoothListTriangle(list,strippedXs=False,degree=5):  
    weight=[]  
    window=degree*2-1  
    smoothed=[0.0]*(len(list)-window)  
    for x in range(1,2*degree):weight.append(degree-abs(degree-x))  
    w=numpy.array(weight)  
    for i in range(len(smoothed)):  
        smoothed[i]=sum(numpy.array(list[i:i+window])*w)/float(sum(w))  
    return smoothed  

def smoothListGaussian(list,strippedXs=False,degree=5):  
    window=degree*2-1  
    weight=numpy.array([1.0]*window)  
    weightGauss=[]  
    for i in range(window):  
        i=i-degree+1  
        frac=i/float(window)  
        gauss=1/(numpy.exp((4*(frac))**2))  
        weightGauss.append(gauss)  
    weight=numpy.array(weightGauss)*weight  
    smoothed=[0.0]*(len(list)-window)  
    for i in range(len(smoothed)):  
        smoothed[i]=sum(numpy.array(list[i:i+window])*weight)/sum(weight)  
    return smoothed 

def smoothGaussian(array,sigma):
    return pylab.ifft(pylab.fft(array)*1/pylab.sqrt(2*pylab.pi)/sigma*pylab.exp(pylab.arange(0,len(array),1.0)**2/2.0/sigma**2))

def smoothGaussian1dMirror(array,sigma):
    array_reversed = list(array.copy())
    array_reversed.reverse()
    array_reversed = pylab.array(array_reversed)
    array_mirrored = pylab.zeros(2*len(array))
    array_mirrored[:len(array)] = array[:]
    array_mirrored[len(array):] = array_reversed[:]
    array_smoothed = pylab.ifft(pylab.fft(array_mirrored)*1/pylab.sqrt(2*pylab.pi)/sigma*pylab.exp(pylab.arange(0,len(array_mirrored),1.0)**2/2.0/sigma**2) )
    array_smoothed = array_smoothed[:len(array)]
    print array_smoothed
    return array_smoothed

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
         
    input:
    x: the input signal 
    window_len: the dimension of the smoothing window; should be an odd integer
    window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
    flat window will produce a moving average smoothing.
   
    output:
    the smoothed signal
           
    example:
    
    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
    
    TODO: the window parameter could be the window itself if an array instead of a string   
    """
    
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    
    
    if window_len<3:
        return x
    
    
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    
    
    s=numpy.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
       #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')
        
    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1] 

def downsample_position(position,downsampling):
    return (position-(downsampling-1)/2.)/(1.*downsampling)
