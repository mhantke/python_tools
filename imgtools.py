import os,re,sys,h5py,pylab,numpy,time


def get_phase(x):
    return abs(pylab.log(x/abs(x))/1.0j)

def get_R_and_Theta_map(Nx,Ny,cx=None,cy=None):
    if not cx:
        cx = (Nx-1)/2.0
    if not cy:
        cy = (Ny-1)/2.0
    x = pylab.arange(0,Nx,1.0)-cx
    y = pylab.arange(0,Ny,1.0)-cy
    X,Y = pylab.meshgrid(x,y)
    R = pylab.sqrt(X**2+Y**2)
    R = R.round()
    Theta = pylab.arctan(-Y/X)
    Theta[X<0] += pylab.pi
    Theta += pylab.pi/2.0
    #pylab.imsave("Theta.png" , Theta)
    #pylab.imsave("X.png" , X)
    #pylab.imsave("Y.png" , Y)
    return [R,Theta]

def cone_pixel_average(image,N_theta,cx=None,cy=None):
    [R,Theta] = get_R_and_Theta_map(image.shape[1],image.shape[0],cx,cy)
    R[pylab.isfinite(image) == False] = -1
    radii = pylab.arange(R.min(),R.max()+1,1)
    if radii[0] == -1:
        radii = radii[1:]
    values = pylab.zeros(shape=(N_theta,len(radii)))
    for j in range(0,N_theta):
        theta_min = j/(1.0*N_theta)*2.0*pylab.pi
        theta_max = (j+1)/(1.0*N_theta)*2.0*pylab.pi
        theta_center = (theta_max+theta_min)/2.0
        theta_interval = theta_max-theta_min
        theta_image = image[abs(Theta-theta_center)<=theta_interval/2.0]
        theta_R = R[abs(Theta-theta_center)<=theta_interval/2.0]
        for i in range(0,len(radii)):
            temp = theta_image[theta_R==radii[i]].copy()
            temp = temp[pylab.isfinite(temp)]
            values[j,i] = temp.mean()
    return [radii,values]

def radial_pixel_sum(image,cx=None,cy=None):
    if not cx:
        cx = (image.shape[1]-1)/2.0
    if not cy:
        cy = (image.shape[0]-1)/2.0
    x = pylab.arange(0,image.shape[1],1.0)-cx
    y = pylab.arange(0,image.shape[1],1.0)-cy
    X,Y = pylab.meshgrid(x,y)
    R = pylab.sqrt(X**2+Y**2)
    R = R.round()
    radii = pylab.arange(R.min(),R.max()+1,1)
    values = pylab.zeros_like(radii)
    for i in range(0,len(radii)):
        values[i] = image[R==radii[i]].sum()
    return pylab.array([radii,values])


def cartesian_to_polar(cartesian_pattern,N_theta,x_center=None,y_center=None,R=None):
    Nx = cartesian_pattern.shape[1]
    Ny = cartesian_pattern.shape[0]
    if R == None:
        R = int(min([Nx,Ny])/2.0-1)
    polar_pattern = pylab.zeros(shape=(R,N_theta))
    if x_center == None:
        x_center = Nx/2.0-0.5
    if y_center == None:
        y_center = Ny/2.0-0.5
    for i_theta in range(0,N_theta):
        for r in range(0,R):
            theta = 2*pylab.pi*i_theta/float(N_theta)
            x = x_center + r * pylab.sin(theta)
            y = y_center + r * pylab.cos(theta)
            # bilinear interpolation
            x1 = int(pylab.floor(x))
            x2 = x1+1
            y1 = int(pylab.floor(y))
            y2 = y1+1
            if x1 < Nx and x2 < Nx and y1 < Ny and y2 < Ny and  x1 >= 0 and x2 >= 0 and y1 >= 0 and y2 >= 0:
                    V11 = cartesian_pattern[int(pylab.floor(y)),int(pylab.floor(x))]
                    V12 = cartesian_pattern[int(pylab.floor(y)),int(pylab.floor(x))+1]
                    V21 = cartesian_pattern[int(pylab.floor(y))+1,int(pylab.floor(x))]
                    V22 = cartesian_pattern[int(pylab.floor(y))+1,int(pylab.floor(x))+1]
                    polar_pattern[r,i_theta] = V11*(x2-x)*(y2-y) + V12*(x-x1)*(y2-y) + V21*(x2-x)*(y-y1) + V22*(x-x1)*(y-y1)
            else:
                polar_pattern[r,i_theta] = pylab.infty
            
    return polar_pattern

def cartesian_to_radial(cartesian,N_theta):
    return pylab.mean(cartesian_to_polar(cartesian,N_theta),1)

def draw_circle(Nx,Ny,diameter):
    X,Y = pylab.meshgrid(pylab.arange(-Nx/2.0+0.5,Nx/2.0+0.5,1),pylab.arange(-Ny/2.0+0.5,Ny/2.0+0.5,1))
    circle = pylab.sqrt(X**2+Y**2)    
    circle[circle>diameter/2.0] = 0
    circle[circle!=0] = 1
    return circle

def gaussian_smooth(I,sm):
    N = 2*sm
    if len(I.shape) == 2:
        import scipy.signal
        kernel = pylab.zeros(shape=(1+2*sm,1+2*sm))
        X,Y = pylab.meshgrid(pylab.arange(0,N,1),pylab.arange(0,N,1))
        X -= (N-1)/2.
        Y -= (N-1)/2.
        R = pylab.sqrt(X**2 + Y**2)
        kernel = pylab.exp(R**2/(1.0*sm**2))
        Ism = scipy.signal.convolve2d(I,pylab.fftshift(kernel),mode='same',boundary='fill')
        return Ism
    elif len(I.shape) == 1:
        kernel = pylab.zeros(1+2*sm)
        X = pylab.arange(0,N,1)
        X -= (N-1)/2.
        kernel = pylab.exp(X**2/(1.0*sm**2))
        Ism = pylab.convolve(I,pylab.fftshift(kernel),mode='full')
        return Ism

def gaussian_smooth_2d1d(I,sm):
    N = 2*sm
    if len(I.shape) == 2:
        import scipy.signal
        kernel = pylab.zeros(shape=(1+2*sm,1+2*sm))
        X,Y = pylab.meshgrid(pylab.arange(0,N,1),pylab.arange(0,N,1))
        X -= (N-1)/2.
        #Y -= (N-1)/2.
        #R = pylab.sqrt(X**2 + Y**2)
        kernel = pylab.exp(X**2/(1.0*sm**2))
        Ism = scipy.signal.convolve2d(I,kernel,mode='same',boundary='wrap')
        return Ism
    elif len(I.shape) == 1:
        print "Error input"
        return []

def gaussian_sharpening(image,sigma):
    imagefourier = pylab.fft2(image)
    Ny = image.shape[0]
    Nx = image.shape[1]
    X,Y = pylab.meshgrid(pylab.arange(-Nx/2.0+0.5,Nx/2.0+0.5,1.0),pylab.arange(-Ny/2.0+0.5,Ny/2.0+0.5,1.0))
    gauss = 1/pylab.sqrt(2*pylab.pi*sigma**2)*pylab.exp(-(X**2+Y**2)/(2*sigma**2))
    gauss = shift(gauss)
    imagefourier *= (1.0-gauss)
    return pylab.ifft2(imagefourier)
    
def downsample(array2d_raw,factor,mode="pick"):
    array2d = pylab.array(array2d_raw,dtype=array2d_raw.dtype)
    factor = int(factor)
    if factor == 1:
        return array2d
    available_modes = ["pick","integrate"]#,"interpolate"]
    if not mode in available_modes:
        print "ERROR: %s is not a valid mode." % mode
        return 0
    Ny = array2d.shape[0]
    Nx = array2d.shape[1]
    if mode == "pick": 
        Ny_new = int(pylab.ceil(1.0*Ny/factor))
        Nx_new = int(pylab.ceil(1.0*Nx/factor))  
        array2d_new = pylab.zeros(Nx_new*Ny_new,dtype=array2d.dtype)  
        array2d_flat = array2d.flatten()
        for i in pylab.arange(0,Nx_new*Ny_new,1):
            ind = i%Nx_new*factor+(i/Nx_new)*Nx*factor
            array2d_new[i] = array2d_flat[ind]
        return pylab.reshape(array2d_new,(Ny_new,Nx_new))
    elif mode == "integrate":
        Ny_new = int(pylab.ceil(1.0*Ny/factor))
        Nx_new = int(pylab.ceil(1.0*Nx/factor))  
        array2d_new = pylab.zeros(shape=(Ny_new,Nx_new),dtype=array2d.dtype)
        for y_new in pylab.arange(0,Ny_new,1):
            for x_new in pylab.arange(0,Nx_new,1):
                y_min = y_new*factor
                y_max = (y_new+1)*factor
                x_min = x_new*factor
                x_max = (x_new+1)*factor
                if y_max < Ny and x_max < Nx:
                    array2d_new[y_new,x_new] = array2d[y_min:y_max,x_min:x_max].mean()
                else:
                    if y_max >= Ny and x_max >= Nx:
                        array2d_new[y_new,x_new] = array2d[y_min:,x_min:].mean()
                    elif y_max >= Ny:
                        array2d_new[y_new,x_new] = array2d[y_min:,x_min:x_max].mean()
                    elif x_max >= Nx:
                        array2d_new[y_new,x_new] = array2d[y_min:y_max,x_min:].mean()
        return array2d_new
    elif mode == "interpolate":
        Ny_new = int(pylab.floor(1.0*Ny/factor))
        Nx_new = int(pylab.floor(1.0*Nx/factor))
        array2d_new = pylab.zeros(shape=(Ny_new,Nx_new),dtype=array2d.dtype)
        for y_new in pylab.arange(0,Ny_new,1):
            for x_new in pylab.arange(0,Nx_new,1):
                y_min = y_new*factor
                y_max = (y_new+1)*factor
                x_min = x_new*factor
                x_max = (x_new+1)*factor

def downsample3d(array_raw,factor,mask):
    array_cp = array_raw.copy()
    factor = int(factor)
    if factor == 1:
        return array_cp
    Nz = array_cp.shape[0]
    Ny = array_cp.shape[1]
    Nx = array_cp.shape[2]
    Nz_new = int(pylab.ceil(1.0*Nz/factor))
    Ny_new = int(pylab.ceil(1.0*Ny/factor))
    Nx_new = int(pylab.ceil(1.0*Nx/factor))  
    array_new = pylab.zeros(shape=(Nz_new,Ny_new,Nx_new),dtype=array_cp.dtype)
    for z_new in pylab.arange(0,Nz_new,1):
        for y_new in pylab.arange(0,Ny_new,1):
            for x_new in pylab.arange(0,Nx_new,1):
                z_min = z_new*factor
                z_max = min([(z_new+1)*factor,Nz])
                y_min = y_new*factor
                y_max = min([(y_new+1)*factor,Ny])
                x_min = x_new*factor
                x_max = min([(x_new+1)*factor,Nx])
                array_new[z_new,y_new,x_new] = array_cp[z_min:z_max,y_min:y_max,x_min:x_max].sum()/(1.*mask[z_min:z_max,y_min:y_max,x_min:x_max].sum())
    return array_new

def downsample_spi(img,factor,mode="pick",ds_msk=True):
    import spimage
    img_array_new = downsample(img.image,factor,mode)
    img_new = spimage.sp_image_alloc(img_array_new.shape[1],img_array_new.shape[0],1)
    img_new.image[:,:] = img_array_new[:,:]
    if ds_msk:
        msk_array_new = downsample(img.mask,factor,mode)
        img_new.mask[:,:] = msk_array_new[:,:]
    else:
        img_new.mask[:,:] = 1
    return img_new



def center_of_mass(pattern,masking_threshold=None):
    if masking_threshold:
        mask = pylab.ones(shape=pattern.shape)
        mask[pattern<=masking_threshold] = -1
        mask[pattern>masking_threshold] = 0
        mask += 1
    else:
        mask = abs(pattern)
    #pylab.figure()
    #pylab.imshow(mask)
    #pylab.show()
    X,Y = pylab.meshgrid(pylab.arange(0,pattern.shape[1],1),pylab.arange(0,pattern.shape[0],1))
    x_center = (X*mask).sum()/mask.sum()
    y_center = (Y*mask).sum()/mask.sum()
    return [y_center,x_center]

def crop(pattern,cropLength,center='middle',bg=0,masking_threshold=None):
    if center == 'middle':
        x_center = (pattern.shape[1] - 1)/2.
        y_center = (pattern.shape[0] - 1)/2.
        temp = pattern.copy()
    else:
        if center == "center_of_mass":
            [y_center,x_center] = center_of_mass(pattern,masking_threshold)
        else:
            x_center = center[1]
            y_center = center[0]
        temp = recenter(pattern,x_center,y_center)

    x_start = (pattern.shape[1]-cropLength)/2
    y_start = (pattern.shape[1]-cropLength)/2
    x_stop = x_start+cropLength
    y_stop = y_start+cropLength

    patternCropped = pylab.ones(shape=(cropLength,cropLength),dtype=pattern.dtype)*bg
    patternCropped = temp[y_start:y_stop,x_start:x_stop]
    return patternCropped

# crop around centre to maximal centrosymmetric size
def crop_max_around_center(array2d,c0,c1):
    d0 = array2d.shape[0] - 1 - c0
    d1 = array2d.shape[1] - 1 - c1
    N_new = min([c0,c1,d0,d1])*2+1
    start0 = c0-(N_new-1)/2.0
    stop0 = start0+N_new
    start1 = c1-(N_new-1)/2.0
    stop1 = start1+N_new
    array2d_new = array2d[start0:stop0,start1:stop1]
    return array2d_new

def diameter_extrema(image):
    N = 32
    polimage = cartesian_to_polar(image,N)
    radii = pylab.zeros(N)
    for i in range(0,N):
        radii[i] = polimage[:,i].argmin()
    diameters = radii[:N/2]+radii[N/2:]
    diameter_max = diameters.max()
    diameter_min = diameters.min()
    return [diameter_min,diameter_max]

    
def shift(pattern,center=None):
    Ny = pattern.shape[0]
    Nx = pattern.shape[1]
    if not center:
        center = [Ny/2,Nx/2]
    center = pylab.array([int(pylab.ceil(center[0])),int(pylab.ceil(center[1]))])
    patternShifted = pylab.zeros(shape=pattern.shape)
    patternShifted[Ny-center[0]:,Nx-center[1]:] = pattern[:center[0],:center[1]]
    patternShifted[Ny-center[0]:,:Nx-center[1]] = pattern[:center[0],center[1]:]
    patternShifted[:Ny-center[0],Nx-center[1]:] = pattern[center[0]:,:center[1]]
    patternShifted[:Ny-center[0],:Nx-center[1]] = pattern[center[0]:,center[1]:]
    return patternShifted

def shift_back(pattern,center):
    Ny = pattern.shape[0]
    Nx = pattern.shape[1]
    new_center = [0,0]
    new_center[0] = Ny-center[0]-1
    new_center[1] = Nx-center[1]-1
    return shift(pattern,new_center)

def turncw(array2d):
    array2d_turned = pylab.zeros_like(array2d)
    for x in range(0,len(array2d[0])):
        temp_list=list(array2d[:,x])
        temp_list.reverse()
        array2d_turned[x,:] = pylab.array(temp_list,dtype=array2d.dtype).T
    return array2d_turned

def turnccw(array2d):
    array2d_turned = pylab.zeros(shape=(array2d.shape[1],array2d.shape[0]),dtype=array2d.dtype)
    N = len(array2d_turned)-1
    for x in range(0,len(array2d[0])):
        array2d_turned[N-x,:] = array2d[:,x].T
    return array2d_turned

def turn180(array2d,cx=None,cy=None):
    if cx == None and cy == None:
        array2d_turned = array2d.swapaxes(0,1)
    else:
        dcx = cx-(array2d.shape[1]-1)/2.
        dcy = cy-(array2d.shape[0]-1)/2.
        array2d_turned = array2d.copy()
        if dcx >= 0 and dcy >= 0: toturn = array2d_turned[2*dcy:,2*dcx:]
        elif dcx < 0 and dcy >= 0: toturn = array2d_turned[2*dcy:,:2*dcx]
        elif dcx >= 0 and dcy < 0: toturn = array2d_turned[:2*dcy,2*dcx:]
        elif dcx < 0 and dcy < 0: toturn = array2d_turned[:2*dcy,:2*dcx]
        toturn[:,:] = toturn.swapaxes(0,1)[:,:]
    return array2d_turned

def slide(img,dx=0,dy=0):
    if dx == 0 and dy == 0: imgout = img.copy()
    else:
        imgout=pylab.zeros_like(img)
        if dx > 0 and dy > 0: imgout[dy:,dx:] = img[:-dy,:-dx]
        elif dx < 0 and dy < 0: imgout[:dy,:dx] = img[-dy:,-dx:]
        elif dx < 0 and dy > 0: imgout[dy:,:dx] = img[:-dy,-dx:]
        elif dx > 0 and dy < 0: imgout[:dy,dx:] = img[-dy:,:-dx]
        elif dx > 0: imgout[:,dx:] = img[:,:-dx]
        elif dy > 0: imgout[dy:,:] = img[:-dy,:]
        elif dx < 0: imgout[:,:dx] = img[:,-dx:]
        elif dy < 0: imgout[:dy,:] = img[-dy:,:]
    return imgout
                                                
def turn180_(array2d,cx=None,cy=None):
    array2d_turned = turnccw(turnccw(array2d))
    dcx = cx-(array2d.shape[1]-1)/2.
    dcy = cy-(array2d.shape[0]-1)/2.
    array2d_turned = slide(array2d_turned,2*dcx,2*dcy)
    return array2d_turned

def horizontalmirr(array2d):
    array2d_mirrored = list(array2d.copy())
    array2d_mirrored.reverse()
    array2d_mirrored = pylab.array(array2d_mirrored)
    return array2d_mirrored

# only for Nx=Ny=Nz
def slice(dm3d,phi,theta,psi):
    #voxelcoordmatrix = get_voxel_coord_matrix(Nx,Ny,Nz)
    N = dm3d.shape[0]
    dm2dslice = pylab.zeros(N**2)
    X2,X1 = pylab.meshgrid(pylab.arange(-(N-1)/2.0,(N/2-1)/2.0,1.0),pylab.arange(-(N-1)/2.0,(N/2-1)/2.0,1.0))
    X0 = pylab.zeros_like(X1)
    coord_slice = pylab.array([X0,X1,X2])
    voxelcoordslice = _rotate_voxel_coord_slice(phi,theta,psi)
    for i in range(0,len(dm2dslice)):
        coord = voxelcoordslice[i]
        dm2dslice[i] = interpolate3d(dm3d,coord)
        
#def _get_voxel_coord_slice(coord_slice,phi,theta,psi):
     ### continue programming here


def interpolate3d(dm3d,coord,interpolation="linear"):
    x0 = coord[0]
    N0 = dm3d.shape[0]
    x1 = coord[1]
    N1 = dm3d.shape[1]
    x2 = coord[2]
    N2 = dm3d.shape[2]
    if x0 > N0-1 or x1 > N1-1 or x2 > N2-1 or x0 < 0 or x1 < 0 or x2 < 0:
        value = 0
    else:
        if interpolation == "linear":
            value = 0
            cases = [pylab.floor,lambda x: 1.0 + pylab.floor(x)]
            for x0_func in cases:
                for x1_func in cases:
                    for x2_func in cases:
                        value +=\
                            (1.0-abs(x0_func(x0) - x0))*\
                            (1.0-abs(x1_func(x1) - x1))*\
                            (1.0-abs(x2_func(x2) - x2))*\
                            dm3d[int(x0_func(x0)),int(x1_func(x1)),int(x2_func(x2))]
        elif interpolation == "nn":
            x0_rounded = pylab.round_(x0) 
            x1_rounded = pylab.round_(x1) 
            x2_rounded = pylab.round_(x2) 
            value = dm3d[int(x0_rounded),int(x1_rounded),int(x2_rounded)]
    return value

def pad_zeros(arr_orig,factor,shifted=False):
    arr = arr_orig.copy()
    Ny = arr.shape[0]
    Nx = arr.shape[1]
    Ny_new = int(round(Ny*factor)) 
    Nx_new = int(round(Nx*factor))
    Nx_d = Nx_new-Nx
    Ny_d = Ny_new-Ny
    arr_out = pylab.zeros(shape=(Ny_new,Nx_new),dtype=arr.dtype)
    if shifted:
        arr = pylab.fftshift(arr)
    arr_out[Ny_d/2:Ny_d/2+Ny,Nx_d/2:Nx_d/2+Nx] = arr[:,:]
    if shifted:
        arr_out = pylab.fftshift(arr_out)
    return arr_out
    
def depixelate(arr_orig,factor,shifted=False):
    arr = arr_orig.copy()
    farr = pylab.fft2(arr)
    if shifted:
        pfarr = pad_zeros(farr,factor,False)
    else:
        pfarr = pad_zeros(farr,factor,True)
    parr = pylab.ifft2(pylab.fftshift(pfarr))
    return parr

def splitx(img,x_split=None):
    if x_split == None:
        x_split = img.shape[1]/2-0.5
    p1 = img[:,:x_split+0.5]
    p2 = img[:,x_split+0.5:]
    return [p1,p2]

def splity(img,y_split=None):
    if y_split == None:
        y_split = img.shape[0]/2-0.5
    p1 = img[:y_split+0.5,:]
    p2 = img[y_split+0.5:,:]
    return [p1,p2]

# Assume vertical slit (x = slitposition)
def stitch(img_in,splitposition,cx,cy,dx,dy):
    [p1,p2] = splitx(img_in,splitposition)
    Ny_out = img_in.shape[0]+abs(dy)
    Nx_out = img_in.shape[1]+abs(dx)
    Nx_p1 = p1.shape[1]
    Ny_p = p1.shape[0]
    img_out = pylab.zeros(shape=(Ny_out,Nx_out))
    if dy>=0:
        img_out[:Ny_p,:Nx_p1] = p1[:,:]
        img_out[dy:,Nx_p1+dx:] = p2[:,:]
    else:
        img_out[-dy:,:Nx_p1] = p1[:,:]
        img_out[:Ny_p,Nx_p1+dx:] = p2[:,:]
    return img_out

def resize2d(arr2d,dX_old,dX_new,X_new=None):
    from scipy.interpolate import interp2d
    from numpy import linspace
    if not X_new: X_new = arr2d.shape[1]
    x,y = pylab.meshgrid(pylab.arange(0,arr2d.shape[0])*dX_old,pylab.arange(0,arr2d.shape[1])*dX_old)
    newfunc = interp2d(x,y,arr2d,fill_value=0.0,kind='linear')
    x_new = pylab.linspace(0,X_new*dX_new,dX_old/dX_new)
    y_new = pylab.linspace(0,X_new*dX_new,dX_old/dX_new)
    map2d_resized = newfunc(x_new,y_new)
    return map2d_resized

def interpolate3d(arr3d,factor):
    import enthought.mayavi.mlab as m
    N = arr3d.shape[0]
    arr3d_new = arr3d.copy()
    farr3d = pylab.fftn(arr3d_new)
    farr3d = pylab.fftshift(farr3d)
    farr3d_new = pylab.zeros(shape=(N*factor,N*factor,N*factor),dtype="complex")
    farr3d_new[round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N,
               round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N,
               round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N] = farr3d[:,:,:]
    farr3d = farr3d_new
    farr3d = pylab.fftshift(farr3d)
    arr3d_new = pylab.ifftn(farr3d)
    return arr3d_new

def smooth1d(arr1d,sm):
    N = 1+2*sm
    x = pylab.arange(0,N,1) - sm
    kernel = pylab.zeros(N)
    kernel = pylab.exp(x**2/(1.0*sm**2))
    arr1d_new = pylab.convolve(arr1d,kernel,'same')
    #print len(arr1d_new)
    #print len(arr1d)
    return arr1d_new

def smooth3d(arr3d,factor):
    import enthought.mayavi.mlab as m
    N = arr3d.shape[0]
    arr3d_new = arr3d.copy()
    farr3d = pylab.fftn(arr3d_new)
    farr3d = pylab.fftshift(farr3d)
    X,Y,Z = pylab.mgrid[-N/2:-N/2+N,-N/2:-N/2+N,-N/2:-N/2+N]
    #R = pylab.sqrt(X**2+Y**2+Z**2)
    farr3d[abs(X)>N/2] = 0
    farr3d[abs(Y)>N/2] = 0
    farr3d[abs(Y)>N/2] = 0
    farr3d = pylab.fftshift(farr3d)
    arr3d_new = pylab.ifftn(farr3d)
    return arr3d_new

def downsample3d_fourier(arr3d,factor):
    import enthought.mayavi.mlab as m
    N = arr3d.shape[0]
    N_new = round(N*factor/2.0)*2
    arr3d_new = arr3d.copy()
    farr3d = pylab.fftn(arr3d_new)
    farr3d = pylab.fftshift(farr3d)
    A = farr3d.sum()
    farr3d = farr3d[(N-N*factor)/2:(N-N*factor)/2+N_new,(N-N*factor)/2:(N-N*factor)/2+N_new,(N-N*factor)/2:(N-N*factor)/2+N_new]
    B = farr3d.sum()
    farr3d /= (N/(1.0*N_new))**3.0
    farr3d = pylab.fftshift(farr3d)
    arr3d_new = pylab.ifftn(farr3d)
    return arr3d_new

def downsample2d_fourier(arr2d,factor):
    N = arr2d.shape[0]
    N_new = round(N*factor/2.0)*2
    arr2d_new = arr2d.copy()
    farr2d = pylab.fftn(arr2d_new)
    farr2d = pylab.fftshift(farr2d)
    A = farr2d.sum()
    farr2d = farr2d[(N-N*factor)/2:(N-N*factor)/2+N_new,(N-N*factor)/2:(N-N*factor)/2+N_new]
    B = farr2d.sum()
    farr2d /= (N/(1.0*N_new))**2.0
    farr2d = pylab.fftshift(farr2d)
    arr2d_new = pylab.ifftn(farr2d)
    return arr2d_new

def interpolate2d(arr2d,factor):
    N = arr2d.shape[0]
    arr2d_new = arr2d.copy()
    farr2d = pylab.fftn(arr2d_new)
    farr2d = pylab.fftshift(farr2d)
    farr2d_new = pylab.zeros(shape=(N*factor,N*factor),dtype="complex")
    farr2d_new[round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N,
               round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N] = farr2d[:,:]
    farr2d = farr2d_new
    farr2d = pylab.fftshift(farr2d)
    arr2d_new = pylab.ifftn(farr2d)
    #pylab.figure()
    #pylab.imshow(pylab.log10(farr2d.real))
    #pylab.imshow(arr2d_new.real)
    #pylab.show()
    return arr2d_new

def interpolate1d(arr1d,factor):
    N = len(arr1d)
    arr1d_new = arr1d.copy()
    farr1d = pylab.fftn(arr1d_new)
    farr1d = pylab.fftshift(farr1d)
    farr1d_new = pylab.zeros(N*factor,dtype="complex")
    farr1d_new[round((N*factor-N)/2.0):round((N*factor-N)/2.0)+N] = farr1d[:]
    farr1d = farr1d_new
    farr1d = pylab.fftshift(farr1d)
    arr1d_new = pylab.ifftn(farr1d)*factor
    #pylab.figure()
    #pylab.imshow(pylab.log10(farr2d.real))
    #pylab.imshow(arr2d_new.real)
    #pylab.show()
    return arr1d_new

def put_besides(img1,img2):
    img3 = pylab.zeros(shape=(max([img1.shape[0],img2.shape[0]]),img1.shape[1]+img2.shape[1]))
    img3[:img1.shape[0],:img1.shape[1]] = img1[:,:]
    img3[:img2.shape[0]:,img1.shape[1]:] = img2[:,:]
    return img3

def radial_pixel_average(image,**kargs):
    if 'cx' in kargs: cx = kargs['cx']
    else: cx = (image.shape[1]-1)/2.0
    if 'cy' in kargs: cy = kargs['cy'] 
    else: cy = (image.shape[0]-1)/2.0
    R = get_R_and_Theta_map(image.shape[1],image.shape[0],cx,cy)[0]
    R[pylab.isfinite(image)==False] = -1
    radii = pylab.arange(R.min(),R.max()+1,1)
    if radii[0] == -1:
        radii = radii[1:]
    values = pylab.zeros_like(radii)
    for i in range(0,len(radii)):
        values[i] = image[R==radii[i]].mean()
    if 'rout' in kargs: return pylab.array([radii,values])
    else:return values

def recenter(I,cx,cy):
    dx = int(pylab.ceil(cx-(I.shape[1]-1)/2.))
    dy = int(pylab.ceil(cy-(I.shape[0]-1)/2.))
    I_recentered = pylab.zeros_like(I)
    if   dx > 0 and dy > 0:
        I_recentered[:-dy,:-dx] = I[dy:,dx:]
    elif dx > 0 and dy < 0:
        I_recentered[-dy:,:-dx] = I[:dy,dx:]
    elif dx < 0 and dy > 0:
        I_recentered[:-dy,-dx:] = I[dy:,:dx]
    elif dx < 0 and dy < 0:
        I_recentered[-dy:,-dx:] = I[:dy,:dx]
    elif dx == 0 and dy < 0:
        I_recentered[-dy:,:] = I[:dy,:]
    elif dx == 0 and dy > 0:
        I_recentered[:-dy,:] = I[dy:,:]
    elif dx > 0 and dy == 0:
        I_recentered[:,:-dx] = I[:,dx:]
    elif dx < 0 and dy == 0:
        I_recentered[:,-dx:] = I[:,:dx]
    elif dx == 0 and dy == 0:
        I_recentered[:,:] = I[:,:]
    return I_recentered

def downsample_position(position,downsampling):
    return (position-(downsampling-1)/2.)/(1.*downsampling)
