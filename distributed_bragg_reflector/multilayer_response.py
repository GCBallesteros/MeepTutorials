import numpy as np

def _T_interface(ni, nj, cos_i, cos_j, pol):
    if pol == 's':
        rij = (ni*cos_i - nj*cos_j)/(ni*cos_i + nj*cos_j)
        tij = (1 + rij)
    elif pol == 'p':
        rij = (nj*cos_i - ni*cos_j)/(nj*cos_i + ni*cos_j)
        tij = (1 + rij)*(ni/nj)
    
    return np.array([[1, rij],
                     [rij, 1]])/tij


def _T_propagation(phase):
    return np.array([[np.exp(-1j * phase), 0],
                     [0, np.exp(1j * phase)]])


def multilayer_response(pol, t_layers, index_layers, freq, theta_in):
    """Response of multilayer structure to normal incident plane wave.
    
    Wave impinges from the left.
    
    Arguments:
      pol (str): 's' or 'p'
      t_layers: Thickness of layers including exterior infinite media.
                The exterior layer thickness is irrelevant.
      index_layers: Index of refraction of all layers. Includes exterior
                    semiinfinite media.
      freq: Frequency in Meep units, i.e. f = 1/lambda
      theta_in: Angle of incidence
    
    Returns:
      Tuple (R, T) with the power Transmitted and Reflected.
    """
    n_layers = len(t_layers)
    if n_layers != len(index_layers):
        print('Number of t_layers and index_layers is different!')
        raise
        
    # Exterior layers correspond to infinite media. So we want
    # the propagation matrix to be the identity
    t_layers[0] = 0
    t_layers[-1] = 0
    
    index_layers = np.asarray(index_layers)
    t_layers = np.asarray(t_layers)
    
    # Compute cos(theta_i) which need for Fresnel Coeff
    # and interface T matrix.
    sin_thetas = (index_layers[0]/index_layers)*np.sin(theta_in)
    cos_thetas = np.sqrt(1-sin_thetas**2)
    
    # Power factor to convert from t coefficient to transmittance
    power_factor = (index_layers[-1]/index_layers[0]) *\
                   (cos_thetas[-1]/cos_thetas[0])
    
    T = np.eye(2)
    for ii in np.arange(0, n_layers-1):
        # Compute propagation matrix
        phase_prop = 2 * np.pi * freq * t_layers[ii] * index_layers[ii] * cos_thetas[ii]
        T_prop = _T_propagation(phase_prop)

        # Compute interface matrix
        n1 = index_layers[ii]
        n2 = index_layers[ii+1]
        cos_1 = cos_thetas[ii]
        cos_2 = cos_thetas[ii+1]
        T_ij = _T_interface(n1, n2, cos_1, cos_2, pol)
        
        # Update system transfer matrix
        T = np.dot(T, np.dot(T_prop, T_ij))
    
    R = np.abs(T[0, 1]/T[0, 0])**2
    T =  power_factor * np.abs(1./T[0, 0])**2

    return (R,T)
    

if __name__ == '__main__':
    # Silicon Dioxide on Silicon Wafer Example
    import matplotlib
    matplotlib.rcParams["backend"] = "TkAgg"
    import matplotlib.pyplot as plt
    n_Si = 3.5
    n_SiO2 = 1.5

    index_layers = [1, n_SiO2, n_Si]
    t_layers = [0, 220, 0]

    lambda_max = 400
    lambda_min = 1000
    freq = 1./np.linspace(lambda_min, lambda_max,400)

    (theory_reflection, theory_transmission) = multilayer_response('s',
                                                                   t_layers,
                                                                   index_layers,
                                                                   freq, 0*np.pi/180)

    plt.figure()
    plt.plot(1/freq, 10*np.log10(theory_reflection), label='Reflection')
    plt.plot(1/freq, 10*np.log10(theory_transmission), label='Tranmission')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Power (dB)')
    plt.legend()
    plt.show()