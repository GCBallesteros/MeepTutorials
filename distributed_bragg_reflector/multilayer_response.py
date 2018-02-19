import numpy as np

def multilayer_response(t_layers, index_layers, index_left, index_right, freq):
    """Response of multilayer structure to normal incident plane wave.
    
    Wave impinges from the left.
    
    Uses meep frequency units.
    """ 
    n_layers = len(t_layers)
    if n_layers != len(index_layers):
        print('Number of t_layers and index_layers is different!')
        raise
        
    T = np.eye(2)
    
    # 1) Compute T-matrix for medium on left and first layer
    n1 = index_layers[0]
    T_left = np.array([[n1 + index_left, n1 - index_left],
                       [n1 - index_left, n1 + index_left]]) / (2 * n1)
    T = np.dot(T, T_left)
    
    for ii in np.arange(0, n_layers):
        # Compute propagation matrix
        phase_prop = 2 * np.pi * freq * t_layers[ii] * index_layers[ii] ;
        T_prop = np.array([[np.exp(-1j * phase_prop), 0],
                           [0, np.exp(1j * phase_prop)]])
        
        # For last layer we only propagate
        if ii == (n_layers-1):
            n1 = index_layers[-1]
            T_interface = np.array([[n1 + index_left, n1 - index_left],
                                    [n1 - index_left, n1 + index_left]]) / (2 * index_right)
        else:
            n1 = index_layers[ii]
            n2 = index_layers[ii+1]
            T_interface = np.array([[n1 + n2, n2 - n1],
                                    [n2 - n1, n1 + n2]]) / (2 * n2)
            
        T = np.dot(T, np.dot(T_prop, T_interface))
    
    R = np.abs(T[0, 1]/T[0, 0])**2
    T = np.abs(1./T[0, 0])**2
    return (R, T)