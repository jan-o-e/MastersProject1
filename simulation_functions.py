import numpy as np
from qutip import *
import matplotlib.pyplot as plt

E_LEVEL_FREQ = 10 * 2*np.pi
U_LEVEL_FREQ = 5 * 2*np.pi
G_LEVEL_FREQ = 2 * 2*np.pi

def run_reduced_2lvl_simulation(initial_state, run_time, num_steps, g0, detuning=0, gamma=0, kappa=0):

    # define characteristics of the atom-cavity system
    e_freq = E_LEVEL_FREQ
    g_freq = G_LEVEL_FREQ
    cav_freq = e_freq - g_freq - detuning * 2*np.pi
    num_photons = 1
    
    # define collapse operator associated with the cavity decay
    photon_c_op = Qobj([
        [0,1,0],
        [0,0,0],
        [0,0,0]
        ])

    # define the collapse operator associated with spontaneous decay of the excited state
    atom_c_op = Qobj([
        [0,0,1],
        [0,0,0],
        [0,0,0]
        ])

    # compile collapse operators into a list with appropriate coefficients
    c_ops = [np.sqrt(kappa)*photon_c_op, np.sqrt(gamma)*atom_c_op]

    # define Hamiltonian
    H = Qobj([
        [0,0,0],
        [0, g_freq + num_photons*cav_freq , -g0],
        [0, -g0, e_freq]])

    # define time grid
    t_steps = np.linspace(0, run_time, num_steps)

    # check normalisation
    norm = sum([abs(coeff)**2 for coeff in initial_state])
    if norm != 1:
        print('WARNING: unnormalised state')
    
    # calculate states that result from cavity dynamics
    result = mesolve(H, initial_state, t_steps, c_ops, [])

    # plot expectations of the measurements of each state
    fig, axes = plt.subplots(1,1)

    g_measure = Qobj([
        [0,0,0],
        [0,1,0],
        [0,0,0]
        ])
    e_measure = Qobj([
        [0,0,0],
        [0,0,0],
        [0,0,1]
        ])

    axes.plot(t_steps, expect(g_measure, result.states))
    axes.plot(t_steps, expect(e_measure, result.states))
    axes.set_xlabel('t')
    axes.set_ylabel('Probability')

    fig.legend(['g','e'])
    return result.states

def run_reduced_3lvl_simulation(initial_state, run_time, num_steps, g0, laser_profile, cav_detuning=0, gamma=0, kappa=0):

    # define characteristic frequencies of the atom-cavity system
    g_freq = G_LEVEL_FREQ
    e_freq = E_LEVEL_FREQ
    mid_freq = U_LEVEL_FREQ
    
    # define time grid
    t_steps = np.linspace(0, run_time, num_steps)

    # define time-independent part of tghe Hamiltonian
    H0 = Qobj([
        [0,0,0,0],
        [0,cav_detuning * 2*np.pi, 0, -g0],
        [0,0, e_freq - mid_freq, 0],
        [0,-g0, 0, 0]
    ])

    H1 = Qobj([
        [0,0,0,0],
        [0,0,0,0],
        [0,0,-1,-1/2],
        [0,0,-1/2,0]
    ])

    H = [H0, [H1, laser_profile]]

    cavity_c_op = Qobj([
        [0,1,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0]
    ])

    spont_c_op = Qobj([
        [0,0,0,1],
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0]
    ])

    c_ops = [np.sqrt(kappa)*cavity_c_op, np.sqrt(gamma)*spont_c_op]
    
    result = mesolve(H, initial_state, t_steps, c_ops, [])

    e_measure = fock_dm(4,3)
    u_measure = fock_dm(4,2)
    g_measure = fock_dm(4,1)

    # plot expectations of the measurements of each state
    fig, axes = plt.subplots(1,1)
    axes.plot(t_steps, expect(g_measure, result.states))
    axes.plot(t_steps, expect(u_measure, result.states))
    axes.plot(t_steps, expect(e_measure, result.states))
    
    axes.set_xlabel('t')
    axes.set_ylabel('Probability')

    fig.legend(['g','u','e'])
    return result.states
    