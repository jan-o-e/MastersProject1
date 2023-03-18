#!python
#cython: language_level=3
# This file is generated automatically by QuTiP.

import numpy as np
cimport numpy as np
import scipy.special as spe
cimport cython
np.import_array()
cdef extern from "numpy/arrayobject.h" nogil:
    void PyDataMem_NEW_ZEROED(size_t size, size_t elsize)
    void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)
from qutip.cy.spmatfuncs cimport spmvpy
from qutip.cy.inter cimport _spline_complex_t_second, _spline_complex_cte_second
from qutip.cy.inter cimport _spline_float_t_second, _spline_float_cte_second
from qutip.cy.inter cimport _step_float_cte, _step_complex_cte
from qutip.cy.inter cimport _step_float_t, _step_complex_t
from qutip.cy.interpolate cimport (interp, zinterp)
from qutip.cy.cqobjevo_factor cimport StrCoeff
from qutip.cy.cqobjevo cimport CQobjEvo
from qutip.cy.math cimport erf, zerf
from qutip.qobj import Qobj
cdef double pi = 3.14159265358979323

include '/opt/anaconda3/envs/qutip3.9/lib/python3.9/site-packages/qutip/cy/complex_math.pxi'

cdef class CompiledStrCoeff(StrCoeff):
    cdef double g0
    cdef double omega
    cdef double delta
    cdef double wSTIRAP
    cdef double T

    def set_args(self, args):
        self.g0=args['g0']
        self.omega=args['omega']
        self.delta=args['delta']
        self.wSTIRAP=args['wSTIRAP']
        self.T=args['T']

    cdef void _call_core(self, double t, complex * coeff):
        cdef double g0 = self.g0
        cdef double omega = self.omega
        cdef double delta = self.delta
        cdef double wSTIRAP = self.wSTIRAP
        cdef double T = self.T

        coeff[0] = g0
        coeff[1] = omega*np.sin(wSTIRAP*t)**2
        coeff[2] = delta
