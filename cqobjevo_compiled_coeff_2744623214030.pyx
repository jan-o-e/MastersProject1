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

include 'c:/Users/ernst/Anaconda3/envs/qutip3.9/lib/site-packages/qutip/cy/complex_math.pxi'

cdef class CompiledStrCoeff(StrCoeff):
    cdef double wSTIRAP
    cdef double a0
    cdef double a1
    cdef double a2
    cdef double T
    cdef double omega

    def set_args(self, args):
        self.wSTIRAP=args['wSTIRAP']
        self.a0=args['a0']
        self.a1=args['a1']
        self.a2=args['a2']
        self.T=args['T']
        self.omega=args['omega']

    cdef void _call_core(self, double t, complex * coeff):
        cdef double wSTIRAP = self.wSTIRAP
        cdef double a0 = self.a0
        cdef double a1 = self.a1
        cdef double a2 = self.a2
        cdef double T = self.T
        cdef double omega = self.omega

        coeff[0] = omega*np.sin(wSTIRAP*t)**2
