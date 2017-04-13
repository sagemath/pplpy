cdef extern from "Python.h":
    ctypedef struct PyObject:
        pass

cdef extern from "gmp.h":
    # gmp's integer numbers structs
    ctypedef struct __mpz_struct:
        pass
    ctypedef __mpz_struct mpz_t[1]
    ctypedef __mpz_struct *mpz_ptr
    ctypedef const __mpz_struct *mpz_srcptr
    
    # gmp's rationnal numbers structs
    ctypedef struct __mpq_struct:
        pass    
    ctypedef __mpq_struct mpq_t[1]
    ctypedef const __mpq_struct *mpq_srcptr

cdef extern from "gmpy2.h":
    cdef (PyObject *)GMPy_MPZ_New(void *)
    cdef int import_gmpy2()

cdef extern from "gmpy2_mpz.h":
    cdef mpz_t MPZ(PyObject *)

cdef extern from "gmpy2_types.h":
    ctypedef struct MPZ_Object:
        pass
    ctypedef struct CTXT_Object:
        pass

cdef extern from "gmpy2_convert_gmp.h":
    cdef (MPZ_Object *)GMPy_MPZ_From_PyIntOrLong(PyObject *obj, CTXT_Object *context)

cdef GMPy_MPZ_From_mpz(mpz_srcptr z)
cdef GMPy_MPQ_From_mpq(mpq_srcptr q)