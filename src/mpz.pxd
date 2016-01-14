cdef extern from "gmp.h":
    ctypedef long mp_limb_t
    ctypedef mp_limb_t * mp_ptr
    ctypedef struct __mpz_struct:
        int _mp_alloc
        int _mp_size
        mp_ptr _mp_d
    ctypedef __mpz_struct mpz_t[1]
    void mpz_set(mpz_t rop, mpz_t op)

    ctypedef __mpz_struct *mpz_ptr
    ctypedef __mpz_struct *mpz_srcptr

    void mpz_import (mpz_t rop, size_t count, int order, int size, int endian, size_t nails, void *op)
    void * mpz_export (void *rop, size_t *countp, int order, int size, int endian, size_t nails, mpz_t op)

    void mpz_init (mpz_t integer)
    signed long int mpz_get_si (mpz_t op)

    bint mpz_fits_slong_p (mpz_t op)
    size_t mpz_sizeinbase (mpz_t op, int base)

    void mpz_neg (mpz_t rop, mpz_t op)
    int mpz_sgn (mpz_t op)


