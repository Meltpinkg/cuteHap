# your_cython_file.pyx
import pysam
cimport pysam
from libc.stdlib cimport malloc, free
cimport cython

# 定义操作列表和参考变更操作
OPLIST=[
    pysam.CBACK,
    pysam.CDEL,
    pysam.CDIFF,
    pysam.CEQUAL,
    pysam.CHARD_CLIP,
    pysam.CINS,
    pysam.CMATCH,
    pysam.CPAD,
    pysam.CREF_SKIP,
    pysam.CSOFT_CLIP
]

# 定义改变表
CHANGETABLE={
    pysam.CMATCH:     (True,True),
    pysam.CINS:       (True,False),
    pysam.CDEL:       (False,True),
    pysam.CREF_SKIP:  (False,True),
    pysam.CPAD:       (False,False),
    pysam.CEQUAL:     (True,True),
    pysam.CDIFF:      (True,True)
}

cdef int MAX_OP = max(OPLIST)

# 定义静态数组
cdef bint* indelop = <bint*>malloc((MAX_OP + 1) * sizeof(bint))
cdef bint* refchangeop = <bint*>malloc((MAX_OP + 1) * sizeof(bint))


# 初始化数组
cdef void init():
    cdef int i
    for i in range(MAX_OP + 1):
        indelop[i] = (i == pysam.CDEL or i == pysam.CINS)
        refchangeop[i] = (i in CHANGETABLE) and CHANGETABLE[i][1]

cpdef bint is_indel(int op):
    return indelop[op]

cpdef bint is_ref_change(int op):
    return refchangeop[op]

def __dealloc__():
    free(indelop)
    free(refchangeop)

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef pre_mal(read, list Combine_sig_in_same_read_ins, list Combine_sig_in_same_read_del, list Combine_sig_in_same_read_tot, int sig_start, int min_siglength, int shift_ins_read):
    cdef int op
    cdef int oplen
    for op, oplen in read.cigartuples:
        if op != 2:#might be fixed later
            shift_ins_read += oplen
        if oplen >= min_siglength and indelop[op]:
            if op==2:
                Combine_sig_in_same_read_del.append([sig_start, oplen])
                Combine_sig_in_same_read_tot.append([sig_start, oplen, "DEL"])
                sig_start += oplen
            else:
                Combine_sig_in_same_read_ins.append([sig_start, oplen,
                    str(read.query_sequence[shift_ins_read-oplen:shift_ins_read])])
                Combine_sig_in_same_read_tot.append([sig_start, oplen, "INS"])
        else:
            # if op in RefChangeOp:
            if refchangeop[op]:
                sig_start += oplen


# 在模块初始化时调用
init()