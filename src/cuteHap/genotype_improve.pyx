from libc.stdlib cimport malloc, free
cimport cython


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef improve_overlap_cover(list svs_list, list reads_list, list sort_list):
    cdef int idx = 0
    for i in reads_list:
        sort_list.append([i[0], 1, idx, i[2], i[3]])
        sort_list.append([i[1], 2, idx, i[2], i[3]])
        idx += 1
    idx = 0
    for i in svs_list:
        sort_list.append([i[0], 3, idx])
        sort_list.append([i[1], 0, idx])
        idx += 1
    


