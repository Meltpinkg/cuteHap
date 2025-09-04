import numpy as np
import logging
from cuteHap.cuteHap_genotype import overlap_cover, assign_gt, bench_phase
import pickle

def resolution_INV(path, chr, svtype, read_count, max_cluster_bias, sv_size, 
    bam_path, action, MaxSize, gt_round, sigs_index):
    '''
    ************************************************************************
    Input file format
    ------------------------------------------------------------------------
    column	#1	#2	#3	    #4	#5  #6  #7
            INV	CHR	STRAND  BP1	BP2	ID  HP
    #1	inversion type
    #2	chromosome number
    #3	breakpoint_1 in each read
    #4	breakpoint_2 in each read
    #5	read ID
    ************************************************************************
    '''

    if chr not in sigs_index["INV"].keys():
        return (chr,[],[])
    # Initialization of some temporary variables
    semi_clusters_list = list()
    pre_pos = list()
    semi_inv_cluster = list()
    semi_inv_cluster.append([0,0,'',0])
    
    # Load inputs & cluster breakpoint from each signature read 

    with open("%s%s.pickle"%(path, "INV"), 'rb') as f:
        f.seek(sigs_index["INV"][chr])
        seqs=pickle.load(f)
    for seq in seqs:
        strand = seq[0]
        breakpoint_1_in_read = int(seq[1])
        breakpoint_2_in_read = int(seq[2])
        read_id = seq[3]
        hp = int(seq[6])

        if (breakpoint_1_in_read - semi_inv_cluster[-1][0] > max_cluster_bias and abs(breakpoint_2_in_read - semi_inv_cluster[-1][1]) > max_cluster_bias): # or strand != semi_inv_cluster[-1][3]:
            if len(semi_inv_cluster) >= read_count:
                if semi_inv_cluster[-1][0] == semi_inv_cluster[-1][1] == 0:
                    pass
                else:
                    semi_clusters_list.append(semi_inv_cluster)
                    pre_pos.append([max(semi_inv_cluster[0][0] - max_cluster_bias/2, 0), semi_inv_cluster[-1][0] + max_cluster_bias/2])
            semi_inv_cluster = []
            semi_inv_cluster.append([breakpoint_1_in_read, breakpoint_2_in_read, read_id, strand, hp])
        else:
            if semi_inv_cluster[-1][0] == semi_inv_cluster[-1][1] == 0:
                semi_inv_cluster = []
                semi_inv_cluster.append([breakpoint_1_in_read, breakpoint_2_in_read, read_id, strand, hp])
            else:
                semi_inv_cluster.append([breakpoint_1_in_read, breakpoint_2_in_read, read_id, strand, hp])

    if len(semi_inv_cluster) >= read_count:
        if semi_inv_cluster[-1][0] == semi_inv_cluster[-1][1] == 0:
            pass
        else:
            semi_clusters_list.append(semi_inv_cluster)
            pre_pos.append([max(semi_inv_cluster[0][0] - max_cluster_bias/2, 0), semi_inv_cluster[-1][0] + max_cluster_bias/2])

    hp_dist_dict, overlap_reads_dict = bench_phase(path, chr, sigs_index, pre_pos)
    unphased_single_SV = list()
    phased_single_SV = list()
    for i in range(len(pre_pos)):
        hp_count = hp_dist_dict[i]
        overlap_reads = overlap_reads_dict[i]
        if hp_count[0]+abs(hp_count[1]-hp_count[2]) < min(hp_count[1],hp_count[2]) or abs(hp_count[1]-hp_count[2]) < hp_count[0] < max(hp_count[1],hp_count[2]):
            a_phased_single_SV = generate_balance(semi_clusters_list[i], 
                            chr, 
                            svtype, 
                            read_count,
                            gt_round,
                            max_cluster_bias,
                            path,
                            sv_size,
                            MaxSize,
                            overlap_reads,
                            sigs_index)
            phased_single_SV.extend(a_phased_single_SV)
        else:
            a_unphased_single_SV = generate_semi_inv_cluster(semi_clusters_list[i], 
                                chr, 
                                svtype, 
                                read_count, 
                                sv_size, 
                                max_cluster_bias,
                                MaxSize,
                                gt_round)
            unphased_single_SV.extend(a_unphased_single_SV)

    gted_unphased_single_SV = call_gt(path, chr, unphased_single_SV, max_cluster_bias, sigs_index)
    logging.info("Finished %s:%s."%(chr, "INV"))
    return (chr, phased_single_SV, gted_unphased_single_SV)

def generate_balance(semi_inv_cluster, chr, svtype, read_count, gt_round, max_cluster_bias, path, sv_size, MaxSize, overlap_reads, sigs_index):
    semi_hp_cluster = {0: [], 1: [], 2:[]}
    clustered_SV = dict()
    a_phased_single_SV = {0: [], 1: [], 2:[]}
    for item in semi_inv_cluster: # [breakpoint_1_in_read, breakpoint_2_in_read, read_id, strand, hp]
        semi_hp_cluster[item[4]].append(item)
    for i in [0, 1, 2]:
        clustered_SV[i] = generate_semi_inv_cluster(semi_hp_cluster[i], chr, svtype, read_count, sv_size, max_cluster_bias, MaxSize, gt_round) # read_count/2 if i>0 else read_count
        for origin_SV in clustered_SV[i]:
            if i > 0:
                revised_single_SV = revise_single_allele(origin_SV, max_cluster_bias, i, overlap_reads)
            else:
                revised_single_SV = [origin_SV[0],origin_SV[1],origin_SV[2],origin_SV[3],origin_SV[4],'.','.|.',origin_SV[5],'.','.','.',origin_SV[6]]
            a_phased_single_SV[i].append(revised_single_SV)
    
    combine_list = combine_alleles(a_phased_single_SV, read_count)
    return combine_list
            

def generate_semi_inv_cluster(semi_inv_cluster, chr, svtype, read_count, sv_size, max_cluster_bias, MaxSize, gt_round):
    # semi_inv_cluster: [bp1, bp2, read_id, strand, hp]
    # return: [chr, svtype, bp1, inv_len, DV, strand, read_id, bp2, read_id_hp]
    if len(semi_inv_cluster) == 0:
        return []
    candidate_single_SV = list()
    strand = semi_inv_cluster[0][-1]

    read_id = [i[2] for i in semi_inv_cluster]
    support_read = len(list(set(read_id)))
    if support_read < read_count:
        return []

    inv_cluster_b2 = sorted(semi_inv_cluster, key = lambda x:x[1])

    last_bp = inv_cluster_b2[0][1]
    temp_count = 1
    temp_sum_b1 = inv_cluster_b2[0][0]
    temp_sum_b2 = last_bp

    # max_sum = 0
    temp_id = dict()
    temp_id_hp = [0, 0, 0]
    temp_id_hp[inv_cluster_b2[0][4]] += 1
    temp_id[inv_cluster_b2[0][2]] = 1 # 0??

    for i in inv_cluster_b2[1:]:
        if i[1] - last_bp > max_cluster_bias:
            if temp_count >= read_count:
                max_count_id = len(temp_id)

                breakpoint_1 = round(temp_sum_b1 / temp_count)
                breakpoint_2 = round(temp_sum_b2 / temp_count)
                inv_len = breakpoint_2 - breakpoint_1
                if inv_len >= sv_size and max_count_id >= read_count:
                    if inv_len <= MaxSize or MaxSize == -1:
                        candidate_single_SV.append([chr, 
                                                    svtype, 
                                                    breakpoint_1, 
                                                    inv_len, 
                                                    max_count_id,
                                                    strand,
                                                    set(temp_id.keys()),
                                                    breakpoint_2,
                                                    temp_id_hp])

            temp_id = dict()
            temp_id_hp = [0, 0, 0]
            temp_count = 1
            temp_sum_b1 = i[0]
            temp_sum_b2 = i[1]
            temp_id[i[2]] = 1
            temp_id_hp[i[4]] += 1
        else:
            if i[2] not in temp_id:
                temp_id[i[2]] = 1
                temp_id_hp[i[4]] += 1
            else:
                temp_id[i[2]] += 1
            temp_count += 1
            temp_sum_b1 += i[0]
            temp_sum_b2 += i[1]
        last_bp = i[1]
    if temp_count >= read_count:
        max_count_id = len(temp_id)
        breakpoint_1 = round(temp_sum_b1 / temp_count)
        breakpoint_2 = round(temp_sum_b2 / temp_count)
        inv_len = breakpoint_2 - breakpoint_1
        if inv_len >= sv_size and max_count_id >= read_count:
            if inv_len <= MaxSize or MaxSize == -1:
                candidate_single_SV.append([chr, 
                                            svtype, 
                                            breakpoint_1, 
                                            inv_len, 
                                            max_count_id,
                                            strand,
                                            set(temp_id.keys()),
                                            breakpoint_2,
                                            temp_id_hp])
    return candidate_single_SV

def run_inv(args):
    return resolution_INV(*args)

def call_gt(temporary_dir, chr, candidate_single_SV, max_cluster_bias, sigs_index):
    # reads_list = list() # [(start, end, 0, 'name'), ...]
    # candidate_single_SV: [chr, svtype, bp1, inv_len, DV, strand, read_id, bp2, read_id_hp]
    # return [chr, svtype, bp1, inv_len, DV, DR, GT, strand, GL, GQ, QUAL, read_id]
    
    if chr not in sigs_index["reads"].keys():
        return []
    readsfile = open("%sreads.pickle"%(temporary_dir), 'rb')
    readsfile.seek(sigs_index["reads"][chr])
    reads_list=pickle.load(readsfile)
    readsfile.close()
    svs_list = list()
    for item in candidate_single_SV:
        svs_list.append((max(item[2] - max_cluster_bias/2, 0), item[2] + max_cluster_bias/2))
        svs_list.append((max(item[7] - max_cluster_bias/2, 0), item[7] + max_cluster_bias/2))
    iteration_dict, primary_num_dict, cover_dict1, read_to_phase_dict = overlap_cover(svs_list, reads_list, 'C') # both key(sv idx), value(set(read id))
    assert len(cover_dict1) == 2 * len(candidate_single_SV), "overlap length error"
    cover_dict = dict()
    read_id_dict = dict()
    phased_id_dict = dict()
    for i in range(len(candidate_single_SV)):
        cover_dict[i] = cover_dict1[i*2]
        for item in cover_dict1[i*2+1]:
            cover_dict[i].add(item)

        read_id_dict[i] = candidate_single_SV[i][6]
        phased_id_dict[i] = candidate_single_SV[i][8]
    assert len(cover_dict) == len(candidate_single_SV), "overlap length error"

    assign_list = assign_gt(iteration_dict, primary_num_dict, cover_dict, read_id_dict, phased_id_dict, read_to_phase_dict, candidate_single_SV, [])
    # [[DV, DR, GT, GL, GQ, QUAL] ...]
    assert len(candidate_single_SV) == len(assign_list), "assign error"
    candidate_single_SV_gt = list()
    for i in range(len(candidate_single_SV)):
        candidate_single_SV_gt.append([candidate_single_SV[i][0], # chr
                                    candidate_single_SV[i][1], # svtype
                                    str(int(candidate_single_SV[i][2])), # bp1
                                    str(int(candidate_single_SV[i][3])), # inv_len
                                    str(candidate_single_SV[i][4]), # read_count DV
                                    str(assign_list[i][1]), # DR
                                    str(assign_list[i][2]), # GT
                                    candidate_single_SV[i][5], # strand
                                    str(assign_list[i][3]), # GL
                                    str(assign_list[i][4]), # GQ
                                    str(assign_list[i][5]), # QUAL
                                    ','.join(candidate_single_SV[i][6])]) # read_id
    return candidate_single_SV_gt

def revise_single_allele(this_SV, max_cluster_bias, hp_idx, overlap_reads):
    # this_SV: [chr, svtype, bp1, inv_len, DV, strand, read_id, bp2, read_id_hp]
    # return_SV: [chr, svtype, bp1, inv_len, DV, DR, GT, strand, GL, GQ, QUAL, read_id]
    def cal_bayes(kdv, ndr):
        err = 0.1
        prior = 0.7022
        ht1 = np.float64(pow((1-err), kdv)*pow(err, ndr)*prior)
        ht2 = np.float64(pow((1-err), ndr)*pow(err, kdv)*(1-prior))
        return True if ht1 > ht2 else False
    
    search_start = max(this_SV[2] - max_cluster_bias, 0)
    search_end = this_SV[2] + max_cluster_bias
    
    read_count = set()
    for read in overlap_reads:
        if read[2] == 1:
            if read[0] <= search_start and read[1] >= search_end and read[5] == hp_idx:
                # cover
                read_count.add(read[3])
    support_read_count = this_SV[4]
    ref_read_count = 0
    for query in read_count:
        if query not in this_SV[6]: # support_read_set
            ref_read_count += 1
    flag_hap = cal_bayes(support_read_count, ref_read_count)
    # this_SV: [chr, svtype, bp1, inv_len, DV, strand, read_id, bp2, read_id_hp]
    # return_SV: [chr, svtype, bp1, inv_len, DV, DR, GT, strand, GL, GQ, QUAL, read_id]
    if hp_idx == 1:
        gt = '1|0' if flag_hap else '.|0'
    if hp_idx == 2:
        gt = '0|1' if flag_hap else '0|.'
    return [this_SV[0],this_SV[1],this_SV[2],this_SV[3],this_SV[4],ref_read_count,gt,this_SV[5],'.','.','.',this_SV[6]]
    
def combine_alleles(a_phased_single_SV, read_count):
    result_list = list()
    used_hp2_id = set()
    used_hp0_read = set()
    for allele_hp1 in a_phased_single_SV[1]:
        used_hp1_id = False
        for allele_hp2_idx in range(len(a_phased_single_SV[2])):
            if allele_hp2_idx in used_hp2_id:
                continue
            allele_hp2 = a_phased_single_SV[2][allele_hp2_idx]
            mean_len_1 = allele_hp1[3]
            mean_len_2 = allele_hp2[3]
            if min(mean_len_1, mean_len_2) / max(mean_len_1, mean_len_2) > 0.9: # a hom SV
                breakpointStart1 = (allele_hp1[2] + allele_hp2[2]) / 2
                breakpointStart2 = (allele_hp1[7] + allele_hp2[7]) / 2
                signalLen = int((mean_len_1 + mean_len_2) / 2)
                genotype = '1|1'
                this_SV = allele_hp1
                this_SV[2] = breakpointStart1
                this_SV[3] = signalLen
                this_SV[4] = allele_hp1[4] + allele_hp2[4]
                this_SV[6] = genotype
                this_SV[11] = allele_hp1[11].union(allele_hp2[11])
                used_hp1_id = True
                used_hp2_id.add(allele_hp2_idx)
                result_list.append(this_SV)
                break
        ### ADD HP=0 IN FUTURE
        if not used_hp1_id: # has not been used for all allele_hp2
            if allele_hp1[4] >= int(read_count / 2):
                result_list.append(allele_hp1)
    for allele_hp2_idx in range(len(a_phased_single_SV[2])):
        if allele_hp2_idx not in used_hp2_id:
            allele_hp2 = a_phased_single_SV[2][allele_hp2_idx]
            if allele_hp2[4] >= int(read_count / 2):
                result_list.append(allele_hp2)
    for allele_hp0 in a_phased_single_SV[0]:
        if allele_hp0[4] >= read_count:
            result_list.append(allele_hp0)
    return result_list
