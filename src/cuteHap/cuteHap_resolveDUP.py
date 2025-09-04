import numpy as np
import logging
from cuteHap.cuteHap_genotype import overlap_cover, assign_gt, bench_phase
import pickle

'''
*******************************************
                TO DO LIST
*******************************************
    1. Identify DP with samfile pointer;
    2. Add CIPOS, CILEN and/or CIEND;
    3. Determine (IM)PRECISE type.
    4. Filter DUP to improve INS FN rate.
*******************************************
'''

def resolution_DUP(path, chr, read_count, max_cluster_bias, sv_size, 
    bam_path, action, MaxSize, gt_round, sigs_index):
    if chr not in sigs_index["DUP"].keys():
        return (chr,[],[])
    # Initialization of some temporary variables
    semi_clusters_list = list()
    pre_pos = list()
    semi_dup_cluster = list()
    semi_dup_cluster.append([0, 0, '', 0])
    candidate_single_SV = list()

    with open("%s%s.pickle"%(path, "DUP"), 'rb') as f:
        f.seek(sigs_index["DUP"][chr])
        seqs=pickle.load(f)
    for seq in seqs:
        # (45461393, 45463990, 'S1S1/167429/ccs', 'DUP', 'chr15', 1)
        pos_1 = int(seq[0])
        pos_2 = int(seq[1])
        read_id = seq[2]
        hp = seq[5]
        
        if pos_1 - semi_dup_cluster[-1][0] > max_cluster_bias and abs(pos_2 - semi_dup_cluster[-1][1]) > max_cluster_bias:
            if len(semi_dup_cluster) >= read_count:
                if semi_dup_cluster[-1][0] == semi_dup_cluster[-1][1] == 0:
                    pass
                else:
                    semi_clusters_list.append(semi_dup_cluster)
                    pre_pos.append([max(semi_dup_cluster[0][0] - max_cluster_bias, 0), semi_dup_cluster[-1][0] + max_cluster_bias])
            semi_dup_cluster = []
            semi_dup_cluster.append([pos_1, pos_2, read_id, hp])
        else:
            if semi_dup_cluster[-1][0] == semi_dup_cluster[-1][1] == 0:
                semi_dup_cluster = []
                semi_dup_cluster.append([pos_1, pos_2, read_id, hp])
            else:
                semi_dup_cluster.append([pos_1, pos_2, read_id, hp])

    if len(semi_dup_cluster) >= read_count:
        if semi_dup_cluster[-1][0] == semi_dup_cluster[-1][1] == 0:
            pass
        else:
            semi_clusters_list.append(semi_dup_cluster)
            pre_pos.append([max(semi_dup_cluster[0][0] - max_cluster_bias, 0), semi_dup_cluster[-1][0] + max_cluster_bias])
    
    hp_dist_dict, overlap_reads_dict = bench_phase(path, chr, sigs_index, pre_pos)

    unphased_single_SV = list()
    phased_single_SV = list()
    for i in range(len(pre_pos)):
        hp_count = hp_dist_dict[i]
        overlap_reads = overlap_reads_dict[i]
        if hp_count[0]+abs(hp_count[1]-hp_count[2]) < 4*min(hp_count[1],hp_count[2]) or abs(hp_count[1]-hp_count[2]) < hp_count[0] < hp_count[1]+hp_count[2]:
        # if 1 == 1:
            a_phased_single_SV = generate_balance(semi_clusters_list[i], 
                            chr, 
                            "DUP", 
                            read_count,
                            gt_round,
                            max_cluster_bias,
                            sv_size,
                            MaxSize,
                            overlap_reads)
            phased_single_SV.extend(a_phased_single_SV)
        else:
            a_unphased_single_SV = generate_dup_cluster(semi_clusters_list[i], 
                                chr, 
                                read_count, 
                                sv_size, 
                                max_cluster_bias,
                                MaxSize,
                                gt_round)
            unphased_single_SV.extend(a_unphased_single_SV)

    gted_unphased_single_SV = call_gt(path, chr, unphased_single_SV, max_cluster_bias, sigs_index)
    
    logging.info("Finished %s:%s."%(chr, "DUP"))
    return (chr, phased_single_SV, gted_unphased_single_SV)

def generate_balance(semi_dup_cluster, chr, svtype, read_count, gt_round, max_cluster_bias, sv_size, MaxSize, overlap_reads):
    semi_hp_cluster = {0: [], 1: [], 2:[]}
    clustered_SV = dict()
    a_phased_single_SV = {0: [], 1: [], 2:[]}
    for item in semi_dup_cluster: # [pos_1, pos_2, read_id, hp]
        semi_hp_cluster[item[3]].append(item)
    for i in [0, 1, 2]:
        clustered_SV[i] = generate_dup_cluster(semi_hp_cluster[i], chr, read_count, sv_size, max_cluster_bias, MaxSize, gt_round) # read_count/2 if i>0 else read_count
        for origin_SV in clustered_SV[i]:
            # origin_SV: [chr, svtype, bp1, bp2, support_read, support_read_hp]
            # revised_single_SV: # [chr, svtype, bp1, dup_len(bp2-bp1), DV, DR, GT, GL, GQ, QUAL, support_read]
            if i > 0:
                revised_single_SV = revise_single_allele(origin_SV, max_cluster_bias, i, overlap_reads)
            else:
                revised_single_SV = [origin_SV[0],origin_SV[1],origin_SV[2],origin_SV[3]-origin_SV[2],len(origin_SV[4]),'.','.|.','.','.','.',origin_SV[4]]
            a_phased_single_SV[i].append(revised_single_SV)
    combine_list = combine_alleles(a_phased_single_SV, read_count)
    return combine_list

def generate_dup_cluster(semi_dup_cluster, chr, read_count, sv_size, max_cluster_bias, MaxSize, gt_round):
    # return [chr, svtype, bp1, bp2, support_read, support_read_hp]
    if len(semi_dup_cluster) == 0:
        return []
    # calculate support reads
    support_read = list(set([i[2] for i in semi_dup_cluster]))
    if len(support_read) < read_count:
        return []

    candidate_single_SV = list()
    semi_dup_cluster.sort(key = lambda x:x[1])
    allele_collect = []
    allele_collect.append([semi_dup_cluster[0]])
    last_len = semi_dup_cluster[0][1]
    for i in semi_dup_cluster[1:]:
        if i[1] - last_len > max_cluster_bias:
            allele_collect.append([])
        allele_collect[-1].append(i)
        last_len = i[1]
    for i in allele_collect:
        # i: [pos_1, pos_2, read_id, hp]
        support_read_hp = [0, 0, 0]
        support_read = set() # set([j[2] for j in i])
        for j in i:
            if j[2] not in support_read:
                support_read.add(j[2])
                support_read_hp[j[3]] += 1
        if len(support_read) < read_count:
            continue
        low_b = int(len(i)*0.4)
        up_b = int(len(i)*0.6)

        if low_b == up_b:
            breakpoint_1 = i[low_b][0]
            breakpoint_2 = i[low_b][1]
        else:
            breakpoint_1 = [i[0] for i in i[low_b:up_b]]
            breakpoint_2 = [i[1] for i in i[low_b:up_b]]
            breakpoint_1 = int(sum(breakpoint_1)/len(i[low_b:up_b]))
            breakpoint_2 = int(sum(breakpoint_2)/len(i[low_b:up_b]))


        if sv_size <= breakpoint_2 - breakpoint_1 <= MaxSize or (sv_size <= breakpoint_2 - breakpoint_1 and MaxSize == -1):
            candidate_single_SV.append([chr,
                                        'DUP', 
                                        breakpoint_1, 
                                        breakpoint_2,
                                        support_read,
                                        support_read_hp])
    return candidate_single_SV

def run_dup(args):
    return resolution_DUP(*args)

def call_gt(temporary_dir, chr, candidate_single_SV, max_cluster_bias, sigs_index):
    # reads_list = list() # [(start, end, 0, 'name'), ...]
    # candidate_single_SV: [chr, svtype, bp1, bp2, support_read, support_read_hp]
    # return [chr, svtype, bp1, dup_len(bp2-bp1), DV, DR, GT, GL, GQ, QUAL, support_read]
    if chr not in sigs_index["reads"].keys():
        return []
    readsfile = open("%sreads.pickle"%(temporary_dir), 'rb')
    readsfile.seek(sigs_index["reads"][chr])
    reads_list=pickle.load(readsfile)
    readsfile.close()
    svs_list = list()
    for item in candidate_single_SV:
        new_cluster_bias = min(max_cluster_bias, item[3] - item[2])
        svs_list.append((item[2], item[3]))

    iteration_dict, primary_num_dict, cover_dict, read_to_phase_dict = overlap_cover(svs_list, reads_list, 'C') # both key(sv idx), value(set(read id))
    read_id_dict = dict()
    phased_id_dict = dict()
    for i in range(len(candidate_single_SV)):

        read_id_dict[i] = candidate_single_SV[i][4]
        phased_id_dict[i] =  candidate_single_SV[i][5]
    assert len(cover_dict) == len(candidate_single_SV), "overlap length error"
    assign_list = assign_gt(iteration_dict, primary_num_dict, cover_dict, read_id_dict, phased_id_dict, read_to_phase_dict, candidate_single_SV, [])
    # [[DV, DR, GT, GL, GQ, QUAL] ...]
    assert len(candidate_single_SV) == len(assign_list), "assign error"
    candidate_single_SV_gt = list()
    # [chr, svtype, bp1, dup_len(bp2-bp1), DV, DR, GT, GL, GQ, QUAL, support_read]
    for i in range(len(candidate_single_SV)):
        candidate_single_SV_gt.append([candidate_single_SV[i][0], 
                                    candidate_single_SV[i][1], 
                                    str(candidate_single_SV[i][2]), 
                                    str(candidate_single_SV[i][3] - candidate_single_SV[i][2]), 
                                    str(len(candidate_single_SV[i][4])),
                                    str(assign_list[i][1]),
                                    str(assign_list[i][2]),
                                    str(assign_list[i][3]),
                                    str(assign_list[i][4]),
                                    str(assign_list[i][5]),
                                    ','.join(candidate_single_SV[i][4])])
    return candidate_single_SV_gt

def revise_single_allele(this_SV, max_cluster_bias, hp_idx, overlap_reads):
    # this_SV: [chr, svtype, bp1, bp2, support_read, support_read_hp]
    # return_SV: [chr, svtype, bp1, dup_len(bp2-bp1), DV, DR, GT, GL, GQ, QUAL, support_read]
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
    support_read_count = len(this_SV[4])
    ref_read_count = 0
    for query in read_count:
        if query not in this_SV[4]: # support_read_set
            ref_read_count += 1
    flag_hap = cal_bayes(support_read_count, ref_read_count)
    if hp_idx == 1:
        gt = '1|0' if flag_hap else '.|0'
    if hp_idx == 2:
        gt = '0|1' if flag_hap else '0|.'
    return [this_SV[0],this_SV[1],this_SV[2],this_SV[3]-this_SV[2],support_read_count,ref_read_count,gt,'.','.','.',this_SV[4]]

def combine_alleles(a_phased_single_SV, read_count):
    # a_phased_single_SV: {key:hp, value:list of [chr, svtype, bp1, dup_len(bp2-bp1), DV, DR, GT, GL, GQ, QUAL, support_read]}
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
                signalLen = int((mean_len_1 + mean_len_2) / 2)
                genotype = '1|1'
                this_SV = allele_hp1
                this_SV[2] = breakpointStart1
                this_SV[3] = signalLen
                this_SV[4] = allele_hp1[4] + allele_hp2[4]
                this_SV[6] = genotype
                this_SV[10] = allele_hp1[10].union(allele_hp2[10])
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