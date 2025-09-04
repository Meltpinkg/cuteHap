import numpy as np
from cuteHap.cuteHap_genotype import cal_CIPOS, overlap_cover, assign_gt_mosaic, bench_phase
import logging
import pickle
import time

def resolution_INV(path, chr, svtype, read_count, max_cluster_bias, sv_size, 
    bam_path, action, MaxSize, sigs_index):
    '''
    cluster INV
    ************************************************************************
    path:	INV.sigs
    chr:	chromosome id
    svtype:	<INV>
    
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
        # seq: ['++', 16975073, 246818027, 'S2S1/176547/ccs', 'INV', 'chr1', 2]
        strand = seq[0]
        breakpoint_1_in_read = int(seq[1])
        breakpoint_2_in_read = int(seq[2])
        read_id = seq[3]
        hp = int(seq[6])

        if breakpoint_1_in_read - semi_inv_cluster[-1][0] > max_cluster_bias and abs(breakpoint_2_in_read - semi_inv_cluster[-1][1]) > max_cluster_bias or strand != semi_inv_cluster[-1][3]:
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
        read_count = max(int(sum(hp_count) * 0.05), 4)
        if hp_count[0]+abs(hp_count[1]-hp_count[2]) < min(hp_count[1],hp_count[2]) or abs(hp_count[1]-hp_count[2]) < hp_count[0] < max(hp_count[1],hp_count[2]):
        # if 1 == 0:
            a_phased_single_SV = generate_inv_balance(semi_clusters_list[i], 
                            chr, 
                            svtype, 
                            read_count,
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
                                MaxSize)
            unphased_single_SV.extend(a_unphased_single_SV)

    gted_unphased_single_SV = call_gt_inv(path, chr, unphased_single_SV, max_cluster_bias, sigs_index)
    logging.info("Finished %s:%s."%(chr, "INV"))
    return (chr, phased_single_SV, gted_unphased_single_SV)

def generate_inv_balance(semi_inv_cluster, chr, svtype, read_count, max_cluster_bias, path, sv_size, MaxSize, overlap_reads, sigs_index):
    semi_hp_cluster = {0: [], 1: [], 2:[]}
    clustered_SV = dict()
    a_phased_single_SV = {0: [], 1: [], 2:[]}
    for item in semi_inv_cluster: # [breakpoint_1_in_read, breakpoint_2_in_read, read_id, strand, hp]
        semi_hp_cluster[item[4]].append(item)
    for i in [0, 1, 2]:
        clustered_SV[i] = generate_semi_inv_cluster(semi_hp_cluster[i], chr, svtype, read_count, sv_size, max_cluster_bias, MaxSize) # read_count/2 if i>0 else read_count
        for origin_SV in clustered_SV[i]:
            # origin_SV: [chr, svtype, bp1, inv_len, DV, strand, read_id, bp2, read_id_hp]
            # revised_single_SV: [chr, svtype, bp1, inv_len, DV, DR, GT, strand, GL, GQ, QUAL, read_id]
            if i > 0:
                revised_single_SV = revise_single_allele(origin_SV, max_cluster_bias, i, overlap_reads)
            else:
                revised_single_SV = [origin_SV[0],origin_SV[1],origin_SV[2],origin_SV[3],origin_SV[4],'.','.|.',origin_SV[5],'.','.','.',origin_SV[6]]
            if '1' not in revised_single_SV[6]:
                a_phased_single_SV[i].append(revised_single_SV)
    combine_list = combine_alleles(a_phased_single_SV, read_count)
    return combine_list

def generate_semi_inv_cluster(semi_inv_cluster, chr, svtype, read_count, sv_size, max_cluster_bias, MaxSize):
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

    # breakpoint_1 = np.mean(breakpoint_1_candidate)
    last_bp = inv_cluster_b2[0][1]
    temp_count = 1
    # max_count = 0
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

def call_gt_inv(temporary_dir, chr, candidate_single_SV, max_cluster_bias, sigs_index):
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

    assign_list = assign_gt_mosaic(iteration_dict, primary_num_dict, cover_dict, read_id_dict, phased_id_dict, read_to_phase_dict)
    # [[DV, DR, GT, GL, GQ, QUAL] ...]
    assert len(candidate_single_SV) == len(assign_list), "assign error"
    candidate_single_SV_gt = list()
    for i in range(len(candidate_single_SV)):
        if '1' in str(assign_list[i][2]):
            continue
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
    # search_start = max(this_SV[2] - max_cluster_bias, 0)
    # search_end = this_SV[2] + max_cluster_bias
    read_count = set()
    for read in overlap_reads:
        if read[2] == 1:
            # if read[0] <= search_start and read[1] >= search_end and read[5] == hp_idx:
            if max(read[0], search_start) < min(read[1], search_end) and read[5] == hp_idx:
                # cover
                read_count.add(read[3])
    support_read_count = this_SV[4]
    ref_read_count = 0
    for query in read_count:
        if query not in this_SV[6]: # support_read_set
            ref_read_count += 1
    flag_hap = True if support_read_count / (support_read_count + ref_read_count) > 0.7 else False
    # this_SV: [chr, svtype, bp1, inv_len, DV, strand, read_id, bp2, read_id_hp]
    # return_SV: [chr, svtype, bp1, inv_len, DV, DR, GT, strand, GL, GQ, QUAL, read_id]
    if hp_idx == 1:
        gt = '1|0' if flag_hap else '.|0'
    if hp_idx == 2:
        gt = '0|1' if flag_hap else '0|.'
    return [this_SV[0],this_SV[1],this_SV[2],this_SV[3],this_SV[4],ref_read_count,gt,this_SV[5],'.','.','.',this_SV[6]]
    
def combine_alleles(a_phased_single_SV, read_count):
    # a_phased_single_SV: {key:hp, value:list of [chr, svtype, bp1, inv_len, DV, DR, GT, strand, GL, GQ, QUAL, read_id]}
    result_list = list()
    used_hp2_id = set()
    used_hp0_id = set()
    for allele_hp1 in a_phased_single_SV[1]:
        used_hp1_id = False
        for allele_hp2_idx in range(len(a_phased_single_SV[2])):
            if allele_hp2_idx in used_hp2_id:
                continue
            allele_hp2 = a_phased_single_SV[2][allele_hp2_idx]
            mean_len_1 = allele_hp1[3]
            mean_len_2 = allele_hp2[3]
            if min(mean_len_1, mean_len_2) / max(mean_len_1, mean_len_2) > 0.9: # a hom SV
                breakpointStart1 = int((allele_hp1[2] + allele_hp2[2]) / 2)
                # breakpointStart2 = (allele_hp1[3] + allele_hp2[3]) / 2
                signalLen = (mean_len_1 + mean_len_2) / 2
                genotype = allele_hp1[6][0] + '|' + allele_hp2[6][2]
                this_SV = allele_hp1
                this_SV[2] = breakpointStart1
                this_SV[3] = int(signalLen)
                this_SV[4] = allele_hp1[4] + allele_hp2[4]
                this_SV[6] = genotype
                this_SV[11] = allele_hp1[11].union(allele_hp2[11])
                used_hp1_id = True
                used_hp2_id.add(allele_hp2_idx)
                result_list.append(this_SV)
                break
        for allele_hp0_idx in range(len(a_phased_single_SV[0])):
            if allele_hp0_idx in used_hp0_id:
                continue
            allele_hp0 = a_phased_single_SV[0][allele_hp0_idx]
            mean_len_1 = allele_hp1[3]
            mean_len_2 = allele_hp0[3]
            if min(mean_len_1, mean_len_2) / max(mean_len_1, mean_len_2) > 0.9: # a hom SV
                breakpointStart1 = int((allele_hp1[2] + allele_hp0[2]) / 2)
                signalLen = (mean_len_1 + mean_len_2) / 2
                genotype = allele_hp1[6][0] + '|.'
                this_SV = allele_hp1
                this_SV[2] = breakpointStart1
                this_SV[3] = int(signalLen)
                this_SV[4] = allele_hp1[4] + allele_hp0[4]
                this_SV[6] = genotype
                this_SV[11] = allele_hp1[11].union(allele_hp0[11])
                used_hp1_id = True
                used_hp0_id.add(allele_hp0_idx)
                result_list.append(this_SV)
                break
        ### ADD HP=0 IN FUTURE
        if not used_hp1_id: # has not been used for all allele_hp2
            if allele_hp1[4] >= read_count:
                result_list.append(allele_hp1)
    for allele_hp2_idx in range(len(a_phased_single_SV[2])):
        allele_hp2 = a_phased_single_SV[2][allele_hp2_idx]
        for allele_hp0_idx in range(len(a_phased_single_SV[0])):
            if allele_hp0_idx in used_hp0_id:
                continue
            allele_hp0 = a_phased_single_SV[0][allele_hp0_idx]
            mean_len_1 = allele_hp2[3]
            mean_len_2 = allele_hp0[3]
            if min(mean_len_1, mean_len_2) / max(mean_len_1, mean_len_2) > 0.9: # a hom SV
                breakpointStart1 = int((allele_hp2[2] + allele_hp0[2]) / 2)
                signalLen = (mean_len_1 + mean_len_2) / 2
                genotype = '.|' + allele_hp0[6][2]
                this_SV = allele_hp2
                this_SV[2] = breakpointStart1
                this_SV[3] = int(signalLen)
                this_SV[4] = allele_hp2[4] + allele_hp0[4]
                this_SV[6] = genotype
                this_SV[11] = allele_hp2[11].union(allele_hp0[11])
                used_hp2_id.add(allele_hp2_idx)
                used_hp0_id.add(allele_hp0_idx)
                result_list.append(this_SV)
                break
        if allele_hp2_idx not in used_hp2_id:
            allele_hp2 = a_phased_single_SV[2][allele_hp2_idx]
            if allele_hp2[4] >= read_count:
                result_list.append(allele_hp2)
    for allele_hp0_idx in range(len(a_phased_single_SV[0])):
        if allele_hp0_idx not in used_hp0_id:
            allele_hp0 = a_phased_single_SV[0][allele_hp0_idx]
            if allele_hp0[4] >= read_count:
                result_list.append(allele_hp0)
    
    return result_list

def run_invm(args):
    return resolution_INV(*args)

def resolution_TRA(path, chr_1, read_count, overlap_size, max_cluster_bias, bam_path, action, sigs_index):
    if chr_1 not in sigs_index["TRA"].keys():
        return (chr_1,[])
    semi_tra_cluster = list()
    semi_tra_cluster.append([0,0,'','N'])
    candidate_single_SV = list()
    candidate_cluster = list()
    with open("%s%s.pickle"%(path, "TRA"), 'rb') as f:
        f.seek(sigs_index["TRA"][chr_1])
        seqs=pickle.load(f)
    chr_2=seqs[0][2]
    for seq in seqs:
        if seq[2]!=chr_2:
            if len(semi_tra_cluster) >= read_count:
                if semi_tra_cluster[-1][0] == semi_tra_cluster[-1][1] == 0:
                    pass
                else:
                    generate_semi_tra_cluster(semi_tra_cluster, 
                                            chr_1, 
                                            chr_2, 
                                            read_count, 
                                            overlap_size, 
                                            max_cluster_bias, 
                                            candidate_cluster,
                                            bam_path,
                                            action)
            
            logging.info("Finished %s-%s:%s."%(chr_1, chr_2, "TRA/BND"))
            semi_tra_cluster = list()
            semi_tra_cluster.append([0,0,'','N'])
            chr_2=seq[2]
        pos_1 = int(seq[1])
        pos_2 = int(seq[3])
        read_id = seq[4]
        BND_type = seq[0]
        hp = int(seq[-1])
    
        if pos_1 - semi_tra_cluster[-1][0] > max_cluster_bias or BND_type != semi_tra_cluster[-1][3]:
            if len(semi_tra_cluster) >= 2:
                if semi_tra_cluster[-1][0] == semi_tra_cluster[-1][1] == 0:
                    pass
                else:
                    generate_semi_tra_cluster(semi_tra_cluster, 
                                            chr_1, 
                                            chr_2, 
                                            read_count, 
                                            overlap_size, 
                                            max_cluster_bias, 
                                            candidate_cluster,
                                            bam_path,
                                            action)
            semi_tra_cluster = []
            semi_tra_cluster.append([pos_1, pos_2, read_id, BND_type, hp])
        else:
            if semi_tra_cluster[-1][0] == semi_tra_cluster[-1][1] == 0:
                semi_tra_cluster = []
                semi_tra_cluster.append([pos_1, pos_2, read_id, BND_type, hp])
            else:
                semi_tra_cluster.append([pos_1, pos_2, read_id, BND_type, hp])

    if len(semi_tra_cluster) >= 2:
        if semi_tra_cluster[-1][0] == semi_tra_cluster[-1][1] == 0:
            pass
        else:
            generate_semi_tra_cluster(semi_tra_cluster, 
                                    chr_1, 
                                    chr_2, 
                                    read_count, 
                                    overlap_size, 
                                    max_cluster_bias, 
                                    candidate_cluster,
                                    bam_path,
                                    action)
    if len(candidate_cluster) == 0:
        logging.info("Finished %s-%s:%s."%(chr_1, chr_2, "TRA/BND"))
        return (chr_1, [], [])

    add_tra_clip(candidate_cluster, path, chr_1, chr_2, sigs_index, max_cluster_bias)
    candidate_single_SV = call_gt_tra(candidate_cluster, path, chr_1, chr_2, sigs_index, max_cluster_bias)
    
    logging.info("Finished %s-%s:%s."%(chr_1, chr_2, "TRA/BND"))
    return (chr_1,candidate_single_SV, [])

def generate_semi_tra_cluster(semi_tra_cluster, chr_1, chr_2, read_count, overlap_size, 
    max_cluster_bias, candidate_cluster, bam_path, action):
    BND_type = semi_tra_cluster[0][3]
    semi_tra_cluster = sorted(semi_tra_cluster, key = lambda x:x[1])
    read_tag = dict()
    temp = list()
    # p1, p2, count
    last_len = semi_tra_cluster[0][1]
    temp.append([semi_tra_cluster[0][0], semi_tra_cluster[0][1], {semi_tra_cluster[0][2]}, [semi_tra_cluster[0][4]]])
    read_tag[semi_tra_cluster[0][2]] = 0
    for element in semi_tra_cluster[1:]:
        if element[1] - last_len > max_cluster_bias:
            temp[-1][0] = temp[-1][0] / len(temp[-1][3])
            temp[-1][1] = temp[-1][1] / len(temp[-1][3])
            temp.append([element[0],element[1],{element[2]},[element[4]]])
            last_len = element[1]
        else:
            if element[2] not in temp[-1][2]:
                temp[-1][0] += element[0]
                temp[-1][1] += element[1]
                temp[-1][2].add(element[2])
                temp[-1][3].append(element[4])
            last_len = element[1]

        if element[2] not in read_tag:
            read_tag[element[2]] = 0
    temp[-1][0] = temp[-1][0] / len(temp[-1][3])
    temp[-1][1] = temp[-1][1] / len(temp[-1][3])
    if len(read_tag) < 2:
        return
    temp = sorted(temp, key = lambda x:-len(set(x[2])))
    
    for item in temp[:2]:
        if len(item[2]) > 1:
            BND_pos = "%s:%s"%(chr_2, int(item[1]))
            if BND_type == 'A':
                TRA = "N[%s["%(BND_pos)
            elif BND_type == 'B':
                TRA = "N]%s]"%(BND_pos)
            elif BND_type == 'C':
                TRA = "[%s[N"%(BND_pos)
            elif BND_type == 'D':
                TRA = "]%s]N"%(BND_pos)
            else:
                return
            
            candidate_cluster.append([int(item[0]), int(item[1]), chr_1, chr_2, item[2], TRA, item[3]])


def add_tra_clip(candidate_cluster, path, chr_1, chr_2, sigs_index, max_cluster_bias):
    # candidate_cluster: [pos_1, pos_2, chr_1, chr_2, read_id_list, TRA_type, hp_list]
    max_cluster_bias = max_cluster_bias / 2
    clipfile = open("%sclip.pickle"%(path), 'rb')
    if chr_1 in sigs_index["clip"]:
        clipfile.seek(sigs_index["clip"][chr_1])
        clip1_list=pickle.load(clipfile)
    else:
        clip1_list = []
    if chr_2 in sigs_index["clip"]:
        clipfile.seek(sigs_index["clip"][chr_2])
        clip2_list=pickle.load(clipfile)
    else:
        clip2_list = []
    clipfile.close()
    sort1_list = list()
    sort2_list = list()
    idx = 0
    for i in clip1_list:
        # ('m84011_220902_175841_s1/153752542/ccs', 0, 10000, 'chr1', 1)
        if i[1] == 2:
            sort1_list.append([i[2], 1, idx, i[0], i[4]]) # [pos, 1, idx, read_name, hp]
        idx += 1
    idx = 0
    for i in clip2_list:
        if i[1] == 2:
            sort2_list.append([i[2], 1, idx, i[0], i[4]]) # [pos, 1, idx, read_name, hp]
        idx += 1
    idx = 0
    for i in candidate_cluster:
        sort1_list.append([max(0, i[0]-max_cluster_bias), 2, idx])
        sort1_list.append([i[0]+max_cluster_bias, 0, idx])
        sort2_list.append([max(0, i[1]-max_cluster_bias), 2, idx])
        sort2_list.append([i[1]+max_cluster_bias, 0, idx])
        idx += 1
    sort1_list = sorted(sort1_list, key = lambda x:(x[0], x[1]))
    sort2_list = sorted(sort2_list, key = lambda x:(x[0], x[1]))
    svs_set = set()
    overlap_dict = dict()
    for node in sort1_list:
        if node[1] == 1: # set2(clip)
            for x in svs_set:
                if candidate_cluster[x][0]+max_cluster_bias == node[0]:
                    continue
                if x not in overlap_dict:
                    overlap_dict[x] = set()
                overlap_dict[x].add(node[2])
        elif node[1] == 2: # set1(pre_pos) left
            svs_set.add(node[2])
        elif node[1] == 0: # set1(pre_pos) right
            svs_set.remove(node[2])
    for i in range(len(candidate_cluster)):
        if i in overlap_dict:
            for j in overlap_dict[i]:
                if clip1_list[j][0] not in candidate_cluster[i][4]:
                    candidate_cluster[i][4].add(clip1_list[j][0])
                    candidate_cluster[i][6].append(clip1_list[j][4])


def run_tram(args):
    return resolution_TRA(*args)

def call_gt_tra(candidate_single_SV, path, chr_1, chr_2, sigs_index, max_cluster_bias):
    # candidate_single_SV: [pos_1, pos_2, chr_1, chr_2, read_id_list, TRA_type, hp_list]
    if chr_1 not in sigs_index["reads"].keys():
        return []
    readsfile = open("%sreads.pickle"%(path), 'rb')
    readsfile.seek(sigs_index["reads"][chr_1])
    reads_list=pickle.load(readsfile)
    readsfile.close()
    svs_list = list()
    for item_idx in range(len(candidate_single_SV)):
        search_start = max(candidate_single_SV[item_idx][0] - max_cluster_bias, 0)
        search_end = candidate_single_SV[item_idx][0] + max_cluster_bias
        svs_list.append((search_start, search_end))

    iteration_dict, primary_num_dict, cover_dict, read_to_phase_dict = overlap_cover(svs_list, reads_list, 'O') # both key(sv idx), value(set(read id))
    assert len(cover_dict) == len(candidate_single_SV), "overlap length error"
    read_id_dict = dict()
    phased_id_dict = dict()
    for i in range(len(candidate_single_SV)):
        read_id_dict[i] = candidate_single_SV[i][4]
        phased_id_dict[i] = [0, 0, 0]
        for hp in candidate_single_SV[i][6]:
            phased_id_dict[i][hp] += 1
        
    assign_list = assign_gt_mosaic(iteration_dict, primary_num_dict, cover_dict, read_id_dict, phased_id_dict, read_to_phase_dict)
    # [[DV, DR, GT, GL, GQ, QUAL] ...]
    assert len(candidate_single_SV) == len(assign_list), "assign error"
    candidate_single_SV_gt = list()
    # candidate_single_SV: [pos_1, pos_2, chr_1, chr_2, read_id_list, TRA_type, hp_list]
    # candidate_single_SV_gt: [chr_1, TRA, pos_1, chr_2, pos_2, DV, DR, GT, GL, GQ, QUAL, read_list]
    for i in range(len(candidate_single_SV)):
        assert len(assign_list[i]) == 6, "assign genotype error"
        assert len(candidate_single_SV[i][4]) == assign_list[i][0], "assign DV error"
        if '1' in assign_list[i][2]:
            continue
        DV = assign_list[i][0]
        DR = assign_list[i][1]
        read_count = max(int((DV + DR) * 0.05), 4)
        if DV < read_count:
            continue
        candidate_single_SV_gt.append([candidate_single_SV[i][2], 
                                    candidate_single_SV[i][5], 
                                    candidate_single_SV[i][0], 
                                    candidate_single_SV[i][3], 
                                    candidate_single_SV[i][1], 
                                    len(candidate_single_SV[i][4]),
                                    assign_list[i][1],
                                    assign_list[i][2],
                                    assign_list[i][3],
                                    assign_list[i][4],
                                    assign_list[i][5],
                                    ','.join(candidate_single_SV[i][4])])
    return candidate_single_SV_gt

def run_tram(args):
    return resolution_TRA(*args)