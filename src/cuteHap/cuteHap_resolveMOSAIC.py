import numpy as np
from cuteHap.cuteHap_genotype import cal_CIPOS, overlap_cover, assign_gt_mosaic, bench_phase
from cuteHap.dp_v12 import solve_dp
import logging
import pickle
import time

def resolution_INDEL(path, chr, svtype, read_count, threshold_gloab, max_cluster_bias,
                 minimum_support_reads, action, min_size, remain_reads_ratio, sigs_index):
    
    semi_clusters_list = list()
    pre_pos = list()
    mean_length = 0
    semi_cluster = list()
    semi_cluster.append(['',0,0,'',0])

    with open("%s%s.pickle"%(path, svtype), 'rb') as f:
        f.seek(sigs_index[svtype][chr])
        seqs=pickle.load(f)
    for seq in seqs:
        if seq[-2] != chr:
            continue
        pos = int(seq[0])
        indel_len = int(seq[1])
        read_id = seq[2]
        hp = int(seq[-1])
        if svtype == "INS":
            alt_seq = seq[3]
        else:
            alt_seq = "<DEL>"

        if pos - semi_cluster[-1][1] > max_cluster_bias:
            if len(semi_cluster) >= read_count or mean_length / len(semi_cluster) >= 3000:
                if semi_cluster[-1][1] == semi_cluster[-1][2] == 0:
                    pass
                else:
                    semi_clusters_list.append(semi_cluster)
                    pre_pos.append([max(semi_cluster[0][1] - max_cluster_bias, 0), semi_cluster[-1][1] + max_cluster_bias])
                    
            semi_cluster = []
            mean_length = 0
            semi_cluster.append([read_id, pos, indel_len, alt_seq, hp])
            mean_length += indel_len
        else:
            if semi_cluster[-1][1] == semi_cluster[-1][2] == 0:
                semi_cluster = []
                mean_length = 0
                semi_cluster.append([read_id, pos, indel_len, alt_seq, hp])
                mean_length += indel_len
            else:
                semi_cluster.append([read_id, pos, indel_len, alt_seq, hp])
                mean_length += indel_len

    if len(semi_cluster) >= read_count or mean_length / len(semi_cluster) >= 3000:
        if semi_cluster[-1][1] == semi_cluster[-1][2] == 0:
            pass
        else:
            semi_clusters_list.append(semi_cluster)
            pre_pos.append([max(semi_cluster[0][1] - max_cluster_bias, 0), semi_cluster[-1][1] + max_cluster_bias])
    
    hp_dist_dict, overlap_reads_dict = bench_phase(path, chr, sigs_index, pre_pos)
    pos2clip_dict = add_clip(path, chr, sigs_index, pre_pos)

    candidate_single_SV = list()
    phased_single_SV = list()
    search_bad_list = list()
    search_flag_list = list()
    for i in range(len(pre_pos)):
        hp_count = hp_dist_dict[i]
        overlap_reads = overlap_reads_dict[i]
        if hp_count[0]+abs(hp_count[1]-hp_count[2]) < (hp_count[1]+hp_count[2])/2 or abs(hp_count[1]-hp_count[2]) < hp_count[0] < max(hp_count[1],hp_count[2]):
        # if 1 == 0:
            a_phased_single_SV = generate_balance(semi_clusters_list[i], 
                                chr, 
                                svtype, 
                                read_count, 
                                threshold_gloab, 
                                minimum_support_reads, 
                                action,
                                min_size,
                                max_cluster_bias,
                                path,
                                overlap_reads,
                                sigs_index,
                                pos2clip_dict[i] if i in pos2clip_dict else [])
            phased_single_SV.extend(a_phased_single_SV)
        else:
            # bad cases   
            if len(semi_clusters_list[i]) < 100:
                generate_bad_ultra(semi_clusters_list[i], 
                                    chr, 
                                    svtype, 
                                    read_count, 
                                    threshold_gloab, 
                                    minimum_support_reads, 
                                    candidate_single_SV,
                                    action,
                                    remain_reads_ratio,
                                    pos2clip_dict[i] if i in pos2clip_dict else [],
                                    search_bad_list,
                                    search_flag_list,
                                    max_cluster_bias)
            else:
                generate_bad_fuzzy(semi_clusters_list[i], 
                                    chr, 
                                    svtype, 
                                    read_count, 
                                    threshold_gloab, 
                                    minimum_support_reads, 
                                    candidate_single_SV,
                                    action,
                                    remain_reads_ratio,
                                    pos2clip_dict[i] if i in pos2clip_dict else [],
                                    search_bad_list,
                                    search_flag_list,
                                    max_cluster_bias)
    candidate_single_SV_gt = call_gt(path, chr, candidate_single_SV, search_bad_list, search_flag_list, max_cluster_bias, svtype, sigs_index, min_size)
    logging.info("Finished %s:%s."%(chr, svtype))
    return (chr,candidate_single_SV_gt, phased_single_SV)
    
def add_clip(path, chr, sigs_index, pre_pos):
    clipfile = open("%sclip.pickle"%(path), 'rb')
    if chr in sigs_index["clip"]:
        clipfile.seek(sigs_index["clip"][chr])
        clip_list=pickle.load(clipfile)
    else:
        clip_list = []
    clipfile.close()
    sort_list = list()
    idx = 0
    for i in clip_list:
        # ('m84011_220902_175841_s1/153752542/ccs', 0, 10000, 'chr1', 1)
        if i[1] == 0 or i[1] == 1:
            sort_list.append([i[2], 1, idx, i[0], i[4]]) # [pos, 1, idx, read_name, hp]
        idx += 1
    idx = 0
    for i in pre_pos:
        sort_list.append([i[0], 2, idx])
        sort_list.append([i[1], 0, idx])
        idx += 1
    sort_list = sorted(sort_list, key = lambda x:(x[0], x[1]))
    svs_set = set()
    overlap_dict = dict()
    for node in sort_list:
        if node[1] == 1: # set2(clip)
            for x in svs_set:
                if pre_pos[x][1] == node[0]:
                    continue
                if x not in overlap_dict:
                    overlap_dict[x] = set()
                overlap_dict[x].add(node[2])
        elif node[1] == 2: # set1(pre_pos) left
            svs_set.add(node[2])
        elif node[1] == 0: # set1(pre_pos) right
            svs_set.remove(node[2])
    pos2clip_dict = dict()
    for idx in overlap_dict:
        pos2clip_dict[idx] = list()
        for x in overlap_dict[idx]:
            pos2clip_dict[idx].append(clip_list[x])
    return pos2clip_dict

def generate_bad_ultra(semi_cluster, chr, svtype, read_count, 
    threshold_gloab, minimum_support_reads, candidate_single_SV, 
    action, remain_reads_ratio, clips, search_bad_list, search_flag_list, max_cluster_bias):

    # [read_id, pos, indel_len, alt_seq, hp]
    real_ans = solve_dp(semi_cluster)
    
    for allele in real_ans:
        breakpointStart = allele[0]
        signalLen = allele[1]
        
    real_ans.sort(key=lambda x: len(x[2]))

    a_candidate_single_SV = list()
    used_read_id = set()
    used_clip_idx = [set(), set()]
    for allele in real_ans:
        breakpointStart = allele[0]
        signalLen = allele[1]
        modified_max_cluster_bias = min(max_cluster_bias, signalLen+50)
        support_read_set = set()
        phase_id_list = list()
        for _ in allele[2]:
            support_read_set.add(_[0])
            phase_id_list.append(_[-1])
        #'''
        used_read_id = used_read_id.union(support_read_set)
        if svtype == 'INS' and signalLen > 2000:
            for _clip_id in range(len(clips)):
                clip = clips[_clip_id]
                if _clip_id in used_clip_idx[clip[1]]:
                    continue
                if clip[0] in used_read_id:
                    continue
                # [read_name, left/right clip, pos, chrom, hp]
                if breakpointStart-modified_max_cluster_bias < clip[2] < breakpointStart+modified_max_cluster_bias and clip[0] not in support_read_set:
                    support_read_set.add(clip[0])
                    used_clip_idx[clip[1]].add(_clip_id)
        if svtype == 'DEL' and signalLen > 2000:
            for _clip_id in range(len(clips)):
                clip = clips[_clip_id]
                if _clip_id in used_clip_idx[clip[1]]:
                    continue
                if clip[0] in used_read_id:
                    continue
                # [read_name, left/right clip, pos, chrom, hp]
                if clip[1] == 1 and breakpointStart-modified_max_cluster_bias < clip[2] < breakpointStart+modified_max_cluster_bias and clip[0] not in support_read_set:
                    support_read_set.add(clip[0])
                    used_clip_idx[clip[1]].add(_clip_id)
        
        if len(support_read_set) >= minimum_support_reads:
            if svtype == 'DEL':
                ideal_seq = '<DEL>'
            else:
                ideal_seq = '<INS>'
            candidate_single_SV.append([chr, 
                                    svtype, 
                                    int(breakpointStart), 
                                    int(signalLen), 
                                    len(allele[2]), 
                                    './.',
                                    './.',
                                    int(breakpointStart), 
                                    support_read_set,
                                    phase_id_list,
                                    ideal_seq])
            search_bad_list.append([int(breakpointStart)-max_cluster_bias, int(breakpointStart)+max_cluster_bias])
            search_flag_list.append(False)

def generate_balance(semi_del_cluster, chrom, svtype, read_count, 
    threshold_gloab, minimum_support_reads, 
    action, min_size, max_cluster_bias, temporary_dir, overlap_reads, sigs_index, clips):

    '''
    generate deletion
    *************************************************************
    threshold_gloab 	threshold_local 	minimum_support_reads
    -------------------------------------------------------------
        0.3					0.7 					5		CLR
        0.4					0.5 				  <=5		CCS
    *************************************************************
    '''
    # [read_id, pos, indel_len, alt_seq, hp]
    # Remove duplicates
    # read_tag [hp] [read_id] = [pos, indel_len, read_id, alt_seq, hp]
    read_tag = [dict(), dict(), dict()]
    for element in semi_del_cluster:
        if element[0] not in read_tag[element[-1]]:
            read_tag[element[-1]][element[0]] = [[element[1], element[2], element[0], element[3], element[4]]]
        else:
            read_tag[element[-1]][element[0]].append([element[1], element[2], element[0], element[3], element[4]])
    read_list = [[], [], []]
    for element in semi_del_cluster:
        read_list[element[-1]].append([element[1], element[2], element[0], element[3], element[4]])
    
    clustered_list = {0: list(), 1: list(), 2: list()}
    remained_alleles = [dict(), dict(), dict()]
    pre_gt = {0: ".|.", 1: "1|0", 2: "0|1"}
    used_read_id = set()
    used_clip_idx = [set(), set()]
    for hp_idx in [1, 2, 0]:
        if len(read_list[hp_idx]) == 0:
            continue
        read_tag2SortedList = sorted(read_list[hp_idx], key = lambda x:x[1])
        global_len = [i[1] for i in read_tag2SortedList]
        DISCRETE_THRESHOLD_LEN_CLUSTER_DEL_TEMP = threshold_gloab * np.mean(global_len)

        last_len = read_tag2SortedList[0][1]
        cur_bias = last_len * threshold_gloab

        allele_collect = list()
        '''
        *********************************************************************************
            #0   			#1			#2  		#3          #4           #5
        ---------------------------------------------------------------------------------
            [breakpoint]	[len]		[#support] 	[read-id]   [alt_seq]    [hp]
        *********************************************************************************
        '''
        allele_collect.append([[read_tag2SortedList[0][0]],[read_tag2SortedList[0][1]],[],
            [read_tag2SortedList[0][2]], [read_tag2SortedList[0][3]], [read_tag2SortedList[0][4]]])

        for i in read_tag2SortedList[1:]:
            if i[1] - last_len > cur_bias:
                allele_collect[-1][2].append(len(allele_collect[-1][0]))
                allele_collect.append([[],[],[],list(),[],[]])

            allele_collect[-1][0].append(i[0])
            allele_collect[-1][1].append(i[1])
            allele_collect[-1][3].append(i[2]) # add
            allele_collect[-1][4].append(i[3])
            allele_collect[-1][5].append(i[4])
            last_len = (last_len * (len(allele_collect[-1][0]) - 1) + i[1]) / len(allele_collect[-1][0])
            cur_bias = last_len * threshold_gloab
        allele_collect[-1][2].append(len(allele_collect[-1][0]))
        allele_sort = sorted(allele_collect, key = lambda x:x[2], reverse=True)
        for allele_hp in allele_sort:
            if (hp_idx > 0 and allele_hp[2][0] < read_count) or (hp_idx == 0 and allele_hp[2][0] < read_count*2):
                continue
            allele_hp[3] = set(allele_hp[3])
            breakpointStart = np.mean(allele_hp[0])
            min_len_1 = allele_hp[1][0]
            max_len_1 = allele_hp[1][-1]
            signalLen = int(np.mean(allele_hp[1]))
            genotype = pre_gt[hp_idx]
            support_read_count = allele_hp[2][0]
            support_read_set = allele_hp[3]
            CIPOS = cal_CIPOS(np.std(allele_hp[0]), len(allele_hp[0]))
            CILEN = cal_CIPOS(np.std(allele_hp[1]), len(allele_hp[1]))
            this_SV = [chrom, svtype, [np.min(allele_hp[0]), int(breakpointStart), np.max(allele_hp[0])], [min_len_1, int(signalLen), max_len_1], support_read_count, str(CIPOS),str(CILEN), '.', genotype, '.', '.', '.', support_read_set, '<%s>'%(svtype)]
            clustered_list[hp_idx].append(this_SV)
    
    combine_list = combine_alleles(clustered_list, read_tag, read_count, read_tag)
    combine_list.sort(key=lambda x: sum(x[4]), reverse=True)
    
    for allele_hp in combine_list:
        breakpointStart = allele_hp[2][1]
        signalLen = allele_hp[3]
        modified_max_cluster_bias = min(max_cluster_bias, signalLen+50)
        used_read_id = used_read_id.union(allele_hp[12])
        allele_hp_info = {0}
        if allele_hp[8][0] == '1': allele_hp_info.add(1)
        if allele_hp[8][2] == '1': allele_hp_info.add(2)
        if svtype == 'INS' and signalLen > 2000:
            for _clip_id in range(len(clips)):
                clip = clips[_clip_id]
                if _clip_id in used_clip_idx[clip[1]]:
                    continue
                if clip[0] in used_read_id:
                    continue
                # [read_name, left/right clip, pos, chrom, hp]
                if min(allele_hp[2][0], breakpointStart-modified_max_cluster_bias) < clip[2] < max(allele_hp[2][2], breakpointStart+modified_max_cluster_bias) and clip[4] in allele_hp_info:
                    if clip[0] not in allele_hp[12]:
                        allele_hp[4][clip[4]] += 1
                        allele_hp[12].add(clip[0])
                        used_clip_idx[clip[1]].add(_clip_id)
        if svtype == 'DEL' and signalLen > 2000:
            for _clip_id in range(len(clips)):
                clip = clips[_clip_id]
                if _clip_id in used_clip_idx[clip[1]]:
                    continue
                if clip[0] in used_read_id:
                    continue
                # [read_name, left/right clip, pos, chrom, hp]
                if clip[1] == 1 and allele_hp[2][0] < clip[2] < allele_hp[2][2] and clip[4] in allele_hp_info:
                    if clip[0] not in allele_hp[12]:
                        allele_hp[4][clip[4]] += 1
                        allele_hp[12].add(clip[0])
                        used_clip_idx[clip[1]].add(_clip_id)
        if '1' in allele_hp[8]:
            if allele_hp[8][0] == '1':
                revise_single_allele(allele_hp, max_cluster_bias, overlap_reads, 1)
            if allele_hp[8][2] == '1':
                revise_single_allele(allele_hp, max_cluster_bias, overlap_reads, 2)
        else:
            revise_single_allele_up(allele_hp, max_cluster_bias, overlap_reads)
        allele_hp[4] = sum(allele_hp[4])
        allele_hp[2] = allele_hp[2][1]
    result_list = list()
    for _ in combine_list:
        if '1' not in _[8] and abs(_[3]) >= min_size:
            result_list.append(_)    
    return result_list

def generate_bad_fuzzy(semi_cluster, chr, svtype, read_count, 
    threshold_gloab, minimum_support_reads, candidate_single_SV, 
    action, remain_reads_ratio, clips, search_bad_list, search_flag_list, max_cluster_bias):

    '''
    generate deletion
    *************************************************************
    threshold_gloab 	threshold_local 	minimum_support_reads
    -------------------------------------------------------------
        0.3					0.7 					5		CLR
        0.4					0.5 				  <=5		CCS
    *************************************************************
    '''
    
    if len(semi_cluster) < read_count:
        return

    # [read_id, pos, indel_len, alt_seq, hp]
    read_tag2SortedList = sorted(semi_cluster, key = lambda x:x[2])

    last_len = read_tag2SortedList[0][2]
    cur_bias = last_len * threshold_gloab

    allele_collect = list()
    '''
    *********************************************************************************
        #0   				#1			    #2  		#3          #4          #5
    ---------------------------------------------------------------------------------
        [breakpoint]	[len]		[#support] 	[read-id]   [alt_seq]    [hp]
    *********************************************************************************
    '''
    allele_collect.append([[read_tag2SortedList[0][1]],[read_tag2SortedList[0][2]],[],
        {read_tag2SortedList[0][0]}, [read_tag2SortedList[0][3]], [read_tag2SortedList[0][4]]])

    for i in read_tag2SortedList[1:]:
        if i[2] - last_len > cur_bias:
            allele_collect[-1][2].append(len(allele_collect[-1][0]))
            allele_collect.append([[],[],[],set(),[],[]])

        allele_collect[-1][0].append(i[1])
        allele_collect[-1][1].append(i[2])
        allele_collect[-1][3].add(i[0])
        allele_collect[-1][4].append(i[3])
        allele_collect[-1][5].append(i[4])
        last_len = (last_len * (len(allele_collect[-1][0]) - 1) + i[2]) / len(allele_collect[-1][0])
        cur_bias = last_len * threshold_gloab
    allele_collect[-1][2].append(len(allele_collect[-1][0]))
    allele_sort = sorted(allele_collect, key = lambda x:x[2], reverse=True)
    a_candidate_single_SV = list()
    used_read_id = set()
    used_clip_idx = [set(), set()]
    for allele in allele_sort:
        breakpointStart = np.mean(allele[0])
        used_read_id = used_read_id.union(allele[3])
        if svtype == 'INS' and allele[1][0] > 2000: 
            for _clip_id in range(len(clips)):
                clip = clips[_clip_id]
                if _clip_id in used_clip_idx[clip[1]]:
                    continue
                if clip[0] in used_read_id:
                    continue
                # [read_name, left/right clip, pos, chrom, hp]
                if min(allele[0][0], breakpointStart-max_cluster_bias) < clip[2] < max(allele[0][-1], breakpointStart+max_cluster_bias) and clip[0] not in allele[3]:
                    allele[2][0] += 1
                    allele[3].add(clip[0])
                    used_clip_idx[clip[1]].add(_clip_id)
        if svtype == 'DEL' and allele[1][0] > 2000:
            for _clip_id in range(len(clips)):
                clip = clips[_clip_id]
                if _clip_id in used_clip_idx[clip[1]]:
                    continue
                if clip[0] in used_read_id:
                    continue
                # [read_name, left/right clip, pos, chrom, hp]
                if clip[1] == 1 and allele[0][0]< clip[2] < allele[0][-1] and clip[0] not in allele[3]:
                    allele[2][0] += 1
                    allele[3].add(clip[0])
                    used_clip_idx[clip[1]].add(_clip_id)

        if allele[2][0] >= minimum_support_reads:
            breakpointStart = np.mean(allele[0])
            search_start = np.min(allele[0])
            search_end = np.max(allele[0])
            CIPOS = cal_CIPOS(np.std(allele[0]), len(allele[0]))
            signalLen = np.mean(allele[1])
            signalLen_STD = np.std(allele[1])
            CILEN = cal_CIPOS(np.std(allele[1]), len(allele[1]))
            if svtype == 'DEL':
                ideal_seq = '<DEL>'
            else:
                ideal_seq = '<INS>'
                for pos,i in zip(allele[0],allele[4]):
                    if len(i) >= int(signalLen):
                        # breakpointStart = pos
                        ideal_seq = i[0:int(signalLen)]
                        break

            candidate_single_SV.append([chr, 
                                        svtype, 
                                        int(breakpointStart), 
                                        int(signalLen), 
                                        allele[2][0], 
                                        str(CIPOS),
                                        str(CILEN),
                                        int(breakpointStart), 
                                        allele[3],
                                        allele[5],
                                        ideal_seq])
            search_bad_list.append([min(search_start, search_end-max_cluster_bias), max(search_end, search_start+max_cluster_bias)])
            search_flag_list.append(True if np.std(allele[0]) < 1 and signalLen_STD < 1 else False)

def run_delm(args):
    return resolution_INDEL(*args)

def run_insm(args):
    return resolution_INDEL(*args)

def call_gt(temporary_dir, chr, candidate_single_SV, search_bad_list, search_flag_list, max_cluster_bias, svtype, sigs_index, min_size):
    # reads_list = list() # [(start, end, 0, 'name'), ...]
    # candidate_single_SV： [chr, svtype, pos, len, DV, CIPOS, CILEN, search_pos, read_id_set, hp_set]
    if chr not in sigs_index["reads"].keys():
        return []
    readsfile = open("%sreads.pickle"%(temporary_dir), 'rb')
    readsfile.seek(sigs_index["reads"][chr])
    reads_list=pickle.load(readsfile)
    readsfile.close()
    assert len(candidate_single_SV) == len(search_bad_list), "search list length error"
    svs_list = list()
    for item_idx in range(len(candidate_single_SV)):
        search_start = max(candidate_single_SV[item_idx][7] - max_cluster_bias, search_bad_list[item_idx][0])
        search_end = min(candidate_single_SV[item_idx][7] + max_cluster_bias, search_bad_list[item_idx][1])
        if search_start == search_end:
            search_end += 1
        svs_list.append((search_start, search_end))
    
    iteration_dict, primary_num_dict, cover_dict, read_to_phase_dict = overlap_cover(svs_list, reads_list, 'O') # both key(sv idx), value(set(read id))
    assert len(cover_dict) == len(candidate_single_SV), "overlap length error"
    read_id_dict = dict()
    phased_id_dict = dict()
    for i in range(len(candidate_single_SV)):
        read_id_dict[i] = candidate_single_SV[i][8]
        phased_id_dict[i] = [0, 0, 0]
        for hp in candidate_single_SV[i][9]:
            phased_id_dict[i][hp] += 1
    assert len(search_flag_list) == len(candidate_single_SV), "search_flag_list length error"
    assign_list = assign_gt_mosaic(iteration_dict, primary_num_dict, cover_dict, read_id_dict, phased_id_dict, read_to_phase_dict)
    # [[DV, DR, GT, GL, GQ, QUAL] ...]
    assert len(candidate_single_SV) == len(assign_list), "assign error"
    candidate_single_SV_gt = list()
    # candidate_single_SV_gt: [chr, svtype, pos, len, DV, CIPOS, CILEN, DR, GT, GL, GQ, QUAL, read_id_set]
    # assign_list: [DV, DR, GT, GL, GQ, QUAL]
    for i in range(len(candidate_single_SV)):
        assert len(assign_list[i]) == 6, "assign genotype error"
        if '1' in assign_list[i][2] or abs(int(candidate_single_SV[i][3])) < min_size:
            continue
        DV = assign_list[i][0]
        DR = assign_list[i][1]
        if DV < 2 or DV / (DV + DR) < 0.05 or DV / (DV + DR) > 0.4:
            continue
        candidate_single_SV_gt.append([candidate_single_SV[i][0], 
                                    candidate_single_SV[i][1], 
                                    candidate_single_SV[i][2], 
                                    candidate_single_SV[i][3], 
                                    candidate_single_SV[i][4], 
                                    candidate_single_SV[i][5],
                                    candidate_single_SV[i][6],
                                    assign_list[i][1],
                                    assign_list[i][2],
                                    assign_list[i][3],
                                    assign_list[i][4],
                                    assign_list[i][5],
                                    ','.join(candidate_single_SV[i][8])])
        if svtype == 'INS':
            candidate_single_SV_gt[-1].append('<INS>')
    return candidate_single_SV_gt

def revise_single_allele(this_SV, max_cluster_bias, overlap_reads, hp_idx):
    def cal_bayes(kdv, ndr):
        if kdv / (kdv+ndr) >= 0.9 or ndr < 2:
            return True
        return False
        err = 0.1
        prior = 0.7022
        ht1 = np.float64(pow((1-err), kdv)*pow(err, ndr)*prior)
        ht2 = np.float64(pow((1-err), ndr)*pow(err, kdv)*(1-prior))
        return True if ht1 > ht2 else False
    # reads_list = list() # [(start, end, 0, 'name'), ...]
    # phased_single_SV: [chr, svtype, pos, len, DV, CIPOS, CILEN, '.', genotype, '.', '.', '.', support_read_set, '<INS>' / None]
    
    search_start = max(this_SV[2][1] - max_cluster_bias, this_SV[2][0])
    search_end = min(this_SV[2][1] + max_cluster_bias, this_SV[2][2])
    read_count = set()
    for read in overlap_reads:
        if read[2] == 1:
            if read[0] <= search_start and read[1] >= search_end and read[5] == hp_idx:
                # cover
                read_count.add(read[3])

    support_read_count = this_SV[4][hp_idx] # sum(this_SV[4])
    ref_read_count = 0
    for query in read_count:
        if query not in this_SV[12]: # support_read_set
            ref_read_count += 1
    flag_hap = cal_bayes(support_read_count, ref_read_count)
    if hp_idx == 1:
        hp_idx = 0 # s = s[:i] + '1' + s[i+1:]
    if not flag_hap and this_SV[8][hp_idx] == '1':
        # this_SV[8][hp_idx] = '.'
        this_SV[8] = this_SV[8][:hp_idx] + '.' + this_SV[8][hp_idx+1:]
        
def revise_single_allele_up(this_SV, max_cluster_bias, overlap_reads):
    def cal_bayes(kdv, ndr):
        if kdv / (kdv+ndr) >= 0.9 or ndr < 2:
            pass
        else:
            return 0
        err = 0.1
        prior = 0.2
        ori_GL00 = np.float64(pow((1-err), ndr)*pow(err, kdv)*(1-prior)/2)
        ori_GL11 = np.float64(pow(err, ndr)*pow((1-err), kdv)*(1-prior)/2)
        ori_GL01 = np.float64(pow(0.5, ndr+kdv)*prior)
        if ori_GL11 > max(ori_GL00, ori_GL01):
            return 2
        elif ori_GL01 > max(ori_GL00, ori_GL11):
            return 1
        else:
            return 0
    # reads_list = list() # [(start, end, 0, 'name'), ...]
    # phased_single_SV: [chr, svtype, pos, len, DV, CIPOS, CILEN, '.', genotype, '.', '.', '.', support_read_set, '<INS>' / None]
        
    search_start = max(this_SV[2][1] - max_cluster_bias, this_SV[2][0])
    search_end = min(this_SV[2][1] + max_cluster_bias, this_SV[2][2])
    read_count = set()
    hp_count = [0, 0, 0]
    for read in overlap_reads:
        if read[2] == 1:
            if read[0] <= search_start and read[1] >= search_end: # cover
                read_count.add(read[3])
                hp_count[read[5]] += 1
    support_read_count = sum(this_SV[4])
    ref_read_count = 0
    for query in read_count:
        if query not in this_SV[12]: # support_read_set
            ref_read_count += 1
    flag_hap = cal_bayes(support_read_count, ref_read_count)    
    if flag_hap == 2:
        this_SV[8] = '1|1'
    elif flag_hap == 1:
        if this_SV[4][1] > this_SV[4][2]:
            this_SV[8] = '1|0'
        else:
            this_SV[8] = '0|1'

def combine_alleles(a_phased_single_SV, read_tag, read_count, remained_alleles):
    # [chr, svtype, [min_pos, int(breakpointStart), max_pos], [min_len_1, int(signalLen), max_len_1], support_read_count, str(CIPOS),str(CILEN), '.', genotype, '.', '.', '.', support_read_set, '<%s>'%(svtype)]

    result_list = list()
    used_hp2_id = set()
    used_hp0_id = set()
    used_hp0_read = set()
    for allele_hp1 in a_phased_single_SV[1]:
        used_hp1_id = False
        for allele_hp2_idx in range(len(a_phased_single_SV[2])):
            if allele_hp2_idx in used_hp2_id:
                continue
            allele_hp2 = a_phased_single_SV[2][allele_hp2_idx]
            min_len_1 = allele_hp1[3][0]
            mean_len_1 = allele_hp1[3][1]
            max_len_1 = allele_hp1[3][2]
            min_len_2 = allele_hp2[3][0]
            mean_len_2 = allele_hp2[3][1]
            max_len_2 = allele_hp2[3][2]
            if min(mean_len_1, mean_len_2) / max(mean_len_1, mean_len_2) > 0.9 and max(min_len_1, min_len_2) <= min(max_len_1, max_len_2): # a hom SV and min_len1 overlap max_len2
                breakpointStart = (allele_hp1[2][1] + allele_hp2[2][1]) / 2
                signalLen = (mean_len_1 + mean_len_2) / 2
                genotype = allele_hp1[8][0] + '|' + allele_hp2[8][2]
                support_read_count = [0, allele_hp1[4], allele_hp2[4]]  # contains sigs, no clips
                support_read_set = allele_hp1[12].union(allele_hp2[12])
                used_hp1_id = True
                used_hp2_id.add(allele_hp2_idx)
                for _read_id in remained_alleles[0]: # if abs(item[0] - breakpointStart) < 1000 and min(item[1], signalLen) / max(item[1], signalLen) > 0.7: [pos, indel_len, read_id, alt_seq, hp]
                    if _read_id in support_read_set:
                        continue
                    for _read_sig in remained_alleles[0][_read_id]: # _read_sig: [853473, 99, 'S1S1/175584/ccs', '<DEL>', 0]
                        if min(min_len_1, min_len_2) <= _read_sig[1] <= max(max_len_1, max_len_2) and min(signalLen, _read_sig[1]) / max(signalLen, _read_sig[1]) > 0.7: # need consider pos?
                            support_read_set.add(_read_id)
                            support_read_count[0] += 1
                            used_hp0_read.add(_read_id)
                result_list.append([allele_hp1[0], allele_hp1[1], [min(allele_hp1[2][0], allele_hp2[2][0]), int(breakpointStart), max(allele_hp1[2][2], allele_hp2[2][2])], int(signalLen), support_read_count, allele_hp1[5], allele_hp1[6], '.', genotype, '.', '.', '.', support_read_set, allele_hp1[13]])
                break
        # has not been used for all allele_hp2
        if not used_hp1_id:
            for allele_hp0_idx in range(len(a_phased_single_SV[0])):
                if allele_hp0_idx in used_hp0_id:
                    continue
                allele_hp0 = a_phased_single_SV[0][allele_hp0_idx]
                min_len_1 = allele_hp1[3][0]
                mean_len_1 = allele_hp1[3][1]
                max_len_1 = allele_hp1[3][2]
                min_len_2 = allele_hp0[3][0]
                mean_len_2 = allele_hp0[3][1]
                max_len_2 = allele_hp0[3][2]
                if min(mean_len_1, mean_len_2) / max(mean_len_1, mean_len_2) > 0.9 and max(min_len_1, min_len_2) <= min(max_len_1, max_len_2): # a hom SV
                    breakpointStart = (allele_hp1[2][1] + allele_hp0[2][1]) / 2
                    signalLen = (mean_len_1 + mean_len_2) / 2
                    genotype = allele_hp1[8]  # allele_hp1[8][0] + '|' + allele_hp0[8][2]
                    support_read_count = [allele_hp0[4], allele_hp1[4], 0]
                    support_read_set = allele_hp1[12].union(allele_hp0[12])
                    used_hp1_id = True
                    used_hp0_id.add(allele_hp0_idx)
                    for _read_id in remained_alleles[2]: # if abs(item[0] - breakpointStart) < 1000 and min(item[1], signalLen) / max(item[1], signalLen) > 0.7:
                        if _read_id in support_read_set:
                            continue
                        for _read_sig in remained_alleles[2][_read_id]: # _read_sig: [853473, 99, 'S1S1/175584/ccs', '<DEL>', 0]
                            if min(min_len_1, min_len_2) <= _read_sig[1] <= max(max_len_1, max_len_2) and min(signalLen, _read_sig[1]) / max(signalLen, _read_sig[1]) > 0.7:
                                support_read_set.add(_read_id)
                                support_read_count[2] += 1
                    result_list.append([allele_hp1[0], allele_hp1[1], [min(allele_hp1[2][0], allele_hp0[2][0]), int(breakpointStart), max(allele_hp1[2][2], allele_hp0[2][2])], int(signalLen), support_read_count, allele_hp1[5], allele_hp1[6], '.', genotype, '.', '.', '.', support_read_set, allele_hp1[13]])
                    break
        if not used_hp1_id: # has not been used for all allele_hp2 or allele_hp0
            if allele_hp1[4] >= int(read_count / 2):
                this_SV = allele_hp1
                this_SV[4] = [0, allele_hp1[4], 0]
                for _read_id in remained_alleles[0]:
                    if _read_id in this_SV[12]:
                        continue
                    for _read_sig in remained_alleles[0][_read_id]:
                        if allele_hp1[3][0] <= _read_sig[1] <= allele_hp1[3][2] and min(allele_hp1[3][1], _read_sig[1]) / max(allele_hp1[3][1], _read_sig[1]) > 0.7:
                            this_SV[12].add(_read_id)
                            this_SV[4][0] += 1
                            used_hp0_read.add(_read_id)
                for _read_id in remained_alleles[2]:
                    if _read_id in this_SV[12]:
                        continue
                    for _read_sig in remained_alleles[2][_read_id]:
                        if allele_hp1[3][0] <= _read_sig[1] <= allele_hp1[3][2] and min(allele_hp1[3][1], _read_sig[1]) / max(allele_hp1[3][1], _read_sig[1]) > 0.7:
                            this_SV[12].add(_read_id)
                            this_SV[4][2] += 1
                this_SV[3] = allele_hp1[3][1]
                result_list.append(this_SV)
    for allele_hp2_idx in range(len(a_phased_single_SV[2])):
        if allele_hp2_idx not in used_hp2_id:
            allele_hp2 = a_phased_single_SV[2][allele_hp2_idx]
            # HP2+HP0
            for allele_hp0_idx in range(len(a_phased_single_SV[0])):
                if allele_hp0_idx in used_hp0_id:
                    continue
                allele_hp0 = a_phased_single_SV[0][allele_hp0_idx]
                min_len_1 = allele_hp2[3][0]
                mean_len_1 = allele_hp2[3][1]
                max_len_1 = allele_hp2[3][2]
                min_len_2 = allele_hp0[3][0]
                mean_len_2 = allele_hp0[3][1]
                max_len_2 = allele_hp0[3][2]
                if min(mean_len_1, mean_len_2) / max(mean_len_1, mean_len_2) > 0.9 and max(min_len_1, min_len_2) <= min(max_len_1, max_len_2): # a hom SV
                    breakpointStart = (allele_hp2[2][1] + allele_hp0[2][1]) / 2
                    signalLen = (mean_len_1 + mean_len_2) / 2
                    genotype = allele_hp2[8]  # allele_hp0[8][0] + '|' + allele_hp2[8][2]
                    support_read_count = [allele_hp0[4], 0, allele_hp2[4]]
                    support_read_set = allele_hp2[12].union(allele_hp0[12])
                    used_hp2_id.add(allele_hp2_idx)
                    used_hp0_id.add(allele_hp0_idx)
                    for _read_id in remained_alleles[1]: # if abs(item[0] - breakpointStart) < 1000 and min(item[1], signalLen) / max(item[1], signalLen) > 0.7:
                        if _read_id in support_read_set:
                            continue
                        for _read_sig in remained_alleles[1][_read_id]:
                            if min(min_len_1, min_len_2) <= _read_sig[1] <= max(max_len_1, max_len_2) and min(signalLen, _read_sig[1]) / max(signalLen, _read_sig[1]) > 0.7:
                                support_read_set.add(_read_id)
                                support_read_count[1] += 1
                    result_list.append([allele_hp2[0], allele_hp2[1], [min(allele_hp0[2][0], allele_hp2[2][0]), int(breakpointStart), max(allele_hp0[2][2], allele_hp2[2][2])], int(signalLen), support_read_count, allele_hp2[5], allele_hp2[6], '.', genotype, '.', '.', '.', support_read_set, allele_hp2[13]])
                    break
        if allele_hp2_idx not in used_hp2_id:
            allele_hp2 = a_phased_single_SV[2][allele_hp2_idx]
            if allele_hp2[4] >= int(read_count / 2):
                this_SV = allele_hp2
                this_SV[4] = [0, 0, allele_hp2[4]]
                for _read_id in remained_alleles[0]:
                    if _read_id in this_SV[12]:
                        continue
                    for _read_sig in remained_alleles[0][_read_id]:
                        if allele_hp2[3][0] <= _read_sig[1] <= allele_hp2[3][2] and min(allele_hp2[3][1], _read_sig[1]) / max(allele_hp2[3][1], _read_sig[1]) > 0.7:
                            this_SV[12].add(_read_id)
                            this_SV[4][0] += 1
                            used_hp0_read.add(_read_id)
                for _read_id in remained_alleles[1]:
                    if _read_id in this_SV[12]:
                        continue
                    for _read_sig in remained_alleles[1][_read_id]:
                        if allele_hp2[3][0] <= _read_sig[1] <= allele_hp2[3][2] and min(allele_hp2[3][1], _read_sig[1]) / max(allele_hp2[3][1], _read_sig[1]) > 0.7:
                            this_SV[12].add(_read_id)
                            this_SV[4][1] += 1
                this_SV[3] = allele_hp2[3][1]
                result_list.append(this_SV)
    for allele_hp0_idx in range(len(a_phased_single_SV[0])):
        if allele_hp0_idx not in used_hp0_id:
            allele_hp0 = a_phased_single_SV[0][allele_hp0_idx]
            if allele_hp0[4] < read_count:
                continue
            support_read_set = set()
            for _read_id in allele_hp0[12]:
                if _read_id not in used_hp0_read:
                    support_read_set.add(_read_id)
            if len(support_read_set) >= read_count:
                this_SV = allele_hp0
                this_SV[3] = allele_hp0[3][1]
                this_SV[4] = [len(support_read_set), 0, 0]
                this_SV[12] = support_read_set
                result_list.append(this_SV)
    return result_list
