import numpy as np
import logging
from cuteHap.cuteHap_genotype import cal_GL, threshold_ref_count, count_coverage, overlap_cover, assign_gt
import pickle

'''
        *********Description*********
        *	TYPE A:		N[chr:pos[	*
        *	TYPE B:		N]chr:pos]	*
        *	TYPE C:		[chr:pos[N	*
        *	TYPE D:		]chr:pos]N	*
        *****************************
'''

def resolution_TRA(path, chr_1, read_count, overlap_size, max_cluster_bias, bam_path, action, gt_round, sigs_index):
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
                                            action,
                                            gt_round)
            
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
                                            action,
                                            gt_round)
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
                                    action,
                                    gt_round)
    if len(candidate_cluster) == 0:
        logging.info("Finished %s-%s:%s."%(chr_1, chr_2, "TRA/BND"))
        return (chr_1, [], [])

    add_tra_clip(candidate_cluster, path, chr_1, chr_2, sigs_index, max_cluster_bias)
    candidate_single_SV = call_gt(candidate_cluster, path, chr_1, chr_2, sigs_index, max_cluster_bias)
    
    logging.info("Finished %s-%s:%s."%(chr_1, chr_2, "TRA/BND"))
    return (chr_1,candidate_single_SV, [])

def generate_semi_tra_cluster(semi_tra_cluster, chr_1, chr_2, read_count, overlap_size, 
    max_cluster_bias, candidate_cluster, bam_path, action, gt_round):
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


def run_tra(args):
    return resolution_TRA(*args)

def call_gt(candidate_single_SV, path, chr_1, chr_2, sigs_index, max_cluster_bias):
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
    assign_list = assign_gt(iteration_dict, primary_num_dict, cover_dict, read_id_dict, phased_id_dict, read_to_phase_dict, candidate_single_SV, [])
    # [[DV, DR, GT, GL, GQ, QUAL] ...]
    assert len(candidate_single_SV) == len(assign_list), "assign error"
    candidate_single_SV_gt = list()
    # candidate_single_SV: [pos_1, pos_2, chr_1, chr_2, read_id_list, TRA_type, hp_list]
    # candidate_single_SV_gt: [chr_1, TRA, pos_1, chr_2, pos_2, DV, DR, GT, GL, GQ, QUAL, read_list]
    for i in range(len(candidate_single_SV)):
        assert len(assign_list[i]) == 6, "assign genotype error"
        assert len(candidate_single_SV[i][4]) == assign_list[i][0], "assign DV error"
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


def call_gt_old(bam_path, pos_1, pos_2, chr_1, chr_2, read_id_list, max_cluster_bias, gt_round):
    import pysam
    bamfile = pysam.AlignmentFile(bam_path)
    querydata = set()
    search_start = max(int(pos_1) - max_cluster_bias, 0)
    search_end = min(int(pos_1) + max_cluster_bias, bamfile.get_reference_length(chr_1))

    up_bound = threshold_ref_count(len(read_id_list))

    status = count_coverage(chr_1, 
                            search_start, 
                            search_end, 
                            bamfile, 
                            querydata, 
                            up_bound, 
                            gt_round)

    if status == -1:
        DR = '.'
        GT = "./."
        GL = ".,.,."
        GQ = "."
        QUAL = "."

    elif status == 1:
        DR = 0
        for query in querydata:
            if query not in read_id_list:
                DR += 1
        GT, GL, GQ, QUAL = cal_GL(DR, len(read_id_list))

    else:
        search_start = max(int(pos_2) - max_cluster_bias, 0)
        search_end = min(int(pos_2) + max_cluster_bias, bamfile.get_reference_length(chr_2))
        status_2 = count_coverage(chr_2, 
                                    search_start, 
                                    search_end, 
                                    bamfile, 
                                    querydata, 
                                    up_bound, 
                                    gt_round)
        # status_2 judgement
        DR = 0
        for query in querydata:
            if query not in read_id_list:
                DR += 1
        GT, GL, GQ, QUAL = cal_GL(DR, len(read_id_list))

    bamfile.close()
    return len(read_id_list), DR, GT, GL, GQ, QUAL
