
from cuteHap.cuteHap_Description import Generation_VCF_header
from math import log10
import numpy as np
import pysam
import pickle
from cuteHap.genotype_improve import improve_overlap_cover

err = 0.1
prior = float(1/3)
Genotype = ["0|0", "0/1", "1|1"]

def log10sumexp(log10_probs):
    # Normalization of Genotype likelihoods
    m = max(log10_probs)
    return m + log10(sum(pow(10.0, x-m) for x in log10_probs))

def normalize_log10_probs(log10_probs):
    # Adjust the Genotype likelihoods
    log10_probs = np.array(log10_probs)
    lse = log10sumexp(log10_probs)
    return np.minimum(log10_probs - lse, 0.0)

def rescale_read_counts(c0, c1, max_allowed_reads=100):
    """Ensures that n_total <= max_allowed_reads, rescaling if necessary."""
    Total = c0 + c1
    if Total > max_allowed_reads:
        c0 = int(max_allowed_reads * float(c0/Total))
        c1 = max_allowed_reads - c0
    return c0, c1

def cal_GL(c0, c1):
    if c0==3 and c1==1:
        return '0/1', '3,3,24', 3, 3.0
    if c0==6 and c1==2:
        return '0/1', '3,3,45', 3, 3.0
    # Approximate adjustment of events with larger read depth
    c0, c1 = rescale_read_counts(c0, c1) # DR, DV

    ori_GL00 = np.float64(pow((1-err), c0)*pow(err, c1)*(1-prior)/2)
    ori_GL11 = np.float64(pow(err, c0)*pow((1-err), c1)*(1-prior)/2)
    ori_GL01 = np.float64(pow(0.5, c0+c1)*prior)

    # normalized genotype likelihood
    prob = list(normalize_log10_probs([log10(ori_GL00), log10(ori_GL01), log10(ori_GL11)]))
    GL_P = [pow(10, i) for i in prob]
    PL = [int(np.around(-10*log10(i))) for i in GL_P]
    GQ = [int(-10*log10(GL_P[1] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[2])), int(-10*log10(GL_P[0] + GL_P[1]))]
    QUAL = abs(np.around(-10*log10(GL_P[0]), 1))

    return Genotype[prob.index(max(prob))], "%d,%d,%d"%(PL[0], PL[1], PL[2]), max(GQ), QUAL

def cal_CIPOS(std, num):
    pos = int(1.96 * std / num ** 0.5)
    return "-%d,%d"%(pos,pos)

def threshold_ref_count(num):
    if num <= 2:
        return 20*num
    elif 3 <= num <= 5:
        return 9*num 
    elif 6 <= num <= 15:
        return 7*num
    else:
        return 5*num

def count_coverage(chr, s, e, f, read_count, up_bound, itround):
    status = 0
    iteration = 0
    primary_num = 0
    for i in f.fetch(chr, s, e):
        iteration += 1
        if i.flag not in [0,16]:
            continue
        primary_num += 1
        if i.reference_start < s and i.reference_end > e:
            read_count.add(i.query_name)
            if len(read_count) >= up_bound:
                status = 1
                break
        if iteration >= itround:
            if float(primary_num/iteration) <= 0.2:
                status = 1
            else:
                status = -1
            break

    return status

def bench_phase(path, chr, sigs_index, pre_pos):
    readsfile = open("%sreads.pickle"%(path), 'rb')
    readsfile.seek(sigs_index["reads"][chr])
    reads_list=pickle.load(readsfile)
    readsfile.close()

    sort_list = list()
    idx = 0
    for i in reads_list:
        sort_list.append([i[0], 1, idx, i[5], i[3]])
        sort_list.append([i[1], 2, idx, i[5], i[3]])
        idx += 1
    idx = 0
    for i in pre_pos:
        sort_list.append([i[0], 3, idx])
        sort_list.append([i[1], 0, idx])
        idx += 1
    sort_list = sorted(sort_list, key = lambda x:(x[0], x[1]))
    svs_set = set()
    read_set = set()
    overlap_dict = dict()
    for node in sort_list:
        if node[1] == 1: # set2(read) left
            read_set.add(node[2])
            for x in svs_set:
                if pre_pos[x][1] == node[0]:
                    continue
                if x not in overlap_dict:
                    overlap_dict[x] = set()
                overlap_dict[x].add(node[2])
        elif node[1] == 2: # set2(read) right
            read_set.remove(node[2])
        elif node[1] == 3: # set1(sv) left
            svs_set.add(node[2])
            overlap_dict[node[2]] = set()
            for x in read_set:
                overlap_dict[node[2]].add(x)
        elif node[1] == 0: # set1(sv) right
            svs_set.remove(node[2])
    hp_dist_dict = dict()
    overlap_reads_dict = dict()
    for idx in overlap_dict:
        hp_count = [0, 0, 0]
        read_id_set = set()
        read_id2read = dict()
        overlap_reads_dict[idx] = list()
        for x in overlap_dict[idx]:
            hp = reads_list[x][5]
            read_id = reads_list[x][3]
            if read_id not in read_id2read:
                read_id2read[read_id] = reads_list[x]
            else:
                if reads_list[x][2] == 1:
                    read_id2read[read_id] = reads_list[x]
        for read_id in read_id2read:
            hp_count[read_id2read[read_id][5]] += 1
            overlap_reads_dict[idx].append(read_id2read[read_id])
        hp_dist_dict[idx] = hp_count
    return hp_dist_dict, overlap_reads_dict

def overlap_cover(svs_list, reads_list, CorO):
    sort_list = list()
    improve_overlap_cover(svs_list, reads_list, sort_list)
    sort_list = sorted(sort_list, key = lambda x:(x[0], x[1]))
    svs_set = set()
    read_set = set()
    overlap_dict = dict()
    cover_dict = dict()
    for node in sort_list:
        if node[1] == 1: # set2(read) left
            read_set.add(node[2])
            for x in svs_set:
                if svs_list[x][1] == node[0]:
                    continue
                if x not in overlap_dict:
                    overlap_dict[x] = set()
                overlap_dict[x].add(node[2])
        elif node[1] == 2: # set2(read) right
            read_set.remove(node[2])
        elif node[1] == 3: # set1(sv) left
            svs_set.add(node[2])
            overlap_dict[node[2]] = set()
            for x in read_set:
                overlap_dict[node[2]].add(x)
            cover_dict[node[2]] = set()
            for x in read_set:
                cover_dict[node[2]].add(x)
        elif node[1] == 0: # set1(sv) right
            svs_set.remove(node[2])
            temp_set = set()
            for x in read_set:
                temp_set.add(x)
            cover_dict[node[2]] = cover_dict[node[2]] & temp_set
    cover2_dict = dict()
    overlap2_dict = dict()
    phaseid_dict = dict()
    read_to_phase_dict = dict()
    iteration_dict = dict()
    primary_num_dict = dict()
    for idx in cover_dict:
        iteration_dict[idx] = len(overlap_dict[idx])
        primary_num_dict[idx] = 0
        for x in overlap_dict[idx]:
            if reads_list[x][2] == 1:
                primary_num_dict[idx] += 1
        overlap2_dict[idx] = set()
        for x in overlap_dict[idx]:
            if reads_list[x][2] == 1:
                overlap2_dict[idx].add(reads_list[x][3])
                read_to_phase_dict[reads_list[x][3]] = reads_list[x][5]
        cover2_dict[idx] = set()
        for x in cover_dict[idx]:
            if reads_list[x][2] == 1:
                if reads_list[x][3] not in cover2_dict[idx]:
                    cover2_dict[idx].add(reads_list[x][3])
                    assert reads_list[x][3] in read_to_phase_dict, "covered reads not in overlap?"
    if CorO == 'C':
        return iteration_dict, primary_num_dict, cover2_dict, read_to_phase_dict
    else:
        return iteration_dict, primary_num_dict, overlap2_dict, read_to_phase_dict

def assign_gt(iteration_dict, primary_num_dict, cover_dict, read_id_dict, phased_id_dict, read_to_phase_dict, candidate_single_SV, search_flag_list):
    def cal_bayes(kdv, ndr):
        err = 0.1
        prior = 0.7022
        ht1 = np.float64(pow((1-err), kdv)*pow(err, ndr)*prior)
        ht2 = np.float64(pow((1-err), ndr)*pow(err, kdv)*(1-prior))
        return True if ht1 > ht2 else False
    assign_list = list()
    for idx in read_id_dict:
        iteration = iteration_dict[idx]
        primary_num = primary_num_dict[idx]
        read_count = cover_dict[idx]
        support_read_count = phased_id_dict[idx]
        ref_read_count = [0, 0, 0]
        
        DV = len(read_id_dict[idx])        
        DR = 0
        for query in read_count:
            if query not in read_id_dict[idx]: # support_read_set
                DR += 1
                query_phase = read_to_phase_dict[query]
                ref_read_count[query_phase] += 1
        bayes_h1 = cal_bayes(support_read_count[1], ref_read_count[1])
        bayes_h2 = cal_bayes(support_read_count[2], ref_read_count[2])
        GT, GL, GQ, QUAL = cal_GL(DR, DV)
        if GT == '0/1':
            GT = '1|0' if support_read_count[1] > support_read_count[2] else '0|1'
        if DV / (DV + DR) < 0.25 and QUAL < 20 and (len(search_flag_list) != 0 and not search_flag_list[idx]):
            QUAL = 0
        hp_pro = [0, 0]
        for _ in [0, 1]:
            __ = _ + 1
            hp_pro[_] = 0 if support_read_count[__] + ref_read_count[__] == 0 else support_read_count[__] / (support_read_count[__] + ref_read_count[__])
        if min(support_read_count[1] + ref_read_count[1], support_read_count[2] + ref_read_count[2]) > 5:
            if hp_pro[0] > 0.95 and 0 <= hp_pro[1] < 0.05 and GT != '1|0':
                GT = '1|0'
            if hp_pro[1] > 0.95 and 0 <= hp_pro[0] < 0.05 and GT != '0|1':
                GT = '0|1'
        assign_list.append([DV, DR, GT, GL, GQ, QUAL])
    return assign_list

def assign_gt_mosaic(iteration_dict, primary_num_dict, cover_dict, read_id_dict, phased_id_dict, read_to_phase_dict):
    assign_list = list()
    for idx in read_id_dict:
        iteration = iteration_dict[idx]
        primary_num = primary_num_dict[idx]
        read_count = cover_dict[idx]
        support_read_count = phased_id_dict[idx]
        ref_read_count = [0, 0, 0]
        DR = 0
        DV = len(read_id_dict[idx])
        for query in read_count:
            if query not in read_id_dict[idx]:
                ref_read_count[read_to_phase_dict[query]] += 1
                DR += 1
        
        GT, GL, GQ, QUAL = cal_GL(DR, len(read_id_dict[idx]))
        assign_list.append([len(read_id_dict[idx]), DR, GT, GL, GQ, QUAL])
    return assign_list

def generate_output(args, semi_result, reference, chrom, temporary_dir):
    
    '''
    Generation of VCF format file.
    VCF version: 4.2
    '''

    def cal_filter(gt, qual):
        if gt == '.|0' or gt == '0|.' or gt == '.|.':
            return "q5"
        if qual == "." or qual == None:
            return "PASS"
        if float(qual) >= 5.0:
            return "PASS"
        return "q5"

    # genotype_trigger = TriggerGT[args.genotype]

    f=open("%sresults/%s.pickle"%(temporary_dir,chrom), "wb")

    semi_result.sort(key = lambda x:int(x[2]))
    action = args.genotype
    fa_file = pysam.FastaFile(reference)
    try:
        ref_chrom=fa_file.fetch(chrom)
    except:
        raise Exception("No corresponding contig in reference with %s."%(chrom))
    fa_file.close()
    lines=[]
    BATCH_SIZE=1000
    # [chr, svtype, pos, len, DV, CIPOS, CILEN, DR, GT, GL, GQ, QUAL, read_id_set]
    for i in semi_result:
        if i[1] in ["DEL", "INS"]:
            if abs(int(i[3])) > args.max_size and args.max_size != -1:
                continue
            if i[1] == "INS":
                cal_end = int(i[2])
            else:
                cal_end = int(i[2]) + abs(int(i[3]))
            info_list = "{PRECISION};SVTYPE={SVTYPE};SVLEN={SVLEN};END={END};CIPOS={CIPOS};CILEN={CILEN};RE={RE};RNAMES={RNAMES}".format(
                PRECISION = "IMPRECISE" if i[8] == "0/0" else "PRECISE", 
                SVTYPE = i[1], 
                SVLEN = i[3] if i[1] == 'INS' else '-' + str(int(i[3])), 
                END = str(cal_end), 
                CIPOS = i[5], 
                CILEN = i[6], 
                RE = i[4],
                RNAMES = i[12] if i[12] == "FROMBAD" else "NULL")
            if action:
                try:
                    info_list += ";AF=" + str(round(int(i[4]) / (int(i[4]) + int(i[7])), 4))
                except:
                    info_list += ";AF=."
            if i[1] =="DEL":
                info_list += ";STRAND=+-"
            filter_lable = cal_filter(i[8], i[11])
            lines.append((i[1],"{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}\t{GT}:{DR}:{RE}:{PL}:{GQ}\n".format(
                CHR = i[0], 
                POS = str(int(i[2])), 
                ID = "cuteSV.%s.<SVID>"%(i[1]),
                REF = ref_chrom[max(int(i[2])-1, 0)] if i[1] == 'INS' else ref_chrom[max(int(i[2])-1, 0):int(i[2])+int(i[3])],
                ALT = "%s"%(ref_chrom[max(int(i[2])-1, 0)]+i[13] if i[1] == 'INS' else ref_chrom[max(int(i[2])-1, 0)]), 
                INFO = info_list, 
                FORMAT = "GT:DR:DV:PL:GQ", 
                GT = i[8],
                DR = i[7],
                RE = i[4],
                PL = i[9],
                GQ = i[10],
                QUAL = i[11],
                PASS = filter_lable)))
        elif i[1] == "DUP":
            if abs(int(float(i[3]))) > args.max_size and args.max_size != -1:
                continue
            cal_end = int(i[2]) + 1 + abs(int(i[3]))
            info_list = "{PRECISION};SVTYPE={SVTYPE};SVLEN={SVLEN};END={END};RE={RE};STRAND=-+;RNAMES={RNAMES}".format(
                PRECISION = "IMPRECISE" if i[6] == "0/0" else "PRECISE", 
                SVTYPE = i[1], 
                SVLEN = i[3], 
                END = str(cal_end), 
                RE = i[4],
                RNAMES = i[10] if args.report_readid else "NULL")
            if action:
                try:
                    info_list += ";AF=" + str(round(int(i[4]) / (int(i[4]) + int(i[5])), 4))
                except:
                    info_list += ";AF=."
            filter_lable = cal_filter(i[6], i[9])
            lines.append((i[1],"{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}\t{GT}:{DR}:{RE}:{PL}:{GQ}\n".format(
                CHR = i[0], 
                POS = str(int(i[2]) + 1), 
                ID = "cuteSV.%s.<SVID>"%(i[1]),
                REF = ref_chrom[int(i[2])],
                ALT = "<%s>"%(i[1]), 
                INFO = info_list, 
                FORMAT = "GT:DR:DV:PL:GQ", 
                GT = i[6],
                DR = i[5],
                RE = i[4],
                PL = i[7],
                GQ = i[8],
                QUAL = i[9],
                PASS = filter_lable)))
        elif i[1] == "INV":
            if args.large_inv or args.max_size == -1 or abs(int(float(i[3]))) <= args.max_size:
                cal_end = int(i[2]) + 1 + abs(int(i[3]))
                info_list = "{PRECISION};SVTYPE={SVTYPE};SVLEN={SVLEN};END={END};RE={RE};STRAND={STRAND};RNAMES={RNAMES}".format(
                    PRECISION = "IMPRECISE" if i[6] == "0/0" else "PRECISE", 
                    SVTYPE = i[1], 
                    SVLEN = i[3], 
                    END = str(cal_end), 
                    RE = i[4],
                    STRAND = i[7],
                    RNAMES = i[11] if args.report_readid else "NULL")
                if action:
                    try:
                        info_list += ";AF=" + str(round(int(i[4]) / (int(i[4]) + int(i[5])), 4))
                    except:
                        info_list += ";AF=."
                filter_lable = cal_filter(i[6], i[10])
                lines.append((i[1],"{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}\t{GT}:{DR}:{RE}:{PL}:{GQ}\n".format(
                    CHR = i[0], 
                    POS = str(int(i[2]) + 1), 
                    ID = "cuteSV.%s.<SVID>"%(i[1]),
                    REF = ref_chrom[int(i[2])],
                    ALT = "<%s>"%(i[1]), 
                    INFO = info_list, 
                    FORMAT = "GT:DR:DV:PL:GQ", 
                    GT = i[6],
                    DR = i[5],
                    RE = i[4],
                    PL = i[8],
                    GQ = i[9],
                    QUAL = i[10],
                    PASS = filter_lable)))
        else:
            # BND
            # info_list = "{PRECISION};SVTYPE={SVTYPE};CHR2={CHR2};END={END};RE={RE};RNAMES={RNAMES}".format(
            info_list = "{PRECISION};SVTYPE={SVTYPE};RE={RE};RNAMES={RNAMES}".format(
                PRECISION = "IMPRECISE" if i[7] == "0/0" else "PRECISE", 
                SVTYPE = "BND", 
                # CHR2 = i[3], 
                # END = str(int(i[4]) + 1), 
                RE = i[5],
                RNAMES = i[11] if args.report_readid else "NULL")
            
            filter_lable = cal_filter(i[7], i[10])
            try:
                ref_bnd = ref_chrom[int(i[2])]
            except:
                ref_bnd = 'N'
            lines.append(("BND","{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{PASS}\t{INFO}\t{FORMAT}\t{GT}:{DR}:{RE}:{PL}:{GQ}\n".format(
                CHR = i[0], 
                POS = str(int(i[2]) + 1), 
                ID = "cuteSV.%s.<SVID>"%("BND"), 
                REF = ref_bnd,
                ALT = i[1], 
                INFO = info_list, 
                FORMAT = "GT:DR:DV:PL:GQ", 
                GT = i[7],
                DR = i[6],
                RE = i[5],
                PL = i[8],
                GQ = i[9],
                QUAL = i[10],
                PASS = filter_lable)))
        if len(lines)>BATCH_SIZE:
            # print(len(lines))
            pickle.dump(lines,f)
            lines=[]
    if len(lines)!=0:
        pickle.dump(lines,f)

          
def load_valuable_chr(path):
    valuable_chr = dict()
    valuable_chr["DEL"] = list()
    valuable_chr["DUP"] = list()
    valuable_chr["INS"] = list()
    valuable_chr["INV"] = list()
    valuable_chr["TRA"] = dict()

    for svtype in ["DEL", "DUP", "INS", "INV"]:
        file = open("%s%s.sigs"%(path, svtype), 'r')
        for line in file:
            chr = line.strip('\n').split('\t')[1]
            if chr not in valuable_chr[svtype]:
                valuable_chr[svtype].append(chr)
        file.close()
        valuable_chr[svtype].sort()

    file = open("%s%s.sigs"%(path, "TRA"), 'r')
    for line in file:
        chr1 = line.strip('\n').split('\t')[1]
        chr2 = line.strip('\n').split('\t')[4]
        
        if chr1 not in valuable_chr["TRA"]:
            valuable_chr["TRA"][chr1] = list()
        if chr2 not in valuable_chr["TRA"][chr1]:
            valuable_chr["TRA"][chr1].append(chr2)

    file.close()
    for chr1 in valuable_chr["TRA"]:
        valuable_chr["TRA"][chr1].sort()

    return valuable_chr

def load_bed(bed_file, Task_list):
    # Task_list: [[chr, start, end], ...]
    bed_regions = dict()
    if bed_file != None:
        # only consider regions in BED file
        with open(bed_file, 'r') as f:
            for line in f:
                seq = line.strip().split('\t')
                if seq[0] not in bed_regions:
                    bed_regions[seq[0]] = list()
                bed_regions[seq[0]].append((int(seq[1]) - 1000, int(seq[2]) + 1000))
        region_list = [[] for i in range(len(Task_list))]
        for chrom in bed_regions:
            bed_regions[chrom].sort(key = lambda x:(x[0], x[1]))
            for item in bed_regions[chrom]:
                for i in range(len(Task_list)):
                    if chrom == Task_list[i][0]:
                        if (Task_list[i][1] <= item[0] and Task_list[i][2] > item[0]) or item[0] <= Task_list[i][1] < item[1]:
                            region_list[i].append(item)
        assert len(region_list) == len(Task_list), "parse bed file error"
        return region_list
    else:
        return None