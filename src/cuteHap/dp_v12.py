from scipy.stats import norm
import math

def solve_dp(input_triple):
    def overlap(l0, r0, l1, r1):
        l = max(l0, l1)
        r = min(r0, r1)
        return max(0, r - l + 1)
    
    def get_sum_info(key, data): # [data, sum(length), sum(center * len), key]
        ans = [data, [0 for _ in range(len(data) + 1)], [0 for _ in range(len(data) + 1)], key]
        for i in range(len(data)):
            length = data[i][1] - data[i][0] + 1; center = (data[i][0] + data[i][1]) / 2
            ans[1][i] = ans[1][i - 1] + length
            ans[2][i] = ans[2][i - 1] + length * center
        return ans

    def sim(l0, l1, o, c0, c1):
        sum = l0 + l1
        mul = l0 * l1
        a = (sum - abs(l0 - l1)) / mul
        b = max(mul / sum, o)
        c = min(l0, l1) / max(l0, l1)
        d = 1 - norm.cdf(6 * abs(c0 - c1) / sum)
        return a * b * c * d

    def sim_set(l0, r0, l1, r1, c0, c1):
        return sim(r0 - l0 + 1, r1 - l1 + 1, overlap(l0, r0, l1, r1), c0, c1)
    
    def sim_list(L, R, data, i, j):
        l = data[0][i][0]; r = data[0][j][1]
        sumL = data[1][j] - data[1][i - 1]
        sumC = data[2][j] - data[2][i - 1]
        return sim_set(L, R, l, r, (L + R) / 2, sumC / sumL)

    def get_set(infos, id_data, id_hp):
        ans = []
        for info in infos:
            key, l, r = info
            x = id_data[key]
            for i in range(l, r + 1):
                ans.append([key, x[i][0], x[i][1] - x[i][0] + 1, id_hp[key]])
        return ans

    sort_input = sorted(input_triple, key = lambda x:x[1])
    id_data = {}
    id_hp = {}
    L = set(); R = set()
    for _ in sort_input:
        id_data.setdefault(_[0], [])
        id_data[_[0]].append((_[1], _[1] + _[2] - 1))
        L.add(_[1]); R.add(_[1] + _[2] - 1)
        id_hp[_[0]] = _[4]
    data = [get_sum_info(key, _) for key, _ in id_data.items()]
    R.add(0); L = sorted(list(L)); R = sorted(list(R))

    score = [[[0 for z in R] for y in R] for x in L]
    for lid in range(len(L)):
        l = L[lid]
        start = []; end = []
        for _ in data:
            id = 0
            while id < len(_[0]) and _[0][id][0] < l:
                id += 1
            start.append(id); end.append(-2)
        for rid in range(len(R)):
            r = R[rid]
            if r <= l:
                continue
            ans_tmp = []
            for id in range(len(data)):
                _ = data[id][0]
                if end[id] == -2 and start[id] < len(_) and _[start[id]][1] <= r:
                    end[id] = start[id] - 1
                while end[id] != -2 and end[id] + 1 < len(_) and _[end[id] + 1][1] <= r:
                    end[id] += 1
                if end[id] != -2:
                    s = sim_list(l, r, data[id], start[id], end[id])
                    if s > 0.3:
                        if 0.5 < s < 0.8:
                            s = 2 * s * s - 0.48
                        elif 0.3 < s <= 0.5:
                            s = 0.1 * s - 0.03
                        sumL = data[id][1][end[id]] - data[id][1][start[id] - 1]
                        ans_tmp.append([_[end[id]][1], s * sumL / (r - l + 1)])
            ans_tmp.sort()
            now = len(ans_tmp) - 1
            ans = 0; cnt = 1
            for rid_ in reversed(range(rid + 1)):
                while now >= 0 and ans_tmp[now][0] > R[rid_]:
                    log_cnt = 1 if cnt == 1 else math.log(cnt)
                    ans /= cnt * log_cnt; cnt += 1
                    ans = (ans + ans_tmp[now][1]) * cnt * log_cnt; now -= 1
                score[lid][rid][rid_] = ans

    f = [[0 for z in R] for y in R]
    fk = [[0 for z in R] for y in R]
    fl = [[0 for z in R] for y in R]
    ans = 0; ansi = 0; ansj = 0
    for i in range(len(R)):
        for j in range(i + 1):
            fmax = f[j][0]; nowk = 0; k = 0
            for l in range(len(L)):
                while k + 1 < len(R) and R[k + 1] < L[l]:
                    k += 1
                    fnow = f[max(k, j)][min(k, j)]
                    if fnow > fmax:
                        fmax = fnow; nowk = k
                if f[i][j] < fmax + score[l][i][j]:
                    f[i][j] = fmax + score[l][i][j]
                    fk[i][j] = nowk; fl[i][j] = l
            if ans < f[i][j]:
                ans = f[i][j]; ansi = i; ansj = j
    
    ans = []
    while f[ansi][ansj] > 0:
        ansk = fk[ansi][ansj]; ansl = fl[ansi][ansj]
        ans.append([ansl, ansi, score[ansl][ansi][max(ansj, ansk)]])
        ansi = max(ansj, ansk)
        ansj = min(ansj, ansk)
    ans = [[L[i], R[j], s] for i, j, s in reversed(ans)]
    ans = sorted(ans, key = lambda x:x[1])
    
    score = [[0, [], []] for x in ans]
    for _ in data:
        f = [[0, -1, -1] for x in range(len(_[0]) + 1)]; nowi = 0
        for j in range(len(ans)):
            while nowi + 1 < len(_[0]) and _[0][nowi][1] < ans[j][1]:
                nowi += 1
                if f[nowi][0] < f[nowi - 1][0]:
                    f[nowi] = [f[nowi - 1][0], nowi - 1, -1]
            for i in range(max(0, nowi - 2), 5):
                if i < len(_[0]) and f[i][0] < f[i - 1][0]:
                    f[i] = [f[i - 1][0], i - 1, -1]
            nowk = nowi
            while nowk >= 0 and _[0][nowk][0] > ans[j][0]:
                nowk -= 1
            for i in reversed(range(max(0, nowi - 2), 5)):
                for k in range(max(0, nowk - 2), 5):
                    if k > i or i >= len(_[0]) or k >= len(_[0]):
                        continue
                    s = sim_list(ans[j][0], ans[j][1], _, k, i)
                    if s > 0.5 and f[i][0] < f[k - 1][0] + s * ans[j][2]:
                        f[i] = [f[k - 1][0] + s * ans[j][2], k - 1, j, (_[0][k][0], _[0][i][1]), (_[3], k, i)]
            for i in range(max(0, nowi - 2), 5):
                if i < len(_[0]) and f[i][0] < f[i - 1][0]:
                    f[i] = [f[i - 1][0], i - 1, -1]
        while nowi + 1 < len(_[0]):
                nowi += 1
                if f[nowi][0] < f[nowi - 1][0]:
                    f[nowi] = [f[nowi - 1][0], nowi - 1, -1]
        while nowi != -1:
            if f[nowi][2] != -1:
                score[f[nowi][2]][0] += 1
                score[f[nowi][2]][1].append(f[nowi][3])
                score[f[nowi][2]][2].append(f[nowi][4])
            nowi = f[nowi][1]
            
    real_ans = []
    for i in range(len(ans)):
        if score[i][0] > 1:
            L = 0; R = 0
            for l, r in score[i][1]:
                L += l; R += r
            L //= score[i][0]; R //= score[i][0]
            real_ans.append([L, R - L + 1, get_set(score[i][2], id_data, id_hp)])
    if len(real_ans) == 0:
        for i in range(len(ans)):
            if score[i][0] == 1:
                real_ans.append([score[i][1][0][0], score[i][1][0][1] - score[i][1][0][0] + 1, get_set(score[i][2], id_data, id_hp)])
    return real_ans
