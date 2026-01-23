from scipy.stats import norm
import bisect
import math
import copy

# interval: [pos, length, length]
# point: [pos, 1, length]

def get_intervals(input_triples, is_interval):
    return sorted([(_[1], _[2] if is_interval else 1, _[2], _) for _ in input_triples], key=lambda x:x[0])

def intervals_splited(intervals):
    id_map = dict()
    for _ in range(len(intervals)):
        i = intervals[_][3]
        id_map.setdefault(i[0], [])
        id_map[i[0]].append((i[1], _))
    return [[_[1] for _ in lst] for lst in id_map.values()]

def sim(l0, l1, p0, p1):
    len_ratio = max(min(l0, l1) / max(l0, l1) * 2.5 - 1.5, 0)
    pos_bias = max(1 - math.log10(1 + abs(p0 - p1)) * 0.25, 0)
    return len_ratio * pos_bias

def iscore(i0, i1):
    s = 1.0 if len(i1) < 4 else max(i0[2] - i1[3], 0) / i0[2]
    return sim(i0[2], i1[2], i0[0], i1[0]) * s

def combine(i0, i1):
    return (i0[0], i1[0] + i1[1] - i0[0], i0[2] + i1[2], i1[0] - i0[0] - i0[1] + (0 if len(i0) < 4 else i0[3]))

def combined_map(splited, intervals):
    ans = dict()
    for lst in splited:
        for i in range(len(lst)):
            s = intervals[lst[i]]
            s = (s[0], s[1], s[2])
            ans[(lst[i], lst[i])] = s
            for j in range(i + 1, len(lst)):
                s = combine(s, intervals[lst[j]])
                ans[(lst[i], lst[j])] = s
    return ans

class BeamParams:
    MinViableSimilarity = 0.4
    MinClusterCenterSimilarity = 0.95
    MaxDiffClusterSimilarity = 0.5

# f[i, j] = max{f[k, l] + iscore(I[j], combine(p, i)) * L} (k < p <= i, l < j)
# g[i, j] = max{f[k, l]} (k <= i, l <= j)
def match_intervals(aims, intervals, splited, combined, params):
    dp_ans = [[] for _ in aims]
    for id in range(len(splited)):
        lst = splited[id]
        valid_sums = [intervals[lst[0]][2]]
        for i in range(1, len(lst)):
            valid_sums.append(intervals[lst[i]][2] + valid_sums[-1])
        dp = [[0 for _ in range(len(aims) + 1)] for i in range(len(lst) + 1)]
        v = copy.deepcopy(dp)
        for i in range(len(lst)):
            for j in range(len(aims)):
                dp[i + 1][j + 1] = max(dp[i + 1][j], dp[i][j + 1])
                v[i + 1][j + 1] = v[i + 1][j] if dp[i + 1][j] >= dp[i][j + 1] else v[i][j + 1]
                _k = bisect.bisect_left(valid_sums, valid_sums[i] - aims[j][2])
                k, score = 0, 0
                for t in range(max(_k - 5, 0), min(_k + 5, i + 1)):
                    s = iscore(aims[j], combined[(lst[t], lst[i])])
                    if s > score:
                        k, score = t, s
                if score < params.MinViableSimilarity:
                    continue
                va = combined[(lst[k], lst[i])][2]
                if dp[i + 1][j + 1] < score * va + dp[k][j]:
                    dp[i + 1][j + 1] = score * va + dp[k][j]
                    v[i + 1][j + 1] = (k, score, i, j)
        i, j = len(lst), len(aims)
        while v[i][j] != 0:
            k, s, i0, j0 = v[i][j]
            dp_ans[j0].append((s, id, k, i0))
            i, j = k, j0
    score = 0.0
    for i in range(len(aims)):
        cnt = 0
        p, n = 1.0, 1.0
        for x in sorted([_[0] for _ in dp_ans[i]], reverse=True)[:10]:
            if x >= params.MinClusterCenterSimilarity:
                cnt += 1
            p *= x
            n *= 1.0 - x
        if cnt < 2:
            return -10000, []
        score += math.log(p) - math.log(p + n)
    return score, dp_ans

def output(answers, dp_ans, intervals, splited, combined):
    ans = []
    for i in range(len(answers)):
        now = [answers[i][0], answers[i][2]]
        s = []
        for _ in dp_ans[i]:
            x = combined[splited[_[1]][_[2]], splited[_[1]][_[3]]]
            ratio = min(x[2], answers[i][2]) / max(x[2], answers[i][2])
            if abs(x[0] - answers[i][0]) > 1000 or ratio < 0.7:
                continue
            for j in range(_[2], _[3] + 1):
                s.append(intervals[splited[_[1]][j]][3])
        if len(s) > 0:
            now.append(s)
            ans.append(now)
    return ans

def pre_generate(intervals, splited, params):
    combined = combined_map(splited, intervals)
    a, ans = [], []
    for lst in splited:
        for i in range(len(lst)):
            s = intervals[lst[i]]
            s = (s[0], s[1], s[2])
            a.append(s)
            for j in range(i + 1, len(lst)):
                s = combine(s, intervals[lst[j]])
                if s[3] > s[2] * 3:
                    break
                a.append(s)
    a = set(a)
    for i in a:
        score = match_intervals((i,), intervals, splited, combined, params)[0]
        if score > math.log(0.5):
            ans.append((score, i))
    ans = [_[1] for _ in sorted(ans, reverse=True, key=lambda x:x[0])]
    return combined, ans

def two_rows(inputs):
    tasks = dict()
    for _ in inputs:
        tasks.setdefault(_[0], [0, 0])
        tasks[_[0]][0] += 1
        tasks.setdefault(_[0] + _[1], [0, 0])
        tasks[_[0] + _[1]][1] += 1
    now = 0
    for _ in sorted(tasks.keys()):
        a, b = tasks[_]
        now += a
        if now > 2:
            return False
        now -= b
    return True

def beam_generate(answers, intervals, splited, combined, chooses, beam_size, params):
    vis = set(answers)
    ans = []
    for choose in chooses:
        if choose in vis:
            continue
        check = True
        for _ in answers:
            if iscore(choose, _) > params.MaxDiffClusterSimilarity:
                check = False
                break
        if not check:
            continue
        if not two_rows(answers):
            continue
        new_ans = (*answers, choose)
        s, _ = match_intervals(new_ans, intervals, splited, combined, params)
        if s <= -10000:
            continue
        ans.append((s, new_ans))
        if len(ans) == beam_size >> 1:
            break
    return ans

def solve_beam(input_triples, is_interval = True, beam_size = 20, params=BeamParams()):
    intervals = get_intervals(input_triples, is_interval)
    splited = intervals_splited(intervals)
    combined, chooses = pre_generate(intervals, splited, params)
    ans = [(0, ())]
    while len(ans) > 0:
        final = ans[0][1]
        new_ans = []
        for _ in ans:
            new_ans.extend(beam_generate(_[1], intervals, splited, combined, chooses, beam_size, params))
        ans = []
        s = set()
        for _ in new_ans:
            x = tuple(sorted(_[1]))
            if x in s:
                continue
            s.add(x)
            ans.append(_)
        ans = sorted(ans, reverse=True, key=lambda x:x[0])[:beam_size]
    s, d = match_intervals(final, intervals, splited, combined, params)
    return output(final, d, intervals, splited, combined)
