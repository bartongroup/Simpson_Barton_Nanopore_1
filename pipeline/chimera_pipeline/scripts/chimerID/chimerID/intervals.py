import itertools as it


def flatten(bundle):
    flattened = []
    all_invs = iter(sorted(it.chain(*bundle)))
    inv_start, inv_end = next(all_invs)
    for start, end in all_invs:
        if start <= inv_end:
            inv_end = max(inv_end, end)
        else:
            flattened.append([inv_start, inv_end])
            inv_start, inv_end = start, end
    if not flattened or flattened[-1] != [inv_start, inv_end]:
        flattened.append([inv_start, inv_end])
    return flattened


def intersect(inv_a, inv_b):
    a_start, a_end = inv_a
    b_start, b_end = inv_b
    if a_end < b_start or a_start > b_end:
        return 0
    else:
        s = max(a_start, b_start)
        e = min(a_end, b_end)
        return e - s


def fraction_of_shorter(inv_a, inv_b):
    min_ln = min(inv_a[1] - inv_a[0], inv_b[1] - inv_b[0])
    return intersect(inv_a, inv_b) / min_ln
    

def intersect_spliced_invs(invs_a, invs_b):
    score = 0
    invs_a = iter(invs_a)
    invs_b = iter(invs_b)
    a_start, a_end = next(invs_a)
    b_start, b_end = next(invs_b)
    while True:
        if a_end < b_start:
            try:
                a_start, a_end = next(invs_a)
            except StopIteration:
                break
        elif a_start > b_end:
            try:
                b_start, b_end = next(invs_b)
            except StopIteration:
                break
        else:
            score += intersect([a_start, a_end], [b_start, b_end])
            if a_end > b_end:
                try:
                    b_start, b_end = next(invs_b)
                except StopIteration:
                    break
            else:
                try:
                    a_start, a_end = next(invs_a)
                except StopIteration:
                    break
    return score