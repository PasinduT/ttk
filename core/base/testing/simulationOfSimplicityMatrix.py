from functools import cmp_to_key

def smaller(I1, I2):
    diff1 = [x for x in I1 if x not in I2]
    if (len(diff1) == 0):
        return True
    i, j = max(diff1)

    diff2 = [x for x in I2 if x not in I1]
    if (len(diff2) == 0):
        return False
    k, l = max(diff2)

    if (i == k):
        return j > l
    return i < k

def mult(mainTerm, parity, det):
    terms = []
    P1, E1 = mainTerm
    if (len(P1) == 1):
        return [(P, E, parity * sign) for P, E, sign in det]
    assert(type(P1) is tuple)
    assert(type(E1) is tuple)
    # P1 is a tuple with a coefficient
    # E1 is a tuple with a epsilon term

    for term in det:
        P, E, sign = term
        # P is a list of coefficients
        # E is a list of epsilon terms
        assert (type(P) is list)
        assert (type(E) is list)
        assert (type(sign) is int)

        if (len(P) == 1 and len(P[0]) == 1):
            terms.append(([P1], [E1], sign * parity))
            continue

        terms.append(([P1] + P, [], sign * parity))

        terms.append(([P1], E, sign * parity))
        terms.append((P, [E1], sign * parity))
        terms.append(([], [E1] + E, sign * parity))
    
    for term in terms:
        P, E, sign = term
        # P is a list of coefficients
        # E is a list of epsilon terms
        assert (type(P) is list)
        assert (type(E) is list)
        assert (type(sign) is int)

    # assert(len(terms) == 4 * len(det))
    return terms

def minor(matrix, m, k):
    n = len(matrix)
    result = []
    for i in range(n):
        if (i == m):
            continue
        row = []
        for j in range(n):
            if (j == k):
                continue
            row.append(matrix[i][j])
        result.append(row)

    return result


def determinant(matrix):
    n = len(matrix)

    if (n == 1):
        P, E = matrix[0][0]
        return [([P], [E], 1)]

    allTerms = []
    parity = 1
    for i in range(n):
        mainTerm = matrix[0][i]
        assert(type(mainTerm) is tuple)

        mMatrix = minor(matrix, 0, i)
        det = determinant(mMatrix)
        # det is a list of tuples of (list, list, int)
        terms = mult(mainTerm, parity, det)
        # terms is a list of tuples of (list, list, int)
        allTerms = allTerms + terms

        parity *= -1

    for term in allTerms:
        P, E, sign = term
        # P is a list of coefficients
        # E is a list of epsilon terms
        assert (type(P) is list)
        assert (type(E) is list)
        assert (type(sign) is int)

    return allTerms

        

def main(d = 2):
    # matrix = []
    # for i in range(1, (d+1) + 1):
    #     lis = [((i, j), (i, j)) for j in range(1, d + 1)]
    #     lis.append(((1,), (1,)))
    #     matrix.append(lis)

    matrix = []
    for i in range(1, (d+1) + 1):
        lis = [((i, j), (i, j)) for j in range(1, d + 2)]
        # lis.append(((1,), (1,)))
        matrix.append(lis)

    for row in matrix:
        print('\t'.join(map(lambda x: f'P{x[0]} + E{x[1]}' if x[0] != 1 else '1', row)))
    # print(matrix)
    print()

    # computes the determinant of a matrix
    allTerms = determinant(matrix)
    # print(allTerms)
    # print(len(allTerms))

    sortedTerms = sorted(allTerms, key=cmp_to_key(lambda x, y: -1 if smaller(x[1], y[1]) else 1))

    print('{:30s}{}'.format('COEFF', 'EPSILON'))
    for s in sortedTerms:
        P = '+' if s[2] > 0 else '-'
        for p in s[0]:
            P += f'P{p} '
        E = ''
        for p in s[1]:
            E += f'E{p} '

        print('{:30s}{}'.format(P, E))

if __name__ == '__main__':
    # print(smaller([(1, 2)], [(1, 1)]))
    # print(smaller([(1, 1)], [(1, 2)]))
    # print(smaller([(1, 1), (2, 2)], [(1, 2), (2, 1)]))
    # print(smaller([(1, 1), (2, 2)], [(1, 2), (2, 1)]))
    # print(smaller([(2, 1)], [(2, 1)]))
    # print(smaller([(1, 1), (2, 2)], [(1, 1)]))

    main(d = 1)