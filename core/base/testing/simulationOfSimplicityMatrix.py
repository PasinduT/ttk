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

def mult(lis, d):
    terms = []
    maximum = 2 ** d
    count = 0
    counters = [0 for i in range(d)]

    while count < maximum:
        P = []
        E = []
        for i in range(d):
            temp = (1 << i)
            c = 1 if (count & temp) else 0
            item = lis[i][c]
            if (c == 0):
                P.append(item)
            else:
                E.append(item)

        # increase the count
        count += 1

        terms.append((P, E, lis[d]))

    return terms

        

def main(d = 2):
    matrix = []
    for i in range(1, (d+1) + 1):
        lis = [((i, j), (i, j)) for j in range(1, d + 1)]
        lis.append((1, 1))
        matrix.append(lis)

    for row in matrix:
        print('\t'.join(map(lambda x: f'P{x[0]} + E{x[1]}' if x[0] != 1 else '1', row)))
    # print(matrix)
    print()

    mainItems = []
    for i in range(d+1):
        items = []
        for j in range(d+1):
            row = j
            col = (i + j) % (d + 1)
            item = matrix[row][col]
            if (item[0] != 1):
                items.append(item)

        items.append('+')
        print(items)
        mainItems.append(items)
    for i in range(d+1):
        items = []
        for j in range(d+1):
            row = d+1 - j - 1
            col = (i + j) % (d + 1)
            item = matrix[row][col]
            if (item[0] != 1):
                items.append(item)

        items.append('-')
        print(items)
        mainItems.append(items)
    print()

    allTerms = []
    for items in mainItems:
        terms = mult(items, d)
        allTerms += terms
    print(allTerms)
    print()

    sortedTerms = sorted(allTerms, key=cmp_to_key(lambda x, y: -1 if smaller(x[1], y[1]) else 1))

    print('{:30s}{}'.format('COEFF', 'EPSILON'))
    for s in sortedTerms:
        P = s[2]
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

    main(d = 2)