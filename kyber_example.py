import math


# 다항식 NTT 변환
def ntt_transform(poly):
    # x^2 + 5 로 나눈 나머지
    a1 = [(poly[0] + poly[2] * -5) % 13, (poly[1] + poly[3] * -5) % 13]

    # x^2 + 8 로 나눈 나머지
    a2 = [(poly[0] + poly[2] * -8) % 13, (poly[1] + poly[3] * -8) % 13]

    return [a1, a2]


# 다항식 역 NTT 변환
def inverse_ntt_transform(ntt_result):
    # NTT 변환 결과
    a, b = ntt_result[0]
    c, d = ntt_result[1]

    # 역변환 계산
    f0 = ((a + c) * 7) % 13
    f1 = ((b + d) * 7) % 13
    f2 = ((a - c) * 9) % 13
    f3 = ((b - d) * 9) % 13

    # 결과 다항식
    return [f0, f1, f2, f3]


# NTT 곱셈
def mod_poly_multiply(a, b, x2_coeff):
    result = [0, 0, 0]  # 2차 다항식
    result[0] = (a[0] * b[0]) % 13
    result[1] = (a[0] * b[1] + a[1] * b[0]) % 13
    result[2] = (a[1] * b[1]) % 13

    # x^2에 대한 계수 대입
    result[0] = (result[0] + (result[2] * x2_coeff) % 13) % 13
    return result[:2]  # 1차 다항식 반환


# NTT 덧셈
def ntt_add(ntt1, ntt2):
    # NTT 변환된 다항식의 덧셈
    result = []
    for i in range(len(ntt1)):
        row = []
        for j in range(len(ntt1[i])):
            # ntt1[i][j]가 리스트이므로, 각 요소에 접근해야 함
            row.append([(ntt1[i][j][0] + ntt2[i][j][0]) % 13,
                         (ntt1[i][j][1] + ntt2[i][j][1]) % 13])
        result.append(row)
    return result


# 행렬과 벡터의 곱
def ntt_multiply_1(ntt_m, ntt_v):
    results = []

    p1 = ntt_m[0][0]
    p2 = ntt_m[0][1]
    p3 = ntt_m[1][0]
    p4 = ntt_m[1][1]
    p5 = ntt_v[0]
    p6 = ntt_v[1]

    a = mod_poly_multiply(p1[0], p5[0], -5)
    b = mod_poly_multiply(p2[0], p6[0], -5)
    first_result = [(a[0] + b[0]) % 13, (a[1] + b[1]) % 13]

    a = mod_poly_multiply(p1[1], p5[1], -8)
    b = mod_poly_multiply(p2[1], p6[1], -8)
    first_result_2 = [(a[0] + b[0]) % 13, (a[1] + b[1]) % 13]

    c = mod_poly_multiply(p3[0], p5[0], -5)
    d = mod_poly_multiply(p4[0], p6[0], -5)
    second_result = [(c[0] + d[0]) % 13, (c[1] + d[1]) % 13]

    c = mod_poly_multiply(p3[1], p5[1], -8)
    d = mod_poly_multiply(p4[1], p6[1], -8)
    second_result_2 = [(c[0] + d[0]) % 13, (c[1] + d[1]) % 13]

    results.append([first_result, first_result_2])
    results.append([second_result, second_result_2])

    return results


# 벡터와 벡터의 곱
def ntt_multiply_2(ntt_v1, ntt_v2):
    p1 = ntt_v1[0]
    p2 = ntt_v1[1]
    p3 = ntt_v2[0]
    p4 = ntt_v2[1]

    a = mod_poly_multiply(p1[0], p3[0], -5)
    b = mod_poly_multiply(p2[0], p4[0], -5)
    first_result = [(a[0] + b[0]) % 13, (a[1] + b[1]) % 13]

    c = mod_poly_multiply(p1[1], p3[1], -8)
    d = mod_poly_multiply(p2[1], p4[1], -8)
    first_result_2 = [(c[0] + d[0]) % 13, (c[1] + d[1]) % 13]

    return [first_result, first_result_2]


# 다항식 덧셈 함수
def poly_add(p1, p2):
    result = [0] * 4
    for i in range(4):
        result[i] = (p1[i] + p2[i]) % 13
    return result


# 다항식 뺄셈 함수
def poly_sub(p1, p2):
    result = [0] * 4
    for i in range(4):
        result[i] = (p1[i] - p2[i]) % 13
    return result


# 행렬 전치 함수
def transpose_matrix(matrix_a):
    return [[matrix_a[j][i] for j in range(2)] for i in range(2)]


# mu 계산 함수
def calculate_mu(poly_m):
    return [7 if coeff == 1 else 0 for coeff in poly_m]


# decompress, compress를 설명해야 한다면 활용할 수 있는 코드
def custom_round(n):
    if n - math.floor(n) == 0.5:
        return math.ceil(n)
    else:
        return round(n)


def decompress(p):
    return [custom_round(13/2 * coeff) for coeff in p]


def compress(p):
    return [custom_round(2/13 * coeff) for coeff in p]


# t 계산 함수
def calculate_t(ntt_m_A, ntt_v1_s, ntt_v2_e):
    A_s = ntt_multiply_1(ntt_m_A, ntt_v1_s)
    ntt_t = ntt_add(A_s, ntt_v2_e)
    return ntt_t


# u 계산 함수
def calculate_u(ntt_m_A_T, ntt_v1_r, vector_e1):
    A_r = ntt_multiply_1(ntt_m_A_T, ntt_v1_r)
    inverse_ntt_A_r = []
    for i in range(2):
        inverse_ntt_A_r.append(inverse_ntt_transform(A_r[i]))
    vector_u = [poly_add(inverse_ntt_A_r[i], vector_e1[i]) for i in range(2)]
    return vector_u


# v 계산 함수
def calculate_v(ntt_t_trans, ntt_v1_r, poly_e2, poly_mu):
    t_r = ntt_multiply_2(ntt_t_trans, ntt_v1_r)
    inverse_ntt_t_r = inverse_ntt_transform(t_r)
    poly_v = poly_add(poly_add(inverse_ntt_t_r, poly_e2), poly_mu)
    return poly_v


# w 계산 함수
def calculate_w(poly_v, ntt_v1_s, ntt_v2_r):
    s_r = ntt_multiply_2(ntt_v1_s, ntt_v2_r)
    inverse_ntt_s_r = inverse_ntt_transform(s_r)
    poly_w = poly_sub(poly_v, inverse_ntt_s_r)
    return poly_w


A = [[[4,0,3,6],[1,7,4,2]],[[0,8,2,1],[6,3,5,0]]]
s = [[1,0,1,1],[0,1,0,1]]
e = [[1,0,0,1],[0,1,0,1]]
e1 = [[0,0,1,1],[1,0,1,0]]
r = [[1,0,1,1],[1,1,0,0]]
e2 = [1,1,0,1]
m = [1,0,0,1]

ntt_A = [[ntt_transform(poly) for poly in row] for row in A]
ntt_s = [ntt_transform(poly) for poly in s]
ntt_e = [ntt_transform(poly) for poly in e]
ntt_r = [ntt_transform(poly) for poly in r]

t = calculate_t(ntt_A, ntt_s, ntt_e)
A_T = transpose_matrix(A)
ntt_A_T = [[ntt_transform(poly) for poly in row] for row in A_T]
u = calculate_u(ntt_A_T, ntt_r, e1)
ntt_u = [ntt_transform(poly) for poly in u]
mu = decompress(m)
v = calculate_v(t, ntt_r, e2, mu)
w = calculate_w(v, ntt_s, ntt_u)

print("\nResult t:")
print(t)

print("\nResult u:")
print(u)

print("\nResult v:")
print(v)

print("\nResult w:")
print(w)

print("\nResult compress(w):")
print(compress(w))




