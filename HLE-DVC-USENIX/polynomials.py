import copy
import numpy as np


def univariate_multiply_polynomials(p1, p2, modulus):
    """将两个一维多项式相乘，返回结果的系数列表"""
    result = [0] * (len(p1) + len(p2) - 1)
    for i in range(len(p1)):
        for j in range(len(p2)):
            result[i + j] += p1[i] * p2[j]
    result = [x % modulus for x in result]
    return result

def univariate_minus_polynomials(p1, p2, modulus):
    """将两个一维多项式相减，返回结果的系数列表，前面的函数减后面的函数。
    并且前面的函数的度数大
    """
    result = p1
    for i in range(len(p2)):
        result[-i-1] -= p2[-i-1]
    result = [x % modulus for x in result]
    return result


def evaluate_univariate_polynomial_at_x(coeff_list, x, modulus):
    """计算一元多项式在x处的取值"""
    result = 0
    n = len(coeff_list)
    for i in range(n):
        result += coeff_list[i] * (x ** (n - i - 1))
    return result % modulus


def univariate_polynomial_division(dividend, divisor, modulus):
    """
    对一元多项式进行有限域除法，返回商多项式和余数多项式。
    dividend 和 divisor 都是系数列表，list[0]代表最高次幂的系数。
    """
    # 初始化商和余数
    quotient = []
    remainder = dividend[:]
    # print("divisor:",divisor)
    # 被除多项式的最高次幂
    while len(remainder) >= len(divisor):
        # 当前项的系数和指数
        coeff = (remainder[0] * mod_inverse(divisor[0], modulus)) % modulus
        power_diff = len(remainder) - len(divisor)

        # 当前项的商
        # quotient_term = [0] * power_diff + [coeff]
        quotient.append(coeff)

        # 从余数中减去当前商项乘以除多项式
        subtract_term = [coeff * d for d in divisor] + [0] * power_diff
        # print("subtract_term",subtract_term)
        remainder = [(a - b) % modulus for a, b in zip(remainder, subtract_term)]

        # 去掉余数的前导零
        while remainder and remainder[0] == 0:
            remainder.pop(0)

    # 处理商多项式
    while len(quotient) < len(dividend) - len(divisor) + 1:
        quotient.insert(0, 0)

    return quotient, remainder



def mod_inverse(a, p):
    """计算模 p 下 a 的乘法逆元，使用扩展欧几里得算法。"""
    return pow(a, p-2, p)

if __name__ == "__main__":

    # 示例
    # dividend = [1, -3, 0, 2]  # 表示 x^3 - 3x^2 + 0x + 2
    # divisor = [1, -1]         # 表示 x - 1

    # quotient, remainder = univariate_polynomial_division(dividend, divisor, 23)
    # print(f"商多项式: {quotient}")
    # print(f"余数多项式: {remainder}")

    minuend = [1, 3, 0, 2]  # 表示 x^3 - 3x^2 + 0x + 2
    subtrahend = [1, 1]         # 表示 x - 1

    difference = univariate_minus_polynomials(minuend, subtrahend, 23)
    print(f"减法结果多项式: {difference}")



    # for row in result_coeffs:
    #     print(row)
    # 示例使用
    # 例如二元函数 M(x, y) = 2*x^2*y^2 + x*y^2 +3*x*y + 4y+4
    # coeff_matrix = [
    #     [2, 0, 0],  #
    #     [1, 3, 0],
    #     [0, 4, 4]
    # ]
    # b = 2
    # # 计算一元函数 f(x) = M(x, b)
    # f_x_coeffs = evaluate_bivariate_polynomial_at_y(coeff_matrix, b)

    # print("f(x) 的系数:", f_x_coeffs)  # 输出: [12, 6, 8]

    # 示例使用
    # bivariate_polynomial_m_xy = [[-1.1666666666666667, 3.0, -0.8333333333333339, -1.0], [0.16666666666666674, 0.0, 0.8333333333333335, 1.0]]
    # univariate_polynomial_p_y =
    # output = gen_bivariate_polynomial_Mxing_xy(bivariate_polynomial_m_xy, univariate_polynomial_p_y, l= 2, n = 4, modulus = 23)
    # 示例
    # coeff_list = [3, 2, 1]  # 3x^2 + 2x + 1

    # x = 2

    # value = evaluate_univariate_polynomial_at_x(coeff_list, x)
    # print(f"一元多项式在x={x}处的值为: {value}")
        # 示例  xy + 2x +3y+4


    # bivariate_poly = [
    #     [1, 2],
    #     [3, 4],
    #     [5, 6]
    # ]
    # # 多项式除以 (x - 1)
    # # xy + 2x +3y+4 = (x - 1)(y+2)+ 4y + 6
    # linear_x = [1, -1]
    # quotient_x, remainder_x = bivariate_polynomial_division_x(bivariate_poly, linear_x)
    # print("Quotient for x:", quotient_x)
    # print("Remainder for x:", remainder_x)


# def bivariate_multiply_polynomials(M_coeffs, N_coeffs):
#     """将两个二元多项式相乘，返回结果的系数矩阵"""

#     # 确定结果矩阵的大小：行数是 M_coeffs 和 N_coeffs 的行数之和减 1，列数同理
#     num_rows_M = len(M_coeffs)
#     num_cols_M = len(M_coeffs[0])
#     num_rows_N = len(N_coeffs)
#     num_cols_N = len(N_coeffs[0])

#     result_rows = num_rows_M + num_rows_N - 1
#     result_cols = num_cols_M + num_cols_N - 1

#     # 初始化结果矩阵
#     result = [[0] * result_cols for _ in range(result_rows)]

#     # 进行多项式乘法
#     for i in range(num_rows_M):
#         for j in range(num_cols_M):
#             for k in range(num_rows_N):
#                 for l in range(num_cols_N):
#                     result[i + k][j + l] += M_coeffs[i][j] * N_coeffs[k][l]

#     return result


def gen_bivariate_polynomial_u_xy(bivariate_polynomial_m_xy,univariate_polynomial_f_x, l, n, modulus):
    # 明面上，应该计算u(x,y)=(m(x,y)-1/n f(x))/y，但是有更简单的方案（正确性验证读者可以尝试下列print函数）
    # print("bivariate_polynomial_m_xy:", bivariate_polynomial_m_xy)
    # print("univariate_polynomial_f_x:", (univariate_polynomial_f_x[0] * mod_inverse(n, modulus))%modulus)
    # 初始化一个空系数矩阵，行列分别代表 x 和 y 的幂次
    bivariate_polynomial_u_xy = [[0 for _ in range(n-1)] for _ in range(l)]
    for i in range(l):
        for j in range(n-1):
            bivariate_polynomial_u_xy[i][j] = bivariate_polynomial_m_xy[i][j]
    # print(len(bivariate_polynomial_u_xy[0]))  # 3
    # # 计算偏移多项式
    # bivariate_polynomial_upie_xy = [[0 for _ in range(n)] for _ in range(l)]
    # for i in range(l):
    #     for j in range(n):
    #         bivariate_polynomial_upie_xy[i][j] = bivariate_polynomial_m_xy[i][j]
    # for i in range(l-1):
    #     bivariate_polynomial_upie_xy[i][n-1] = 0
    # print("bivariate_polynomial_m_xy:", bivariate_polynomial_m_xy)
    # print("bivariate_polynomial_u_xy:", bivariate_polynomial_u_xy)
    # print("bivariate_polynomial_upie_xy:", bivariate_polynomial_upie_xy)
    return bivariate_polynomial_u_xy

# def phase_three(bivariate_polynomial_m_xy, univariate_polynomial_p_y,xs, ys, tau_r, tau_x, tau_y, l, n, modulus):
#     # 先计算二元多项式M(x,y) = m(x,y)(m(x,y)-p(y))
#     # tmp_polynoimal = m()
#     tmp_polynoimal = copy.deepcopy(bivariate_polynomial_m_xy)
#     for j in range(n):
#         tmp_polynoimal[-1][j] = bivariate_polynomial_m_xy[-1][j]-univariate_polynomial_p_y[j]
#     bivariate_polynomial_M_xy = bivariate_multiply_polynomials(bivariate_polynomial_m_xy, tmp_polynoimal) #3*7

#     # print(evaluate_bivariate_polynomial_at_x_y(bivariate_polynomial_M_xy, xs[1], ys[0])% modulus) 全都是0
#     # 当r确定，M*是一元函数
#     univariate_polynomial_Mxing_x = [0] * len(bivariate_polynomial_M_xy) # len = 3
#     # print("len(bivariate_polynomial_M_xy):",len(bivariate_polynomial_M_xy))  3*7
#     # print("len(ys):",len(ys)) 4
#     for j in range(n):
#         for i in range(len(bivariate_polynomial_M_xy)):
#             tmp_polynomial = evaluate_bivariate_polynomial_at_y(bivariate_polynomial_M_xy, ys[j])
#             # print(evaluate_univariate_polynomial_at_x(tmp_polynomial, xs[1])% modulus)
#             # print("len(tmp_polynomial):",len(tmp_polynomial)) 3
#             # 二元随机多项式R(r,y)=(r^n-y^n)/(r-y)
#             value_R = (pow(tau_r, n, modulus)- pow(ys[j], n, modulus)) * mod_inverse(tau_r-ys[j], modulus)
#             # value_R = 1
#             univariate_polynomial_Mxing_x[i] += tmp_polynomial[i]* value_R
#     # vanishing_polynomial_vHl_x = x^l - 1
#     vanishing_polynomial_vHl_x = [0] * (l+1)
#     vanishing_polynomial_vHl_x[0] = 1
#     vanishing_polynomial_vHl_x[-1] = -1 % modulus
#     print(evaluate_univariate_polynomial_at_x(univariate_polynomial_Mxing_x, xs[1])% modulus) # 全是0
#     print(vanishing_polynomial_vHl_x)
#     # print(len(univariate_polynomial_Mxing_x)) 3
#     univariate_polynomial_e_x, remainder_x = univariate_polynomial_division(univariate_polynomial_Mxing_x, vanishing_polynomial_vHl_x, modulus)
#     value_Mxing = evaluate_univariate_polynomial_at_x(univariate_polynomial_Mxing_x, tau_x)
#     assert remainder_x == [], "多项式不能整除"

#     # 这步很难理解，就是把系数转换，通过转置矩阵把X换成Y了，这样就不用再写另外的函数了。
#     # 使用列表推导进行转置
#     transposed_bivariate_polynomial_M_xy = [[bivariate_polynomial_M_xy[j][i] for j in range(len(bivariate_polynomial_M_xy))] for i in range(len(bivariate_polynomial_M_xy[0]))]
#     bivariate_polynomial_r_ry = [[0 for _ in range(n)] for _ in range(n)]
#     # 根据次方差公式，写出r(r,y)
#     for i in range(n):
#         bivariate_polynomial_r_ry[i][n-1-i] =1
#     univariate_polynomial_Mr_y = \
#         univariate_multiply_polynomials( evaluate_bivariate_polynomial_at_y(transposed_bivariate_polynomial_M_xy, tau_x), \
#                                         evaluate_bivariate_polynomial_at_y(bivariate_polynomial_r_ry, tau_r))
#     # vanishing_polynomial_vHl_x = x^n - 1
#     vanishing_polynomial_vHn_y = [0] * (n+1)
#     vanishing_polynomial_vHn_y[0] = 1
#     vanishing_polynomial_vHn_y[-1] = -1 % modulus
#     univariate_polynomial_q_y, remainder_y = univariate_polynomial_division(univariate_polynomial_Mr_y, vanishing_polynomial_vHn_y, modulus)
#     univariate_polynomial_g_y, remainder= univariate_polynomial_division(remainder_y, [1,0], modulus)
#     assert remainder[0] % modulus== (value_Mxing * mod_inverse(n, modulus))% modulus

#     value_m = evaluate_bivariate_polynomial_at_x_y(bivariate_polynomial_m_xy, tau_x,tau_y)

#     value_p = evaluate_univariate_polynomial_at_x(univariate_polynomial_p_y, tau_y)



#     assert (value_m * (value_m - value_p) * (pow(tau_r,n, modulus)-pow(tau_y, n, modulus))* mod_inverse(tau_r - tau_y, modulus)) % modulus == \
#         evaluate_univariate_polynomial_at_x(univariate_polynomial_Mr_y, tau_y) % modulus, "等式左边计算错误"

#     # 验证者验证




#     return univariate_polynomial_e_x, value_Mxing, vanishing_polynomial_vHl_x, vanishing_polynomial_vHn_y, univariate_polynomial_q_y, univariate_polynomial_g_y
