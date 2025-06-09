def univariate_multiply_polynomials(p1, p2):
    """将两个一维多项式相乘，返回结果的系数列表"""
    result = [0] * (len(p1) + len(p2) - 1)
    for i in range(len(p1)):
        for j in range(len(p2)):
            result[i + j] += p1[i] * p2[j]
    return result


def lagrange_basis_polynomial(points, k, modulus):
    """
    计算第k个插值点对应的拉格朗日基函数的多项式形式
    """
    n = len(points)
    basis_poly = [1]  # 初始化为常数1的多项式
    x_k = points[k]

    for i in range(n):
        if i != k:
            factor = [1, -points[i]]  # (x - x_i)
            basis_poly = univariate_multiply_polynomials(basis_poly, factor)
            denominator = x_k - points[i]
            # basis_poly = [coef / denominator for coef in basis_poly]
            basis_poly = [coef * mod_inverse(denominator, modulus)  for coef in basis_poly]


    return basis_poly

def bivariate_lagrange_interpolation_coefficients(x_points, y_points, z_points, modulus):
    """
    计算二元函数的拉格朗日插值的系数矩阵
    特别地,返回的list中coefficients[0]为最高位，对应x^{n-1}
    """
    n = len(x_points)
    m = len(y_points)

    # 初始化一个空系数矩阵，行列分别代表 x 和 y 的幂次
    coefficients = [[0 for _ in range(m)] for _ in range(n)]

    for i in range(n):
        for j in range(m):
            l_x = lagrange_basis_polynomial(x_points, i, modulus)
            m_y = lagrange_basis_polynomial(y_points, j, modulus)

            # 对应的z值乘以两个基多项式的乘积
            for p in range(len(l_x)):
                for q in range(len(m_y)):
                    coefficients[p][q] += z_points[i * m + j] * l_x[p] * m_y[q] % modulus
            for p in range(len(l_x)):
                for q in range(len(m_y)):
                    coefficients[p][q] = coefficients[p][q] % modulus
    return coefficients







def mod_inverse(a, p):
    """计算模 p 下 a 的乘法逆元，使用扩展欧几里得算法。"""
    return pow(a, p-2, p)

def univariate_lagrange_interpolation_coefficients(x_points, y_points, modulus):
    n = len(x_points)
    # 初始化多项式系数列表，最多有 n 个系数，从 x^(n-1) 到 x^0
    coefficients = [0] * n

    for i in range(n):
        # 初始化当前基多项式的系数 (l_i(x))
        basis_poly = [1]

        for j in range(n):
            if i != j:
                # (x - xj) / (xi - xj) 部分的多项式展开
                factor = [1, -x_points[j]]  # x - xj 对应的多项式 [1, -xj]
                basis_poly = univariate_multiply_polynomials(basis_poly, factor)
                # 除以分母 (xi - xj)
                denominator = (x_points[i] - x_points[j]) % modulus
                basis_poly = [coef * mod_inverse(denominator, modulus)  for coef in basis_poly]

        # 将基多项式乘以 yi，并累加到最终结果
        for k in range(n):
            coefficients[k] += basis_poly[k] * y_points[i]

    return [ coef % modulus for coef in coefficients ]


if __name__ == "__main__":
    # 示例使用
    x_points = [0, 1, 2]
    y_points = [5, 3, 7]

    coefficients = univariate_lagrange_interpolation_coefficients(x_points, y_points, 23)
    print("Lagrange Interpolation Polynomial Coefficients:", coefficients)

    #     # 示例使用
    # x_points = [1, 2, 3]
    # y_points = [1, 2, 3]
    # z_points = [1, 4, 9, 16, 25, 36, 49, 64, 81]  # 假设这是某二元函数的z值

    # # 获取二元函数的插值系数矩阵
    # coefficients = bivariate_lagrange_interpolation_coefficients(x_points, y_points, z_points, 23)
    # print("Lagrange Interpolation Coefficients Matrix:")
    # # print(coefficients)
    # for row in coefficients:
    #     print(row)




