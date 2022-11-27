# Example for LU Factrization, QR Factrization, Givens Reduction, Housholder Reduction and URV Factrization. 
# 3x3 非奇异矩阵
m1 = [[1.0, 2.0, 3.0], [2.0, 2.0, 8.0], [-3.0, -10.0, -2.0]]
b1 = [0.0, -4.0, -11.0]
# 3x3 非奇异矩阵
m2 = [[2.0, -3.0, 1.0], [1.0, -1.0, -2.0], [3.0, 1.0, -1.0]]
b2 = [7.0, -2.0, 0.0]
# 3x3 非奇异矩阵
m3 = [[1.0, 2.0, 5.0], [4.0, 20.0, 6.0], [7.0, 1.0, 9.0]]
# 3x3 非奇异矩阵
m4 = [[10.0, 2.0, 4.0], [6.0, 5.0, 2.0], [7.0, 15.0, 9.0]]
# 3x3 奇异矩阵
m5 = [[10.0, 2.0, 4.0], [10.0, 2.0, 4.0], [7.0, 15.0, 9.0]]
# 3x4 矩阵
m6 = [[10.0, 2.0, 4.0, 8.0], [6.0, 5.0, 2.0, 1.0], [7.0, 15.0, 9.0, 6.0]]
# 4x3 矩阵
m7 = [[10.0, 6.0, 7.0], [2.0, 5.0, 15.0], [4.0, 2.0, 9.0], [8.0, 1.0, 6.0]] 

matrix_option = {'1': (m1, b1),\
                 '2': (m2, b2),\
                 '3': (m3, []),\
                 '4': (m4, []),\
                 '5': (m5, []),\
                 '6': (m6, []),\
                 '7': (m7, []),}

solve_option = {'1': (m1, b1),\
                '2': (m2, b2)}