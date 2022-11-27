import numpy as np
from numpy.linalg import norm
from numpy.linalg import matrix_rank


def check_unique_solution(A, b):
    '''
    Check if the Ax = b has unique solution

    Parameters:
    A, b

    Returns
    True, False
    '''
    Ab = np.c_[A, b.T]
    rank_A = matrix_rank(A)
    rank_Ab = matrix_rank(Ab)
    if rank_A == rank_Ab:
        return True
    else:
        return False

    
def Householder_Reduction(matrix, print_result=True):
    '''
    Housholder Reduction

    Parameters:
    matrix: the input matrix

    Returns:
    P, T
    '''
    np.set_printoptions(precision=3, suppress=True)
    matrix = np.array(matrix)

    row_len = matrix.shape[0]
    col_len = matrix.shape[1]

    # Cleck the input matrix is a square matrix
    if row_len == col_len:
        # The matrix is a square matrix
        P = np.identity(row_len)
        T = np.copy(matrix)

        for index in range(row_len - 1):
            cur_col = T[index:, index]
            e = np.zeros_like(cur_col)
            e[0] = 1

            norm_factor = norm(cur_col)
            u = cur_col - norm_factor * e

            R = np.identity(row_len)
            R[index:, index:] -= (2.0 * 1 / (u.T.dot(u))) * np.dot(u[:, None], u.T[None, :])

            T = np.dot(R, T)
            P = np.dot(R, P)

    else:
        # The matrix is not a square matrix
        P = np.identity(row_len)
        T = np.copy(matrix)
        
        num = min(row_len, col_len) - 1
        if row_len > col_len :
            num += 1

        for i in range(num) :
            cur_col = np.reshape(T[i:, i], (row_len - i, 1))
            e = np.zeros_like(cur_col)
            e[0,] = 1
            
            norm_factor = norm(cur_col)
            u = cur_col - norm_factor * e

            R = np.identity(row_len)

            R[i:, i:] -= 2.0 * 1 * np.dot(u, u.T) / np.dot(u.T, u)
            T = np.dot(R, T)
            P = np.dot(R, P)

    if print_result == True:
        print("The input matrix is: ")
        print(matrix)
        print("The result of the Housholder Reduction is shown as follows: ")
        print("-----The P matrix is: ")
        print(P)
        print("-----The T matrix is: ")
        print(T)

    return (P, T)


def Householder_solve(A, b):
    '''
    Use Housholder Reduction to solve Ax = b

    Parameters:
    A, b

    Returns:
    x
    '''
    np.set_printoptions(precision=3, suppress=True)
    A = np.array(A, dtype=np.float32)
    b = np.array(b, dtype=np.float32)
    print("The input matrix A is: ")
    print(A)
    print("The input matrix b is: ")
    print(b)
    if check_unique_solution(A, b) == False:
        print("Error! The Ax = b does not have unique solution!")
    else:
        P, T = Householder_Reduction(A)
        n = len(A)
        # Solve Tx = Pb
        Pb = np.dot(P, b.T)

        x = np.zeros((n, 1))
        for i in range(len(A)-1, -1, -1):
            t = 0
            for j in range(i+1,len(A)):
                t += T[i][j] * x[j][0]
            t = Pb[i] - t
            if t != 0 and T[i][i] == 0:
                return 0
            x[i] = t / T[i][i]

        print("The solution of Tx = Pb is:")
        print(x)

        return x

        
if __name__ == "__main__":
    matrix = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
    Householder_Reduction(matrix)
    print("----------------------------------------------")
    matrix = [[10.0, 2.0, 4.0], [6.0, 5.0, 2.0], [7.0, 15.0, 9.0]]
    Householder_Reduction(matrix)
    print("----------------------------------------------")
    matrix = [[10.0, 2.0, 4.0], [10.0, 2.0, 4.0], [7.0, 15.0, 9.0]]
    Householder_Reduction(matrix)
    print("----------------------------------------------")
    matrix = [[10.0, 2.0, 4.0, 8.0], [6.0, 5.0, 2.0, 1.0], [7.0, 15.0, 9.0, 6.0]]
    Householder_Reduction(matrix)
    print("----------------------------------------------")
    matrix = [[10.0, 6.0, 7.0], [2.0, 5.0, 15.0], [4.0, 2.0, 9.0], [8.0, 1.0, 6.0], [10.0, 2.0, 4.0]]
    Householder_Reduction(matrix)
    print("----------------------------------------------")
    A = [[1.0, 2.0, 3.0], [2.0, 2.0, 8.0], [-3.0, -10.0, -2.0]]
    b = [0.0, -4.0, -11.0]  
    Householder_solve(A, b)
    print("----------------------------------------------")
    A = [[2.0, -3.0, 1.0], [1.0, -1.0, -2.0], [3.0, 1.0, -1.0]]
    b = [7.0, -2.0, 0.0]
    Householder_solve(A, b)


