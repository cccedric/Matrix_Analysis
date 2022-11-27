import numpy as np
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

    
def Givens_Reduction(matrix):
    '''
    Givens Reduction

    Parameters:
    matrix: the input matrix

    Returns:
    P, T
    '''
    np.set_printoptions(precision=3, suppress=True)
    matrix = np.array(matrix)
    print("The input matrix is: ")
    print(matrix)

    row_len = np.shape(matrix)[0]
    col_len = np.shape(matrix)[1]

    P = np.identity(row_len)
    T = np.copy(matrix)

    for i in range(col_len):
        for j in range(row_len - 1, i, -1):
            x1 = T[i, i]
            x2 = T[j, i]
            norm_factor = np.sqrt((x1 ** 2 + x2 ** 2))
            
            c = x1 / norm_factor
            s = x2 / norm_factor

            Pi = np.eye(row_len)
            Pi[i, i] = c
            Pi[j, j] = c
            Pi[i, j] = s
            Pi[j, i] = -s
            
            P = np.dot(Pi,P)
            T = np.dot(Pi,T)

    print("The result of the Givens Reduction is shown as follows: ")
    print("-----The P matrix is: ")
    print(P)
    print("-----The T matrix is: ")
    print(T)

    return (P, T)


def Givens_solve(A, b):
    '''
    Use Givens Reduction to solve Ax = b

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
        P, T = Givens_Reduction(A)
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
    Givens_Reduction(matrix)
    print("----------------------------------------------")
    matrix = [[10.0, 2.0, 4.0], [6.0, 5.0, 2.0], [7.0, 15.0, 9.0]]
    Givens_Reduction(matrix)
    print("----------------------------------------------")
    matrix = [[10.0, 2.0, 4.0], [10.0, 2.0, 4.0], [7.0, 15.0, 9.0]]
    Givens_Reduction(matrix)
    print("----------------------------------------------")
    matrix = [[10.0, 2.0, 4.0, 8.0], [6.0, 5.0, 2.0, 1.0], [7.0, 15.0, 9.0, 6.0]]
    Givens_Reduction(matrix)
    print("----------------------------------------------")
    matrix = [[10.0, 6.0, 7.0], [2.0, 5.0, 15.0], [4.0, 2.0, 9.0], [8.0, 1.0, 6.0]]
    Givens_Reduction(matrix)
    print("----------------------------------------------")
    A = [[1.0, 2.0, 3.0], [2.0, 2.0, 8.0], [-3.0, -10.0, -2.0]]
    b = [0.0, -4.0, -11.0]  
    Givens_solve(A, b)
    print("----------------------------------------------")
    A = [[2.0, -3.0, 1.0], [1.0, -1.0, -2.0], [3.0, 1.0, -1.0]]
    b = [7.0, -2.0, 0.0]
    Givens_solve(A, b)