import numpy as np
from HousholderReduction import Householder_Reduction
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

    
def URV_Factorization(matrix):
    '''
    URV Factorization

    Parameters:
    matrix: the input matrix

    Returns:
    '''
    np.set_printoptions(precision=3, suppress=True)
    matrix = np.array(matrix
    )

    P, B = Householder_Reduction(matrix.tolist(), print_result = False)
    Q, T = Householder_Reduction(B.T.tolist(), print_result = False)

    U = P.T
    V = Q

    R = np.zeros_like(matrix, dtype=float)
    n = matrix_rank(T)
    R[ :(n+1), :(n+1)] = T.T

    print("The input matrix is: ")
    print(matrix)
    print("The result of the URV Factorization is shown as follows: ")
    print("-----The U matrix is: ")
    print(U)
    print("-----The R matrix is: ")
    print(R)
    print("-----The V matrix is: ")
    print(V)
    return U, R, V


def URV_solve(A, b):
    '''
    Use URV Factorization to solve Ax = b

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
        U, R, V = URV_Factorization(A)
        n = len(A)
        # Solve the Uy = b
        y = np.dot(U.T, b.T)

        # Solve RVx = y
        x = np.dot(V.T, np.linalg.inv(R))
        x = np.dot(x, y)

        print("The solution of Ax = b is:")
        print(x)

        return x


if __name__ == "__main__":
    matrix = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
    URV_Factorization(matrix)
    print("----------------------------------------------")
    matrix = [[10.0, 2.0, 4.0], [10.0, 2.0, 4.0], [7.0, 15.0, 9.0]]
    URV_Factorization(matrix)
    print("----------------------------------------------")
    A = [[1.0, 2.0, 3.0], [2.0, 2.0, 8.0], [-3.0, -10.0, -2.0]]
    b = [0.0, -4.0, -11.0]  
    URV_solve(A, b)
    print("----------------------------------------------")
    A = [[2.0, -3.0, 1.0], [1.0, -1.0, -2.0], [3.0, 1.0, -1.0]]
    b = [7.0, -2.0, 0.0]
    URV_solve(A, b)
