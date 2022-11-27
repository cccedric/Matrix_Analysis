import numpy as np
from numpy.linalg import matrix_rank

def check(matrix):
    '''
    Check if the number of the columns is equal to the rank

    Parameters:
    matrix: the input matrix

    Returns
    True, False
    '''
    rank = matrix_rank(matrix)
    n = matrix.shape[-1]
    if rank == n:
        return True
    else:
        return False
    

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

    
def row_trans(matrix,i,j):
    '''
    Elementary tow transmation
    Parameters:
    matrix: the original matrix
    i, j: the index of the 2 interchange rows 

    Return:
    matrix: the final matrix
    '''
    temp = np.copy(matrix[i])
    matrix[i] = matrix[j]
    matrix[j] = temp

    return matrix


def LU_Factorization(matrix):
    '''
    LU Factorization
    
    Parameters:
    matrix: the input matrix

    Returns:
    L, U, P
    '''
    np.set_printoptions(precision=3, suppress=True)
    matrix = np.array(matrix, dtype=np.float32)
    print("The input matrix is: ")
    print(matrix)

    row_len = matrix.shape[0]
    col_len = matrix.shape[1]

    # Cleck the input matrix is a square invertible matrix, or otherwise report it and end the function
    if row_len != col_len:
        print("ERROR: The input matrix is not a square matrix.")
        return
    elif check(matrix) == False:
        print("ERROR: The matrix is not invertible.")
        return
    else:
        # Expand it to a (matrix|b) matrix, the b illustrates the order of the rows
        row_order = []
        for i in range(row_len):
            row_order.append(i)
        
        row_order = np.expand_dims(np.array(row_order), axis=0)
        matrix = np.c_[matrix,row_order.T]
        col_len += 1

        # Implement the Gaussion Elimination
        L_matrix = np.zeros(matrix.shape)
        col_index = 0
        for row_index in range(row_len - 1):
            pivot_row = row_index + np.argmax(np.abs(matrix[row_index:row_len, col_index]))
            while(matrix[pivot_row, col_index] == 0):
                if(col_index <= row_len):
                    col_index += 1
                else:
                    break

            if(pivot_row != row_index):
                matrix = row_trans(matrix, row_index, pivot_row)
                L_matrix = row_trans(L_matrix, row_index, pivot_row)

            for j in range(row_index+1, row_len):
                if(matrix[j, col_index] != 0):
                    L_matrix[j, col_index] = 1.0 * matrix[j, col_index] / matrix[row_index, col_index]
                    matrix[j, :row_len] = matrix[j, :row_len] - 1.0 * matrix[row_index, :row_len] * matrix[j, col_index] / matrix[row_index, col_index]
            
            col_index += 1

        # If the matrix is invertible, find the L, U, P
        U = matrix[:row_len, :row_len]

        L = L_matrix[:row_len, :row_len]
        for i in range(row_len):
            L[i, i]=1

        P = np.zeros((row_len, row_len))
        for i in range(row_len):
            P[i][int(matrix[i, -1])] = 1

        print("The result of the LU factorization is shown as follows: ")
        print("-----The L matrix is: ")
        print(L)
        print("-----The U matrix is: ")
        print(U)
        print("-----The P matrix is: ")
        print(P)
        return (L, U, P)

    
def LU_solve(A, b):
    '''
    Use LU Factorization to solve Ax = b
    
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
        L, U, P = LU_Factorization(A) 
        n = len(A)

        # Solve Ly = Pb
        y = np.zeros((n, 1))
        Pb = np.dot(P, b.T)

        Pb = np.array(Pb).reshape(n,1) # 把b列表格式变成向量格式
        for i in range(len(A)):
            t = 0
            for j in range(i):
                t += L[i][j] * y[j][0]
            y[i][0] = Pb[i][0] - t

        print("The solution of Ly = Pb is:")
        print(y)

        # Solve Ux = y
        x = np.zeros((n, 1))
        for i in range(len(A)-1, -1, -1):
            t = 0
            for j in range(i+1, len(A)):
                t += U[i][j] * x[j][0]
            t = y[i][0] - t
            if t != 0 and U[i][i] == 0:
                return 0
            x[i] = t / U[i][i]
            
        print("The solution of Ax = b is:")
        print(x)

        return x


if __name__ == "__main__":
    matrix = [[1.0, 2.0, 5.0], [4.0, 20.0, 6.0], [7.0, 1.0, 9.0]]
    LU_Factorization(matrix)
    print("----------------------------------------------")
    matrix = [[10.0, 2.0, 4.0], [6.0, 5.0, 2.0], [7.0, 15.0, 9.0]]
    LU_Factorization(matrix)
    print("----------------------------------------------")
    matrix = [[10.0, 2.0, 4.0, 8.0], [6.0, 5.0, 2.0, 1.0], [7.0, 15.0, 9.0, 6.0]]
    LU_Factorization(matrix)
    print("----------------------------------------------")
    matrix = [[10.0, 2.0, 4.0], [10.0, 2.0, 4.0], [7.0, 15.0, 9.0]]
    LU_Factorization(matrix)
    print("----------------------------------------------")
    A = [[1.0, 2.0, 3.0], [2.0, 2.0, 8.0], [-3.0, -10.0, -2.0]]
    b = [0.0, -4.0, -11.0]  
    LU_solve(A, b)
    print("----------------------------------------------")
    A = [[2.0, -3.0, 1.0], [1.0, -1.0, -2.0], [3.0, 1.0, -1.0]]
    b = [7.0, -2.0, 0.0]
    LU_solve(A, b)

