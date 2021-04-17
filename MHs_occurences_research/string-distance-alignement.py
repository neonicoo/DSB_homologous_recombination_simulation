import numpy as np


def hamming_distance(S,T):
    if len(S) != len(T):
        return None
    d = 0
    for i,s in enumerate(S):
        d += s != T[i]
    return d
    
    
    
def Levenstein_distance(S,T,return_matrix=False):
    L = np.zeros((len(S)+1,len(T)+1))
    L[0,:] = range(len(T)+1)
    L[:,0] = range(len(S)+1)
    for i in range(1,len(S)+1):
        for j in range(1,len(T)+1):
            #print(L)
            comparison = int(S[i-1]!=T[j-1])
            L[i,j] = min([L[i-1,j] + 1, L[i,j-1] + 1, L[i-1,j-1] + comparison])
    if return_matrix:
        return L
    return(L[-1,-1])


def NWmatrix2strings(S,T,D):
    S_WITH_GAP = -1
    T_WITH_GAP = 1
    S_T_ALIGN = 0
    i = D.shape[0]-1
    j = D.shape[1]-1
    S_aligned = ""
    T_aligned = ""
    while i > 0 or j > 0:
        #print(i,j,T[j-1],S[i-1],D[i,j])
        if D[i,j] == S_WITH_GAP:
            S_aligned = S[i-1]+S_aligned
            T_aligned = "_"+T_aligned
            #print(S_aligned,T_aligned)
            i -= 1
        elif D[i,j] == T_WITH_GAP:
            T_aligned = T[j-1]+ T_aligned
            S_aligned = "_"+S_aligned
            #print(S_aligned,T_aligned)
            j -= 1
        elif D[i,j] == S_T_ALIGN:
            S_aligned = S[i-1] + S_aligned
            T_aligned = T[j-1] + T_aligned
            #print(S_aligned,T_aligned)
            j -= 1
            i -= 1
    return(S_aligned,T_aligned)


def Needleman_Wunch(S,T,return_matrix=False):
    T_WITH_GAP = 1
    S_WITH_GAP = -1
    S_T_ALIGN = 0
    decision = [T_WITH_GAP,S_WITH_GAP,S_T_ALIGN]
    F = np.zeros((len(S)+1,len(T)+1))
    F[0,:] = range(len(T)+1)
    F[:,0] = range(len(S)+1)
    F *= -1
    D = np.zeros((len(S)+1,len(T)+1))
    D[0,1:] = T_WITH_GAP 
    D[1:,0] = S_WITH_GAP 
    for i in range(1,len(S)+1):
        for j in range(1,len(T)+1):
            comparison = int(S[i-1]==T[j-1]) # * 10
            options = [F[i,j-1] - 1, F[i-1,j] - 1, F[i-1,j-1] + comparison]
            F[i,j] = options[0]
            D[i,j] = T_WITH_GAP
            for o,val in enumerate(options):
                if val > F[i,j]:
                    F[i,j] = val
                    D[i,j] = decision[o]
    if return_matrix:
        return F,D
        

def SWmatrix2strings(S,T,D,F):
    S_WITH_GAP = -1
    T_WITH_GAP = 1
    S_T_ALIGN = 0
    i,j = np.unravel_index(F.argmax(),F.shape)
    print(i,j)
    S_aligned = ""
    T_aligned = ""
    while i > 0 and j > 0:
        #print(i,j,T[j-1],S[i-1],D[i,j])
        if D[i,j] == S_WITH_GAP:
            S_aligned = S[i-1]+S_aligned
            T_aligned = "_"+T_aligned
            #print(S_aligned,T_aligned)
            i -= 1
        elif D[i,j] == T_WITH_GAP:
            T_aligned = T[j-1]+ T_aligned
            S_aligned = "_"+S_aligned
            #print(S_aligned,T_aligned)
            j -= 1
        elif D[i,j] == S_T_ALIGN:
            S_aligned = S[i-1] + S_aligned
            T_aligned = T[j-1] + T_aligned
            #print(S_aligned,T_aligned)
            j -= 1
            i -= 1
    return(S_aligned,T_aligned)

def Smith_Waterman(S,T,return_matrix=False):
    T_WITH_GAP = 1
    S_WITH_GAP = -1
    S_T_ALIGN = 0
    decision = [T_WITH_GAP,S_WITH_GAP,S_T_ALIGN]
    F = np.zeros((len(S)+1,len(T)+1))
    D = np.zeros((len(S)+1,len(T)+1))
    D[0,1:] = T_WITH_GAP 
    D[1:,0] = S_WITH_GAP 
    for i in range(1,len(S)+1):
        for j in range(1,len(T)+1):
            comparison = int(S[i-1]==T[j-1])
            options = [F[i,j-1] - 1, F[i-1,j] - 1, F[i-1,j-1] + comparison]
            F[i,j] = options[0]
            D[i,j] = T_WITH_GAP
            for o,val in enumerate(options):
                if val > F[i,j]:
                    F[i,j] = val
                    D[i,j] = decision[o]
    if return_matrix:
        return F,D
        
############################################################
############################################################
if __name__ == "__main__":   
	pass

