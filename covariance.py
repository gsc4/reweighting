import numpy as np
import pymbar

if __name__ == "__main__":

    def covariance(u_kn, N_k, A_kn, B_kn)
        
        mbar = pymbar.MBAR(u_kn, N_k)
        
        A_avg = mbar.computeExpectations(A_kn)
        B_avg = mbar.computeExpectations(B_kn)
        AB_avg = mbar.computeExpectations(A_kn*B_kn)

        cov_AB = AB_avg - A_avg*B_avg

        return cov_AB

