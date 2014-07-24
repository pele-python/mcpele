#include <cmath>

#include "../../source/mcpele/progress.h"

int main()
{
    const unsigned int N(100000000);
    mcpele::progress stat(N);
    unsigned int niter(0);
    double tmp(N);
    while (niter < N) {
        tmp = sqrt(N);
        tmp *= tmp;
        tmp = pow(tmp,4.2);
        tmp = pow(tmp,1.0/4.2);
        ++niter;
        stat.next(niter);
        //std::cout << "\r";
        //std::cout << niter;
    }
    std::cout << tmp << "\n";
    return 0;   
}

