#include <iostream>
#include <iomanip>
#include <cmath>

using std::cout ;
using std::setw ;
using std::pow ;

const double varepsilon_1 = -0.5782 ;
const double varepsilon_2 = 0.6703 ;
const double j_11 = 0.6746 ;
const double j_12 = 0.6636 ;
const double j_22 = 0.6975 ;
const double k_12 = 0.1813 ;
const double delta = varepsilon_2 - varepsilon_1 + 0.5 * ( j_11 + j_22 ) - 2 * j_12 + k_12 ;

int main(){

    cout << std::fixed << "delta = " << delta << "\n" ;
    cout << setw( 3 ) << "  N\t" << "energy_dci\t"  << "err_1\t\t" ;
    cout << "energy_total\t" << "err_2" << "\t\ttenergy_exact\n" ;
    for( unsigned N = 1 ; N <= 100 ; ++N ){

        double energy_dci = delta - std::sqrt( std::pow( delta , 2 ) + N * pow( k_12 , 2 ) ) ;        // from (4.67)
        double energy_exact = N * ( delta - std::sqrt( std::pow( delta , 2 ) + pow( k_12 , 2 ) ) ) ;  // from (4.68)
        double energy_davidson = pow( energy_dci , 3 ) / ( N * pow( k_12 , 2 ) + pow( energy_dci , 2 ) ) ;
        double err_1 = ( energy_dci - energy_exact ) / energy_exact ;
        double err_2 = ( energy_dci + energy_davidson - energy_exact ) / energy_exact ;

        std::setprecision(7) ;
        cout << setw(3) << N << "\t" << energy_dci << "\t" << 100 * err_1 << "%\t" ;
        cout << ( energy_dci + energy_davidson ) << "\t" << 100 * err_2 << "%\t" ;
        cout << energy_exact << '\n' ;
        if( N == 20 )
            N += 79 ;
    }

    return 0 ;
}