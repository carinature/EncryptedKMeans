
#include <properties.h>
#include "Point.h"
//#include "utils/aux.h" // for including KeysServer.h

using std::cout;
using std::endl;
Point::Point(const helib::PubKey &public_key, const long coordinates[])  :
        public_key(public_key),
        cCoordinates(DIM, std::vector(bitSize, helib::Ctxt(public_key))) {
    cout << " Point Init" << endl;
    if (coordinates)
        for (int dim = 0; dim < DIM; ++dim)
            for (long bit = 0; bit < bitSize; ++bit) // Extract the i'th bit of coordinates[dim]
                this->public_key.Encrypt(cCoordinates[dim][bit],
                                         NTL::to_ZZX((coordinates[dim] >> bit) & 1));
}