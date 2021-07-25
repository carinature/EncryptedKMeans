

#ifndef ENCKMEAN_CLIENT_H
#define ENCKMEAN_CLIENT_H


#include "utils/aux.h"


class Client {
    const helib::SecKey encryptionKey;
    const helib::EncryptedArray ea;

public:
    Client(KeysServer & keysServer, const long *coordinates);

protected:
    NTL::Vec<helib::Ctxt> cCoordinates{};
    std::vector<helib::Ctxt> cCoordinatesStd{};

private:
#if DBG
    long *pCoordinates; // for dbg
#endif
};


#endif //ENCKMEAN_CLIENT_H
