

#ifndef ENCKMEAN_DATASERVER_H
#define ENCKMEAN_DATASERVER_H


#include "Client.h"

class DataServer {

//    [[maybe_unused]] const helib::SecKey encryptionKey;
//    [[maybe_unused]] const helib::PubKey encryptionKey;
//    [[maybe_unused]] const helib::EncryptedArray ea;
     const helib::PubKey encryptionKey;
     const helib::EncryptedArray ea;
    const helib::Ctxt scratch;
protected:
//    const helib::PubKey &public_key;// = encryptionKey;

public:
    /**
     * Constructor for \class{Client},
     * @param keysServer binds to the \class{KeysServer} responsible for the distributing the appropriate key
     * @brief simulates a mini-protocol between the keys server and the Data Server.
     * ks and ds will create a scratch that will allow to encrypt the results:
     * result bit in case of compare, and result num in case of add/multiplication.
     * */
    explicit DataServer(KeysServer &keysServer);

    //void BGV_binary_arithmetic();
    //void BGV_packed_arithmetic();
    //void KT_packed_arithmetic();
    //void add_numbers();
    //void mult_numbers();
    //void cmp_numbers();

    //    helib::Ctxt compare(Client &client1, Client &client2) {
    //#if VERBOSE
    //        cout << "compare(Client &client) " << endl;
    //#endif
    //        Ctxt &encryptionKey;
    //        helib::Ctxt isBiggerFlag(encryptionKey), isSmallerFlag(encryptionKey);
    //        cout << "odin" << endl;
    //        compareTwoNumbers(isBiggerFlag, isSmallerFlag,
    //                          helib::CtPtrs_VecCt(this->cCoordinatesNTL),
    //                          helib::CtPtrs_VecCt(client.cCoordinatesNTL),
    //                          false, nullptr);
    //        cout << "onaes" << endl;
    //        return isBiggerFlag;
    //    }
    Point scratchPoint();
};


#endif //ENCKMEAN_DATASERVER_H
