
#include "utils/aux.h"

#include "TestAux.h"
#include "src/Client.h"

void TestAux::testGenerateDataClients() {
    cout << " ------ testGenerateDataClients ------ " << endl << endl;
    KeysServer keysServer;
    auto t0_main = CLOCK::now();

    const std::vector<Client> clients = generateDataClients(keysServer);

    for (const Client &client:clients)
        for (const Point &point:client.getPoints())
            printPoint(point, keysServer);

    cout << endl << printDuration(t0_main, "testGenerateDataClients");
    cout << " ------ testGenerateDataClients finished ------ " << endl << endl;

}

void TestAux::minitest() {
    cout << " ------ minitest ------ " << endl;
    KeysServer keysServer(5);
    int n = NUMBER_OF_POINTS;

    helib::PubKey &pubKey = keysServer.getPublicKey();
    helib::Ctxt ctxt(pubKey), ctxt2(pubKey), ctxt3(pubKey), ctxt4(pubKey);
    long arr[] = {0, 1, 2, 3};
    pubKey.Encrypt(ctxt3, NTL::ZZX(3));
    for (int i = 0; i < n; ++i) {
        cout << "================";
        printNameVal(i);
        //        cout << "================"<<endl;
        pubKey.Encrypt(ctxt, NTL::ZZX(i));
        pubKey.Encrypt(ctxt2, NTL::ZZX(i));
        printNameVal(keysServer.decryptCtxt(ctxt));
        printNameVal(keysServer.decryptCtxt(ctxt2));
        printNameVal((ctxt) == (ctxt2));
        printNameVal(ctxt.equalsTo(ctxt2));
        printNameVal(&(ctxt) == &(ctxt2));
        cout << "-------------" << endl;
        ctxt = ctxt3;
        ctxt2 = ctxt3;
        printNameVal((ctxt) ==
                     (ctxt2));      //  true when both initialized to the same ctxt, false when initialize to the same (p)
        printNameVal(ctxt.equalsTo(
                ctxt2)); //  true when both initialized to the same ctxt, false when initialize to the same (p)value
        printNameVal(&(ctxt) == &(ctxt2));    //  never true

        cout << "-------1------" << endl;
        ctxt *= ctxt3;
        ctxt2 *= ctxt3;
        printNameVal((ctxt) ==
                     (ctxt2));      //  true when both initialized to the same ctxt, false when initialize to the same (p)
        printNameVal(ctxt.equalsTo(
                ctxt2)); //  true when both initialized to the same ctxt, false when initialize to the same (p)value
        printNameVal(&(ctxt) == &(ctxt2));    //  never true

        cout << "-------2------" << endl;
        ctxt *= ctxt4;
        ctxt2 *= ctxt4;
        printNameVal((ctxt) ==
                     (ctxt2));      //  true when both initialized to the same ctxt, false when initialize to the same (p)
        printNameVal(ctxt.equalsTo(
                ctxt2)); //  true when both initialized to the same ctxt, false when initialize to the same (p)value
        printNameVal(&(ctxt) == &(ctxt2));    //  never true

        cout << "-------3------" << endl;
        ctxt = ctxt3;
        ctxt2 = ctxt4;
        ctxt *= ctxt4;
        ctxt2 *= ctxt3;
        printNameVal((ctxt) ==
                     (ctxt2));      //  true when both initialized to the same ctxt, false when initialize to the same (p)
        printNameVal(ctxt.equalsTo(
                ctxt2)); //  true when both initialized to the same ctxt, false when initialize to the same (p)value
        printNameVal(&(ctxt) == &(ctxt2));    //  never true
        cout << "================" << endl;

    }
    pubKey.Encrypt(ctxt, NTL::ZZX(3));
    pubKey.Encrypt(ctxt2, NTL::ZZX(4));
    pubKey.Encrypt(ctxt3, NTL::ZZX(3));
    pubKey.Encrypt(ctxt4, NTL::ZZX(4));
    ctxt *= ctxt4;
    ctxt2 *= ctxt3;
    printNameVal((ctxt) ==
                 (ctxt2));      //  true when both initialized to the same ctxt, false when initialize to the same (p)
    printNameVal(ctxt.equalsTo(
            ctxt2)); //  true when both initialized to the same ctxt, false when initialize to the same (p)value
    printNameVal(&(ctxt) == &(ctxt2));    //  never true

    cout << "flalallsdf" << endl;
    // Create a vector of long with nslots elements
    helib::Ptxt<helib::BGV> ptxt(pubKey);
    // Set it with numbers 0..nslots - 1
    // ptxt = [0] [1] [2] ... [nslots-2] [nslots-1]
    for (int i = 0; i < ptxt.size(); ++i) {
        ptxt[i] = i;
        printNameVal(i);
        printNameVal(ptxt[i]);
        cout << ptxt << endl;
    }
    pubKey.Encrypt(ctxt, ptxt);
    keysServer.decryptCtxt(ctxt);


    cout << " ------ minitest finished ------ " << endl << endl;

}

struct Bla {
    Bla *bla;

    Bla() {
        bla = this;
    }

    void printBla() {
        printNameVal(bla);
        printNameVal(this);
    }
};

void TestAux::minitest2() {
    Bla bla;
    bla.printBla();
    printNameVal(bla.bla);
    printNameVal(&bla);
    cout << "fin" << endl;
}

void TestAux::testBGVPackedArithmetics_Original() {
    /*  Example of BGV scheme  */

    // Plaintext prime modulus
    unsigned long p = PLAINTEXT_PRIME_MODULUS;
    //    unsigned long p = 4999;
    // Cyclotomic polynomial - defines phi(m)
    unsigned long m = 32109;
    // Hensel lifting (default = 1)
    unsigned long r = 1;
    // Number of bits of the modulus chain
    unsigned long bits = 500;
    // Number of columns of Key-Switching matrix (default = 2 or 3)
    unsigned long c = 2;

    cout << "\n*********************************************************";
    cout << "\n*         Basic Mathematical Operations Example         *";
    cout << "\n*         =====================================         *";
    cout << "\n*                                                       *";
    cout << "\n* This is a sample program for education purposes only. *";
    cout << "\n* It attempts to show the various basic mathematical    *";
    cout << "\n* operations that can be performed on both ciphertexts  *";
    cout << "\n* and plaintexts.                                       *";
    cout << "\n*                                                       *";
    cout << "\n*********************************************************";
    cout << endl;

    cout << "Initialising context object..." << endl;
    // Initialize context
    // This object will hold information about the algebra created from the
    // previously set parameters
    helib::Context context = helib::ContextBuilder<helib::BGV>()
            .m(m)
            .p(p)
            .r(r)
            .bits(bits)
            .c(c)
            .build();

    // Print the context
    context.printout();
    cout << endl;

    // Print the security level
    cout << "Security: " << context.securityLevel() << endl;

    // Secret key management
    cout << "Creating secret key..." << endl;
    // Create a secret key associated with the context
    helib::SecKey secret_key(context);
    // Generate the secret key
    secret_key.GenSecKey();
    cout << "Generating key-switching matrices..." << endl;
    // Compute key-switching matrices that we need
    helib::addSome1DMatrices(secret_key);

    // Public key management
    // Set the secret key (upcast: SecKey is a subclass of PubKey)
    const helib::PubKey &public_key = secret_key;

    // Get the EncryptedArray of the context
    const helib::EncryptedArray &ea = context.getEA();

    // Get the number of slot (phi(m))
    long nslots = ea.size();
    cout << "Number of slots: " << nslots << endl;

    // Create a vector of long with nslots elements
    helib::Ptxt<helib::BGV> ptxt(context);
    // Set it with numbers 0..nslots - 1
    // ptxt = [0] [1] [2] ... [nslots-2] [nslots-1]
    for (int i = 0; i < ptxt.size(); ++i) {
        ptxt[i] = i;
    }

    // Print the plaintext
    cout << "Initial Plaintext: " << ptxt << endl;

    // Create a ciphertext object
    helib::Ctxt ctxt(public_key);
    // Encrypt the plaintext using the public_key
    public_key.Encrypt(ctxt, ptxt);

    /********** Operations **********/
    // Ciphertext and plaintext operations are performed
    // "entry-wise".

    // Square the ciphertext
    // [0] [1] [2] [3] [4] ... [nslots-1]
    // -> [0] [1] [4] [9] [16] ... [(nslots-1)*(nslots-1)]
    ctxt.multiplyBy(ctxt);
    // Plaintext version
    ptxt.multiplyBy(ptxt);

    // Divide the ciphertext by itself
    // To do this we must calculate the multiplicative inverse using Fermat's
    // Little Theorem.  We calculate a^{-1} = a^{p-2} mod p, where a is non-zero
    // and p is our plaintext prime.
    // First make a copy of the ctxt using copy constructor
    helib::Ctxt ctxt_divisor(ctxt);
    // Raise the copy to the exponent p-2
    // [0] [1] [4] ... [16] -> [0] [1] [1] ... [1]
    // Note: 0 is a special case because 0^n = 0 for any power n
    ctxt_divisor.power(p - 2);
    // a^{p-2}*a = a^{-1}*a = a / a = 1;
    ctxt.multiplyBy(ctxt_divisor);

    // Plaintext version
    helib::Ptxt<helib::BGV> ptxt_divisor(ptxt);
    ptxt_divisor.power(p - 2);
    ptxt.multiplyBy(ptxt_divisor);

    // Double it (using additions)
    // [0] [1] [1] ... [1] [1] -> [0] [2] [2] ... [2] [2]
    ctxt += ctxt;
    // Plaintext version
    ptxt += ptxt;

    // Subtract it from itself (result should be 0)
    // i.e. [0] [0] [0] [0] ... [0] [0]
    ctxt -= ctxt;
    // Plaintext version
    ptxt -= ptxt;

    // Create a plaintext for decryption
    helib::Ptxt<helib::BGV> plaintext_result(context);
    // Decrypt the modified ciphertext
    secret_key.Decrypt(plaintext_result, ctxt);

    cout << "Operation: 2(a*a)/(a*a) - 2(a*a)/(a*a) = 0" << endl;
    // Print the decrypted plaintext
    // Should be [0] [0] [0] ... [0] [0]
    cout << "Decrypted Result: " << plaintext_result << endl;
    // Print the plaintext version result, should be the same as the ctxt version
    cout << "Plaintext Result: " << ptxt << endl;

    // We can also add constants
    // [0] [0] [0] ... [0] [0] -> [1] [1] [1] ... [1] [1]
    ctxt.addConstant(NTL::ZZX(1l));
    // Plaintext version
    ptxt.addConstant(NTL::ZZX(1l));

    // And multiply by constants
    // [1] [1] [1] ... [1] [1]
    // -> [1*1] [1*1] [1*1] ... [1*1] [1*1] = [1] [1] [1] ... [1] [1]
    ctxt *= 1l;
    // Plaintext version
    ptxt *= 1l;

    // We can also perform ciphertext-plaintext operations
    // ctxt = [1] [1] [1] ... [1] [1], ptxt = [1] [1] [1] ... [1] [1]
    // ctxt + ptxt = [2] [2] [2] ... [2] [2]
    // Note: the output of this is also a ciphertext
    ctxt += ptxt;

    // Decrypt the modified ciphertext into a new plaintext
    helib::Ptxt<helib::BGV> new_plaintext_result(context);
    secret_key.Decrypt(new_plaintext_result, ctxt);

    cout << "Operation: Enc{(0 + 1)*1} + (0 + 1)*1" << endl;
    // Print the decrypted plaintext
    // Should be [2] [2] [2] ... [2] [2]
    cout << "Decrypted Result: " << new_plaintext_result << endl;

}

class PointMini {
public:
    std::vector<helib::Ctxt> coors;
    const helib::PubKey &pubKey;
    std::vector<long> pCoors;
    std::vector<helib::Ptxt<helib::BGV> > ptxtCoors;

    PointMini(const helib::PubKey &pubKey, const std::vector<long> &pCoors) :
            pubKey(pubKey),
            pCoors(pCoors) {
        cout << "PointMini" << endl;
        for (int dim = 0; dim < DIM; ++dim) {
            //            ptxtCoors[dim]
            helib::Ptxt<helib::BGV> tempPtxt(pubKey);
            //            tempPtxt.convertToSlot()
        }
    }
};

void TestAux::testBGVPackedArithmetics__Comparison() {
    /*  Example of BGV scheme  */

    // Plaintext prime modulus
    unsigned long p = 17;
    //    unsigned long p = 4999;
    // Cyclotomic polynomial - defines phi(m)
    unsigned long m = 105;
    // Hensel lifting (default = 1)
    unsigned long r = 1;
    // Number of bits of the modulus chain
    unsigned long bits = 500;
    // Number of columns of Key-Switching matrix (default = 2 or 3)
    unsigned long c = 2;
    //
    //    // Plaintext prime modulus
    //    unsigned long p = PLAINTEXT_PRIME_MODULUS;
    //    //    unsigned long p = 4999;
    //    // Cyclotomic polynomial - defines phi(m)
    //    unsigned long m = 32109;
    //    // Hensel lifting (default = 1)
    //    unsigned long r = 1;
    //    // Number of bits of the modulus chain
    //    unsigned long bits = 500;
    //    // Number of columns of Key-Switching matrix (default = 2 or 3)
    //    unsigned long c = 2;

    cout << "\n*********************************************************";
    cout << "\n*         Basic Mathematical Operations Example         *";
    cout << "\n*         =====================================         *";
    cout << "\n*                                                       *";
    cout << "\n* This is a sample program for education purposes only. *";
    cout << "\n* It attempts to show the various basic mathematical    *";
    cout << "\n* operations that can be performed on both ciphertexts  *";
    cout << "\n* and plaintexts.                                       *";
    cout << "\n*                                                       *";
    cout << "\n*********************************************************";
    cout << endl;

    cout << "Initialising context object..." << endl;
    // Initialize context
    // This object will hold information about the algebra created from the
    // previously set parameters
    helib::Context context = helib::ContextBuilder<helib::BGV>()
            .m(m)
            .p(p)
            .r(r)
            .bits(bits)
            .c(c)
            .build();

    // Print the context
    context.printout();
    cout << endl;

    // Print the security level
    cout << "Security: " << context.securityLevel() << endl;

    // Secret key management
    cout << "Creating secret key..." << endl;
    // Create a secret key associated with the context
    helib::SecKey secret_key(context);
    // Generate the secret key
    secret_key.GenSecKey();
    cout << "Generating key-switching matrices..." << endl;
    // Compute key-switching matrices that we need
    helib::addSome1DMatrices(secret_key);

    // Public key management
    // Set the secret key (upcast: SecKey is a subclass of PubKey)
    const helib::PubKey &public_key = secret_key;

    // Get the EncryptedArray of the context
    const helib::EncryptedArray &ea = context.getEA();

    // Get the number of slot (phi(m))
    long nslots = ea.size();
    cout << "Number of slots: " << nslots << endl;

    // Create a vector of long with nslots elements
    helib::Ptxt<helib::BGV> ptxt(public_key);
    helib::Ptxt<helib::BGV> ptxt0(public_key);
    helib::Ptxt<helib::BGV> ptxt1(public_key);
    helib::Ptxt<helib::BGV> ptxt2(public_key);

    // Print the plaintext
    cout << "Initial Plaintext: " << ptxt << endl;
    cout << "Initial Plaintext: " << ptxt0 << endl;
    cout << "Initial Plaintext: " << ptxt1 << endl;
    cout << "Initial Plaintext: " << ptxt2 << endl;

    printNameVal(ptxt.size());
    printNameVal(ptxt0.size());
    printNameVal(ptxt1.size());
    printNameVal(ptxt2.size());
    // Set it with numbers 0..nslots - 1
    // ptxt = [0] [1] [2] ... [nslots-2] [nslots-1]
    for (int i = 0; i < ptxt.size(); ++i) {
        //                ptxt[i] = 2;
        //                ptxt0[i] = 3;
        //                ptxt1[i] = 5;
        //                ptxt2[i] = 9;
        ptxt[i] = i;
        ptxt0[i] = 0;
        ptxt1[i] = 1;
        ptxt2[i] = 2 * i;
    }
    //    ptxt[0] = 2;
    //    ptxt0[0] = 3;
    //    ptxt1[0] = 5;
    //    ptxt2[0] = 9;

    // Print the plaintext
    cout << "Initial Plaintext: " << ptxt << endl;
    cout << "Initial Plaintext: " << ptxt1 << endl;
    cout << "Initial Plaintext: " << ptxt2 << endl;

    // Create a ciphertext object
    helib::Ctxt ctxt(public_key);
    helib::Ctxt ctxt0(public_key);
    helib::Ctxt ctxt1(public_key);
    helib::Ctxt ctxt2(public_key);
    // Encrypt the plaintext using the public_key
    public_key.Encrypt(ctxt, ptxt);
    public_key.Encrypt(ctxt0, ptxt0);
    public_key.Encrypt(ctxt1, ptxt1);
    public_key.Encrypt(ctxt2, ptxt2);

    // Create a plaintext for decryption
    helib::Ptxt<helib::BGV> plaintext_result(public_key);
    helib::Ptxt<helib::BGV> plaintext_result0(public_key);
    helib::Ptxt<helib::BGV> plaintext_result1(public_key);
    helib::Ptxt<helib::BGV> plaintext_result2(public_key);
    // Decrypt the modified ciphertext
    secret_key.Decrypt(plaintext_result, ctxt);
    secret_key.Decrypt(plaintext_result0, ctxt0);
    secret_key.Decrypt(plaintext_result1, ctxt1);
    secret_key.Decrypt(plaintext_result2, ctxt2);
    cout << "Decrypted ctxt: " << plaintext_result << endl;
    cout << "Decrypted ctxt0: " << plaintext_result0 << endl;
    cout << "Decrypted ctxt1: " << plaintext_result1 << endl;
    cout << "Decrypted ctxt2: " << plaintext_result2 << endl;
    // Print the plaintext version result, should be the same as the ctxt version
    cout << "Plaintext ptxt: " << ptxt << endl;
    cout << "Plaintext ptxt0: " << ptxt0 << endl;
    cout << "Plaintext ptxt1: " << ptxt1 << endl;
    cout << "Plaintext ptxt2: " << ptxt2 << endl;

    /********** Operations **********/
    // Ciphertext and plaintext operations are performed
    // "entry-wise".

    helib::Ctxt tempCtxt(public_key);
    helib::Ptxt<helib::BGV> tempPtxt(public_key);
    cout << "Plaintext Result: " << tempPtxt << endl;
    secret_key.Decrypt(plaintext_result0, tempCtxt);
    cout << "Decrypted Result 0: " << plaintext_result0 << endl;

    //    tempCtxt = ctxt;
    //    secret_key.Decrypt(plaintext_result, tempCtxt);
    //    cout << "Decrypted Result 1: " << plaintext_result << endl;
    //
    //    tempCtxt += ctxt1;
    //    secret_key.Decrypt(plaintext_result, tempCtxt);
    //    cout << "Decrypted Result 2: " << plaintext_result << endl;
    //
    //    tempCtxt += ptxt1;
    //    secret_key.Decrypt(plaintext_result, tempCtxt);
    //    cout << "Decrypted Result 3: " << plaintext_result << endl;
    //
    //    tempCtxt += ptxt1;
    //    secret_key.Decrypt(plaintext_result, tempCtxt);
    //    cout << "Decrypted Result 4: " << plaintext_result << endl;

    helib::Ctxt mu(public_key), ni(public_key);

    std::vector<Ctxt> ctxtVec{ctxt};
    std::vector<Ctxt> ctxtVec0{ctxt0};
    std::vector<Ctxt> ctxtVec1{ctxt1};
    std::vector<Ctxt> ctxtVec2{ctxt2};
    std::vector<Ctxt> tempVec{tempCtxt};
    std::vector<Ctxt> tempVec2{ctxt};
    helib::CtPtrs_vectorCt bla{tempVec};

    cout << "------------------------------------" << endl;
    cout << " ctxt > ctxt " << endl;
    helib::compareTwoNumbers(mu, ni, helib::CtPtrs_vectorCt(ctxtVec),
                             helib::CtPtrs_vectorCt(ctxtVec));
    helib::Ptxt<helib::BGV> dMu(public_key);
    secret_key.Decrypt(dMu, mu);
    cout << "Decrypted dMu: " << dMu << endl;
    helib::Ptxt<helib::BGV> dNi(public_key);
    secret_key.Decrypt(dNi, ni);
    cout << "Decrypted dNi: " << dNi << endl;
    cout << "------------------------------------" << endl;

    cout << "------------------------------------" << endl;
    cout << " ctxt0 > ctxt0 " << endl;
    helib::compareTwoNumbers(mu, ni, helib::CtPtrs_vectorCt(ctxtVec0),
                             helib::CtPtrs_vectorCt(ctxtVec0));
    secret_key.Decrypt(dMu, mu);
    cout << "Decrypted dMu: " << dMu << endl;
    secret_key.Decrypt(dNi, ni);
    cout << "Decrypted dNi: " << dNi << endl;
    cout << "------------------------------------" << endl;

    cout << "------------------------------------" << endl;
    cout << " ctxt1 > ctxt1 " << endl;
    helib::compareTwoNumbers(mu, ni, helib::CtPtrs_vectorCt(ctxtVec1),
                             helib::CtPtrs_vectorCt(ctxtVec1));
    secret_key.Decrypt(dMu, mu);
    cout << "Decrypted dMu: " << dMu << endl;
    secret_key.Decrypt(dNi, ni);
    cout << "Decrypted dNi: " << dNi << endl;
    cout << "------------------------------------" << endl;

    cout << "------------------------------------" << endl;
    cout << " ctxt2 > ctxt2 " << endl;
    helib::compareTwoNumbers(mu, ni, helib::CtPtrs_vectorCt(ctxtVec2),
                             helib::CtPtrs_vectorCt(ctxtVec2));
    secret_key.Decrypt(dMu, mu);
    cout << "Decrypted dMu: " << dMu << endl;
    secret_key.Decrypt(dNi, ni);
    cout << "Decrypted dNi: " << dNi << endl;
    cout << "------------------------------------" << endl;

    cout << "------------------------------------" << endl;
    cout << " ctxt1 > ctxt0 " << endl;
    helib::compareTwoNumbers(mu, ni, helib::CtPtrs_vectorCt(ctxtVec1),
                             helib::CtPtrs_vectorCt(ctxtVec0));
    secret_key.Decrypt(dMu, mu);
    cout << "Decrypted dMu: " << dMu << endl;
    secret_key.Decrypt(dNi, ni);
    cout << "Decrypted dNi: " << dNi << endl;
    cout << "------------------------------------" << endl;

    cout << "------------------------------------" << endl;
    cout << " ctxt0 > ctxt1 " << endl;
    helib::compareTwoNumbers(mu, ni, helib::CtPtrs_vectorCt(ctxtVec0),
                             helib::CtPtrs_vectorCt(ctxtVec1));
    secret_key.Decrypt(dMu, mu);
    cout << "Decrypted dMu: " << dMu << endl;
    secret_key.Decrypt(dNi, ni);
    cout << "Decrypted dNi: " << dNi << endl;
    cout << "------------------------------------" << endl;

    cout << "------------------------------------" << endl;
    cout << " ctxt > ctxt2 " << endl;
    helib::compareTwoNumbers(mu, ni, helib::CtPtrs_vectorCt(ctxtVec),
                             helib::CtPtrs_vectorCt(ctxtVec2));
    secret_key.Decrypt(dMu, mu);
    cout << "Decrypted dMu: " << dMu << endl;
    secret_key.Decrypt(dNi, ni);
    cout << "Decrypted dNi: " << dNi << endl;
    cout << "------------------------------------" << endl;

    cout << "------------------------------------" << endl;
    cout << " ctxt2 > ctxt " << endl;
    helib::compareTwoNumbers(mu, ni, helib::CtPtrs_vectorCt(ctxtVec2),
                             helib::CtPtrs_vectorCt(ctxtVec));
    secret_key.Decrypt(dMu, mu);
    cout << "Decrypted dMu: " << dMu << endl;
    secret_key.Decrypt(dNi, ni);
    cout << "Decrypted dNi: " << dNi << endl;
    cout << "------------------------------------" << endl;

    cout << "------------------------------------" << endl;
    cout << " ctxt1 > ctxt " << endl;
    helib::compareTwoNumbers(mu, ni, helib::CtPtrs_vectorCt(ctxtVec1),
                             helib::CtPtrs_vectorCt(ctxtVec));
    secret_key.Decrypt(dMu, mu);
    cout << "Decrypted dMu: " << dMu << endl;
    secret_key.Decrypt(dNi, ni);
    cout << "Decrypted dNi: " << dNi << endl;
    cout << "------------------------------------" << endl;

    cout << "------------------------------------" << endl;
    cout << " ctxt > ctxt1 " << endl;
    helib::compareTwoNumbers(mu, ni, helib::CtPtrs_vectorCt(ctxtVec),
                             helib::CtPtrs_vectorCt(ctxtVec1));
    secret_key.Decrypt(dMu, mu);
    cout << "Decrypted dMu: " << dMu << endl;
    secret_key.Decrypt(dNi, ni);
    cout << "Decrypted dNi: " << dNi << endl;
    cout << "------------------------------------" << endl;


    cout << "------------------------------------" << endl;
    cout << " ctxt2 > ctxt   ---   With Extract Bits\n";
    std::vector<Ctxt> bitsVec;
    ctxt.extractBits(bitsVec);
    std::vector<Ctxt> bitsVec2;
    ctxt2.extractBits(bitsVec2);

    helib::compareTwoNumbers(mu, ni,
                             helib::CtPtrs_vectorCt(bitsVec2),
                             helib::CtPtrs_vectorCt(bitsVec));

    secret_key.Decrypt(dMu, mu);
    cout << "Decrypted dMu: " << dMu << endl;
    secret_key.Decrypt(dNi, ni);
    cout << "Decrypted dNi: " << dNi << endl;
    cout << "------------------------------------" << endl;


    cout << "------------------------------------" << endl;
    cout << " ctxt2 > ctxt   ---   With Extract Bits\n";
    helib::compareTwoNumbers(mu, ni,
                             helib::CtPtrs_vectorCt(bitsVec2),
                             helib::CtPtrs_vectorCt(bitsVec),
                             true
    );

    secret_key.Decrypt(dMu, mu);
    cout << "Decrypted dMu: " << dMu << endl;
    secret_key.Decrypt(dNi, ni);
    cout << "Decrypted dNi: " << dNi << endl;
    cout << "------------------------------------" << endl;

    cout << "------------------------------------" << endl;
    cout << " ctxt2 > ctxt2   ---   With Extract Bits\n";
    helib::compareTwoNumbers(mu, ni,
                             helib::CtPtrs_vectorCt(bitsVec2),
                             helib::CtPtrs_vectorCt(bitsVec2)
    );

    secret_key.Decrypt(dMu, mu);
    cout << "Decrypted dMu: " << dMu << endl;
    secret_key.Decrypt(dNi, ni);
    cout << "Decrypted dNi: " << dNi << endl;
    cout << "------------------------------------" << endl;

    cout << "------------------------------------" << endl;
    cout << " ctxt2 > ctxt2   ---   With Extract Bits\n";
    helib::compareTwoNumbers(mu, ni,
                             helib::CtPtrs_vectorCt(bitsVec2),
                             helib::CtPtrs_vectorCt(bitsVec2),
                             true
    );

    secret_key.Decrypt(dMu, mu);
    cout << "Decrypted dMu: " << dMu << endl;
    secret_key.Decrypt(dNi, ni);
    cout << "Decrypted dNi: " << dNi << endl;
    cout << "------------------------------------" << endl;

    cout << "------------------------------------" << endl;
    std::vector<Ctxt> digits;
    helib::extractDigits(digits, ctxt);

    std::vector<long> pNums;
    helib::decryptBinaryNums(pNums, helib::CtPtrs_vectorCt(digits), secret_key, ea, false, true);

    cout << "==========" << endl;
    for (auto const &n: pNums) cout << n << ",";
    cout << endl;
    cout << "==========" << endl;
    printNameVal(pNums.size());
    printNameVal(digits.size());
    cout << "==========" << endl;
    helib::Ptxt<helib::BGV> plaintxt(public_key);
    for (auto const &d: digits) {
        secret_key.Decrypt(plaintxt, d);
        cout << "[" << plaintxt << "]";
    }
    cout << endl;
    cout << "==========" << endl;
    cout << "------------------------------------" << endl;

    cout << "------------------------------------" << endl;
    cout << "------------------------------------" << endl;
    cout << "------------------------------------" << endl;
    helib::Ctxt copyCtxt(ctxt);
    helib::Ctxt copyCtxt0(ctxt0);
    helib::Ctxt copyCtxt1(ctxt1);
    helib::Ctxt copyCtxt2(ctxt2);
    rotate(copyCtxt, 2);
    rotate(copyCtxt0, 2);
    rotate(copyCtxt1, 2);
    rotate(copyCtxt2, 2);
    helib::Ptxt<helib::BGV> temp(public_key);
    helib::Ptxt<helib::BGV> temp0(public_key);
    helib::Ptxt<helib::BGV> temp1(public_key);
    helib::Ptxt<helib::BGV> temp2(public_key);
    secret_key.Decrypt(temp, copyCtxt);
    secret_key.Decrypt(temp0, copyCtxt0);
    secret_key.Decrypt(temp1, copyCtxt1);
    secret_key.Decrypt(temp2, copyCtxt2);
    printNameVal(temp);
    printNameVal(temp0);
    printNameVal(temp1);
    printNameVal(temp2);
    //    std::vector<helib::CtxtPart> parts = ctxt.parts
    cout << "------------------------------------" << endl;
    cout << "------------------------------------" << endl;
    cout << "------------------------------------" << endl;


    // Double it (using additions)
    // [0] [1] [1] ... [1] [1] -> [0] [2] [2] ... [2] [2]
    ctxt += ctxt;
    // Plaintext version
    ptxt += ptxt;

    // Subtract it from itself (result should be 0)
    // i.e. [0] [0] [0] [0] ... [0] [0]
    ctxt -= ctxt;
    // Plaintext version
    ptxt -= ptxt;

    // Create a plaintext for decryption
    //    helib::Ptxt<helib::BGV> plaintext_result(public_key);
    // Decrypt the modified ciphertext
    secret_key.Decrypt(plaintext_result, ctxt);

    cout << "Operation: 2(a*a)/(a*a) - 2(a*a)/(a*a) = 0" << endl;
    // Print the decrypted plaintext
    // Should be [0] [0] [0] ... [0] [0]
    cout << "Decrypted Result: " << plaintext_result << endl;
    // Print the plaintext version result, should be the same as the ctxt version
    cout << "Plaintext Result: " << ptxt << endl;

    // We can also add constants
    // [0] [0] [0] ... [0] [0] -> [1] [1] [1] ... [1] [1]
    ctxt.addConstant(NTL::ZZX(1l));
    // Plaintext version
    ptxt.addConstant(NTL::ZZX(1l));

    // And multiply by constants
    // [1] [1] [1] ... [1] [1]
    // -> [1*1] [1*1] [1*1] ... [1*1] [1*1] = [1] [1] [1] ... [1] [1]
    ctxt *= 1l;
    // Plaintext version
    ptxt *= 1l;

    // We can also perform ciphertext-plaintext operations
    // ctxt = [1] [1] [1] ... [1] [1], ptxt = [1] [1] [1] ... [1] [1]
    // ctxt + ptxt = [2] [2] [2] ... [2] [2]
    // Note: the output of this is also a ciphertext
    ctxt += ptxt;

    // Decrypt the modified ciphertext into a new plaintext
    helib::Ptxt<helib::BGV> new_plaintext_result(context);
    secret_key.Decrypt(new_plaintext_result, ctxt);

    cout << "Operation: Enc{(0 + 1)*1} + (0 + 1)*1" << endl;
    // Print the decrypted plaintext
    // Should be [2] [2] [2] ... [2] [2]
    cout << "Decrypted Result: " << new_plaintext_result << endl;


    // Square the ciphertext
    // [0] [1] [2] [3] [4] ... [nslots-1]
    // -> [0] [1] [4] [9] [16] ... [(nslots-1)*(nslots-1)]
    ctxt.multiplyBy(ctxt);
    // Plaintext version
    ptxt.multiplyBy(ptxt);

    // Divide the ciphertext by itself
    // To do this we must calculate the multiplicative inverse using Fermat's
    // Little Theorem.  We calculate a^{-1} = a^{p-2} mod p, where a is non-zero
    // and p is our plaintext prime.
    // First make a copy of the ctxt using copy constructor
    helib::Ctxt ctxt_divisor(ctxt);
    // Raise the copy to the exponent p-2
    // [0] [1] [4] ... [16] -> [0] [1] [1] ... [1]
    // Note: 0 is a special case because 0^n = 0 for any power n
    ctxt_divisor.power(p - 2);
    // a^{p-2}*a = a^{-1}*a = a / a = 1;
    ctxt.multiplyBy(ctxt_divisor);

    // Plaintext version
    helib::Ptxt<helib::BGV> ptxt_divisor(ptxt);
    ptxt_divisor.power(p - 2);
    ptxt.multiplyBy(ptxt_divisor);
}

/**
 * learning threads
 * */
#include <thread>

int GetRandom(int max) {
    srand(time(NULL));
    return rand() % max;
}

std::string GetTime() {
    auto nowTime = CLOCK::now();
    std::time_t sleepTime = CLOCK::to_time_t(nowTime);

    tm myLocalTime = *localtime((&sleepTime));
    cout << "Hours: " << myLocalTime.tm_hour << " ";
    cout << "Minutes: " << myLocalTime.tm_min << " ";
    cout << "Seconds: " << myLocalTime.tm_sec << " ";

    return std::ctime(&sleepTime);
}

void ExecuteThread(int id) {
    cout << "Thread " << id << " Sleep Time : " << GetTime() << endl;
    std::this_thread::sleep_for(std::chrono::seconds(GetRandom(5)));
    cout << "Thread " << id << " Awake Time : " << GetTime() << endl;
}

double accBal = 100;
std::mutex accLock;

void GetMoney(int id, double amount) {
    std::lock_guard<std::mutex> lockGuard(accLock);
    std::this_thread::sleep_for(std::chrono::seconds(3));
    cout << "\nid-" << id << " tries to withdraw $" << amount << " on "
         << GetTime();// << endl;
    if (accBal - amount >= 0) {
        accBal -= amount;
        cout << "New Balance is $" << accBal << endl;
    } else {
        cout << "Not Enough Money in Account\nCurrent Balance is $" << accBal << endl;
    }
}

void FindPrimes(unsigned int start, unsigned int end, std::vector<unsigned int> &vect) {
    for (unsigned int x = start; x <= end; x += 2)
        for (unsigned int y = 2; y < x; ++y) {
            if (0 == (x % y)) break;
            else if (x == (y + 1)) vect.push_back(x);
        }
}

std::mutex vecLock;
std::vector<unsigned int> primeVecGlobal;

void FindPrimesWithGlobalVec(unsigned int start, unsigned int end) {
    for (unsigned int x = start; x <= end; x += 2)
        for (unsigned int y = 2; y < x; ++y) {
            if (0 == (x % y)) break;
            else if (x == (y + 1)) {
                vecLock.lock();
                primeVecGlobal.push_back(x);
                vecLock.unlock();
            }
        }
}

void FindPrimesWithThreads(unsigned int start, unsigned int end, unsigned int numThreads) {
    std::vector<std::thread> threadVec;
    unsigned int threadSpread = end / numThreads;
    unsigned int newEnd = start + threadSpread - 1;
    for (unsigned int x = 0; x < numThreads; ++x) {
        threadVec.emplace_back(FindPrimesWithGlobalVec, start, newEnd);
        start += threadSpread;
        newEnd += threadSpread;
    }
    for (auto &t:threadVec) t.join();

}

/**
 * fin learning threads
 * */

void TestAux::testMultithreading() {
    std::thread th1(ExecuteThread, 1);
    th1.join();
    std::thread th2(ExecuteThread, 2);
    th2.join();
    //  -------------------------------------------------------
    int num = 10;
    std::thread threads[num];
    for (int i = 0; i < num; ++i) {
        threads[i] = std::thread(GetMoney, i, 15);
    }
    for (int i = 0; i < num; ++i) {
        threads[i].join();
    }
    //  -------------------------------------------------------
    std::vector<unsigned int> primeVec;
    int startTime = clock();
    FindPrimes(1, 100000, primeVec);
    int endTime = clock();
    for (auto i:primeVec) cout << i << ", ";
    cout << endl;
    cout << "Execution Time : " <<
         (endTime - startTime) / double(CLOCKS_PER_SEC) << endl;
    primeVec.clear();
    //  -------------------------------------------------------
    startTime = clock();
    FindPrimesWithThreads(1, 100000, 3);
    endTime = clock();
    for (auto i:primeVec) cout << i << ", ";
    cout << endl;
    cout << "Execution Time : " <<
         (endTime - startTime) / double(CLOCKS_PER_SEC) << endl;

}

#include <sstream>
#include <string>

#include "src/coreset/run1meancore.h"

void TestAux::testPythonRun() {

    //    string commandLS ="ls -l */*"; system(commandLS.c_str());
    string command = "/home/karina/CLionProjects/EncryptedKMeans/venv/bin/python3.7 -m src.coreset.chkchk";
    system(command.c_str());
    //     command = "/home/karina/CLionProjects/EncryptedKMeans/venv/bin/python3.7 src/coreset/chkchk.py";
    //    system(command.c_str());
    //     command = "/home/karina/CLionProjects/EncryptedKMeans/venv/bin/python3.7 ~/CLionProjects/EncryptedKMeans/src/coreset/chkchk.py";
    //    system(command.c_str());
    //     command = "/home/karina/CLionProjects/EncryptedKMeans/venv/bin/python3.7 -m .src.coreset.chkchk";
    //    system(command.c_str());
    //    command = "/home/karina/CLionProjects/EncryptedKMeans/venv/bin/python3.7 -m ..src.coreset.chkchk";
    //    system(command.c_str());
    //    command = "/home/karina/CLionProjects/EncryptedKMeans/venv/bin/python3.7 -m EncryptedKMeans.src.coreset.chkchk";
    //    system(command.c_str());
    //    command = "/home/karina/CLionProjects/EncryptedKMeans/venv/bin/python3.7 -m cmake-build-debug.src.coreset.chkchk";
    //    system(command.c_str());
    //    command = "/home/karina/CLionProjects/EncryptedKMeans/venv/bin/python3.7 -m .coreset.chkchk";
    //    system(command.c_str());
    //    command = "/home/karina/CLionProjects/EncryptedKMeans/venv/bin/python3.7 -m coreset.chkchk";
    //    system(command.c_str());
    //    command = "/home/karina/CLionProjects/EncryptedKMeans/venv/bin/python3.7 -m chkchk";
    //    system(command.c_str());
    //    string commandl = "ls -l ";
    //    system(commandl.c_str());
    //    string commandl2 = "pwd ";
    //    system(commandl2.c_str());


    printNameVal(" - TIME - ");
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    printNameVal("VERSION 1");
    cout << std::put_time(&tm, "%d-%m-%Y %H-%M-%S") << endl;
    cout << std::put_time(&tm, "%Y-%m-%d-%H-%M-%S") << endl;
    _Put_time<char> bla = std::put_time(&tm, "%Y-%m-%d-%H-%M-%S");
    string blabla = string(bla._M_fmt);
    printNameVal(blabla);

    printNameVal("VERSION 2");
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y %H-%M-%S") << endl;
    oss << std::put_time(&tm, "%Y-%m-%d-%H-%M-%S") << endl;
    cout << oss.str() << endl;
}

std::vector<long> randomLongsVector(short size, short upperLimit = NUMBERS_RANGE) {
    std::vector<long> rlv(size);
    static std::uniform_int_distribution<long>
            randomLongInRangeMiniTest(0, upperLimit);
    for (int i = 0; i < size; ++i) {
        //        long rl = randomLongInRange(mt);
        long rl = randomLongInRangeMiniTest(mt);
        cout << rl << " ";
        rlv[i] = rl;
    }
    //    cout << endl;
    return rlv;
}

void TestAux::testComparison_diffCtxtRepresentation() {
    /*  Example of BGV scheme  */

    KeysServer keysServer;
    cout << "Initialising context object..." << endl;
    // Initialize context
    // This object will hold information about the algebra created from the
    // previously set parameters
    const helib::Context &context = keysServer.getContextDBG();

    // Print the context
    context.printout();
    cout << endl;

    // Print the security level
    cout << "Security: " << context.securityLevel() << endl;

    // Secret key management
    const helib::SecKey &secret_key = keysServer.getSecKeyDBG();

    // Public key management
    // Set the secret key (upcast: SecKey is a subclass of PubKey)
    const helib::PubKey &public_key = secret_key;

    // Get the EncryptedArray of the context
    const helib::EncryptedArray &ea = context.getEA();

    // Get the number of slot (phi(m))
    long nslots = ea.size();
    cout << "Number of slots: " << nslots << endl;

    //    long arrOfArrs[3][5] = {
    //            {0, 1, 2, 3, 5},
    //            {1, 2, 3, 4, 5},
    //            {4, 7, 2, 5, 2}
    //    };
    std::vector<std::vector<long>> arrOfArrs(5);

    for (int i = 0; i < arrOfArrs.size(); ++i)
        arrOfArrs[i] = randomLongsVector(BIT_SIZE);
    cout << "\n";

    helib::Ptxt<helib::BGV> plaintext_result(public_key);
    for (const auto &arrayOfSomeNums : arrOfArrs) {
        cout << "===========================================================\n";
        cout << "Running with array of numbers: \n";
        cout << "{ ";
        for (long n:arrayOfSomeNums) cout << n << " ";
        cout << "}  \n";
        cout << "===========================================================\n";

        //  bitVecI contains all the i'th bits of all the points
        //  bitVecI[J] contains I'th bit of J'th point
        std::vector<std::vector<Ctxt>>
                arr2bits(BIT_SIZE, std::vector(BIT_SIZE, helib::Ctxt(public_key)));
        std::vector<std::vector<Ctxt>>
                num2bits(BIT_SIZE, std::vector(BIT_SIZE, helib::Ctxt(public_key)));
        std::vector<Ctxt> number_encrypted(BIT_SIZE, helib::Ctxt(public_key));

        Ctxt mu(public_key), ni(public_key);
        std::vector<Ctxt>
                resVectorMu(BIT_SIZE, Ctxt(public_key)),
                resVectorNi(BIT_SIZE, Ctxt(public_key)),
                resVectorMu2(BIT_SIZE, Ctxt(public_key)),
                resVectorNi2(BIT_SIZE, Ctxt(public_key));

        for (int number_compared = 0; number_compared < NUMBERS_RANGE; ++number_compared) {
            cout << "array: { ";
            for (long n:arrayOfSomeNums) cout << n << " ";
            cout << " }  compared with number " << number_compared << " \t[ encrypted : ";
            for (int bit = BIT_SIZE - 1; bit >= 0; --bit) cout << ((number_compared >> bit) & 1);
            cout << " ] " << endl;
            //            cout << "---------------\n";
            cout << "\t --- \tcompare with number as bit matrix TRANSPOSED \t --- " << endl;
            for (short bit = 0; bit < BIT_SIZE; ++bit) {
                for (int num = 0; num < BIT_SIZE; ++num) {
                    public_key.Encrypt(
                            arr2bits[bit][num],
                            NTL::to_ZZX((arrayOfSomeNums[num] >> bit) & 1));
                    public_key.Encrypt(
                            num2bits[bit][num],
                            NTL::to_ZZX((number_compared >> bit) & 1));
                }
                public_key.Encrypt(
                        number_encrypted[bit],
                        NTL::to_ZZX((number_compared >> bit) & 1));

                cout << "num2bits[" << bit << "]: ";
                for (int j = BIT_SIZE - 1; j >= 0; --j) cout << ((number_compared >> bit) & 1);
                cout << "  \t  ";
                cout << "arr2bits[" << bit << "]: ";
                for (int j = BIT_SIZE - 1; j >= 0; --j) cout << ((arrayOfSomeNums[j] >> bit) & 1);
                cout << "  = " << keysServer.decryptNum(arr2bits[bit]);
                cout << "   \t  ";

                helib::compareTwoNumbers(mu, ni,
                                         helib::CtPtrs_vectorCt(arr2bits[bit]),
                                         helib::CtPtrs_vectorCt(num2bits[bit]));
                resVectorMu[bit] = mu;
                resVectorNi[bit] = ni;

                secret_key.Decrypt(plaintext_result, resVectorMu[bit]);
                cout << "resVectorMu[" << bit << "]: " << plaintext_result[0]
                     << "   \t  ";
                secret_key.Decrypt(plaintext_result, resVectorNi[bit]);
                cout << "resVectorNi[" << bit << "]: " << plaintext_result[0]
                     << "  \t  ";
                cout << endl;

            }
            cout << " \t \t \t \t \t \t \t \t \t \t \t \t "
                 << "\t \t resVectorMu: " << keysServer.decryptNum(resVectorMu)
                 << "\t \t resVectorNi: " << keysServer.decryptNum(resVectorNi)
                 << endl;

            cout << "\t --- \tcompare with number as bit matrix (original)\t ---" << endl;
            for (short bit = 0; bit < BIT_SIZE; ++bit) {

                cout << "num2bits[" << bit << "]: ";
                for (int j = BIT_SIZE - 1; j >= 0; --j) cout << ((number_compared >> j) & 1);
                cout << "  \t  ";
                cout << "arr2bits[" << bit << "]: ";
                for (int j = BIT_SIZE - 1; j >= 0; --j) cout << ((arrayOfSomeNums[j] >> bit) & 1);
                cout << " = " << keysServer.decryptNum(arr2bits[bit]);
                cout << "   \t  ";

                helib::compareTwoNumbers(mu, ni,
                                         helib::CtPtrs_vectorCt(arr2bits[bit]),
                                         helib::CtPtrs_vectorCt(number_encrypted));
                resVectorMu2[bit] = mu;
                resVectorNi2[bit] = ni;

                secret_key.Decrypt(plaintext_result, resVectorMu2[bit]);
                //                cout << "resVectorMu2[" << bit << "]: " << plaintext_result
                cout << "resVectorMu2[" << bit << "]: " << plaintext_result[0]
                     << "    \t  ";
                secret_key.Decrypt(plaintext_result, resVectorNi2[bit]);
                //                cout << "resVectorNi2[" << bit << "]: " << plaintext_result << endl;
                cout << "resVectorNi2[" << bit << "]: " << plaintext_result[0] << endl;
            }
            cout << " \t \t \t \t \t \t \t \t \t \t \t \t "
                 << "\t \t resVectorMu2: " << keysServer.decryptNum(resVectorMu2)
                 << "\t \t resVectorNi2: " << keysServer.decryptNum(resVectorNi2)
                 << endl;
            cout << "----------------------------------------------------\n";
        }
    }
}

void TestAux::testComparison_diffCtxtRepresentation_BGVPackedArithmetics() {

    KeysServer keysServer;
    cout << "Initialising context object..." << endl;
    // Initialize context
    // This object will hold information about the algebra created from the
    // previously set parameters
    const helib::Context &context = keysServer.getContextDBG();

    // Print the context
    context.printout();
    cout << endl;

    // Print the security level
    cout << "Security: " << context.securityLevel() << endl;

    // Secret key management
    const helib::SecKey &secret_key = keysServer.getSecKeyDBG();

    // Public key management
    // Set the secret key (upcast: SecKey is a subclass of PubKey)
    const helib::PubKey &public_key = secret_key;

    // Get the EncryptedArray of the context
    const helib::EncryptedArray &ea = context.getEA();

    // Get the number of slot (phi(m))
    long nslots = ea.size();
    cout << "Number of slots: " << nslots << endl;

    cout << "==========================\n";
    cout << "==========================\n";
    cout << "==========================\n";

    // Create Data Arrs
    std::vector<std::vector<long>> arrOfArrs(10);
    for (int i = 0; i < arrOfArrs.size(); ++i)
        arrOfArrs[i] = randomLongsVector(nslots, pow(2, nslots) - 1);
    cout << "\n";

    std::vector<helib::Ptxt<helib::BGV> >
            arrsTransposedAsPtxt(nslots, helib::Ptxt<helib::BGV>(public_key));
    std::vector<helib::Ctxt>
            arrAsCtxtTransposed(nslots, helib::Ctxt(public_key));
    std::vector<helib::Ptxt<helib::BGV> >
            someNumberAsArrTransposed(nslots, helib::Ptxt<helib::BGV>(public_key));
    std::vector<helib::Ctxt>
            someNumberAsArrTransposed_Ctxt(nslots, Ctxt(public_key));

    helib::Ptxt<helib::BGV> plaintext_result(public_key);
    helib::Ctxt mu(public_key), ni(public_key);

    // for each array in the list of arrays
    for (int arrIndex = 0; arrIndex < arrOfArrs.size(); ++arrIndex) {
        cout << "==========================\n";

        //  Convert into binary data
        for (int numIndex = 0; numIndex < nslots; ++numIndex) // for each number
            for (int bitIndex = 0; bitIndex < nslots; ++bitIndex) // for each bit
                arrsTransposedAsPtxt[bitIndex][numIndex] =
                        (arrOfArrs[arrIndex][numIndex] >> bitIndex) & 1;

        for (int i = 0; i < nslots; ++i)
            public_key.Encrypt(arrAsCtxtTransposed[i], arrsTransposedAsPtxt[i]);

        //   to get a bulk of results for getting a sense of a pattern
        //      compare current array with all numbers in range
        for (int someNumber = 0; someNumber < pow(2, nslots) - 1; ++someNumber) {

            cout << "__________________________\n";
            printNameVal(someNumber);
            cout << "Compare to Array: { ";
            for (int i = 0; i < nslots; ++i) cout << arrOfArrs[arrIndex][i] << " ";
            cout << "} " << endl;

            //  Convert into binary data
            for (int bit = 0; bit < nslots; ++bit)
                for (int i = 0; i < nslots; ++i)
                    someNumberAsArrTransposed[bit][i] = (someNumber >> bit) & 1;

            for (int i = 0; i < nslots; ++i)
                public_key.Encrypt(someNumberAsArrTransposed_Ctxt[i], someNumberAsArrTransposed[i]);

            //  Compare
            helib::compareTwoNumbers(mu, ni,
                                     helib::CtPtrs_vectorCt(someNumberAsArrTransposed_Ctxt),
                                     helib::CtPtrs_vectorCt(arrAsCtxtTransposed));

            cout << "   ---     Results:    ---- \n";
            secret_key.Decrypt(plaintext_result, mu);
            cout << "mu: " << plaintext_result << endl;
            secret_key.Decrypt(plaintext_result, ni);
            cout << "ni: " << plaintext_result << endl;

            //            printNameVal(mu.bitCapacity());
            /*           std::vector<Ctxt> extractedBits;//(nslots, Ctxt(public_key));
                       Ctxt mush(public_key);

                       mush = mu.operator>>=(1);
                       secret_key.Decrypt(plaintext_result, mush);
                       cout << "mu.operator>>=(1);: " << plaintext_result << endl;

                       //            rotate(mu, 1);
                       extractDigits(extractedBits, mu, 1);
                       for (const auto &eb:extractedBits) {
                           secret_key.Decrypt(plaintext_result, eb);
                           cout << "mu.extractBits: " << plaintext_result << endl;
                       }

                       mush = mu.operator>>=(1);
                       secret_key.Decrypt(plaintext_result, mush);
                       cout << "mu.operator>>=(1);: " << plaintext_result << endl;
                       //            rotate(mu, 1);
                       extractDigits(extractedBits, mu, 2);
                       for (const auto &eb:extractedBits) {
                           secret_key.Decrypt(plaintext_result, eb);
                           cout << "mu.extractBits: " << plaintext_result << endl;
                       }

                       mush = mu.operator>>=(1);
                       secret_key.Decrypt(plaintext_result, mush);
                       cout << "mu.operator>>=(1);: " << plaintext_result << endl;
                       //            rotate(mu, 1);
                       extractDigits(extractedBits, mu, 3);
                       for (const auto &eb:extractedBits) {
                           secret_key.Decrypt(plaintext_result, eb);
                           cout << "mu.extractBits: " << plaintext_result << endl;
                       }

                       mush = mu.operator>>=(1);
                       secret_key.Decrypt(plaintext_result, mush);
                       cout << "mu.operator>>=(1);: " << plaintext_result << endl;
                       //            rotate(mu, 1);
                       extractDigits(extractedBits, mu, 4);
                       for (const auto &eb:extractedBits) {
                           secret_key.Decrypt(plaintext_result, eb);
                           cout << "mu.extractBits: " << plaintext_result << endl;
                       }

           */
        }
        cout << "==========================\n";
    }

    //    TODO
    // todo post a question on stackoverflow:
    //    how to get a single slot data from a Ctxt

}

void TestAux::testIsMatchImplementation() {
    KeysServer keysServer;
    cout << "Initialising context object..." << endl;
    // Initialize context
    // This object will hold information about the algebra created from the
    // previously set parameters
    const helib::Context &context = keysServer.getContextDBG();

    // Print the context
    context.printout();
    cout << endl;

    // Print the security level
    cout << "Security: " << context.securityLevel() << endl;

    // Secret key management
    const helib::SecKey &secret_key = keysServer.getSecKeyDBG();

    // Public key management
    // Set the secret key (upcast: SecKey is a subclass of PubKey)
    const helib::PubKey &public_key = secret_key;

    // Get the EncryptedArray of the context
    const helib::EncryptedArray &ea = context.getEA();

    // Get the number of slot (phi(m))
    long nslots = ea.size();
    cout << "Number of slots: " << nslots << endl;


    std::vector<Ctxt> number_encrypted(BIT_SIZE, helib::Ctxt(public_key));
    std::vector<Ctxt> number_encrypted_2(BIT_SIZE, helib::Ctxt(public_key));
    Ctxt mu(public_key), ni(public_key);
    for (int i = 0; i < 5; ++i) {

        long number_compared = giveMeRandomLong();
        long number_compared_2 = giveMeRandomLong();

        cout << "---------------\n";
        printNameVal(number_compared);
        cout << "In Binary: ";
        for (int bit = BIT_SIZE - 1; bit >= 0; --bit) cout << ((number_compared >> bit) & 1);
        cout << endl;
        printNameVal(number_compared_2);
        cout << "In Binary: ";
        for (int bit = BIT_SIZE - 1; bit >= 0; --bit) cout << ((number_compared_2 >> bit) & 1);
        cout << endl;

        for (short bit = 0; bit < BIT_SIZE; ++bit) {

            public_key.Encrypt(
                    number_encrypted[bit],
                    NTL::to_ZZX((number_compared >> bit) & 1));

            public_key.Encrypt(
                    number_encrypted_2[bit],
                    NTL::to_ZZX((number_compared_2 >> bit) & 1));

        }

        helib::compareTwoNumbers(mu, ni,
                                 helib::CtPtrs_vectorCt(number_encrypted),
                                 helib::CtPtrs_vectorCt(number_encrypted_2));

        helib::Ptxt<helib::BGV> plaintext_result(public_key);
        secret_key.Decrypt(plaintext_result, mu);
        cout << "mu: " << plaintext_result << endl;
        secret_key.Decrypt(plaintext_result, ni);
        cout << "ni: " << plaintext_result << endl;

        printNameVal(mu.bitCapacity());

        std::vector<Ctxt> subs(BIT_SIZE, Ctxt(public_key));
        std::vector<Ctxt> squares(BIT_SIZE, Ctxt(public_key));
        std::vector<Ctxt> nots(BIT_SIZE, Ctxt(public_key));
        Ctxt product(public_key);
        Ctxt result(public_key);

        /*
        for (int bit = 0; bit < BIT_SIZE; ++bit) {
//            subs[i].negate();
//            subs[i].addCtxt(number_encrypted[i]);
            subs[i].addCtxt(number_encrypted_2[i], true);
            squares[i] = subs[i];
            squares[i] *= subs[i];
//            squares[i].power(2)
            nots[i] = squares[i];
            nots[i].negate();
            nots[i].addConstant(1l);
        }
*/

        //      a_i - b_i
        helib::vecCopy(subs, number_encrypted);
        cout << "subs Vector: ";
        for (int bit = BIT_SIZE - 1; bit >= 0; --bit) {
            secret_key.Decrypt(plaintext_result, subs[bit]);
            cout << plaintext_result[0];
        }
        cout << endl;

        for (int bit = 0; bit < BIT_SIZE; ++bit) subs[i] -= number_encrypted_2[i];
        cout << "subs Vector: ";
        for (int bit = BIT_SIZE - 1; bit >= 0; --bit) {
            secret_key.Decrypt(plaintext_result, subs[bit]);
            cout << plaintext_result[0];
        }
        cout << endl;
        //        for (int bit = BIT_SIZE - 1; bit >= 0; --bit) {
        //            secret_key.Decrypt(plaintext_result, subs[bit]);
        //            printNameVal(plaintext_result);
        //        }

        //      (a_i - b_i) ^ 2
        //        for (int bit = 0; bit < BIT_SIZE; ++bit) squares[i] = subs[i];
        helib::vecCopy(squares, subs);
        cout << "squares Vector: ";
        for (int bit = BIT_SIZE - 1; bit >= 0; --bit) {
            secret_key.Decrypt(plaintext_result, squares[bit]);
            cout << plaintext_result[0];
        }
        cout << endl;

        for (int bit = 0; bit < BIT_SIZE; ++bit)
            squares[i].square(); // squares[i] *= subs[i];    //    squares[i].power(2)
        cout << "squares Vector: ";
        for (int bit = BIT_SIZE - 1; bit >= 0; --bit) {
            secret_key.Decrypt(plaintext_result, squares[bit]);
            cout << plaintext_result[0];
        }
        cout << endl;
        //        for (int bit = BIT_SIZE - 1; bit >= 0; --bit) {
        //            secret_key.Decrypt(plaintext_result, subs[bit]);
        //            printNameVal(plaintext_result);
        //        }

        //      1 - [ (a_i - b_i) ^ 2 ]
        for (int bit = 0; bit < BIT_SIZE; ++bit) nots[i] = squares[i];
        cout << "nots Vector: ";
        for (int bit = BIT_SIZE - 1; bit >= 0; --bit) {
            secret_key.Decrypt(plaintext_result, nots[bit]);
            cout << plaintext_result[0];
        }
        cout << endl;

        for (int bit = 0; bit < BIT_SIZE; ++bit) nots[i].negate();
        cout << "nots Vector: ";
        for (int bit = BIT_SIZE - 1; bit >= 0; --bit) {
            secret_key.Decrypt(plaintext_result, nots[bit]);
            cout << plaintext_result[0];
        }
        cout << endl;

        for (int bit = 0; bit < BIT_SIZE; ++bit) nots[i].addConstant(1l);
        cout << "nots Vector: ";
        for (int bit = BIT_SIZE - 1; bit >= 0; --bit) {
            secret_key.Decrypt(plaintext_result, nots[bit]);
            cout << plaintext_result[0];
        }
        cout << endl;

    }
}


void TestAux::testPrefixAndSuffix() {
    cout << " ------ testPrefixAndSuffix ------ " << endl << endl;
    auto t0_main = CLOCK::now();

    KeysServer keysServer;
    const helib::Context &context = keysServer.getContextDBG();
    context.printout();    // Print the context
    cout << endl;
    cout << "Security: " << context.securityLevel() << endl;    // Print the security level
    const helib::SecKey &secret_key = keysServer.getSecKeyDBG();    // Secret key management
    // Public key management
    // Set the secret key (upcast: SecKey is a subclass of PubKey)
    const helib::PubKey &public_key = secret_key;
    const helib::EncryptedArray &ea = context.getEA();    // Get the EncryptedArray of the context
    long nslots = ea.size();    // Get the number of slot (phi(m))
    cout << "Number of slots: " << nslots << endl;
    cout << "---------------\n";
    cout << "---------------\n";
    cout << "---------------\n";
    cout << endl;

    std::vector<Ctxt> number_encrypted(BIT_SIZE, helib::Ctxt(public_key));
    std::vector<Ctxt> prefix_encrypted, suffix_encrypted;
    Ctxt mu(public_key), ni(public_key);
    helib::Ptxt<helib::BGV> plaintext_result(public_key);

    //  Create Data & Encrypt
    long number_compared = giveMeRandomLong();
    printNameVal(number_compared);
    cout << "In Binary: ";
    for (int bit = BIT_SIZE - 1; bit >= 0; --bit) cout << ((number_compared >> bit) & 1);
    cout << endl;

    for (short bit = 0; bit < BIT_SIZE; ++bit)
        public_key.Encrypt(
                number_encrypted[bit],
                NTL::to_ZZX((number_compared >> bit) & 1));

    for (int i = 0; i < BIT_SIZE; ++i) {
        cout << "---------------\n";
        printNameVal(i);

        cout << "Prefix In Binary: ";
        for (int bit = i-1; bit >= 0; --bit) cout << ((number_compared >> bit) & 1);
        cout << endl;

        cout << "Encrypted Prefix\n";
        prefix_encrypted = prefix(number_encrypted, long(i));

        for (int bit = prefix_encrypted.size()-1; bit >= 0; --bit){
            secret_key.Decrypt(plaintext_result, prefix_encrypted[bit]);
            cout << plaintext_result[0];
        }
        cout  << endl;

        cout << "Suffix In Binary: ";
        for (int bit = BIT_SIZE - 1; bit >= i; --bit) cout << ((number_compared >> bit) & 1);
        cout << endl;

        cout << "Encrypted Suffix\n";
        suffix_encrypted = suffix(number_encrypted, long(i));

        for (int bit = suffix_encrypted.size()-1; bit >= 0; --bit){
            secret_key.Decrypt(plaintext_result, suffix_encrypted[bit]);
            cout << plaintext_result[0];
        }
        cout  << endl;

    }

    cout << endl << printDuration(t0_main, "testPrefixAndSuffix");
    cout << " ------ testPrefixAndSuffix finished ------ " << endl << endl;
}

void TestAux::testIsEqualImplementation() {
    cout << " ------ testIsEqualImplementation ------ " << endl << endl;
    auto t0_main = CLOCK::now();

    KeysServer keysServer;
    const helib::Context &context = keysServer.getContextDBG();
    context.printout();    // Print the context
    cout << endl;
    cout << "Security: " << context.securityLevel() << endl;    // Print the security level
    const helib::SecKey &secret_key = keysServer.getSecKeyDBG();    // Secret key management
    // Public key management
    // Set the secret key (upcast: SecKey is a subclass of PubKey)
    const helib::PubKey &public_key = secret_key;
    const helib::EncryptedArray &ea = context.getEA();    // Get the EncryptedArray of the context
    long nslots = ea.size();    // Get the number of slot (phi(m))
    cout << "Number of slots: " << nslots << endl;

    std::vector<Ctxt> number_encrypted(BIT_SIZE, helib::Ctxt(public_key));
    std::vector<Ctxt> number_encrypted_2(BIT_SIZE, helib::Ctxt(public_key));
    Ctxt mu(public_key), ni(public_key);
    for (int i = 0; i < 5; ++i) {

        //  Create Data & Encrypt
        long number_compared = giveMeRandomLong();
        long number_compared_2 = giveMeRandomLong();

        cout << "---------------\n";
        printNameVal(number_compared);
        cout << "In Binary: ";
        for (int bit = BIT_SIZE - 1; bit >= 0; --bit) cout << ((number_compared >> bit) & 1);
        cout << endl;
        printNameVal(number_compared_2);
        cout << "In Binary: ";
        for (int bit = BIT_SIZE - 1; bit >= 0; --bit) cout << ((number_compared_2 >> bit) & 1);
        cout << endl;

        for (short bit = 0; bit < BIT_SIZE; ++bit) {

            public_key.Encrypt(
                    number_encrypted[bit],
                    NTL::to_ZZX((number_compared >> bit) & 1));

            public_key.Encrypt(
                    number_encrypted_2[bit],
                    NTL::to_ZZX((number_compared_2 >> bit) & 1));

        }

        //        Ctxt res = isEqual()

        helib::Ptxt<helib::BGV> plaintext_result(public_key);
        secret_key.Decrypt(plaintext_result, mu);
        cout << "mu: " << plaintext_result << endl;
        secret_key.Decrypt(plaintext_result, ni);
        cout << "ni: " << plaintext_result << endl;

        printNameVal(mu.bitCapacity());

    }

    cout << endl << printDuration(t0_main, "testIsEqualImplementation");
    cout << " ------ testIsEqualImplementation finished ------ " << endl << endl;
}
