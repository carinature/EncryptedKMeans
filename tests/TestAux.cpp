
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

    cout << endl<<printDuration(t0_main, "testGenerateDataClients");
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
    unsigned long p = 4999;
    // Cyclotomic polynomial - defines phi(m)
    unsigned long m = 32109;
    // Hensel lifting (default = 1)
    unsigned long r = 1;
    // Number of bits of the modulus chain
    unsigned long bits = 500;
    // Number of columns of Key-Switching matrix (default = 2 or 3)
    unsigned long c = 2;

    std::cout << "\n*********************************************************";
    std::cout << "\n*         Basic Mathematical Operations Example         *";
    std::cout << "\n*         =====================================         *";
    std::cout << "\n*                                                       *";
    std::cout << "\n* This is a sample program for education purposes only. *";
    std::cout << "\n* It attempts to show the various basic mathematical    *";
    std::cout << "\n* operations that can be performed on both ciphertexts  *";
    std::cout << "\n* and plaintexts.                                       *";
    std::cout << "\n*                                                       *";
    std::cout << "\n*********************************************************";
    std::cout << std::endl;

    std::cout << "Initialising context object..." << std::endl;
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
    std::cout << std::endl;

    // Print the security level
    std::cout << "Security: " << context.securityLevel() << std::endl;

    // Secret key management
    std::cout << "Creating secret key..." << std::endl;
    // Create a secret key associated with the context
    helib::SecKey secret_key(context);
    // Generate the secret key
    secret_key.GenSecKey();
    std::cout << "Generating key-switching matrices..." << std::endl;
    // Compute key-switching matrices that we need
    helib::addSome1DMatrices(secret_key);

    // Public key management
    // Set the secret key (upcast: SecKey is a subclass of PubKey)
    const helib::PubKey& public_key = secret_key;

    // Get the EncryptedArray of the context
    const helib::EncryptedArray& ea = context.getEA();

    // Get the number of slot (phi(m))
    long nslots = ea.size();
    std::cout << "Number of slots: " << nslots << std::endl;

    // Create a vector of long with nslots elements
    helib::Ptxt<helib::BGV> ptxt(context);
    // Set it with numbers 0..nslots - 1
    // ptxt = [0] [1] [2] ... [nslots-2] [nslots-1]
    for (int i = 0; i < ptxt.size(); ++i) {
        ptxt[i] = i;
    }

    // Print the plaintext
    std::cout << "Initial Plaintext: " << ptxt << std::endl;

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

    std::cout << "Operation: 2(a*a)/(a*a) - 2(a*a)/(a*a) = 0" << std::endl;
    // Print the decrypted plaintext
    // Should be [0] [0] [0] ... [0] [0]
    std::cout << "Decrypted Result: " << plaintext_result << std::endl;
    // Print the plaintext version result, should be the same as the ctxt version
    std::cout << "Plaintext Result: " << ptxt << std::endl;

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

    std::cout << "Operation: Enc{(0 + 1)*1} + (0 + 1)*1" << std::endl;
    // Print the decrypted plaintext
    // Should be [2] [2] [2] ... [2] [2]
    std::cout << "Decrypted Result: " << new_plaintext_result << std::endl;

}
