#include "DataServer.h"

#include "utils/aux.h"

using std::cout;
using std::endl;
static Logger dataServerLogger(log_debug, "dataServerLogger");

DataServer::DataServer(KeysServer &keysServer) :
        encryptionKey(keysServer.getPublicKey()),
//        public_key(keysServer.getPublicKey()),
//        public_key(encryptionKey),
//        ea(keysServer.getEA()),
        scratch(keysServer.getPublicKey()){

    dataServerLogger.log("DataServer()");
    cout << "DataServer()" << endl;
}

Point DataServer::scratchPoint() {
    cout << " scratchPoint"<<endl;
    return Point(encryptionKey);//, nullptr);
}

//std::vector<Client> DataServer::generateDataClients(const KeysServer &server) {
//    int uniquePointsNum = 3 + rand() % 10, clientsNum = uniquePointsNum;
//    //  init coordinate arrays
//    long arrs[uniquePointsNum][DIM];
//    for (auto &arr: arrs) for (int dim = 0; dim < DIM; ++dim) arr[dim] = rand() % NUMBERS_RANGE;
//    /*  init clients vector
//     *      client [0] stays empty
//     *      client [1] has 1 point - {point0}
//     *      client [2] has 2 points - {point0, point1}
//     *      ...
//     *      client [n] has n points -  {point0, point1, ... , pointN}
//     */
//    std::vector<Client> clients(clientsNum, Client(server));
//    for (int i = 1; i < clients.size(); ++i)
//        for (int j = 0; j < i; ++j)
//            clients[i].encryptPoint(arrs[j]);
//    return clients;
//}




int main2() {
    cout << "Hello, World!" << endl;
//    BGV_binary_arithmetic();
//    BGV_packed_arithmetic();
//    KT_packed_arithmetic();
//    add_numbers();
//    mult_numbers();
//    cmp_numbers();


    return 0;
}

constexpr static long mValues[8][15] = {
        // clang-format off
        // {p,phi(m),     m,  d,   m1,  m2,  m3,   g1,   g2,   g3,ord1,ord2,ord3,  B, c}
        {2,    48,   105, 12,    3,  35,   0,   71,   76,    0,   2,   2,   0, 25, 2},
        {2,   600,  1023, 10,   11,  93,   0,  838,  584,    0,  10,   6,   0, 25, 2},
        {2,  2304,  4641, 24,    7,   3, 221, 3979, 3095, 3760,   6,   2,  -8, 25, 3},
        {2,  5460,  8193, 26, 8193,   0,   0,   46,    0,    0, 210,   0,   0, 25, 3},
        {2,  8190,  8191, 13, 8191,   0,   0,   39,    0,    0, 630,   0,   0, 25, 3},
        {2, 10752, 11441, 48,   17, 673,   0, 4712, 2024,    0,  16, -14,   0, 25, 3},
        {2, 15004, 15709, 22,   23, 683,   0, 4099,13663,    0,  22,  31,   0, 25, 3},
        {2, 27000, 32767, 15,   31,   7, 151,11628,28087,25824,  30,   6, -10, 28, 4}
        // clang-format on
};

void add_numbers() {

}



void KT_packed_arithmetic() {
    /*  Example of BGV scheme  */

    unsigned long p = 4999;    // Plaintext prime modulus
    unsigned long m = 32109;    // Cyclotomic polynomial - defines phi(m)
    unsigned long r = 1;    // Hensel lifting (default = 1)
    unsigned long bits = 500;    // Number of bits of the modulus chain
    unsigned long c = 2;    // Number of columns of Key-Switching matrix (default = 2 or 3)

    cout << "Initialising context object..." << endl;
    // Initialize context
    // This object will hold information about the algebra created from the previously set parameters
    helib::Context context = helib::ContextBuilder<helib::BGV>().m(m).p(p).r(r).bits(bits).c(c).build();

    // Print the context
    context.printout();
    cout << endl;
    // Print the security level
    cout << "Security: " << context.securityLevel() << endl;

    // Secret key management
    cout << "Creating secret key..." << endl;
    helib::SecKey secret_key(context);    // Create a secret key associated with the context
    secret_key.GenSecKey();    // Generate the secret key
    cout << "Generating key-switching matrices..." << endl;
    helib::addSome1DMatrices(secret_key);    // Compute key-switching matrices that we need


    // Public key management
    const helib::PubKey &public_key = secret_key;    // Set the secret key (upcast: SecKey is a subclass of PubKey)

    const helib::EncryptedArray &ea = context.getEA();    // Get the EncryptedArray of the context

    long nslots = ea.size();    // Get the number of slot (phi(m))
    cout << "Number of slots: " << nslots << endl;

    // Create a vector of long with nslots elements
    helib::Ptxt<helib::BGV> ptxt(context);
    // Set it with numbers 0..nslots - 1
    // ptxt = [13] [0] [0] [0] ... [0]
    ptxt[0] = 3;
    for (int i = 1; i < ptxt.size(); ++i) {
        ptxt[i] = 0;
    }

    // Print the plaintext
    cout << "Initial Plaintext: " << ptxt << endl;

    // Create a ciphertext object
    helib::Ctxt ctxt(public_key);
    // Encrypt the plaintext using the public_key
    public_key.Encrypt(ctxt, ptxt);

    // Print the Ciphertext
    helib::Ptxt<helib::BGV> decrypted_ctxt(context);    // Create a plaintext for decryption
    secret_key.Decrypt(decrypted_ctxt, ctxt);    // Decrypt the modified ciphertext
    cout << "Initial Ciphertext: " << decrypted_ctxt << endl;

    /********** Operations **********/
    // Ciphertext and plaintext operations are performed
    // "entry-wise".

    // Square the ciphertext
    // [13] [0] [0] [0] ... [0]
    // -> [169] [0] [0] [0] ... [0]
    ctxt.multiplyBy(ctxt);
    // Plaintext version
    ptxt.multiplyBy(ptxt);

    // Print the Ciphertext
    secret_key.Decrypt(decrypted_ctxt, ctxt);    // Decrypt the modified ciphertext
    cout << "ctxt.multiplyBy(ctxt): " << decrypted_ctxt << endl;

    // Divide the ciphertext by itself
    // To do this we must calculate the multiplicative inverse using Fermat's
    // Little Theorem.  We calculate a^{-1} = a^{p-2} mod p, where a is non-zero
    // and p is our plaintext prime.
    // First make a copy of the ctxt using copy constructor
    helib::Ctxt ctxt_divisor(ctxt);
    // Raise the copy to the exponent 2
    // [0] [1] [4] ... [16] -> [0] [1] [1] ... [1]
    // Note: 0 is a special case because 0^n = 0 for any power n
    ctxt_divisor.power(2);
    // Print the Ciphertext
    secret_key.Decrypt(decrypted_ctxt, ctxt);    // Decrypt the modified ciphertext
    cout << "ctxt_divisor.power(2): " << decrypted_ctxt << endl;

    // a^{2}*a = a^{-1}*a = a / a = 1;
    ctxt.multiplyBy(ctxt_divisor);
    // Print the Ciphertext
    secret_key.Decrypt(decrypted_ctxt, ctxt);    // Decrypt the modified ciphertext
    cout << "ctxt.multiplyBy(ctxt_divisor): " << decrypted_ctxt << endl;

    // Plaintext version
    helib::Ptxt<helib::BGV> ptxt_divisor(ptxt);
    ptxt_divisor.power(2);
    ptxt.multiplyBy(ptxt_divisor);
    cout << "ptxt.multiplyBy(ptxt_divisor): " << ptxt << endl;

    cout << endl << endl;

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

void BGV_packed_arithmetic() {

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


void BGV_binary_arithmetic() {/*  Example of binary arithmetic using the BGV scheme  */

    // First set up parameters.
//
// NOTE: The parameters used in this example code are for demonstration only.
// They were chosen to provide the best performance of execution while
// providing the context to demonstrate how to use the "Binary Arithmetic
// APIs". The parameters do not provide the security level that might be
// required by real use/application scenarios.

    // Plaintext prime modulus.
    long p = 2;
    // Cyclotomic polynomial - defines phi(m).
    long m = 4095;
    // Hensel lifting (default = 1).
    long r = 1;
    // Number of bits of the modulus chain.
    long bits = 500;
    // Number of columns of Key-Switching matrix (typically 2 or 3).
    long c = 2;
    // Factorisation of m required for bootstrapping.
    std::vector<long> mvec = {7, 5, 9, 13};
    // Generating set of Zm* group.
    std::vector<long> gens = {2341, 3277, 911};
    // Orders of the previous generators.
    std::vector<long> ords = {6, 4, 6};

    cout << "\n*********************************************************";
    cout << "\n*            Basic Binary Arithmetic Example            *";
    cout << "\n*            ===============================            *";
    cout << "\n*                                                       *";
    cout << "\n* This is a sample program for education purposes only. *";
    cout << "\n* It attempts to demonstrate the use of the API for the *";
    cout << "\n* binary arithmetic operations that can be performed.   *";
    cout << "\n*                                                       *";
    cout << "\n*********************************************************";
    cout << endl;

    cout << "Initialising context object..." << endl;
    // Initialize the context.
// This object will hold information about the algebra created from the
// previously set parameters.
    helib::Context context = helib::ContextBuilder<helib::BGV>()
            .m(m)
            .p(p)
            .r(r)
            .gens(gens)
            .ords(ords)
            .bits(bits)
            .c(c)
            .bootstrappable(true)
            .mvec(mvec)
            .build();

    // Print the context.
    context.printout();
    cout << endl;

    // Print the security level.
    cout << "Security: " << context.securityLevel() << endl;

    // Secret key management.
    cout << "Creating secret key..." << endl;
    // Create a secret key associated with the context.
    helib::SecKey secret_key(context);
    // Generate the secret key.
    secret_key.GenSecKey();

    // Generate bootstrapping data.
    secret_key.genRecryptData();

    // Public key management.
// Set the secret key (upcast: SecKey is a subclass of PubKey).
    const helib::PubKey &public_key = secret_key;

    // Get the EncryptedArray of the context.
    const helib::EncryptedArray &ea = context.getEA();

    // Build the unpack slot encoding.
    std::vector<helib::zzX> unpackSlotEncoding;
    buildUnpackSlotEncoding(unpackSlotEncoding, ea);

    // Get the number of slot (phi(m)).
    long nslots = ea.size();
    cout << "Number of slots: " << nslots << endl;

    // Generate three random binary numbers a, b, c.
// Encrypt them under BGV.
// Calculate a * b + c with HElib's binary arithmetic functions, then decrypt
// the result.
// Next, calculate a + b + c with HElib's binary arithmetic functions, then
// decrypt the result.
// Finally, calculate popcnt(a) with HElib's binary arithmetic functions,
// then decrypt the result.  Note that popcnt, also known as hamming weight
// or bit summation, returns the count of non-zero bits.

    // Each bit of the binary number is encoded into a single ciphertext. Thus
// for a 16 bit binary number, we will represent this as an array of 16
// unique ciphertexts.
// i.e. b0 = [0] [0] [0] ... [0] [0] [0]        ciphertext for bit 0
//      b1 = [1] [1] [1] ... [1] [1] [1]        ciphertext for bit 1
//      b2 = [1] [1] [1] ... [1] [1] [1]        ciphertext for bit 2
// These 3 ciphertexts represent the 3-bit binary number 110b = 6

    // Note: several numbers can be encoded across the slots of each ciphertext
// which would result in several parallel slot-wise operations.
// For simplicity we place the same data into each slot of each ciphertext,
// printing out only the back of each vector.
// NB: fifteenOrLess4Four max is 15 bits. Later in the code we pop the MSB.
    long bitSize = 16;
    long outSize = 2 * bitSize;
    long a_data = NTL::RandomBits_long(bitSize);
    long b_data = NTL::RandomBits_long(bitSize);
    long c_data = NTL::RandomBits_long(bitSize);

    cout << "Pre-encryption data:" << endl;
    cout << "a = " << a_data << endl;
    cout << "b = " << b_data << endl;
    cout << "c = " << c_data << endl;

    // Use a scratch ciphertext to populate vectors.
    helib::Ctxt scratch(public_key);
    std::vector<helib::Ctxt> encrypted_a(bitSize, scratch);
    std::vector<helib::Ctxt> encrypted_b(bitSize, scratch);
    std::vector<helib::Ctxt> encrypted_c(bitSize, scratch);
    // Encrypt the data in binary representation.
    for (long i = 0; i < bitSize; ++i) {
        std::vector<long> a_vec(ea.size());
        std::vector<long> b_vec(ea.size());
        std::vector<long> c_vec(ea.size());
        // Extract the i'th bit of a,b,c.
        for (auto &slot : a_vec)
            slot = (a_data >> i) & 1;
        for (auto &slot : b_vec)
            slot = (b_data >> i) & 1;
        for (auto &slot : c_vec)
            slot = (c_data >> i) & 1;
        ea.encrypt(encrypted_a[i], public_key, a_vec);
        ea.encrypt(encrypted_b[i], public_key, b_vec);
        ea.encrypt(encrypted_c[i], public_key, c_vec);
    }

    // Although in general binary numbers are represented here as
// std::vector<helib::Ctxt> the binaryArith APIs for HElib use the PtrVector
// wrappers instead, e.g. helib::CtPtrs_vectorCt. These are nothing more than
// thin wrapper classes to standardise access to different vector types, such
// as NTL::Vec and std::vector. They do not take ownership of the underlying
// object but merely provide access to it.
//
// helib::CtPtrMat_vectorCt is a wrapper for
// std::vector<std::vector<helib::Ctxt>>, used for representing a list of
// encrypted binary numbers.

    // Perform the multiplication first and put it in encrypted_product.
    std::vector<helib::Ctxt> encrypted_product;
    helib::CtPtrs_vectorCt product_wrapper(encrypted_product);
    helib::multTwoNumbers(
            product_wrapper,
            helib::CtPtrs_vectorCt(encrypted_a),
            helib::CtPtrs_vectorCt(encrypted_b),
            /*rhsTwosComplement=*/false, // This means the rhs is unsigned rather
            // than 2's complement.
            outSize, // Outsize is the limit on the number of bits in the output.
            &unpackSlotEncoding); // Information needed for bootstrapping.

    // Now perform the encrypted sum and put it in encrypted_result.
    std::vector<helib::Ctxt> encrypted_result;
    helib::CtPtrs_vectorCt result_wrapper(encrypted_result);
    helib::addTwoNumbers(
            result_wrapper,
            product_wrapper,
            helib::CtPtrs_vectorCt(encrypted_c),
            /*negative=*/false, // This means the number are unsigned rather than 2's
            // complement.
            &unpackSlotEncoding); // Information needed for bootstrapping.

    // Decrypt and print the result.
    std::vector<long> decrypted_result;
    helib::decryptBinaryNums(decrypted_result, result_wrapper, secret_key, ea);
    cout << "a*b+c = " << decrypted_result.back() << endl;

    // Now calculate the sum of a, b and c using the addManyNumbers function.
    encrypted_result.clear();
    decrypted_result.clear();
    std::vector<std::vector<helib::Ctxt>> summands = {encrypted_a,
                                                      encrypted_b,
                                                      encrypted_c};
    helib::CtPtrMat_vectorCt summands_wrapper(summands);
    helib::addManyNumbers(
            result_wrapper,
            summands_wrapper,
            0,                    // sizeLimit=0 means use as many bits as needed.
            &unpackSlotEncoding); // Information needed for bootstrapping.

    // Decrypt and print the result.
    helib::decryptBinaryNums(decrypted_result, result_wrapper, secret_key, ea);
    cout << "a+b+c = " << decrypted_result.back() << endl;

    // This section calculates popcnt(a) using the fifteenOrLess4Four
// function.
// Note: the output i.e. encrypted_result should be of size 4
// because 4 bits are required to represent numbers in [0,15].
    encrypted_result.resize(4lu, scratch);
    decrypted_result.clear();
    encrypted_a.pop_back(); // drop the MSB since we only support up to 15 bits.
    helib::fifteenOrLess4Four(result_wrapper,
                              helib::CtPtrs_vectorCt(encrypted_a));

    // Decrypt and print the result.
    helib::decryptBinaryNums(decrypted_result, result_wrapper, secret_key, ea);
    cout << "popcnt(a) = " << decrypted_result.back() << endl;
}