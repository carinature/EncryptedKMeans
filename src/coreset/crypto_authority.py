import numpy as np
from phe import paillier


class CA:
    def __init__(self, n_length=1024):
        self.secure = n_length > 0
        if self.secure:
            self.public_key, self.__private_key = paillier.generate_paillier_keypair(n_length=n_length)
        else:
            self.public_key = None

    def decrypt(self, enc_vec):
        if self.secure:
            return np.array([self.__private_key.decrypt(x) for x in enc_vec])
        else:
            return enc_vec
