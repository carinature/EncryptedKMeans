# not sure which one of the two (KT 13apr)
# from multiprocessing.dummy import Pool

from phe import paillier

from src.coreset.utils import *


class Client:
    def __init__(self, p: np.ndarray, param: dict, public_key: paillier.PaillierPublicKey, private: bool = True):
        """
        Initialize Client object

        :param p: data-point the client possesses
        :param param: parameters dictionary provided by the server
        :param public_key: public key for encryption provided by C.A.
        :param private: indicates if noise should be added
        """
        self.p = p
        self.public_key = public_key
        self.n = param['n']
        self.d = param['d']
        self.alpha = param['alpha']
        self.eps = param['eps']
        self.teta = param['teta']
        self.private = private
        self.line_amnt = np.round(np.pi / self.teta).astype(int)  # round is just for safety

    # def encode_original(self) -> np.ndarray:
    #     """
    #     Encodes the client's point to the expanded histogram vector
    #
    #     :return: a one-hot vector
    #     """
    #     normalized_p = self.p - self.collective_mean
    #     angles = cartesian_to_polar(normalized_p)
    #     line = self.teta * np.round(angles / self.teta)
    #
    #     direction = polar_to_cartesian(line)
    #     projection = direction @ normalized_p
    #     interval = self.find_interval(projection)
    #
    #     angles_idx = np.round(angles / self.teta).astype(int)
    #     angles_idx[:-1] %= self.line_amnt + 1
    #     angles_idx[-1] %= 2 * self.line_amnt
    #     flatter = (self.line_amnt+1)**np.arange(self.d - 1)
    #     line_num = flatter @ angles_idx
    #
    #     i = int(line_num * (self.t + 1) + interval)
    #     v = np.zeros((self.t + 1) * 2 * self.line_amnt * (self.line_amnt+1) ** (self.d - 2) + 1, dtype=int)
    #     v[i] = 1
    #     return v

    def encode_line_mean(self) -> np.ndarray:  # TODO CHANGE DOCUMENTATION
        """
        Encodes the client's point to the expanded histogram vector

        :return: a one-hot vector
        """
        normalized_p = self.p - self.collective_mean
        angles = cartesian_to_polar(normalized_p)
        line = (self.teta * np.round(angles / self.teta)) % (2 * np.pi)
        if np.any(line >= np.pi):
            line[0] += np.pi
            line[0] %= 2 * np.pi

        self.direction = polar_to_cartesian(line)
        line = cartesian_to_polar(
            self.direction)  # easiest way to make translate line with angles>180, to the same line with all angles<180

        angles_idx = np.round(line / self.teta).astype(int)
        flatter = self.line_amnt ** np.arange(self.d - 1)
        self.line_num = flatter @ angles_idx
        self.projection = self.direction @ normalized_p
        if self.d == 2:  # TODO PLOTS
            projected_p = self.collective_mean + self.direction * self.projection
            plt.plot(projected_p[0], projected_p[1], 'k*')

        i = 2 * self.line_num
        v = np.zeros(2 * self.line_amnt ** (self.d - 1), dtype=float)
        v[i] = self.projection
        v[i + 1] = 1
        return v

    def encode(self) -> np.ndarray:
        normalized_projected_p = self.projection - self.line_mean
        abs_proj = abs(normalized_projected_p)
        abs_interval = self.find_interval(abs_proj)
        if normalized_projected_p > 0:
            interval = abs_interval
        else:
            interval = self.max_t + 1 + abs_interval
        self.interval_amnt = (2 * self.max_t + 3)
        i = int(self.line_num * self.interval_amnt + interval)
        v = np.zeros(self.interval_amnt * self.line_amnt ** (self.d - 1), dtype=int)
        v[i] = 1
        return v

    def find_interval(self, p: float) -> int:
        """
        Returns the index of the interval which p is on

        :param p: 0 <= p <= 1, distance of a point on the line
        :return: index of interval, between 0 and t
        """
        if abs(p) < self.r0:
            return 0
        elif abs(p) > 1:
            return self.t
        else:
            ci = np.round(np.emath.logn(1 + self.eps, abs(p) / self.r0)).astype(int)
            if ci > self.t:  # TODO CHANGE TO MIN
                return self.t
            else:
                return ci

    def encrypt_no_parallel(self, l: np.ndarray) -> np.ndarray:
        """
        Encrypts the given values in serial

        :param l: array of values to encrypt
        :return: array of encrypted values
        """
        return np.array([self.public_key.encrypt(float(x)) for x in l])

    def encrypt_vector(self, vec: np.ndarray, pool: Pool = None) -> np.ndarray:
        """
        Encrypts the given values in parallel (if possible)

        :param vec: array of values to encrypt
        :param pool: processes pool for paralleling
        :return: array of encrypted values
        """
        if self.public_key is None:
            return vec
        elif pool is None:
            return self.encrypt_no_parallel(vec)
        else:
            return np.concatenate(pool.map(self.encrypt_no_parallel, np.array_split(vec, pool._processes)))

    def get_noise(self, scale: float, size: int) -> np.ndarray:
        """
        Generates noise vector drawn from difference between Gamma distributions

        :param scale: scale of the gamma distribution
        :param size: length of noise vector
        :return:
        """
        if self.private:
            # return np.random.gamma(1/self.n, scale, size) - np.random.gamma(1/self.n, scale, size)
            return np.random.normal(0, scale / (self.n ** 0.5), size)
        else:
            return np.zeros(size)

    def send_point(self, pool: Pool = None) -> np.ndarray:
        """
        Sends the point divided by n, encrypted and noised (used to calculate secure mean)

        :param pool: processes pool for paralleling
        :return: array of encrypted (noised) values
        """
        to_send = self.p / self.n
        noise = self.get_noise(self.d / (self.n * self.alpha), self.d)
        return self.encrypt_vector(to_send + noise, pool)

    def send_line_encoded(self, pool: Pool = None) -> np.ndarray:
        l_encoded = self.encode_line_mean()
        noise = self.get_noise(1 / self.alpha, l_encoded.size)
        return self.encrypt_vector(l_encoded + noise, pool)

    def send_encoded(self, pool: Pool = None) -> np.ndarray:
        """
        Sends the encoding of the point encrypted and noised

        :param pool: processes pool for paralleling
        :return: array of encrypted (noised) values
        """
        encoded = self.encode()
        noise = self.get_noise(1 / self.alpha, encoded.size)
        return self.encrypt_vector(encoded + noise, pool)

    def set_mean(self, noisy_mean: np.ndarray) -> None:
        """
        Sets the mean of all the points calculated in the server

        :param noisy_mean: the (noisy) mean
        :return: None
        """
        self.collective_mean = np.array(noisy_mean)

    def set_line_mean_and_params(self, line_means: np.ndarray, r0: np.ndarray, t: np.ndarray) -> None:
        """
        Sets the mean the class calculated in the server

        :param line_mean: the (noisy) mean of the points on the specific line
        :return: None
        """
        self.line_mean = line_means[self.line_num]
        self.r0 = r0[self.line_num]
        self.t = t[self.line_num]
        self.max_t = t.max()
        if self.d == 2:  # TODO PLOTS
            line_mean_d = self.collective_mean + self.direction * self.line_mean
            plt.plot(line_mean_d[0], line_mean_d[1], 'bx')
