from src.coreset.utils import *


class Server:
    def __init__(self, n: int, d: int, eps: float, alpha: float, delta: float):
        """
        Initialization of a Server object

        :param n: number of clients
        :param d: number of dimensions/features
        :param eps: multiplicative error
        :param alpha: privacy parameter
        :param delta: failure probability
        """
        self.n = n
        self.d = d
        self.eps = eps
        self.alpha = alpha
        self.delta = delta
        self.teta = np.pi / np.ceil(np.pi / (2 * self.eps))
        self.params = {'n': self.n,
                       'd': self.d,
                       'eps': self.eps,
                       'alpha': self.alpha,
                       'delta': self.delta,
                       'teta': self.teta}

    def __init__original(self, n: int, d: int, eps: float, alpha: float, delta: float):
        """
        Initialization of a Server object

        :param n: number of clients
        :param d: number of dimensions/features
        :param eps: multiplicative error
        :param alpha: privacy parameter
        :param delta: failure probability
        """
        self.n = n
        self.d = d
        self.eps = eps
        self.alpha = alpha
        self.delta = delta
        self.context_coreset_original()

        self.teta = self.params['teta']
        self.r0 = self.params['r0']
        self.t = self.params['t']

        self.collective_mean = None

    # def context_coreset_original(self) -> None:
    #     """
    #     Creates sets the parameters dictionary
    #
    #     :return: None
    #     """
    #     teta = np.pi / np.ceil(np.pi / (2 * self.eps))
    #     n0 = np.log(1.0 / (self.delta * self.eps ** (self.d - 1))) / (self.alpha * self.eps ** (self.d + 1))
    #     r0 = np.log(1.0 / self.delta) / (self.alpha * n0)
    #     t = int(np.ceil(np.emath.logn(1 + self.eps, 1 / r0)))  # interval count
    #     self.params = {'n': self.n,
    #                    'd': self.d,
    #                    'eps': self.eps,
    #                    'alpha': self.alpha,
    #                    'delta': self.delta,
    #                    'teta': teta,
    #                    'r0': r0,
    #                    't': t}

    def sum_ciphers(self, vectors_list: list, pool: Pool = None) -> np.ndarray:
        """
        Accumulates the cipher-texts in parallel (if possible)

        :param vectors_list: list of encrypted vectors
        :param pool: processes pool for paralleling
        :return: sum vector of the inputs (encrypted)
        """
        # if pool is not None:
        #     processes = pool._processes
        #     l = np.array_split(vectors_list, processes)
        #     s = np.sum(pool.map(np.sum, ))
        # else:
        return np.sum(vectors_list, axis=0)

    def set_collective_mean(self, mean: np.ndarray) -> None:
        """
        Sets the mean

        :param mean: mean of the points
        :return: None
        """
        self.collective_mean = mean
        return

    def decode_to_coreset(self, weight_vector: np.ndarray) -> list:
        C = []
        interval_amnt = 2 * self.t.max() + 3
        print('---', len(weight_vector), '---')
        for i in range(len(weight_vector)):
            l_idx = i // interval_amnt
            interval = i % interval_amnt
            if interval == interval_amnt - 1:
                q = 0
                w = -np.sum(weight_vector[l_idx * interval_amnt: i]) + weight_vector[i] + self.n0[l_idx]
            elif interval % (self.t.max() + 1) > self.t[l_idx]:
                print('continue', l_idx, interval)
                continue
            elif interval < self.t.max() + 1:
                q = self.unitLineUnFlat(interval, l_idx)
                w = weight_vector[i]
            else:
                q = -self.unitLineUnFlat(interval - (self.t.max() + 1), l_idx)
                w = weight_vector[i]
            v = self.unitBallUnFlat(l_idx)
            p = v * (q + self.means[l_idx]) + self.collective_mean
            C.append((p, w))
            print(C[-1])
        return C

    # def decode_to_coreset_original(self, weight_vector: np.ndarray) -> list:
    #     # TODO Z CENTERED MUST FIX
    #     total_noise = np.sum(weight_vector[:-1]) - self.n
    #     z_weight = weight_vector[-1] - total_noise
    #     C = [(self.collective_mean, z_weight)]
    #     for i in range(len(weight_vector)):
    #         l = i // (self.t + 1)
    #         interval = i % (self.t + 1)
    #         v = self.unitBallUnFlat(l)
    #         q = self.unitLineUnFlat(interval)
    #         p = v * q + self.collective_mean
    #         C.append((p, weight_vector[i]))
    #     return C

    def unitLineUnFlat(self, i: int, line_idx: int) -> float:
        return min(self.r0[line_idx] * (1 + self.eps) ** i, 1)  # r(i)

    def unitBallUnFlat(self, l_idx: int) -> np.ndarray:
        """
        Translates a index of the expanded vector to its represented unit vector

        :param l_idx: the index
        :return: unit vector
        """
        line_amnt = np.round(np.pi / self.teta).astype(int)
        line = np.array(
            [(l_idx // (line_amnt ** k)) % line_amnt for k in range(self.d - 1)])
        line = line * self.teta
        return polar_to_cartesian(line)

    # def unitBallUnFlat_original(self, l_idx: int) -> np.ndarray:
    #     """
    #     Translates a index of the expanded vector to its represented unit vector
    #
    #     :param l_idx: the index
    #     :return: unit vector
    #     """
    #     line_amnt = np.round(np.pi / self.teta).astype(int)
    #     line = np.array(
    #         [(l_idx // line_amnt ** k) % (2 * line_amnt if k == 0 else line_amnt + 1) for k in range(self.d - 1)])
    #     line = line * self.teta
    #     return polar_to_cartesian(line)

    def decode_means_and_params(self, means_vector: np.ndarray) -> (np.ndarray, np.ndarray):
        self.means = np.zeros(len(means_vector) // 2)
        self.n0 = np.zeros(len(means_vector) // 2)
        for i in range(len(self.means)):
            s, n = means_vector[2 * i], means_vector[2 * i + 1]
            self.means[i] = s / n if n != 0 else 0
            self.n0[i] = max(n, np.log(1 / self.delta) / self.alpha)
        self.r0 = np.log(1 / self.delta) / (self.alpha * self.n0)
        self.t = np.ceil(np.emath.logn(1 + self.eps, 1 / self.r0)).astype(int)
        return self.means, self.r0, self.t
