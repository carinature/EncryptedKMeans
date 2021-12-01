# from pyparsing import delimitedList

import argparse

# from pies.yonis.client import Client
# from pies.yonis.crypto_authority import CA
# from pies.yonis.server import Server
# from pies.yonis.utils import *

# from .client import Client
# from .crypto_authority import CA
# from .server import Server
# from .utils import *
from src.yonis.client import Client
from src.yonis.crypto_authority import CA
from src.yonis.server import Server
from src.yonis.utils import *


class Simulator:
    def __init__(self, P: np.ndarray, eps: float, alpha: float, delta: float, private: bool = True,
                 security: int = 1024, processes: int = None):
        """
        Initialization of the simulator

        :param P: input points
        :param eps: epsilon - multiplicative error
        :param alpha: privacy parameter
        :param delta: failure probability
        :param private: indicates if noise should be added
        :param security: security parameter
        :param processes: number of processes in the pool (None will result in the number of CPUs in the PC)
        """
        self.pool = Pool(processes)
        n, d = P.shape
        self.ca = CA(security)
        self.server = Server(n, d, eps, alpha, delta)
        self.clients = [Client(p, self.server.params, self.ca.public_key, private) for p in P]
        self.private = private

    def mean_round(self) -> np.ndarray:
        """
        Simulates the mean round

        :return: the noisy mean
        """
        points_ciphers = [client.send_point(self.pool) for client in self.clients]
        enc_noisy_mean = self.server.sum_ciphers(points_ciphers)
        return self.ca.decrypt(enc_noisy_mean)

    def lines_means_round(self) -> np.ndarray:
        points_ciphers = [client.send_line_encoded(self.pool) for client in self.clients]
        enc_noisy_lines_means = self.server.sum_ciphers(points_ciphers)
        return self.ca.decrypt(enc_noisy_lines_means)

    def sum_round(self) -> np.ndarray:
        """
        Simulates the sum of weights round

        :return: weights vector
        """
        encoded_ciphers = [client.send_encoded(self.pool) for client in self.clients]
        enc_noisy_sum_vector = self.server.sum_ciphers(encoded_ciphers)
        return self.ca.decrypt(enc_noisy_sum_vector)

    def coreset(self, plot: bool = False) -> list:
        """
        Simulates full coreset protocol

        :param plot: Plots the coreset if plot == True
        :return: weighted Coreset
        """
        noisy_mean = self.mean_round()
        # noisy_mean = np.zeros(self.server.d)
        for client in self.clients:
            client.set_mean(noisy_mean)
        self.server.set_collective_mean(noisy_mean)

        noisy_lines_means_vector = self.lines_means_round()
        means, r0, t = self.server.decode_means_and_params(noisy_lines_means_vector)
        for client in self.clients:
            client.set_line_mean_and_params(means, r0, t)

        # plot_coreset([], self.server.params['teta'], self.server.collective_mean)

        noisy_sum_vector = self.sum_round()
        coreset = self.server.decode_to_coreset(noisy_sum_vector)
        if plot:
            plot_coreset(coreset, self.server.params['teta'], self.server.collective_mean)
        return coreset


def main():
    parser = argparse.ArgumentParser(description='Run 1-mean coreset')
    parser.add_argument('points', type=str, help='File to read input points')
    parser.add_argument('eps', type=float, help='Epsilon accuracy value')
    parser.add_argument('alpha', type=float, help='Alpha privacy value')
    parser.add_argument('delta', type=float, help='Delta failure probability')
    parser.add_argument('-n', '--noprivacy', default=False, action='store_true')
    parser.add_argument('-s', '--security', type=int, default=1024, help='Security parameter (0 means no encryption)')
    parser.add_argument('-p', '--plot', default=False, action='store_true')
    parser.add_argument('-f', '--file', type=str, help='File to write coreset')
    args = parser.parse_args()
    print(args)
    P = np.genfromtxt(args.points, delimiter=',')
    sim = Simulator(P, args.eps, args.alpha, args.delta, not args.noprivacy, args.security)
    core = sim.coreset(args.plot)
    core_table = np.array([np.append(p, w) for p, w in core])
    if args.plot:
        plt.show()

        plt.savefig("/home/rbd/workspace/rbd/rbd_helib_with_remote_debugger/yonis_res.png")

    if args.file is not None:
        np.savetxt(args.file, core_table, delimiter=',')
    else:
        print('Coreset:')
        print(core_table)


if __name__ == '__main__':
    main()
