import os

import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Pool

PATH = os.getcwd() + r'\fhe_tests'
RGB = ['#ff0000', '#ff8000', '#ffff00', '#00ff00', '#00ffff', '#0000ff', '#8000ff', '#000000']


def cartesian_to_polar(p: np.ndarray):
    '''
    Translates from cartesian coordinates to spherical coordinates

    :param p: a point
    :return: angles representing p (all between 0 and π, except last one which is between 0 and 2π)
    '''
    d = len(p)
    angles = np.zeros(d - 1)
    for i in range(d - 1):
        if np.allclose(p[i:],
                       np.zeros(d - i)):  # if the remaining coordinates are 0 than the angle is ambiguous - we choose 0
            angles[i] = 0
        else:
            angles[i] = np.arccos(p[i] / np.linalg.norm(p[i:]))
    if np.round(p[-1], 8) < 0:
        angles[-1] = 2 * np.pi - angles[-1]
    return angles


def polar_to_cartesian(angles: np.ndarray):
    """
    Translates from cartesian coordinates to spherical coordinates

    :param angles: vector of angles (spherical coordinates)
    :return: unit vector @angles represents, in a cartesian coordinates system
    """
    # From https://stackoverflow.com/questions/20133318/n-sphere-coordinate-system-to-cartesian-coordinate-system
    a = np.concatenate((np.array([2 * np.pi]), angles))
    si = np.sin(a)
    si[0] = 1
    si = np.cumprod(si)
    co = np.cos(a)
    co = np.roll(co, -1)
    return si * co


def draw_sun(teta, c, r=1):
    ang = 0
    while ang < 2 * np.pi:
        plt.plot([c[0], r * np.cos(ang) + c[0]], [c[1], r * np.sin(ang) + c[1]], 'y')
        plt.plot([c[0], r * np.cos(ang + teta / 2.0) + c[0]], [c[1], r * np.sin(ang + teta / 2.0) + c[1]], 'y:')
        ang += teta


def rgb(i: int):
    if i >= len(RGB):
        return 'k'
    elif i < 0:
        return 'r'
    else:
        return RGB[int(i)]


def plot_coreset(coreset, teta, mean):
    draw_sun(teta, mean, 1)
    for i, c in enumerate(RGB):
        plt.plot(i * 0.05, 1.3, color=c, marker='$' + str(i) + '$')
    for p, w in coreset:
        plt.plot(p[0], p[1], color=rgb(w), marker='o')
    plt.gca().set_aspect('equal', adjustable='box')


def random_in_unit_ball(d):
    return random_in_ball(d, np.zeros(d), 1)


def random_in_ball(d, c, r):
    while True:
        p = np.random.uniform(-r, r, d)
        if np.linalg.norm(p) <= 1:
            return c + p


def random_set_in_unit_ball(n, d):
    return np.array([random_in_unit_ball(d) for _ in range(n)])


def random_ON_unit_ball(d):
    p = np.random.uniform(-1, 1, d)
    p = p / np.linalg.norm(p)
    return p


def random_set_ON_unit_ball(n, d):
    return np.array([random_ON_unit_ball(d) for _ in range(n)])


def ignore():
    print()
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # t = 0.25*np.pi
    # # abc = np.array([(np.pi-t)/2, (np.pi-t)/2])
    # # abcd = np.array([(np.pi)/2, (np.pi-t)/2])
    # # vv = polar_to_cartesian(abc)*2
    # # vvv = polar_to_cartesian(abcd)*2
    # # # ax.plot([0 * -vvv[0], vvv[0]], [0 * -vvv[1], vvv[1]], [0 * -vvv[2], vvv[2]], 'b')
    # #
    # # vv =
    # vv = np.zeros(3)
    # c = 0
    # for ang1 in np.arange(0, np.pi, t):
    #     for ang2 in np.arange(0, np.pi, t):
    #         angs = np.array([ang1, ang2])
    #         v = polar_to_cartesian(angs).round(7)
    #         vv += v
    #         c += 1
    # vv /= c
    # print(vv.round(7))
    # for ang1 in np.arange(0, np.pi, t):
    #     vvv = np.zeros(3)
    #     c = 0
    #     for ang2 in np.arange(0, np.pi, t):
    #         # for ang3 in np.arange(0, np.pi, t):
    #         #     for ang4 in np.arange(0, np.pi, 0.25*np.pi):
    #         angs = np.array([ang1, ang2])
    #         v = polar_to_cartesian(angs).round(7)
    #         ax.plot([0*-v[0], v[0]], [0*-v[1], v[1]], [0*-v[2], v[2]])
    #         # print(np.rad2deg(angs), np.rad2deg(aaa).round(5))
    #         if (v @ vv).round(7) <= 0:
    #             print('err', (v @ vv).round(7), v)
    #         else:
    #             print('ok', (v @ vv).round(7))
    #         c += 1
    #     vvv /= c
    #     ax.plot([0 * -vvv[0], vvv[0]], [0 * -vvv[1], vvv[1]], [0 * -vvv[2], vvv[2]], 'b')
    #
    #     print('---------------')
    #
    # ax.plot([0 * -vv[0], vv[0]], [0 * -vv[1], vv[1]], [0 * -vv[2], vv[2]], 'k')
    # plt.show()
    # # for i in range(100):
    # #     v = random_ON_unit_ball(np.random.randint(2, 16))
    # #     ang = cartesian_to_polar(v)
    # #     # print(ang)
    # #     ang[0] += np.pi
    # #     vr = polar_to_cartesian(ang)
    # #     # print(np.nonzero((v+vr).round(7))[0])
    # #     print(np.abs((v+vr).round(7)))
    # # d = 4
    # # threshold = np.ones(d - 1) * 0.5 * np.pi
    # # print(threshold)


if __name__ == '__main__':
    v = np.array([1, 0, 0, 0, 0])
    print(cartesian_to_polar(v))
    print(cartesian_to_polar(-v))
