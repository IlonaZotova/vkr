import numpy as np
import tqdm
import json
import argparse
import math
from scipy.integrate import odeint
from NeProfile import NeProfile
from scipy.optimize._tstutils import functions

from special_functions import *
import matplotlib.pyplot as plt

profile = NeProfile().linear(0, 1000)

def P_init_states(a, b, r_init, alpha, ):
    n = a + b * r_init
    P_r = n * np.cos(alpha)
    P_tetta = n * np.sin(alpha)
    return [P_r, P_tetta]

# def system(functions, tao, a, model):
   # if model == 'easy':
   #     coeff = 0
   # elif model == 'hard':
   #     coeff = 1
   # n = 1 + (coeff * functions[0] / a)
   # dn_dr = coeff * a ** (-1)
   # dn_dtetta = 0
   # dr_dtau = functions[2]
   # dtetta_dtau = functions[3] / functions[0]
   # dPr_dtau = n * dn_dr + functions[3] ** 2 / functions[0]
   # dPt_dtau = (n * dn_dtetta - functions[2] * functions[3]) / functions[0]
   # return [dr_dtau, dtetta_dtau, dPr_dtau, dPt_dtau]

def system_2(functions, tao, a, b):
    z = functions[0]

    # Linear
    # n = 1 + z / 6371

    # Parabola
    hm = 300
    h0 = 60
    if z < h0:
        n = 1 + z / 6371
        # print("??? z =", z, "n =", n)
    else:
        fc = 6 # Крит частота в мгц
        nm = 1.24 * 10**10 * fc**2 # макс эл. конц (эл. на м^3)

        ym = hm - h0
        x = z - hm
        n = nm * (1 - (x*x) / (ym * ym)) / 1e6
        # print("z =", z, "n =", n)

    # n = profile.get_n(z)
    # print("N(", z, ") =", n)
    dn_dr = a
    dn_dtetta = 0
    dR_dtau = n * np.cos(functions[2])
    dT_dtau = n * np.sin(functions[2]) / functions[0]
    dA_dtau = (dn_dtetta * np.cos(functions[2]) - functions[0] * dn_dr * np.sin(functions[2]) - \
                                                                                n * np.sin(functions[2])) / functions[0]
    return [dR_dtau, dT_dtau, dA_dtau]

def analitycal_relations(tao, a, b, alpha):
    R0 = 6371
    def calc_r():
        top = R0 * (np.sin(alpha)**2 + np.exp(2*tao/R0) * (1 + np.cos(alpha))**2)
        bottom = 2 * np.exp(tao/R0) * (1 + np.cos(alpha))

        res = top/bottom - R0
        return  res

    def calc_tau():
        um = np.log10(1 - np.exp(tao/R0))
        vc = np.log10(1 - np.cos(alpha)-np.exp(tao/R0)*(1+np.cos(alpha)))
        lg = (um-vc)
        res = np.tan(alpha) * lg
        return res


    R = calc_r()
    tau = calc_tau()
    #R = (((np.sin(alpha) ** 2 + np.exp(2*tao) * (a + (a ** 2 - np.sin(alpha) ** 2) ** 1/2) ** 2)/(2 * np.exp(tao) * (a + (a ** 2 - np.sin(alpha) ** 2) ** 1/2))) - a) / b / 1e148
    # B =
    return [R, tau]

# def analitycal_relations(tao, a, alpha):
    # R =
    # B =
    # P_R =
    # P_B =
    # return [R, B, P_R, P_B]

class PropagationOfWaves:

    def __init__(self):
        #Directory
        self.directory_for_saving = path_creating_for_data_saving("traektorii")
        # Main parameters
        self.a = 1.5
        self.b = 0.005
        # Init states
        self.step_of_ALPHA = 0.1 * np.pi
        self.ALFA = np.arange(0.03 * np.pi, 0.5 * np.pi, self.step_of_ALPHA)
        self.r_k = 30000
        self.tao_init = 1
        self.r_init = 1
        self.tetta_init = 0
        # Configuration
        self.accuracy = 5
        self.quantity_of_points = 100
        # Arrays of parameters
        self.main_parameters = []

    def waves_propagation(self, graph_type):
        fig, ax = plt.subplots(figsize=(8, 5))

        for item_ALFA in tqdm.tqdm(self.ALFA, desc='ALFA init'):
            self.Pr = P_init_states(self.a, self.b, self.r_init, item_ALFA)[0]
            self.Pt = P_init_states(self.a, self.b, self.r_init, item_ALFA)[1]
            ALFA_index = str(self.ALFA.tolist().index(item_ALFA))
            self.data_complect["ALFA_" + ALFA_index] = dict()
            index = 0
            func_r = 0
            number_of_graph = 0
            while func_r <= self.r_k:
                tao = np.linspace(0, self.tao_init + index * self.accuracy, self.quantity_of_points)
                init_states = [self.r_init, self.tetta_init, self.Pr, self.Pt]
                init_states_2 = [self.r_init, self.tetta_init, item_ALFA]
                # functions = odeint(system, init_states, tao, args=(self.a, model))
                functions = odeint(system_2, init_states_2, tao, args=(self.a, self.b))
                index += 1
                func_r = functions[:, 0][-1]
                # print(func_r)
            else:
                self.data_complect["ALFA_" + ALFA_index]['r'] = functions[:, 0].tolist()
                self.data_complect["ALFA_" + ALFA_index]["t"] = functions[:, 1].tolist()
                self.data_complect["ALFA_" + ALFA_index]['Pr'] = functions[:, 2].tolist()
                # self.data_complect["ALFA_" + ALFA_index]['Pt'] = functions[:, 3].tolist()
                ax.plot(tao, functions[:, graph_type], label='Численное решение')
                # комент. нижние две строчки для аналитики
                anal_rel = analitycal_relations(tao, self.a, self.b, item_ALFA)
                ax.plot(tao, anal_rel[graph_type], label='Аналитическое решение')
                ax.legend()
        plt.xlabel("tau, sec")
        # plt.ylabel("r, km")
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.grid()
        filename = filename_in_2D(number_of_graph, graph_type)
        fig.savefig(self.directory_for_saving + '/' + filename, dpi=500)
        plt.show()
        number_of_graph += 1
        with open(self.directory_for_saving + '/Complect_of_data.json', "a") as file_for_saving:
            json.dump(self.data_complect, file_for_saving, indent=4)
        file = open(self.directory_for_saving + '/state.txt', "w+")
        file.write("SUCCESS!")
        file.close()

    def parameters_saving(self):
        self.data_complect = dict()
        self.data_complect["main_parameters"] = self.main_parameters
        self.data_complect["Alpha"] = [self.ALFA[0], self.ALFA[-1], self.step_of_ALPHA]
        self.data_complect["configuration"] = [self.accuracy, self.quantity_of_points]


def main():
    parser = argparse.ArgumentParser(description="my program")
    parser.add_argument("type_of_graph", type=int, help="Set type of graph")
    args = parser.parse_args()

    acquire_data = PropagationOfWaves()
    acquire_data.parameters_saving()
    acquire_data.waves_propagation(args.type_of_graph)


if __name__ == '__main__':
    main()
