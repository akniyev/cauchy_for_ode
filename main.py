from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtQml import QQmlApplicationEngine
from PyQt5.QtCore import QUrl, QObject
import matplotlib.pyplot as plt
from matplotlib import style

style.use('fivethirtyeight')

import sys

from CauchySolver import CauchySolver
from laguerre_functions import *
from main_window import Ui_MainWindow


class CauchySolverWindow(Ui_MainWindow):
    def __init__(self):
        self.solver = CauchySolver()
        self.window = QMainWindow()
        self.setupUi(self.window)
        self.setup()

    def setup(self):
        orderValidator = QIntValidator(0, 1000)
        nPartValidator = QIntValidator(0, 1000)
        alphaValidator = QDoubleValidator()
        abValidator = QDoubleValidator()
        self.txtOrder.setValidator(orderValidator)
        self.txtAlpha.setValidator(alphaValidator)
        self.txtA.setValidator(abValidator)
        self.txtB.setValidator(abValidator)
        self.txtNPart.setValidator(nPartValidator)
        self.calculateRootsButton.clicked.connect(self.draw_laguerre_with_roots)
        self.btnResetIteration.clicked.connect(self.reset_iteration)
        self.btnNextIteration.clicked.connect(self.next_iteration)

    def set_parameters(self):
        if not self.txtOrder.hasAcceptableInput() or not self.txtAlpha.hasAcceptableInput() or\
                not self.txtNPart.hasAcceptableInput() or not self.txtA.hasAcceptableInput() or\
                not self.txtB.hasAcceptableInput():
            return False

        n = int(self.txtOrder.text())
        n_part = int(self.txtNPart.text())
        alpha = float(self.txtAlpha.text())
        a = float(self.txtA.text())
        b = float(self.txtB.text())

        self.solver.set_parameters(n, alpha, n_part, a, b)

        return True

    def draw_laguerre_with_roots(self):
        if not self.set_parameters():
            return

        roots_xs = self.solver.roots()
        a = min(roots_xs) - 1
        b = max(roots_xs) + 1
        xs_polynomial, ys_polynomial = self.solver.laguerre_on_grid(a, b, 3000)
        roots_ys = [0 for _ in roots_xs]

        plt.clf()
        plt.plot(roots_xs, roots_ys, 'ro', label="Roots")
        plt.plot(xs_polynomial, ys_polynomial, label="Laguerre Polynomial", linewidth=1)
        plt.show()
        plt.close()

    def reset_iteration(self):
        if not self.set_parameters():
            return

        self.solver.reset_iterations()
        self.btnNextIteration.setEnabled(True)

    def next_iteration(self):
        self.solver.next_iterations()
        xs_solution, ys_solution = self.solver.get_solution(2000)

        plt.clf()
        plt.plot(xs_solution, ys_solution, label="Solution", linewidth=1)
        plt.show()
        plt.close()



def run():
    app = QApplication(sys.argv)
    engine = QQmlApplicationEngine()
    engine.load('qt_gui/main_window.qml')
    if not engine.rootObjects():
        return -1

    win = engine.rootObjects()[0]

    okButton = win.findChild(QObject, "okButton")
    okButton.clicked.connect(action)

    return app.exec_()


def action():
    print("Action!")


if __name__ == "__main__":
    # Setup GUI old school
    app = QApplication([])

    solver_ui = CauchySolverWindow()
    solver_ui.window.show()

    app.exec_()


