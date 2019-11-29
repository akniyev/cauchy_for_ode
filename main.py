from PyQt5.QtWidgets import *
from PyQt5.QtGui import *

from laguerre_functions import *
from main_window import Ui_MainWindow


class CauchySolverWindow(Ui_MainWindow):
    def __init__(self):
        self.window = QMainWindow()
        self.setupUi(self.window)
        self.setup()

    def setup(self):
        validator = QIntValidator(0, 1000)
        self.polynomialOrder.setValidator(validator)
        self.calculateRootsButton.clicked.connect(self.draw_laguerre_with_roots)

    def draw_laguerre_with_roots(self):
        if not self.polynomialOrder.hasAcceptableInput():
            return
        n = int(self.polynomialOrder.text())
        print(n)


if __name__ == "__main__":
    # Setup GUI
    app = QApplication([])

    solver_ui = CauchySolverWindow()
    solver_ui.window.show()

    app.exec_()
