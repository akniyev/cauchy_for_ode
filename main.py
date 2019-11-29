from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtQml import QQmlApplicationEngine
from PyQt5.QtCore import QUrl, QObject
import sys

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
    # app = QApplication([])
    #
    # solver_ui = CauchySolverWindow()
    # solver_ui.window.show()
    #
    # app.exec_()

    # Setup GUI with QML
    sys.exit(run())

