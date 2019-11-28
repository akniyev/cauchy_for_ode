from laguerre_functions import *
from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QPushButton, QVBoxLayout

counter = 1


def button1Pressed():
    global counter
    print(f"Hello {counter}!")
    counter += 1


def button2Pressed():
    print("Bye!")


if __name__ == "__main__":
    app = QApplication([])
    window = QWidget()
    label = QLabel('Press the button')
    button1 = QPushButton("Button1")
    button2 = QPushButton("Button2")

    button1.clicked.connect(button1Pressed)
    button2.clicked.connect(button2Pressed)

    vbox = QVBoxLayout()
    vbox.addWidget(label)
    vbox.addWidget(button1)
    vbox.addWidget(button2)

    window.setLayout(vbox)

    window.show()
    app.exec_()
