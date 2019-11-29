import QtQuick.Window 2.2
import QtQuick.Controls 2.2
import QtQuick.Layouts 1.3


Window {
    visible: true

    RowLayout {
        Button {
            signal okButtonSignal
            objectName: "okButton"
            text: "Ok"
        }
        Button {
            signal cancelButtonSignal
            objectName: "cancelButton"
            text: "Cancel"
        }
    }
}
