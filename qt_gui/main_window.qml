import QtQuick 2.3
import QtQuick.Window 2.2
import QtQuick.Controls 2.2
import QtQuick.Layouts 1.3


Window {
    id: window
    visible: true

    ColumnLayout {
        TextField {
            id: textField
            x: 40
            y: -74
            text: qsTr("Text Field")
            Layout.fillWidth: true
        }

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

}
