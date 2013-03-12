import enaml
from enaml.qt.qt_application import QtApplication
from enaml.stdlib.sessions import SimpleSession
from necofs_model import OceanModel

if __name__ == '__main__':
    model = OceanModel()

    with enaml.imports():
        from ocean_view import Main

    app = QtApplication(
        [SimpleSession.factory('default', 'Ocean Model',
                           lambda: Main(model=model))])
    app.start_session('default')
    app.start()
    app.destroy()
