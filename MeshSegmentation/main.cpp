#include "meshsegmentation.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MeshSegmentation w;
    w.show();

    return a.exec();
}
