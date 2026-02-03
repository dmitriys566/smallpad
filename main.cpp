#include "pad.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    pad w;
    w.show();
    return a.exec();
}
