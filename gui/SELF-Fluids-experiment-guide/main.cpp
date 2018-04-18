#include "self_fluids_guide.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    SELF_Fluids_Guide w;
    w.show();

    return a.exec();
}
