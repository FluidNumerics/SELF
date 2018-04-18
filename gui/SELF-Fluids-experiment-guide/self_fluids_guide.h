#ifndef SELF_FLUIDS_GUIDE_H
#define SELF_FLUIDS_GUIDE_H

#include <QMainWindow>

namespace Ui {
class SELF_Fluids_Guide;
}

class SELF_Fluids_Guide : public QMainWindow
{
    Q_OBJECT

public:
    explicit SELF_Fluids_Guide(QWidget *parent = 0);
    ~SELF_Fluids_Guide();

private:
    Ui::SELF_Fluids_Guide *ui;
};

#endif // SELF_FLUIDS_GUIDE_H
