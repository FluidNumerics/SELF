#include "self_fluids_guide.h"
#include "ui_self_fluids_guide.h"

SELF_Fluids_Guide::SELF_Fluids_Guide(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::SELF_Fluids_Guide)
{
    ui->setupUi(this);
}

SELF_Fluids_Guide::~SELF_Fluids_Guide()
{
    delete ui;
}
