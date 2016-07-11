#include "meshsegmentation.h"

MeshSegmentation::MeshSegmentation(QWidget *parent) : QMainWindow(parent) {
    ui.setupUi(this);

    modelViewerLen = 8;
    widget = new QWidget(this);
    modelViewer = new QVTKModelViewer(this);
    mainLayout = new QGridLayout;
    openFileButton = new QPushButton(tr("Open STL File"));
    
    connect(openFileButton, SIGNAL(released()), this, SLOT(SetModelFileName()));

    mainLayout->addWidget(modelViewer, 0, 0, modelViewerLen, modelViewerLen);
    mainLayout->addWidget(openFileButton, 0, modelViewerLen);

    widget->setLayout(mainLayout);
    this->setCentralWidget(widget);
    this->setFixedSize(900, 800);
}

MeshSegmentation::~MeshSegmentation() {

}

void MeshSegmentation::SetModelFileName() {
    path = QFileDialog::getOpenFileName(this, tr("Open model file"), tr("../../objects/"), tr("STL Model Files(*.stl)"));
    if (path.isEmpty()) {
        return;
    }
    modelViewer->RenderModel(path.toStdString());
}