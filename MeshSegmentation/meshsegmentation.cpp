#include "meshsegmentation.h"

#include <iostream>

using namespace std;

MeshSegmentation::MeshSegmentation(QWidget *parent) : QMainWindow(parent) {
    ui.setupUi(this);

    /* ================================ Initliazation ================================ */

    seedCnt = 32;
    colorNum = 16;
    modelViewerLen = 18;
    widget = new QWidget(this);
    modelViewer = new QVTKModelViewer(this);
    mainLayout = new QGridLayout;
    openFileButton = new QPushButton(tr("Open STL File"));
    colorButtons = new QPushButton*[16];
    segmentButton = new QPushButton(tr("Start Segmentation"));
    resetButton = new QPushButton(tr("Set Reset Mode Open"));
    currentColorLabel = new QLabel(tr(""));
    for (int i = 0; i < 16; ++i) {
        colorButtons[i] = new QPushButton(tr(""));
    }
    clusterNumSlider = new QSlider(Qt::Horizontal);
    
    /* =============================================================================== */

    /* ================================ Connections ================================ */
    
    connect(openFileButton, &QPushButton::released, this, &MeshSegmentation::SetModelFileName);
    for (int i = 0; i < 16; ++i) {
        connect(colorButtons[i], &QPushButton::released, this, [=] {
            SetBrushColor(i);
        });
    }
    connect(segmentButton, &QPushButton::released, this, &MeshSegmentation::StartSegmentation);
    connect(resetButton, &QPushButton::released, this, &MeshSegmentation::SetResetMode);
    connect(clusterNumSlider, SIGNAL(valueChanged(int)), this, SLOT(SetClusterNum(int)));
    connect(clusterNumSlider, &QSlider::sliderReleased, this, &MeshSegmentation::DisplayCluster);

    /* ============================================================================= */

    /* ================================ Other actions ================================ */

    clusterNumSlider->setMinimum(2);
    clusterNumSlider->setMaximum(seedCnt);
    clusterNumSlider->setValue(seedCnt);
    clusterNumSlider->setTickInterval(1);

    currentColorLabel->setStyleSheet("background-color : rgb(0, 0, 0);");

    mainLayout->setMargin(10);
    mainLayout->addWidget(modelViewer, 0, 0, modelViewerLen, modelViewerLen);
    mainLayout->addWidget(openFileButton, 0, modelViewerLen, 1, 4);
    mainLayout->addWidget(segmentButton, 1, modelViewerLen, 1, 4);
    mainLayout->addWidget(resetButton, 2, modelViewerLen, 1, 4);
    mainLayout->addWidget(clusterNumSlider, 3, modelViewerLen, 1, 4);
    for (int i = 0; i < 16; ++i) {
        mainLayout->addWidget(colorButtons[i], modelViewerLen + (i / 8), 6 + i % 8);
    }
    mainLayout->addWidget(currentColorLabel, modelViewerLen, 4, 2, 2);

    widget->setLayout(mainLayout);

    /* =============================================================================== */

    this->setCentralWidget(widget);
    this->setFixedSize(900, 800);
}

MeshSegmentation::~MeshSegmentation() {

}

void MeshSegmentation::SetModelFileName() {
    path = QFileDialog::getOpenFileName(this, tr("Open STL file"), tr("../../objects/"), tr("STL Model Files(*.stl)"));
    if (path.isEmpty()) {
        return;
    }
    uiManager = modelViewer->RenderModel(path.toStdString());
    colors = uiManager->GetColors();

    QString tmpStr;
    for (int i = colorNum - 1; i >= 0; --i) {
        unsigned char *color = colors[i];
        tmpStr = QString("background-color : rgb(%1, %2, %3);").arg(QString::number(color[0]), QString::number(color[1]), QString::number(color[2]));
        colorButtons[i]->setStyleSheet(tmpStr);
    }
    currentColorLabel->setStyleSheet(tmpStr);
}

void MeshSegmentation::SetBrushColor(int k) {
    uiManager->SetCurrentCluster(k);
    QString tmpStr = QString("background-color : rgb(%1, %2, %3);").arg(QString::number(colors[k][0]), QString::number(colors[k][1]), QString::number(colors[k][2]));
    currentColorLabel->setStyleSheet(tmpStr);
}

void MeshSegmentation::StartSegmentation() {
    uiManager->AutomaticSelectSeeds(seedCnt, modelViewer->GetInteractor());
    uiManager->StartSegmentation(modelViewer->GetInteractor());
    uiManager->MergeClusters(seedCnt, modelViewer->GetInteractor());
    clusterNumSlider->setValue(seedCnt);
}

void MeshSegmentation::SetResetMode() {
    bool& tmp = modelViewer->style->isHideButtonDown;

    if (tmp) {
        resetButton->setText(tr("Set Reset Mode Open"));
        tmp = false;
    } else {
        resetButton->setText(tr("Set Reset Mode Close"));
        tmp = true;
    }
}

void MeshSegmentation::SetClusterNum(int k) {
    currentClusterNum = k;
}

void MeshSegmentation::DisplayCluster() {
    uiManager->SetClusterNum(seedCnt, currentClusterNum, modelViewer->GetInteractor());

}