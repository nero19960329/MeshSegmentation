#include "meshsegmentation.h"

#include <iostream>

using namespace std;

MeshSegmentation::MeshSegmentation(QWidget *parent) : QMainWindow(parent) {
    ui.setupUi(this);

    /* ================================ Initliazation ================================ */

    seedCnt = 64;
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
    uiManager = NULL;
    colors = NULL;
    
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

    /* ================================ Other operations ================================ */

    clusterNumSlider->setMinimum(2);
    clusterNumSlider->setMaximum(seedCnt);
    clusterNumSlider->setValue(seedCnt);
    clusterNumSlider->setTickInterval(1);
    clusterNumSlider->setTickPosition(QSlider::TicksBelow);
    clusterNumSlider->setDisabled(true);

    currentColorLabel->setStyleSheet("background-color : rgb(0, 0, 0);");

    mainLayout->setMargin(10);
    mainLayout->addWidget(modelViewer, 0, 0, modelViewerLen, modelViewerLen);
    mainLayout->addWidget(openFileButton, 0, modelViewerLen, 1, 4);
    mainLayout->addWidget(segmentButton, 1, modelViewerLen, 1, 4);
    mainLayout->addWidget(resetButton, 2, modelViewerLen, 1, 4);
    mainLayout->addWidget(clusterNumSlider, modelViewerLen, 0, 1, modelViewerLen);
    /*for (int i = 0; i < 16; ++i) {
        mainLayout->addWidget(colorButtons[i], modelViewerLen + 1 + (i / 8), 6 + i % 8);
    }
    mainLayout->addWidget(currentColorLabel, modelViewerLen + 1, 4, 2, 2)*/;

    widget->setLayout(mainLayout);

    /* ================================================================================== */

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

    if (uiManager) {
        delete uiManager;
    }
    uiManager = modelViewer->RenderModel(path.toStdString());
    colors = uiManager->GetColors();

    resetButton->setText(tr("Set Reset Mode Open"));
    modelViewer->style->isHideButtonDown = false;

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
    double dur[3];
    clock_t begin, end;

    cout << "Step 1 : Automatic selecting seeds . . ." << endl;
    begin = clock();
    uiManager->AutomaticSelectSeeds(seedCnt, modelViewer->GetInteractor());
    end = clock();
    dur[0] = (end - begin) * 1.0 / CLOCKS_PER_SEC;

    cout << "Step 2 : Segmenting . . ." << endl;
    begin = clock();
    uiManager->StartSegmentation(modelViewer->GetInteractor());
    end = clock();
    dur[1] = (end - begin) * 1.0 / CLOCKS_PER_SEC;

    cout << "Step3 : Merging clusters . . ." << endl;
    begin = clock();
    uiManager->MergeClusters(seedCnt, modelViewer->GetInteractor());
    end = clock();
    dur[2] = (end - begin) * 1.0 / CLOCKS_PER_SEC;

    cout << "time 1 : " << dur[0] << "s" << endl;
    cout << "time 2 : " << dur[1] << "s" << endl;
    cout << "time 3 : " << dur[2] << "s" << endl;

    clusterNumSlider->setValue(seedCnt);
    clusterNumSlider->setDisabled(false);
}

void MeshSegmentation::SetResetMode() {
    bool& tmp = modelViewer->style->isHideButtonDown;

    if (tmp) {
        resetButton->setText(tr("Set Reset Mode Open"));
        tmp = false;
    } else {
        clusterNumSlider->setDisabled(true);
        uiManager->ConfirmClusterSegmentation(seedCnt, currentClusterNum);
        resetButton->setText(tr("Set Reset Mode Close"));
        tmp = true;
    }
}

void MeshSegmentation::SetClusterNum(int k) {
    currentClusterNum = k;
}

void MeshSegmentation::DisplayCluster() {
    uiManager->SetClusterStep(seedCnt, currentClusterNum, modelViewer->GetInteractor());

}