#include "meshsegmentation.h"

#include <QFileDialog>

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
    segmentButton = new QPushButton(tr("Start Segmentation"));
    mergeButton = new QPushButton(tr("Open Merge Mode"));
    divideButton = new QPushButton(tr("Open Divide Mode"));
    clusterNumSlider = new QSlider(Qt::Horizontal);
    uiManager = NULL;
    
    /* =============================================================================== */

    /* ================================ Connections ================================ */
    
    connect(openFileButton, &QPushButton::released, this, &MeshSegmentation::SetModelFileName);
    connect(segmentButton, &QPushButton::released, this, &MeshSegmentation::StartSegmentation);
    connect(mergeButton, &QPushButton::released, this, &MeshSegmentation::SetMergeMode);
    connect(divideButton, &QPushButton::released, this, &MeshSegmentation::SetDivideMode);
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

    mainLayout->setMargin(10);
    mainLayout->addWidget(modelViewer, 0, 0, modelViewerLen, modelViewerLen);
    mainLayout->addWidget(openFileButton, 0, modelViewerLen, 1, 4);
    mainLayout->addWidget(segmentButton, 1, modelViewerLen, 1, 4);
    mainLayout->addWidget(mergeButton, 2, modelViewerLen, 1, 4);
    mainLayout->addWidget(divideButton, 3, modelViewerLen, 1, 4);
    mainLayout->addWidget(clusterNumSlider, modelViewerLen, 0, 1, modelViewerLen);

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

    mergeButton->setText(tr("Open Merge Mode"));
    divideButton->setText(tr("Open Divide Mode"));
    modelViewer->style->isMergeButtonDown = false;
    modelViewer->style->isDivideButtonDown = false;
}

void MeshSegmentation::StartSegmentation() {
    double dur[4], *dur_2;
    clock_t begin, end;
    clock_t totalBegin, totalEnd;

    cout << "=============================================" << endl;
    cout << "Step 1 : Converting model to dual graph . . ." << endl;
    begin = clock();
    totalBegin = begin;
    uiManager->ConvertPolydataToDualGraph();
    end = clock();
    dur[0] = (end - begin) * 1.0 / CLOCKS_PER_SEC;

    cout << "Step 2 : Automatic selecting seeds . . ." << endl;
    begin = clock();
    uiManager->AutomaticSelectSeeds(seedCnt, modelViewer->GetInteractor());
    end = clock();
    dur[1] = (end - begin) * 1.0 / CLOCKS_PER_SEC;

    cout << "Step 3 : Segmenting . . ." << endl;
    begin = clock();
    dur_2 = uiManager->StartSegmentation(modelViewer->GetInteractor());
    end = clock();
    dur[2] = (end - begin) * 1.0 / CLOCKS_PER_SEC;

    cout << "Step 4 : Merging clusters . . ." << endl;
    begin = clock();
    uiManager->MergeClusters(seedCnt, modelViewer->GetInteractor());
    end = clock();
    totalEnd = end;
    dur[3] = (end - begin) * 1.0 / CLOCKS_PER_SEC;
    cout << "=============================================" << endl;

    cout << "Runtime analysis : " << endl;
    double totalDur = (totalEnd - totalBegin) * 1.0 / CLOCKS_PER_SEC;
    printf("time 1 : \t%.3lf\t\t%.1lf%%\n", dur[0], dur[0] * 100.0 / totalDur);
    printf("time 2 : \t%.3lf\t\t%.1lf%%\n", dur[1], dur[1] * 100.0 / totalDur);
    printf("time 3 : \t%.3lf\t\t%.1lf%%\n", dur[2], dur[2] * 100.0 / totalDur);

    printf("- time 3.1 : \t%.3lf\t\t%.1lf%%\n", dur_2[0], dur_2[0] * 100.0 / totalDur);
    printf("- time 3.2 : \t%.3lf\t\t%.1lf%%\n", dur_2[1], dur_2[1] * 100.0 / totalDur);
    printf("- time 3.3 : \t%.3lf\t\t%.1lf%%\n", dur_2[2], dur_2[2] * 100.0 / totalDur);
    printf("- time 3.4 : \t%.3lf\t\t%.1lf%%\n", dur_2[3], dur_2[3] * 100.0 / totalDur);
    printf("- time 3.5 : \t%.3lf\t\t%.1lf%%\n", dur_2[4], dur_2[4] * 100.0 / totalDur);

    printf("time 4 : \t%.3lf\t\t%.1lf%%\n", dur[3], dur[3] * 100.0 / totalDur);
    printf("total : \t%.3lf\n", totalDur);
    cout << "=============================================" << endl;

    clusterNumSlider->setValue(seedCnt);
    clusterNumSlider->setDisabled(false);

    delete[] dur_2;
}

void MeshSegmentation::SetMergeMode() {
    bool& tmp = modelViewer->style->isMergeButtonDown;

    if (tmp) {
        mergeButton->setText(tr("Open Merge Mode"));
        tmp = false;
    } else {
        if (clusterNumSlider->isEnabled()) {
            clusterNumSlider->setDisabled(true);
            uiManager->ConfirmClusterSegmentation(seedCnt, currentClusterNum);
        }
        mergeButton->setText(tr("Close Merge Mode"));
        tmp = true;
    }
}

void MeshSegmentation::SetDivideMode() {
    bool& tmp = modelViewer->style->isDivideButtonDown;

    if (tmp) {
        divideButton->setText(tr("Open Divide Mode"));
        tmp = false;
    } else {
        if (clusterNumSlider->isEnabled()) {
            clusterNumSlider->setDisabled(true);
            uiManager->ConfirmClusterSegmentation(seedCnt, currentClusterNum);
        }
        divideButton->setText(tr("Close Divide Mode"));
        tmp = true;
    }
}

void MeshSegmentation::SetClusterNum(int k) {
    currentClusterNum = k;
}

void MeshSegmentation::DisplayCluster() {
    uiManager->SetClusterStep(seedCnt, currentClusterNum, modelViewer->GetInteractor());

}