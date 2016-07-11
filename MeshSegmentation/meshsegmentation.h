#ifndef MESHSEGMENTATION_H
#define MESHSEGMENTATION_H

#include <QtWidgets/QMainWindow>
#include "ui_meshsegmentation.h"

#include <QFileDialog>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QString>
#include <QVBoxLayout>
#include <QVTKWidget.h>
#include <QWidget>

#include "QVTKModelViewer.h"

class MeshSegmentation : public QMainWindow
{
    Q_OBJECT

public:
    MeshSegmentation(QWidget *parent = 0);
    ~MeshSegmentation();

private:
    Ui::MeshSegmentationClass ui;

    int modelViewerLen;
    QString path;
    QWidget *widget;
    QVTKModelViewer *modelViewer;
    QGridLayout *mainLayout;
    QPushButton *openFileButton;

private slots:
    void SetModelFileName();
};

#endif // MESHSEGMENTATION_H
