#include "mainwindow.h"
#include "drawingwidget.h"
#include <QPushButton>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QWidget>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    central = new QWidget(this);
    setCentralWidget(central);

    drawing = new DrawingWidget(this);
    drawing->setMinimumSize(640, 480);

    runButton = new QPushButton("Run Convex Hull", this);
    clearButton = new QPushButton("Clear", this);

    QHBoxLayout *hButtons = new QHBoxLayout;
    hButtons->addWidget(runButton);
    hButtons->addWidget(clearButton);
    hButtons->addStretch();

    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(drawing);
    layout->addLayout(hButtons);

    central->setLayout(layout);

    connect(runButton, &QPushButton::clicked, drawing, &DrawingWidget::runBothAlgorithms);
    connect(clearButton, &QPushButton::clicked, drawing, &DrawingWidget::clearAll);
}
