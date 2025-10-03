#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

class DrawingWidget;
class QPushButton;
class QHBoxLayout;
class QVBoxLayout;

class MainWindow : public QMainWindow
{
    Q_OBJECT
public:
    MainWindow(QWidget *parent = nullptr);

private:
    DrawingWidget *drawing;
    QWidget *central;
    QPushButton *runButton;
    QPushButton *clearButton;
};

#endif // MAINWINDOW_H
