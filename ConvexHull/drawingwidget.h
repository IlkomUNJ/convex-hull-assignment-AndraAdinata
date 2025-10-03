#ifndef DRAWINGWIDGET_H
#define DRAWINGWIDGET_H

#include <QWidget>
#include <QVector>
#include <QPointF>

class DrawingWidget : public QWidget
{
    Q_OBJECT
public:
    explicit DrawingWidget(QWidget *parent = nullptr);

    // called by mainwindow buttons
public slots:
    void runBothAlgorithms();
    void clearAll();

protected:
    void paintEvent(QPaintEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;

private:
    QVector<QPointF> points;

    // hulls
    QVector<int> hullFast; // indices into points (Graham)
    QVector<int> hullSlow; // indices into points (from brute edges then ordered)

    // iteration counts
    qint64 iterationsFast;
    qint64 iterationsSlow;

    // algorithm implementations
    void computeGrahamScan(qint64 &iterations, QVector<int> &outHull);
    void computeSlowConvexHull(qint64 &iterations, QVector<int> &outHull);

    // helpers
    static double cross(const QPointF &o, const QPointF &a, const QPointF &b);
    static double crossVec(const QPointF &a, const QPointF &b);
    static double dist2(const QPointF &a, const QPointF &b);
    void orderHullPoints(const QVector<QPointF> &pts, QVector<int> &indices);
};

#endif // DRAWINGWIDGET_H
