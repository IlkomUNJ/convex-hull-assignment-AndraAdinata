#include "drawingwidget.h"
#include <QPainter>
#include <QMouseEvent>
#include <algorithm>
#include <set>
#include <cmath>

DrawingWidget::DrawingWidget(QWidget *parent)
    : QWidget(parent),
      iterationsFast(0),
      iterationsSlow(0)
{
    setAutoFillBackground(true);
    setBackgroundRole(QPalette::Base);
}

void DrawingWidget::mousePressEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton) {
        QPointF p = event->position(); // Qt 6: use position(); for Qt5 use event->posF()
#ifdef QT_VERSION_MAJOR
#if QT_VERSION_MAJOR < 6
        p = event->posF();
#endif
#endif
        points.append(p);
        // reset hulls until user presses Run again
        hullFast.clear();
        hullSlow.clear();
        iterationsFast = iterationsSlow = 0;
        update();
    }
}

void DrawingWidget::paintEvent(QPaintEvent * /*event*/)
{
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);

    // clear
    p.fillRect(rect(), Qt::white);

    // draw points
    p.setPen(Qt::black);
    for (const QPointF &pt : points) {
        p.drawEllipse(pt, 4, 4);
    }

    // draw slow hull in red (opaque)
    if (!hullSlow.isEmpty()) {
        QPen pen(Qt::red, 2);
        p.setPen(pen);
        QPolygonF poly;
        for (int idx : hullSlow) poly << points[idx];
        if (poly.size() > 1) {
            p.drawPolygon(poly);
            // close polygon
            p.drawLine(poly.last(), poly.first());
        }
    }

    // draw fast hull in blue dashed (overlay)
    if (!hullFast.isEmpty()) {
        QPen pen(Qt::blue, 2, Qt::DashLine);
        p.setPen(pen);
        QPolygonF poly;
        for (int idx : hullFast) poly << points[idx];
        if (poly.size() > 1) {
            p.drawPolygon(poly);
            p.drawLine(poly.last(), poly.first());
        }
    }

    // iteration counts and instructions
    p.setPen(Qt::black);
    QFont f = p.font();
    f.setPointSize(10);
    p.setFont(f);

    QString info = QString("Points: %1\nFast (Graham) iterations: %2\nSlow (brute) iterations: %3\n\nLeft click to add points.")
            .arg(points.size())
            .arg(iterationsFast)
            .arg(iterationsSlow);
    p.drawText(8, 16, info);
}

void DrawingWidget::clearAll()
{
    points.clear();
    hullFast.clear();
    hullSlow.clear();
    iterationsFast = iterationsSlow = 0;
    update();
}

double DrawingWidget::cross(const QPointF &o, const QPointF &a, const QPointF &b)
{
    return (a.x() - o.x()) * (b.y() - o.y()) - (a.y() - o.y()) * (b.x() - o.x());
}

double DrawingWidget::crossVec(const QPointF &a, const QPointF &b)
{
    return a.x()*b.y() - a.y()*b.x();
}

double DrawingWidget::dist2(const QPointF &a, const QPointF &b)
{
    double dx = a.x()-b.x();
    double dy = a.y()-b.y();
    return dx*dx + dy*dy;
}

void DrawingWidget::runBothAlgorithms()
{
    hullFast.clear();
    hullSlow.clear();
    iterationsFast = iterationsSlow = 0;

    if (points.size() < 3) {
        // nothing to do
        update();
        return;
    }

    computeGrahamScan(iterationsFast, hullFast);
    computeSlowConvexHull(iterationsSlow, hullSlow);

    update();
}

// Graham scan (fast) implementation. iterations counts comparisons and cross ops roughly.
void DrawingWidget::computeGrahamScan(qint64 &iterations, QVector<int> &outHull)
{
    iterations = 0;
    int n = points.size();
    // create indices
    QVector<int> idx(n);
    for (int i = 0; i < n; ++i) idx[i] = i;

    // find pivot = lowest y (and lowest x if tie)
    int pivot = 0;
    for (int i = 1; i < n; ++i) {
        ++iterations;
        if (points[i].y() < points[pivot].y() || (points[i].y() == points[pivot].y() && points[i].x() < points[pivot].x()))
            pivot = i;
    }

    QPointF p0 = points[pivot];
    // sort by angle wrt pivot, ties by distance
    std::sort(idx.begin(), idx.end(), [&](int a, int b){
        if (a == pivot) return true;
        if (b == pivot) return false;
        QPointF va = QPointF(points[a].x() - p0.x(), points[a].y() - p0.y());
        QPointF vb = QPointF(points[b].x() - p0.x(), points[b].y() - p0.y());
        double cr = crossVec(va, vb);
        ++iterations;
        if (std::abs(cr) < 1e-9) {
            // collinear: closer one first
            return dist2(p0, points[a]) < dist2(p0, points[b]);
        }
        return cr > 0; // a before b if left of b (i.e. smaller angle)
    });

    // remove duplicates with same angle keeping farthest (typical Graham variant) OR keep as is but handle in stack
    QVector<int> filtered;
    for (int i = 0; i < idx.size(); ++i) {
        if (!filtered.isEmpty() && idx[i] == pivot) continue;
        if (filtered.isEmpty()) { filtered.append(idx[i]); continue; }
        // if same angle as previous, keep the farthest
        QPointF A = points[filtered.last()];
        QPointF B = points[idx[i]];
        QPointF VA = QPointF(A.x()-p0.x(), A.y()-p0.y());
        QPointF VB = QPointF(B.x()-p0.x(), B.y()-p0.y());
        double cr = crossVec(VA, VB);
        ++iterations;
        if (std::abs(cr) < 1e-9) {
            // choose farthest
            if (dist2(p0, A) < dist2(p0, B))
                filtered.last() = idx[i];
            // else keep existing
        } else {
            filtered.append(idx[i]);
        }
    }

    if (filtered.size() < 3) {
        // everything collinear or too small
        outHull = filtered;
        return;
    }

    // stack
    QVector<int> st;
    st.append(filtered[0]);
    st.append(filtered[1]);

    for (int i = 2; i < filtered.size(); ++i) {
        while (st.size() >= 2) {
            int s1 = st[st.size()-2];
            int s2 = st[st.size()-1];
            int s3 = filtered[i];
            ++iterations;
            double cr = cross(points[s1], points[s2], points[s3]);
            if (cr <= 0) { // non-left turn -> pop (use <= to exclude collinear non-left)
                st.removeLast();
            } else {
                break;
            }
        }
        st.append(filtered[i]);
    }

    outHull = st;
}

// Slow brute-force convex hull: check every pair (i,j) if all points are on same side of line (i->j).
// We accumulate endpoints of valid edges and then order unique endpoints into polygon by centroid angle.
void DrawingWidget::computeSlowConvexHull(qint64 &iterations, QVector<int> &outHull)
{
    iterations = 0;
    int n = points.size();
    if (n < 3) return;

    std::set<std::pair<int,int>> edges;
    for (int i = 0; i < n; ++i) {
        for (int j = i+1; j < n; ++j) {
            bool pos = false, neg = false;
            for (int k = 0; k < n; ++k) {
                if (k==i || k==j) continue;
                ++iterations;
                double c = cross(points[i], points[j], points[k]);
                if (c > 1e-9) pos = true;
                else if (c < -1e-9) neg = true;
                if (pos && neg) break; // not an edge
            }
            if (!(pos && neg)) {
                // all points on one side -> (i,j) is an edge of convex hull (or collinear)
                edges.insert({i,j});
                edges.insert({j,i});
            }
        }
    }

    // collect unique vertices used in edges
    std::set<int> verts;
    for (auto &e : edges) {
        verts.insert(e.first);
        verts.insert(e.second);
    }

    if (verts.empty()) return;

    QVector<QPointF> pts;
    QVector<int> idx;
    for (int v : verts) {
        pts.append(points[v]);
        idx.append(v);
    }

    // order vertices by angle around centroid
    orderHullPoints(pts, idx);
    outHull = idx;
}

// orders idx (indices into original points) by computing centroid and sorting by angle
void DrawingWidget::orderHullPoints(const QVector<QPointF> &pts, QVector<int> &indices)
{
    int m = pts.size();
    if (m <= 1) return;
    // compute centroid
    double cx = 0, cy = 0;
    for (const QPointF &p : pts) { cx += p.x(); cy += p.y(); }
    cx /= m; cy /= m;
    QPointF center(cx, cy);

    // build vector of pairs (angle, index)
    std::vector<std::pair<double,int>> arr;
    arr.reserve(m);
    for (int i = 0; i < m; ++i) {
        QPointF v = QPointF(pts[i].x() - center.x(), pts[i].y() - center.y());
        double ang = std::atan2(v.y(), v.x());
        arr.push_back({ang, indices[i]});
    }
    std::sort(arr.begin(), arr.end(), [](const auto &a, const auto &b){ return a.first < b.first; });
    // rewrite indices in sorted order
    for (int i = 0; i < m; ++i) indices[i] = arr[i].second;
}
