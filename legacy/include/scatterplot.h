/*
 ScatterPlot, derived from IncrementalPlot. Includes a plot zoomer.
 */

#ifndef _SCATTERPLOT_H_
#define _SCATTERPLOT_H_ 1

#include "incrementalplot.h"
#include <qwt-qt3/qwt_scale_widget.h>
#include <qwt-qt3/qwt_scale_draw.h>

class QTimer;

class ScatterPlot : public IncrementalPlot
{
Q_OBJECT

public:
  ScatterPlot(QWidget *parent);
  virtual
  ~ScatterPlot()
  {
  }
  ;
  //    virtual QSize sizeHint() const;

  signals:

public slots:
  void
  clear();

private slots:
  void
  appendPoint(double x, double y);
  void
  appendPoint(double x, double y, QwtSymbol::Style s);

};

#endif // _SCATTERPLOT_H_
